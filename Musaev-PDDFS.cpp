/**
 * Implementation of Musaev's Parallel Distributed Depth First Search algorithm
 * The program takes all the edges of the graph as input on STDIN
 * An edge is defined by the two nodes it connects, and the two nodes are separated by a space
 * Directed graphs are supported as input, but they are not supported by the algorithm.
 * For undirected graphs edges need to be specified in both directions.
 * Edges are expected in sorted order, sorted by the first (source) node
 * For a complete graph with two nodes {0,1} (and one edge), input will look as follows:
 * 0 1
 * 1 0
 * When the program is done, each node knows its parent and children in the DFS tree (if successful)
 * Each node prints its children list on termination such that proper execution can be verified
 * 
 * @author Gert Hartzema <g.hartzema@student.vu.nl>
 * Date: August 26, 2020
 * 
 * Compatible with openMPI and MPICH. MPICH was used for testing.
 */

#include <mpi.h>
#include <string>
#include <signal.h>
#include <set>
#include <thread>
#include <chrono>
#include <iostream>

#define DEBUG_PRINT false // toggle debug printing
#define DISCOVER_TYPE 1
#define REJECT_TYPE 2
#define TERMINATE_TYPE 3

/**
 * Array to string for debug printing
 */
std::string to_str(int n, int arr[])
{
    std::string out = "[";
    for (int i = 0; i < n; i++)
        out += std::to_string(arr[i]) + ", ";
    out += "]";
    return out;
}

/**
 * Set to string for debug printing
 */
std::string to_arr(std::set<int> elems)
{
    std::string out;
    out.push_back('[');
    for (auto elem : elems)
        out.append(std::to_string(elem) + ", ");
    out.push_back(']');
    return out;
}

void handle_sigint(int n)
{
    MPI_Finalize();
    MPI_Abort(MPI_COMM_WORLD, 1);
    exit(1);
}

/**
 * Takes edges sorted by source on stdin, returns graph communicator
 * 
 * @param rank The MPI process ID of the current process
 * @param size The amount of processes in the graph
 * @param comm The graph communicator that is written to
 * @return The passed graph communicator is loaded with topology data, each process becomes aware of its neighbours
 */
void load_graph(int rank, int size, MPI_Comm *comm)
{
    if (rank == 0)
    {
        int n = 0;
        int edges = 0;
        int sources[size];
        int degrees[size] = {0};
        int destinations[size * (size - 1)];
        MPI_Info info;
        MPI_Info_create(&info);
        std::string line;
        int source, dest;

        while (getline(std::cin, line))
        {
            sscanf(line.c_str(), "%i %i", &source, &dest);
            if (n == 0 || sources[n - 1] != source)
            {
                sources[n] = source;
                n++;
            }
            degrees[n - 1]++;
            destinations[edges] = dest;
            edges++;
        }
        MPI_Dist_graph_create(MPI_COMM_WORLD, n, sources, degrees, destinations, MPI_UNWEIGHTED, info, false, comm);
    }
    else
    {
        MPI_Info info;
        MPI_Info_create(&info);
        MPI_Dist_graph_create(MPI_COMM_WORLD, 0, NULL, NULL, NULL, MPI_UNWEIGHTED, info, false, comm);
    }
}

/**
 * Send DISCOVER message to set of children
 * 
 * @param dests The destinations to send DISCOVER to
 * @param path The path vector at the current node
 * @param path_length The size of the path vector
 * @param comm The communicator to write on
 * @return DISCOVER messages with the path vector (with destination ID appended) written to each destination channel
 */
void send_discover(std::set<int> dests, int path[], int path_length, MPI_Comm comm)
{
    MPI_Request request[dests.size()];

    path_length++;
    int i = 0;
    for (auto dest : dests)
    {
        path[path_length - 1] = dest;
        MPI_Issend(path, path_length, MPI_INT, dest, DISCOVER_TYPE, comm, &request[i]);
        i++;
    }
}

/**
 * send_discover() overload for ease of use with a single destination
 */
void send_discover(int dest, int path[], int path_length, MPI_Comm comm)
{
    std::set<int> dest_wrapper;
    dest_wrapper.insert(dest);
    send_discover(dest_wrapper, path, path_length, comm);
}

/**
 * Send REJECT message to destination
 * 
 * @param dest The destination to send REJECT to
 * @param comm The communicator to write on
 * @return REJECT message written to the destination channel
 */
void send_reject(int dest, MPI_Comm comm)
{
    MPI_Request request;
    MPI_Issend(NULL, 0, MPI_INT, dest, REJECT_TYPE, comm, &request);
}

/**
 * Send TERMINATE message to destination
 * 
 * @param parent The destination to send TERMINATE to
 * @param comm The communicator to write on
 * @return TERMINATE message written to the destination channel
 */
void send_terminate(int parent, MPI_Comm comm, int rank = -1)
{
    MPI_Request request;
    MPI_Issend(NULL, 0, MPI_INT, parent, TERMINATE_TYPE, comm, &request);
}

/**
 * Calculate path order, which path is 'more depth-first'
 * 
 * @param path_length_1 The length of the first path
 * @param path1 The first path vector
 * @param path_length2 The length of the second path
 * @param path2 The second path vector
 * @return -1 for path1 >_{df} path2, 1 for path1 <_{df} path2, 0 for path1 \subset_{df} path2 
 */
int path_order(int path_length_1, int path1[], int path_length_2, int path2[])
{
    for (int i = 0; i < std::min(path_length_1, path_length_2); i++)
    {
        if (path1[i] < path2[i])
            return -1;
        else if (path1[i] > path2[i])
            return 1;
    }
    return 0;
}

int main()
{

    if (!DEBUG_PRINT)
        freopen("/dev/null", "w", stderr); // send stderr to dev/null
    struct sigaction sigIntHandler;

    sigIntHandler.sa_handler = handle_sigint;
    sigemptyset(&sigIntHandler.sa_mask);
    sigIntHandler.sa_flags = 0;

    sigaction(SIGINT, &sigIntHandler, NULL); // handle SIGINT

    int world_size, world_rank;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    // containers for graph communicator topology info
    int sources[world_size], dests[world_size], sweights[world_size],
        dweights[world_size], outdegree, indegree, weighted;
    MPI_Comm local;

    load_graph(world_rank, world_size, &local);

    MPI_Dist_graph_neighbors(local, world_size, sources, sweights, world_size, dests, dweights);
    MPI_Dist_graph_neighbors_count(local, &indegree, &outdegree, &weighted);

    // containers for algorithm functionality
    std::set<int> children;
    for (int i = 0; i < outdegree; i++)
        children.insert(dests[i]); // all neighbours are added to children list, parent is removed upon first discovery

    std::set<int> terminated_children;
    int graph_path[world_size];
    int path_length = 0;
    int parent = -1;
    int recv_path_length;
    int recv_graph_path[world_size];
    bool mounted = false;
    bool is_parent_rejected = false;
    int msgct = 0;

    if (DEBUG_PRINT)
        freopen(("./debug_log/" + std::to_string(world_rank)).c_str(), "w+", stderr); // send debugprints to files, debug info from different processes is separated

    if (world_rank == 0) // If current process is root, start the algorithm
    {
        mounted = true;
        graph_path[0] = 0;
        path_length = 1;
        send_discover(children, graph_path, path_length, local);
    }

    while (1)
    {
        MPI_Status status;
        MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, local, &status);

        msgct++;

        switch (status.MPI_TAG)
        {
        case DISCOVER_TYPE:
            std::cerr << "[" << world_rank << "]:"
                      << "Got DISCOVER msg FROM: " << status.MPI_SOURCE << "\t\t";
            {
                if (!mounted) // Node is not yet attached to DFS tree
                {
                    std::cerr << "For the first time " << std::endl;

                    mounted = true;
                    parent = status.MPI_SOURCE;
                    children.erase(parent);
                    MPI_Get_count(&status, MPI_INT, &path_length);

                    MPI_Recv(&graph_path, path_length, MPI_INT, parent, DISCOVER_TYPE, local, MPI_STATUS_IGNORE);

                    send_discover(children, graph_path, path_length, local);
                }
                else if (status.MPI_SOURCE == parent)
                { // sometimes you may get the same path you already have, ignore this.
                    MPI_Get_count(&status, MPI_INT, &recv_path_length);
                    MPI_Recv(recv_graph_path, recv_path_length, MPI_INT, status.MPI_SOURCE, DISCOVER_TYPE, local, MPI_STATUS_IGNORE); // take message off the queue
                    std::cerr << "From parent with path: " << to_str(recv_path_length, recv_graph_path) << std::endl;
                    if (path_order(path_length, graph_path, recv_path_length, recv_graph_path) == 1)
                    {
                        std::copy(recv_graph_path, recv_graph_path + recv_path_length, graph_path); // sometimes the path from parent is better, update own path to save some work
                        path_length = recv_path_length;
                    }
                }
                else // Node is already part of DFS tree
                {
                    MPI_Get_count(&status, MPI_INT, &recv_path_length);
                    MPI_Recv(&recv_graph_path, recv_path_length, MPI_INT, status.MPI_SOURCE, DISCOVER_TYPE, local, &status);
                    std::cerr << "WITH PATH: " << to_str(recv_path_length, recv_graph_path) << std::endl;

                    int order = path_order(path_length, graph_path, recv_path_length, recv_graph_path);
                    if (order == 1) // recv path >df curr path: update own path, update parent, send DISCOVER to old parent
                    {
                        std::copy(recv_graph_path, recv_graph_path + recv_path_length, graph_path);
                        path_length = recv_path_length;

                        if (!is_parent_rejected)
                        {
                            children.insert(parent);                               // old parent becomes child
                            send_discover(parent, graph_path, path_length, local); // send updated path to old parent
                        }
                        parent = status.MPI_SOURCE; // change parent
                        is_parent_rejected = false;
                        children.erase(parent); // remove new parent from children
                    }
                    else if (order == 0) // curr path \subsetdf recv path: remove sender or t from children, send reject to sender
                    {
                        int t = recv_graph_path[path_length]; // link t is the other link that connects p to the loop
                        if (t < status.MPI_SOURCE)
                        { // if the path through t is more df, sender needs to be rejected
                            children.erase(status.MPI_SOURCE);
                            send_reject(status.MPI_SOURCE, local);
                        }
                        else
                        { // t is rejected
                            children.erase(t);
                            send_reject(t, local);
                        }
                    }
                    else if (order == -1)
                    { // curr path more df than recv path, send path back to sender
                        send_discover(status.MPI_SOURCE, graph_path, path_length, local);
                    }
                }
                break;
            }
        case REJECT_TYPE:
        {
            std::cerr << "[" << world_rank << "]:"
                      << "Got REJECT msg FROM: " << status.MPI_SOURCE << "\t\t" << msgct << std::endl;
            MPI_Recv(NULL, 0, MPI_INT, status.MPI_SOURCE, REJECT_TYPE, local, MPI_STATUS_IGNORE);

            if (status.MPI_SOURCE == parent)
            {
                is_parent_rejected = true;
            }
            else
            {
                children.erase(status.MPI_SOURCE);
            }
            break;
        }
        case TERMINATE_TYPE:
        {
            std::cerr << "[" << world_rank << "]:"
                      << "Got TERMINATE msg FROM: " << status.MPI_SOURCE << "\t\t" << msgct << std::endl;
            MPI_Recv(NULL, 0, MPI_INT, status.MPI_SOURCE, TERMINATE_TYPE, local, MPI_STATUS_IGNORE);

            terminated_children.insert(status.MPI_SOURCE);
            break;
        }
        }
        std::cerr << "[" << world_rank << "]: "
                  << " parent: " << parent << " curr-path[" << to_str(path_length, graph_path) << "]"
                  << " with len:" << std::to_string(path_length) << " children: " << to_arr(children) << " terminated: " << to_arr(terminated_children) << " parent-rejected?: " << is_parent_rejected << std::endl
                  << std::endl;
        if (terminated_children.size() == children.size()) // all children have trminated
        {
            if (world_rank != 0)
            {
                send_terminate(parent, local, world_rank);
            }

            std::string out = "[" + std::to_string(world_rank) + "]:\t DONE - Children: " + to_arr(children) + "\t\t" + std::to_string(msgct) + "\n";
            std::cout << out;
            MPI_Finalize();
            return 0; // stop the infinite loop and finalise
        }
    }
}

/**
* Generate random graphs as input for the PPDFS algorithm
* Random generation using the erdos renyi algorithm
* The program takes two command line arguments
* The first argument is the number of nodes
* The second argument is the chance that two nodes are connected (between 0 and 1)
*/

#include <iostream>
#include <vector>
#include <bits/stdc++.h>

using namespace std;

int main(int argc, char *argv[])
{
    srand(time(NULL));
    int n = stoi(argv[1]);
    float frac = stof(argv[2]);
    vector<pair<int, int>> edges;
    // cout << rand() << " " << rand() << " " << rand();
    for (int i = 0; i < n; i++)
    {
        for (int j = i + 1; j < n; j++)
        {
            if ((double)rand() / RAND_MAX <= frac)
            {
                edges.push_back(make_pair(i, j));
                edges.push_back(make_pair(j, i));
            }
        }
    }
    sort(edges.begin(), edges.end());
    for (int i = 0; i < edges.size(); i++)
    {
        cout << edges[i].first << " " << edges[i].second << endl;
    }
    return 0;
}

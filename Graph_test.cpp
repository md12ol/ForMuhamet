#include "Graph/Graph.h"

int main() {
    // Create a graph with 5 nodes
    Graph graph(5);

    // Create a vector for the weights of the edges in the nodes
    vector weights = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1};

    // Fill the graph with these weights
    graph.fill(weights, true);

    // Print the graph to console
    graph.print(cout);
    return 0;
}
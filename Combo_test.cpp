#include "Graph/Graph.h"
#include "SDA/SDA.h"

int main() {
    /*
     * TODO: See if you can create an SDA and use its output to create a graph.
     *
     * Consider a graph that is stored in an upper-triangular adjacency matrix
     * That means if a graph has 4 nodes: A, B, C, and D with the following connections
     * A <-> B
     * B <-> C
     * B <-> D
     * Then the matrix would look like this:
     * 0 1 0 0
     * 0 0 1 1
     * 0 0 0 0
     * 0 0 0 0
     *
     * In general, a graph with N nodes needs N*(N-1)/2 weights to store all the possible edges.
     * NUMER OF WEIGHTS == NUMBER OF EDGES? or is it because edge with weight = 0 ist not an existing edge but still counts as weight
    */
    SDA sda_test(3, 3, 2, 6);

    vector<int> output(6);

    sda_test.fillOutput(output, true, std::cout);

    Graph graph_test(4);

    graph_test.fill(output, false);

    graph_test.print(cout);
    return 0;

}
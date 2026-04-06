import os
import ast
import pandas as pd
from graphviz import Graph


RUN_NAME = "Output (NSGA2 - spread vs cost) - 96PS, 20000Mevs, 2.7%MuR, 99%CrR, 2ETA, 1Islands, 30SEpis, 12ST"


BASE_INPUT_DIR = "cmake-build-release-wsl/Output"
BASE_OUTPUT_DIR = "GraphOutput"

'''
For a particular folder we need to:
1. Determine which of the particular lines has the best hypervolume
2. Get a list of all the SDAs
3. Sort this list of SDAs based on a decreasing of one of the objective functions [[(#edges, length), SDA], [(#edges, length), SDA]], etc.]
4. For each SDA:
    4.1. Convert the floating point representation into an SDA
    4.2. Get the output from the SDA for a sufficient number of characters (n*(n-1)/2)
    4.3. Input the characters into an adjacency matrix
    4.4. Generate a network of that adjacency matrix (modify edge_list(...) below)

Note: if we have the folder "Output - 123" then just name the networks "Output - 123 - 1.png", "Output - 123 - 2.png", etc.
'''

'''
Currently this method expects a string formatted as follows:
2 3 4
1 3 4
1 2
1 2

Where each line represents the edges associated with a single node
For example: the line '2 3 4' is the first line so it represents the edges for node one
And that means that node 1 has an edge going to nodes 2, 3, and 4
'''


def edge_list(edge_lists, verts: int):
    # Make an adjacency matrix
    adjM = [[False for _ in range(verts)] for _ in range(verts)]
    for from_node, line in enumerate(edge_lists):
        line = line.rstrip()
        line = line.split(" ")
        for to_node in line:
            if to_node != '':
                if from_node < int(to_node):
                    adjM[from_node][int(to_node) - 1] = True
                    adjM[int(to_node) - 1][from_node] = True
                    pass
                pass
            pass
        pass

    # Get the list of edges: [[1, 3], [2, 7], etc.]
    edge_lists = []
    for row in range(verts):
        for col in range(row + 1, verts):
            if adjM[row][col]:
                edge_lists.append([row, col])
                pass
            pass
        pass
    return edge_lists, adjM


# added new function
def get_edges_from_flat_list(weights_list, verts: int):
    adjM = [[False for _ in range(verts)] for _ in range(verts)]
    edge_lists = []
    k = 0

    for iter in range (1, verts):
        for row in range(verts - iter):
            col = row + iter
            if weights_list[k] == 1:
                adjM[row][col] = True
                adjM[col][row] = True
                edge_lists.append([row, col])
            k += 1

    return edge_lists, adjM



def make_graph(el, out_file: str, verts: int, out_dir: str):
    g = Graph(engine='sfdp')
    e_cout = 0

    g.graph_attr.update(dpi='1000', size="6,6", outputorder='edgesfirst', overlap='false', splines='true')
    g.node_attr.update(color='black', shape='point', width='0.04', height='0.04')
    # g.node_attr.update(color='black', shape='circle', fixedsize='true', width='0.25', fontsize='8')
    g.edge_attr.update(color='black', penwidth='0.5')


    # for i in range(verts):
    #     g.node(str(i))
    #     pass

    for n in range(verts):
        if n == 0:
            g.node(str(n), label=str(n), color='red')
            pass
        else:
            g.node(str(n), label=str(n))
            pass
        pass

    for idx, d in enumerate(el):
        if d[0] < d[1]:
            g.edge(str(d[0]), str(d[1]), color='black')
            e_cout += 1
            pass
        pass


    g.render(filename=out_file, directory=out_dir, cleanup=True, format='png')

    # g.save(filename=out_file, directory=outp)
    # g.clear()
    print("Made network: " + os.path.join(out_dir, out_file) + ".png with " + str(e_cout) + " Edges")
    pass


def main():

    csv_file_path = os.path.join(BASE_INPUT_DIR, RUN_NAME, "best_run_networks.csv")

    target_output_folder = os.path.join(BASE_OUTPUT_DIR, RUN_NAME)

    os.makedirs(target_output_folder, exist_ok=True)

    NUM_VERTS = 256

    try:
        # Read and sort data based on edges (decreasing)
        df = pd.read_csv(csv_file_path, sep=";")
        df_sorted = df.sort_values(by="Edges", ascending=False).reset_index(drop=True)

        # save the sorted dataframe as a new csv overview
        sorted_csv_path = os.path.join(target_output_folder, "best_run_networks_sorted.csv")

        df_sorted.to_csv(sorted_csv_path, sep=";", index=False)
        print(f"Saved sorted CSV overview to: {sorted_csv_path}")


        # Loop through all networks
        for index, row in df_sorted.iterrows():
            weights = ast.literal_eval(row['EdgeWeights'])

            # Use our new function to get the edges
            el, adjM = get_edges_from_flat_list(weights, NUM_VERTS)

            # Create file name
            out_name = f"{RUN_NAME} - {index + 1}"

            # Generate the network
            make_graph(el, out_name, NUM_VERTS, target_output_folder)

    except FileNotFoundError:
        print(f"File {csv_file_path} not found. Please check the path.")


if __name__ == "__main__":
    main()
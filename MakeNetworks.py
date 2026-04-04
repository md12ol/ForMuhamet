import os
from graphviz import Graph


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

    edge_lists = []
    for row in range(verts):
        for col in range(row + 1, verts):
            if adjM[row][col]:
                edge_lists.append([row, col])
                pass
            pass
        pass
    return edge_lists, adjM


def make_graph(el, out_file: str, verts: int):
    g = Graph(engine='sfdp')
    e_cout = 0

    g.graph_attr.update(dpi='1000', size="6,6", outputorder='edgesfirst', overlap='false', splines='true')
    g.node_attr.update(color='black', shape='point', width='0.04', height='0.04')
    # g.node_attr.update(color='black', shape='circle', fixedsize='true', width='0.25', fontsize='8')
    g.edge_attr.update(color='black', penwidth='1.5')

    for i in range(verts):
        g.node(str(i))
        pass

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
    g.render(filename=out_file, directory='./', cleanup=True, format='png')
    # g.save(filename=out_file, directory=outp)
    # g.clear()
    print("Made network: " + out_file + " with " + str(e_cout) + " Edges")
    pass

def main():
    with open("./exampleGraph.txt") as f:
        print(f)
        el, adjM = edge_list(f, 4)
        make_graph(el, "exampleGraph", 4)

main()
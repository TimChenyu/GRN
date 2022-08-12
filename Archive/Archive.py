
def add_nodes(G, parent, num_to_add):
    num_Nodes = G.number_of_nodes()
    nodes_to_add = []
    for n in range(1, 3):
        nodes_to_add.append(num_Nodes+n)
    for node in nodes_to_add:
        G.add_edge(parent, node)


def generate_gene_network(maxNodes, p_split, p_3):
    """
    generate a tree with depth n

    each level can have a maximum of 6n + 1 nodes

    first generate a center node
        if max depth not reached:
            generate 0-> 6n+1 nodes:
                distribute the nodes across the nodes
                nodes remaining = 6n+1
                for nodes in this level:
                    num child = random(0, nodes remaining)
                    nodes remaining = nodes remanining - numchild

    """
    G = nx.path_graph(1)
    num_Nodes = G.number_of_nodes()
    parentNodes = []
    newparents = []
    if random.random() > p_3:
        add_nodes(G, 0, 3)
        num_Nodes = G.number_of_nodes()
        for n in range(1, 3):
            parentNodes.append(num_Nodes-n)
    else:
        add_nodes(G, 0, 2)
        num_Nodes = G.number_of_nodes()
        for n in range(1, 2):
            parentNodes.append(num_Nodes-n)
    count = 0
    newparents = parentNodes
    preNum = 111
    while G.number_of_nodes() < maxNodes:
        count += 1
        #print(count)
        print(G.number_of_nodes())
        if preNum != G.number_of_nodes():
            parentNodes = newparents
            newparents = []
        preNum = G.number_of_nodes()
        for node in parentNodes:
            if random.random() > p_split:
                if random.random() > p_3:
                    add_nodes(G, node, 5)
                    num_Nodes = G.number_of_nodes()
                    for n in range(1, 3):
                        newparents.append(num_Nodes-n)
                else:
                    add_nodes(G, node, 3)
                    num_Nodes = G.number_of_nodes()
                    for n in range(1, 3):
                        newparents.append(num_Nodes-n)

        if preNum != G.number_of_nodes():
            parentNodes = []
    return G

from itertools import combinations, groupby
import networkx as nx
import random
import numpy as np
import matplotlib.pyplot as plt


class BoltzmannGRN():
"""
定义一个玻尔兹曼机基因调控网络的类
Define a Boltzmann Gene Regulatory Network Class
"""
    def __init__(self, size):
        self.time = 0
        self.size = size
        self.relationProbability = 0.1
        #self.graph = gnp_random_connected_graph(self.size, self.relationProbability)
        self.graph = gnp_random_graph(self.size, 4)
        self.matrix = make_wmatrix(self.graph, self.size)

    def update(self):
        self.time += 1

def make_wmatrix(G, size):
    """
    Make a weight matrix for the graph
    The simpliest way of setting weights, no additional rules
        are added in this method, all edges have a weight of 1
    """
    nums = nx.get_node_attributes(G, "number")
    print(nums)
    wmatrix = np.zeros((size, size))
    for node in G.nodes():
        ycount = 0
        for neighbour in G.neighbors(node):
            wmatrix[node][neighbour] = 1
    return wmatrix

def gnp_random_connected_graph(n, p):
    """
    Generates a random undirected graph, similarly to an Erdős-Rényi
    graph, but enforcing that the resulting graph is conneted
    """
    edges = combinations(range(n), 2)
    G = nx.Graph()
    G.add_nodes_from(range(n))
    if p <= 0:
        return G
    if p >= 1:
        return nx.complete_graph(n, create_using=G)
    for _, node_edges in groupby(edges, key=lambda x: x[0]):
        node_edges = list(node_edges)
        random_edge = random.choice(node_edges)
        G.add_edge(*random_edge)
        for e in node_edges:
            if random.random() < p:
                G.add_edge(*e)
    return G

def gnp_random_graph(n, max_nei):
    G = nx.path_graph(n)
    #pos = nx.circular_layout(G)
    """
    Generates a random undirected graph, similarly to an Erdős-Rényi
    graph, but enforcing that the resulting graph is conneted
    """
    #G = nx.Graph()
    #G.add_nodes_from(range(n))
    num_nei_list = []
    for node in G.nodes:
        num_nei = random.randint(0, max_nei)
        num_nei_list.append(num_nei)
    for node in G.nodes:
        list_nei = []
        list_nei = random.sample(G.nodes, num_nei_list[node])
        if check_max(G, node, max_nei) == False:
            for nei in list_nei:
                if check_max(G, nei, max_nei) == False:
                    G.add_edge(node, nei)
                else:
                        #print('hit')
                    for x in G.nodes:
                        if check_max(G, nei, max_nei) == False:
                            G.add_edge(node, x)
        else:
            print('exception caught')
    ooo = []
    for node in G.nodes:
        check = [n for n in G.neighbors(node)]
        print(check)
        checks = len(check)
        ooo.append(checks)
    print(ooo)
    print(num_nei_list)

    return G

def check_max(G, node, max):
    if len([n for n in G.neighbors(node)]) >= max:
        return True
    else:
        return False

numGenes = 100
def getHit(gene_i, time):
    """
    Get the H_{i,t} value for gene i at timestep t
    ##############################################
    1. Gene_i : The gene that we want to target to
        get the H value for
    2. time: The numerical timestep
    ##############################################
    """
    for neighbour in gene:
        phit = weight(gene_i, neighbour) * state(neighbour, time)
        hit += phit
    return hit

def probability(gene, time, hit):
    """
    Get the probability of a specific gene being on
    ##############################################
    1. Gene: The target gene
    2. hit: The H_{i,t} value we obtained in the
        getHit function
    3. time: The numerical timestep
    ##############################################
    """
    p = (1/(1+2.718281828459045**(getHit(gene, time)*(-1))))
    return p

def update():
    print('fml')
size = 100
test = BoltzmannGRN(size)
print(test.matrix)
G = test.graph

#nx.add_path(G, range(10))
#nx.add_star(G, range(9, 13))
#pos = nx.spring_layout(G, seed=225)  # Seed for reproducible layout
nx.draw(G)
plt.show()

from itertools import combinations, groupby
import networkx as nx
import random
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from networkx import is_connected, connected_components


def makeDistanceMatrix(G, size):
    Distance_Matrix = np.zeros((size, size))
    for node in G.nodes():
        p = nx.shortest_path(G, source=node)
        for node1 in G.nodes():
            if node != node1:
                try:
                    distance = len(p[node1])
                except KeyError:
                    G.add_edge(node, node1)
                    p = nx.shortest_path(G, source=node)
                    distance = len(p[node1])
                    print('hit')

                Distance_Matrix[node-1][node1-1] = distance
    return Distance_Matrix
def generate_Gene_tree(n):
    G=nx.Graph()
    G.add_node(0)
    for i in range(1, n):
        num = G.number_of_nodes()
        print(num)
        G.add_node(i)
        toConnect = random.randint(0, num-1)
        G.add_edge(toConnect, i)
    return G
def add_nodes(G, parent, num_to_add):
    num_Nodes = G.number_of_nodes()
    nodes_to_add = []
    for n in range(1, 3):
        nodes_to_add.append(num_Nodes+n)
    for node in nodes_to_add:
        G.add_edge(parent, node)

G = generate_Gene_tree(100)
dm = makeDistanceMatrix(G, G.number_of_nodes())





def get_3_4_neighbours(G, matrix):
    nodeList = []
    rownum = 0
    for eachrow in matrix:
        temp = []
        row = eachrow.tolist()
        #print(row)
        for n in range(0, len(row)):
            if (row[n] ==3 ):
                temp.append(n)
        nodeList.append(temp)
        #print(temp)
        rownum += 1
        connected = False
        if rownum > 1:
            while connected == False:
                #print('hit')
                #print(temp)
                for i in temp:
                    if random.random() > 0.8:
                        if i < rownum:
                            pass
                        #print(type(i))
                        G.add_edge(rownum, i)
                        pass
                        connected = True
    return G
    print(nodeList)



#G = generate_gene_network(100, 0.2, 0.2)
#components = list(connected_components(G))
#biggest_component_size = max(len(c) for c in components)
#problem_components = [c for c in components if len(c) != biggest_component_size]
#for component in problem_components:
#    for node in component:
#        G.remove_node(node)
#print(is_connected(G))
#G.delete_nodes_from(solitary)
#print(nodes)
#dm = makeDistanceMatrix(G, G.number_of_nodes())
#G = get_3_4_neighbours(G, dm)
#print(dm)
#pos = nx.spring_layout(G, seed=225)  # Seed for reproducible layout
#nx.draw(G)
#plt.show()
#exit()












tree = nx.random_tree(n=100, seed=0)










"""
???????????????????????????????????????????????????
Define a Boltzmann Gene Regulatory Network Class
"""
class BoltzmannGRN():

    def __init__(self, size):
        self.time = 0
        self.size = size
        #self.relationProbability = 0.1
        #self.graph = gnp_random_connected_graph(self.size, self.relationProbability)
        #self.graph = gnp_random_graph(self.size, 5)
        self.graph = generate_Gene_tree(size)
        self.matrix = make_wmatrix(self.graph, self.size)
        for node in self.graph.nodes:
            chance = random.random()
            if chance > 1:
                self.graph.nodes[node]["state"] = 1
            else:
                self.graph.nodes[node]["state"] = 0

    def update(self):
        self.time += 1

def make_wmatrix(G, size):
    """
    Make a weight matrix for the graph
    The simpliest way of setting weights, no additional rules
        are added in this method, all edges have a weight of 1
    """
    nums = nx.get_node_attributes(G, "number")
    #print(nums)
    wmatrix = np.zeros((size, size))
    for node in G.nodes():
        ycount = 0
        for neighbour in G.neighbors(node):
            weight = 0.3#random.random()
            wmatrix[node][neighbour] = weight
    return wmatrix

def gnp_random_connected_graph(n, p):
    """
    Generates a random undirected graph, similarly to an Erd??s-R??nyi
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
    Generates a random undirected graph, similarly to an Erd??s-R??nyi
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
        #print(check)
        checks = len(check)
        ooo.append(checks)
    #print(ooo)
    #print(num_nei_list)

    return G

def check_max(G, node, max):
    if len([n for n in G.neighbors(node)]) >= max:
        return True
    else:
        return False

numGenes = 100

def state(G, node):
    state = G.nodes[node]["state"]
    return state

def getHit(G, gene_i, model):
    """
    Get the H_{i,t} value for gene i at timestep t
    ##############################################
    1. Gene_i : The gene that we want to target to
        get the H value for
    2. time: The numerical timestep
    ##############################################
    """
    hit = 0
    for neighbour in G.neighbors(gene_i):
        phit = model.matrix[gene_i][neighbour] * state(G, neighbour)
        hit += phit
    return hit

def probability(G, gene, hit, model):
    """
    Get the probability of a specific gene being on
    ##############################################
    1. Gene: The target gene
    2. hit: The H_{i,t} value we obtained in the
        getHit function
    3. time: The numerical timestep
    ##############################################
    """
    p = (1/(1+2.718281828459045**(getHit(G, gene, model)*(-1))))
    return p

def update(G, model):
    for node in G:
        hit = getHit(G, node, model)
        p = probability(G, node, hit, model)
        chance = random.random()
        if chance <= p:
            G.nodes[node]["state"] = 1
        else:
            G.nodes[node]["state"] = 0

def group_duplicate_index(df):
    a = df.values
    sidx = np.lexsort(a.T)
    b = a[sidx]

    m = np.concatenate(([False], (b[1:] == b[:-1]).all(1), [False] ))
    idx = np.flatnonzero(m[1:] != m[:-1])
    I = df.index[sidx].tolist()
    return [I[i:j] for i,j in zip(idx[::2],idx[1::2]+1)]

size = 100
title = []
listi = []
for i in range(0, size):
    title.append(str(i))
    listi.append('0')
#listi = ['0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0']
#
#print(test.matrix)
#G = test.graph
t = 0


#print(statedf)
test = BoltzmannGRN(size)
list_size = []
list_states = []
state_time = []
time = []
timestep = 0
state_num = 0
for t in range(0, 1000):
    """if t > 0:
        #print(state_time)
        plt.rcParams["figure.figsize"] = [17.50, 3.50]
        plt.rcParams["figure.autolayout"] = True
        plt.plot(time, list_size, 'r*')
        #plt.show()"""
    time.append(timestep)
    timestep += 1
    states = []
    update(G, test)
    for node in G.nodes:
        s = G.nodes[node]["state"]
        #print(s)
        states.append(s)
        sums = sum(states)
        sum_index = sums/size
    """if states not in list_states:
        state_time.append(state_num)
        state_num += 1
        max_state_num = state_num
    if states in list_states:
        state_num = list_states.index(states)
        print('Timestep ' + str(t) + ' : hit state: ' + str(state_num))
        state_time.append(state_num)
        state_num = max_state_num"""
    list_size.append(sum_index)
    list_states.append(states)
    #print(len(states))
    #print(states)
df = pd.DataFrame(list_states, columns = listi)
print(df)

graphy = []
for index, row in df.iterrows():
    graphy.append(row.sum())

print(state_time)
plt.rcParams["figure.figsize"] = [17.50, 3.50]
#plt.rcParams["figure.autolayout"] = True
plt.plot(time, list_size, 'r*')
plt.show()


duplicates = group_duplicate_index(df)
print("==========================")
print("Duplicated Rows")
for element in duplicates:
    print(element)

#print(duplicates)
#nx.add_path(G, range(10))
#nx.add_star(G, range(9, 13))
#pos = nx.spring_layout(G, seed=225)  # Seed for reproducible layout
#nx.draw(G)
#plt.show()

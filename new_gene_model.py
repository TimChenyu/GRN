import random
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
from itertools import combinations, groupby
from networkx import is_connected, connected_components


"""
定义一个玻尔兹曼机基因调控网络的类
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

def generate_Gene_tree(n):
    G=nx.DiGraph()
    G.add_node(0)
    for i in range(1, n):
        num = G.number_of_nodes()
        print(num)
        G.add_node(i)
        toConnect = random.randint(0, num-1)
        G.add_nodes_from([toConnect, i])
        if random.random() >= 0.5:
            G.add_nodes_from([toConnect, i])
        #G.add_edge(toConnect, i)
    return G
def add_nodes(G, parent, num_to_add):
    num_Nodes = G.number_of_nodes()
    nodes_to_add = []
    for n in range(1, 3):
        nodes_to_add.append(num_Nodes+n)
    for node in nodes_to_add:
        G.add_edge(parent, node)

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
            weight = 0.11#random.random()
            wmatrix[node][neighbour] = weight
    return wmatrix

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

def state(G, node):
    state = G.nodes[node]["state"]
    return state
size = 1000
title = []
listi = []
for i in range(0, size):
    title.append(str(i))
    listi.append('0')
#listi = ['0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0','0']
#
#print(test.matrix)

t = 0


#print(statedf)
test = BoltzmannGRN(size)
G = test.graph
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
plt.xlim([0, t+20])
plt.ylim([0, 1])
plt.rcParams["figure.figsize"] = [17.50, 3.50]
#plt.rcParams["figure.autolayout"] = True
plt.plot(time, list_size, 'r*')
plt.show()


duplicates = group_duplicate_index(df)
print("==========================")
print("Duplicated Rows")
for element in duplicates:
    print(element)

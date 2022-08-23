# GRNPY - Gene Regulatory Network .py

This python package aims to be a Gene Regulatory Network Graph generator that is fully customizable parameters such as sparsity, directed graphs, weight of edges, type of system (binary or spin) and etc.

Aug 23: Work in progress
- take a *dict* or *list* of *lists* as input and converting them to graphs
- take a graph as an input to the system


# Usage

In the terminal, type:

    pip install grnpy

In the code, import grnpy and the class *BoltzmannGRN*

    from grnpy import BoltzmannGRN

# Functions
## Generate_Gene_Network(self, n, sparsity, on_p, directed):
 **Generates a graph based on the parameters :**
 - *n* - int - number of nodes
 - *sparsity* - float - define the sparsity of the graph (1 being sparse, and 0 being dense) [0,1]
 - *on_p* - float - the probability of a node being on at the beginning [0,1]
 - *directed* - bool - if the graph is directed

returns graph

## make_wmatrix(self, G, size, weight):
**Generates an edge weight matrix**

 - *size* - int - size of the graph
 - *weight* - defined set weight of the graph (going to be removed or chagned in later versions)

returns np array
## getHit(self, G, gene_i, field)
**Generates the h value at timestep t for node i**

returns H_{i,t}

## probability(self, G, gene, hit, field)
returns probability of a node being at state "1"

returns p
## change_state(self, G, field)
Updates the state of all nodes in G

returns NONE
## state(self, G, node)
returns the state of node n in graph G
returns state(node)
## update(self, size, G, timesteps, field)
For each timestep, iterate through all nodes in G and update their state

 - *timestep* - int - number of timesteps to be updated
 - *size* - int - size of the graph
 - *G* - nx.Graph - Graph of nodes
 - *field* - float - constant that is added as a probablity
 - *time* - list - list of timesteps
 - *states_mean* - list - mean state value across timesteps
 - *list_states* - list - list of list containing the specific states of each node for each timestep

returns time, states_mean, list_states

## show_plot(self, time, states_mean)
shows plot of state transitions


# Define and Use
Define the model

    BoltzmannModel = BoltzmannGRN(size, gene_start_probability, weight, sparsity, directed, 'binary')

Define the graph

    BoltzmannGraph = BoltzmannModel.graph

Drawing the graph

    nx.draw(BoltzmannGraph)
    plt.show()

Run the state transition model

    time, states_mean, states = BoltzmannModel.update(size, BoltzmannGraph, 100, field)

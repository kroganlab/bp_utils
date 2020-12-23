#! /usr/local/bin/python3

# TODO:  build as an actual script with command line arguments.
#   as of now, this is just a copy of code run from a command block in RStudio

import networkx as nx
import numpy as np
import pandas as pd
import time

# For math guidance/explanation
# see https://doi.org/10.1038/nrg.2017.38
#  Network propagation: a universal amplifier of genetic associations
# AU  - Cowen, Lenore
# AU  - Ideker, Trey
# AU  - Raphael, Benjamin J.
# AU  - Sharan, Roded
# PY  - 2017
# DA  - 2017/09/01
# TI  - Network propagation: a universal amplifier of genetic associations
# JO  - Nature Reviews Genetics
# 
# This uses an S matrix (terminology from above paper) which is the transpose of 'inv_denom' used elsewhere in Krogan lab (see Mehdi's code)
# to propagate with an S matrix do: S %*% heat_0
# to propagate with inverse_denom do: heat_0 %*% inv_denom * pr
# Note:
# ! Order around %*% matters and S includes the pr factor !
# col sums of S = 1
# row sums of inv_denom = 1/pr 


mouseEdges = pd.read_csv ("~/Downloads//10090.protein.links.v11.0.txt.gz", sep = " ")
filteredEdges = mouseEdges[mouseEdges['combined_score'] > 600][['protein1','protein2']]
edge_list = filteredEdges.values.tolist()

# Create Graph
G = nx.Graph(edge_list)
G.remove_edges_from(nx.selfloop_edges(G)) # removing selfloops
largest_cc = max(nx.connected_components(G), key=len) # only taking the largest connected component
Gc = G.subgraph(largest_cc)

# Calculate matrix inputs
nodes = Gc.nodes()
A = nx.adjacency_matrix(Gc, nodelist=nodes)
D = np.array([x[1] for x in nx.degree(Gc)]) # degree of each node

# integrity check: node order matches
degreeNodes = np.array([x[0] for x in nx.degree(Gc)])
assert (all(degreeNodes == nodes))

Dinverse = np.diag(1/D)
W = np.matmul(A.todense(), Dinverse)
I = np.identity(W.shape[0])

for pr in (0.2, 0.3):
  print ("starting matrix inverse calculation at {}".format(time.asctime(time.localtime())))
  S_matrix = pr * np.linalg.inv(I-(1-pr)*W)
  print ("finished matrix inverse calculation at {}".format(time.asctime(time.localtime())))
  np.save("S_matrix.mouse.string.gt600.pr%0.2f.npy"%(pr), S_matrix)
  # node names
  nodes_network = pd.DataFrame(nodes) # Getting node names(gene names)
  nodes_network.rename(columns={0: "gene_name"}, inplace=True) # rename column header
  nodes_network.to_csv("S_matrix.mouse.string.gt600.pr%0.2f.nodeNames.csv" %(pr), index=False)
  

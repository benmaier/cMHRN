import numpy as np
import cMHRN
import networkx as nx
import time

N = 1000
k = 3
p = k/(N-1.)

G = nx.fast_gnp_random_graph(N,p)


tic = time.time()
comp_sizes_cMHRN = sorted([len(g) for g in cMHRN.get_components(N,G.edges())])
toc = time.time()

ticnx = time.time()
comp_sizes_nx = sorted([len(g) for g in list(nx.connected_component_subgraphs(G))])
tocnx = time.time()

print toc-tic,comp_sizes_cMHRN
print tocnx-ticnx,comp_sizes_nx


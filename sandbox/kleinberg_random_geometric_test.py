import mhrn
import cMHRN
import networkx as nx
import numpy as np
import time
import scipy.sparse as sprs
import progressbar

N = 10000
k = 7
mu = -1

N_meas = 10

k_nx, k_c = 0.0, 0.0
t_nx, t_c = 0.0, 0.0

seed = 371

bar = progressbar.ProgressBar()
for meas in bar(range(N_meas)):
    tic = time.time()
    N, G2 = cMHRN.kleinberg_network(N,k,mu)
    toc = time.time()
    t_c += (toc-tic) / N_meas
    k_c += len(G2)*2./N / float(N_meas)

    tic = time.time()
    N, G = cMHRN.random_geometric_kleinberg_network(N,k,mu,epsilon=1e-3)
    toc = time.time()
    t_nx += (toc-tic) / N_meas
    k_nx += len(G)*2./N / float(N_meas)

print("lattice  : k = %4.2f, t = %fs" % (k_c,t_c))
print("continuum: k = %4.2f, t = %fs" % (k_nx,t_nx))
print() 



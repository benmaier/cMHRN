import mhrn
import cMHRN
import networkx as nx
import numpy as np
import time


B = 4
L = 3
k = 3
xi = 0.3

N_meas = 1000

k_nx, k_c = 0.0, 0.0
t_nx, t_c = 0.0, 0.0

seed = 1

for meas in xrange(N_meas):
    seed += 1
    tic = time.time()
    G1 = mhrn.fast_mhr_graph(B,L,k,xi)
    toc = time.time()
    t_nx += (toc-tic) / N_meas
    k_nx += G1.number_of_edges()*2./B**L / N_meas

    tic = time.time()
    G2 = cMHRN.fast_mhrn(B,L,k,xi,seed=seed)
    toc = time.time()
    t_c += (toc-tic) / N_meas
    k_c += len(G2)*2./B**L / N_meas

print "python: k = %4.2f, t = %fs" % (k_nx,t_nx)
print "c++   : k = %4.2f, t = %fs" % (k_c,t_c)

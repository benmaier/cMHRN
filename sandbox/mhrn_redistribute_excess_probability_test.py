import mhrn
import cMHRN
import networkx as nx
import numpy as np
import time
import scipy.sparse as sprs
import progressbar


B = 10
L = 3
k = 13
xi = 0.3

N_meas = 20

k_nx, k_c = 0.0, 0.0
t_nx, t_c = 0.0, 0.0

seed = 1

print("measuring only creation")

bar = progressbar.ProgressBar()
for meas in bar(range(N_meas)):
    tic = time.time()
    N, G2 = cMHRN.fast_mhrn(B,L,k,xi,seed=seed+1,allow_probability_redistribution=True,use_giant_component=False)
    toc = time.time()
    t_c += (toc-tic) / N_meas
    k_c += len(G2)*2./B**L / N_meas

#print "python: k = %4.2f, t = %fs" % (k_nx,t_nx)
print("c++   : k = %4.2f, t = %fs" % (k_c,t_c))


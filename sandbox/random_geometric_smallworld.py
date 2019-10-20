import mhrn
import cMHRN
import networkx as nx
import numpy as np
import time
import scipy.sparse as sprs
import progressbar


N = 100
k = 32
beta= 0.3

N_meas = 200

k_mod, k_rgg = 0.0, 0.0
t_nx, t_c = 0.0, 0.0

seed = 1

bar = progressbar.ProgressBar()
for meas in bar(range(N_meas)):
    tic = time.time()
    _N, Grgg = cMHRN.random_geometric_small_world_network(N,k,beta,use_giant_component=False)
    _N, Gmod = cMHRN.modified_small_world_network(N,k,beta,use_giant_component=False)
    toc = time.time()
    t_c += (toc-tic) / N_meas
    k_rgg += len(Grgg)*2./N / N_meas
    k_mod += len(Gmod)*2./N / N_meas

#print "python: k = %4.2f, t = %fs" % (k_nx,t_nx)
print("RGG   : k = %4.2f, t = %fs" % (k_rgg,t_c))
print("MOD   : k = %4.2f, t = %fs" % (k_mod,t_c))


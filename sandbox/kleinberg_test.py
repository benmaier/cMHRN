import mhrn
import cMHRN
import networkx as nx
import numpy as np
import time
import scipy.sparse as sprs
import progressbar

drawing = False
B = 7
L = 3
k = 7
xi = 0.3

N_meas = 10

k_nx, k_c = 0.0, 0.0
t_nx, t_c = 0.0, 0.0

seed = 371

bar = progressbar.ProgressBar()
for meas in bar(xrange(N_meas)):
    tic = time.time()
    N, G2 = cMHRN.kleinberg_network(B**L,k,np.log(xi)/np.log(B))
    toc = time.time()
    t_c += (toc-tic) / N_meas
    k_c += len(G2)*2./B**L / float(N_meas)

    tic = time.time()
    G = mhrn.continuous_hierarchical_graph(B**L,k,np.log(xi)/np.log(B),redistribute_probability=True)
    toc = time.time()
    t_nx += (toc-tic) / N_meas
    k_nx += G.number_of_edges()*2./B**L / float(N_meas)

print "python: k = %4.2f, t = %fs" % (k_nx,t_nx)
print "c++   : k = %4.2f, t = %fs" % (k_c,t_c)
print 


if drawing:
    N, G2 = cMHRN.kleinberg_network(B**L,k,np.log(xi)/np.log(B),seed=seed+meas)
    G = nx.Graph()
    G.add_edges_from(G2)

    import pylab as pl

    nx.draw(G)
    pl.show()


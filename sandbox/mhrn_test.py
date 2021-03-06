import mhrn
import cMHRN
import networkx as nx
import numpy as np
import time
import scipy.sparse as sprs
import progressbar


B = 10
L = 5
k = 6
xi = 0.3

N_meas = 10

k_nx, k_c = 0.0, 0.0
t_nx, t_c = 0.0, 0.0

seed = 1

print("measuring only creation")

bar = progressbar.ProgressBar()
for meas in bar(range(N_meas)):
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

print("python: k = %4.2f, t = %fs" % (k_nx,t_nx))
print("c++   : k = %4.2f, t = %fs" % (k_c,t_c))
print() 

print("measuring with extraction of giant component")

k_nx, k_c = 0.0, 0.0
t_nx, t_c = 0.0, 0.0


bar = progressbar.ProgressBar()
for meas in bar(range(N_meas)):
    seed += 1
    tic = time.time()
    G1 = mhrn.fast_mhr_graph(B,L,k,xi)
    giant = max(nx.connected_component_subgraphs(G1), key=len)
    N_new = giant.number_of_nodes()
    A = nx.to_scipy_sparse_matrix(giant)
    current_k = float(np.array(A.sum(axis=0)).flatten().sum(axis=0))
    toc = time.time()

    t_nx += (toc-tic) / N_meas
    k_nx += current_k / N_new / N_meas

    tic = time.time()
    N,rows,cols = cMHRN.fast_mhrn_coord_lists(B,L,k,xi,seed=seed,use_giant_component=True)
    rows = np.array(rows)
    cols = np.array(cols)
    data = np.ones_like(rows)
    A = sprs.csr_matrix((data,(rows,cols)),shape=(N,N))
    current_k = float(np.array(A.sum(axis=0)).flatten().sum(axis=0))
    toc = time.time()
    t_c += (toc-tic) / N_meas
    k_c += current_k / N_new / N_meas

print("python: k = %4.2f, t = %fs" % (k_nx,t_nx))
print("c++   : k = %4.2f, t = %fs" % (k_c,t_c))

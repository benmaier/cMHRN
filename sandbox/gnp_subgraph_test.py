import cMHRN
import networkx as nx
import numpy as np
import time


N = 10
p = 0.02
seed = int(time.time())

N_meas = 100


k_c, k_nx = 0., 0.


edges = cMHRN.fast_gnp(10,p,start_node=10,seed=seed)
print edges


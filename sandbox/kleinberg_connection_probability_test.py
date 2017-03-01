import numpy as np
import pylab as pl
import mhrn
import cMHRN

N = 9
k = 5
mu = -0.5

p_py = mhrn.continuous_hierarchical_graph_distribution(N,k,mu,True)
p_c = cMHRN.get_kleinberg_connection_probability(N,k,mu)

x = np.arange(1,N)
pl.plot(x,p_py)
pl.plot(x,p_c)

print "sum p_py = ", sum(p_py)
print "sum p_c = ", sum(p_c)

pl.show()

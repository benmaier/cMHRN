import cMHRN 
import numpy as np
import scipy.sparse as sprs
from scipy.special import binom 

import matplotlib.pyplot as pl


def P_theory(N,k,beta,kmax):
    k = int(k)
    p0 = (1-beta+beta*(N-1)/k)**(-1)


    P = np.zeros(kmax+1)

    for c in range(kmax+1):
        s = 0.

        n = np.min([c,k])
        for ks in range(n+1):
            s += binom(k,ks) * p0**ks * (1-p0)**(k-ks) * binom(N-1-k,c-ks) * \
                 (beta*p0)**(c-ks) * (1-beta*p0)**(N-1-k-(c-ks))
        P[c] = s

    return P



N_meas = 1000


betas = np.logspace(-3,0,10)
k = 8.
N = 100


for beta in betas:
    hist = np.zeros((N,))

    for meas in range(N_meas):

        _, row, col = cMHRN.modified_small_world_network_coord_lists(
                N,
                k,
                beta,
                use_giant_component = False,
                )

        A = sprs.csr_matrix((np.ones_like(row),(row,col)), shape=(N,N))

        degree = np.asarray(A.sum(axis=1)).reshape((N,))

        for k_ in degree:
            hist[k_] += 1.0

    hist /= hist.sum()

    kmax = (np.where(hist>0)[0]).max()
    degrees = np.arange(kmax+1,dtype=int)

    this_plot, = pl.step(degrees,hist[:kmax+1],where='mid')
    pl.step(degrees,P_theory(N,k,beta,kmax),c=this_plot.get_color(),lw=3, alpha=0.4,linestyle = '--',where='mid')
pl.show()


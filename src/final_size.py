#!/usr/bin/env python3
import numpy as np

def calc_n_class_S_final(B_mat, S0s, I0s,
                         recoveries=None,
                         maxiter = 10000):
    if recoveries is None:
        inv_recov = np.identity(S0s.size)
    else:
        inv_recov = np.diag(1 / recoveries)

    def S_iter(X):
        BE = B_mat @ inv_recov
        internal = X - S0s - I0s
        term = BE @ internal
        return S0s * np.exp(term)
    
    X = S0s
    k = 0
    notstop = True
    while(notstop and k < maxiter):
        Xprime = S_iter(X)
        if np.all(np.abs(Xprime - X) < 1e-10):
            notstop = False
        X = Xprime
        k += 1
        
    if notstop:
        raise RuntimeError("Failed to converge in {}"
                           "iterations".format(maxiter))
    return X


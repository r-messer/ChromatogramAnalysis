import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
from scipy import sparse
from scipy.sparse.linalg import spsolve

# def drift_correction(x):
#     """ Correct baseline drift """
#
#     L = len(y)
#     D = sparse.csc_matrix(np.diff(np.eye(L), 2))
#     w = np.ones(L)
#     for i in xrange(niter):
#     W = sparse.spdiags(w, 0, L, L)
#     Z = W + lam * D.dot(D.transpose())
#     z = spsolve(Z, w*y)
#     w = p * (y > z) + (1-p) * (y < z)
#     return z


def gen_peaks(x):
    return 100*np.exp(np.negative(((x-300)/15)**2)) + 200*np.exp(np.negative(((x-750)/30)**2)) + 100*np.exp(np.negative(((x-800)/15)**2))

def sig_shift(x, y):
    sig = np.linspace(1, 10, len(x))
    return 40/(.5+np.exp(-sig+5)) + y

def lin_shift(x, y):
    lin = np.linspace(10, 40, len(x))
    return 2*lin+y

def baseline_als(y, lam, p, ratio=50):
    """
    :param y: 1D array of data
    :param lam: smoothness parameter
    :param p: skewness parameter
    :param ratio: number of iterations
    :return: calculated baseline of data
    """
    L = len(y)
    D = sparse.diags([1,-2,1],[0,-1,-2], shape=(L,L-2))
    w = np.ones(L)
    for i in range(ratio):
        W = sparse.spdiags(w, 0, L, L)
        Z = W + lam * D.dot(D.transpose())
        z = spsolve(Z, w*y)
        w = p * (y > z) + (1-p) * (y < z)
    return z
    # while True:
    #     W = scipy.sparse.spdiags(w, 0, L, L)
    #     z = np.linalg.solve(H, w*y)
    #     d = y-z
    #     dn = d(d < 0)
    #     m = np.mean(dn)
    #     s = np.std(dn)
    #     wt = 1/(1+np.exp(2*(d-(2*s-m))/s))
    #     if np.linalg.norm(w-wt)/np.linalg.norm(w) < ratio:
    #         break
    #     else:
    #         w = wt


x = np.arange(1,1000)
y = gen_peaks(x)
y_linshift = lin_shift(x, y)
y_sigshift = sig_shift(x, y)
z = baseline_als(y_linshift, 100000000, .01)

plt.plot(x, y_linshift, x, z)
plt.show()


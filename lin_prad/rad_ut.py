

'''
This module contains funcitons
for caluclating purposes in the reconstruction alogorithm
'''
import math

from constants import M_PROTON_G, ESU, C, V_PER_E

from scipy.fftpack import fftn, ifftn
import numpy as np

dmax = 0
delta = 0

def idx2vec(idx):
    xx = -dmax + (idx[0]+0.5)*delta
    yy = -dmax + (idx[1]+0.5)*delta
    return np.array([xx,yy])

def vec2idx(vec):
    i = int((vec[0] + dmax)/delta)
    j = int((vec[1] + dmax)/delta)
    return (i,j)

def position(flux_ref, bin_um):
    num_bins = flux_ref.shape[0]
    delta = bin_um/10000.0
    dmax = (delta * num_bins)/2.0
    x = np.zeros((num_bins, num_bins))
    y = np.zeros((num_bins, num_bins))
    for i in range(num_bins):
        for j in range(num_bins):
            xx = -dmax + ((i + 0.5) * delta)
            yy = -dmax + ((j + 0.5) * delta)
            x[i,j] = xx
            y[i,j] = yy

    return (x,y)

def gradient(fn, idx):
    N0 = fn.shape[0]
    N1 = fn.shape[1]
    dfdx = ( fn[(idx[0]+1)%N0, idx[1]] - fn[(idx[0]-1)%N0, idx[1]] ) / (2*delta)
    dfdy = ( fn[idx[0], (idx[1]+1)%N1] - fn[idx[0], (idx[1]-1)%N1] ) / (2*delta)
    return np.array([dfdx, dfdy])

def fnorm(fn):
    buf = fn.flatten()
    nrm = math.sqrt(buf.dot(buf))
    return nrm

def fconvolve(farr):
    N0 = farr.shape[0]
    N1 = farr.shape[1]
    buf = np.zeros(farr.shape, dtype=complex)
    for i0 in range(N0):
        a0 = 2.0*math.pi*i0/N0
        c0 = math.cos(a0)
        for i1 in range(N1):
            a1 = 2.0*math.pi*i1/N0
            c1 = math.cos(a1)
            if i0 != 0 or i1 != 0:
                q = 0.5 / (c0 + c1 - 2.0)
                buf[i0,i1] = farr[i0,i1] * q
    return buf

def solve_poisson(src):
    buf = fftn(src)
    buf = fconvolve(buf)
    buf = ifftn(buf)
    return buf.real

def GS_Iteration(x, y, D, O, b):
    for i in range(x.shape[0]):
        for j in range(x.shape[1]):
            x[i,j] = (b[i,j] - O(i,j,x,y)) / D(i,j,y)

def residual(x, y, D, O, b):
    r = np.zeros(x.shape)
    for i in range(x.shape[0]):
        for j in range(x.shape[1]):
            r[i,j] = O(i,j,x,y) + D(i,j,y)*x[i,j] - b[i,j]
    return r

def Gauss_Seidel(x, y, D, O, b, maxiter=2000, tol=1.0E-02, talk=0):

    L2b = fnorm(b)
    for itn in range(maxiter):
        xprev = x.copy()
        GS_Iteration(x, y, D, O, b)
        r = residual(x, y, D, O, b)
        L2r = fnorm(r) / L2b
        if talk > 0 and itn % talk == 0:
            print "Iteration # %d, L2 of residual = %10.3E" % (itn, L2r)
        if L2r <= tol: break

    return (L2r, itn)

def bc_enforce_D(x, i, j):
    if i < 0:
        return -x[0,j]
    elif i >= x.shape[0]:
        return -x[-1,j]
    elif j < 0:
        return -x[i,0]
    elif j >= x.shape[1]:
        return -x[i,-1]
    else:
        return x[i,j]

def bc_enforce_N(x, i, j):
    ii, jj = i, j
    if i < 0: ii = 0
    if i >= x.shape[0]: ii = x.shape[0]-1
    if j < 0 : jj = 0
    if j >= x.shape[1]: jj = x.shape[1]-1

    return x[ii,jj]

def bc_enforce_P(x, i, j):

    return x[i % x.shape[0], j % x.shape[1]]

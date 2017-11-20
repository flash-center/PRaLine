

'''
Provides function that will plot the reconstructed magnetic field caluctlated from
the alogrithm in recconstruct.py
'''

import math
import sys
import os.path

import rad_ut as ru
from constants import M_PROTON_G, ESU, C, V_PER_E

import matplotlib as mpl
mpl.use('Agg') # Headless plotting (avoids python-tk GUI requirement)
from matplotlib.ticker import FixedLocator
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import re


def magnetic_field(Br):
    '''
    The goal is to calculate the magnetic field per bin

    Parameters
    ----------
    Br (2D array of (x,y)): B Field per (x,y)

    Returns
    -------
    BrMag (2D array): Log reconstructed B Field per bin
    '''
    num_bins = Br.shape[0]  # num_bins x num_bins
    BrMag = np.zeros((num_bins, num_bins))
    for i in range(num_bins):
        for j in range(num_bins):
            # Field Strength on a logarithmic scale
            BrMag[i, j] = 0.5 * math.log10(Br[i, j, 0]**2 + Br[i, j, 1]**2)
    return BrMag


def B_plot(B, flux_ref, bin_um, type, title):
    '''
    Genereates the  B perpendicular Projection

    Parameters
    ----------
    B (2D array of (x,y)): B Field per (x,y)
    flux_ref (2D array): Number protons per bin without an interaction region
    bin_um (float): Length of the side of a bin, in cm
    type (string): the type of input
    title (string): the title of the plot

    Returns
    -------
    B_recon.png (image):  B perpendicular Projection Plot
    '''
    print (r"Constructing Log " + title + " Projection Plot")
    font = {'family': 'serif',
            'color':  'black',
            'weight': 'normal',
            'size': 32,
            }
    X, Y = ru.position(flux_ref, bin_um)
    BMag = magnetic_field(B)
    stretch = 13.1 / 10.2
    #plt.rc('text', usetex=True)

    #################################################
    fig = plt.figure()
    fig.set_figwidth(6.0 * stretch)
    fig.set_figheight(6.0)
    ax = fig.add_subplot(1, 1, 1)
    vmin = BMag.T.min()
    vmax = BMag.T.max()
    norm = mpl.colors.Normalize(vmin= vmin,vmax= vmax)

    strm = ax.streamplot(X[:, 0], Y[0, :], B[:, :, 0].T, B[:, :, 1].T, color=BMag.T,
                         linewidth=2, cmap=cm.RdYlGn, density=2.0, arrowsize=2.0, norm = norm)
    fig.colorbar(strm.lines)
    #################################################

    xmin = round(min(X[:, 0]), 1)
    xmax = round(max(X[:, 0]), 1)
    ymin = round(min(Y[0, :]), 1)
    ymax = round(max(Y[0, :]), 1)

    ax.set_xlim(int(xmin) - 0.5, int(xmax) + 0.5)
    ax.set_ylim(int(ymin) - 0.5, int(ymax) + 0.5)

    if type == 'carlo':
        x = "Carlo"
    elif type == 'flash4':
        x = "Flash"
    elif type == 'mitcsv':
        x = 'MITCSV'

    ax.set(title=x + ": Log " + title + r" $B_\perp$ Projection (G cm)",
           ylabel=r"Y (cm)", xlabel=r"X (cm)")
    ax.tick_params(labelsize='large')

    fig.savefig("B_" + title + ".png", format='png')

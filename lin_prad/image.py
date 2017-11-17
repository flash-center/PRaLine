
'''
Provides tools to analyze the inputted data from a proton radiogtaphy experiment
'''

import math
import sys

import rad_ut as ru
from constants import M_PROTON_G, ESU, C, V_PER_E

import matplotlib
matplotlib.use('Agg') # Headless plotting (avoids python-tk GUI requirement)
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np



def hist2D_plot(array, bin_um, type, title):
    '''
    Genereates the 2D histogram plot based on 2D array

    Parameters
    ----------
    array (2D array): The array that will be plotted on a 2D histogram
    bin_um (float): Length of the side of a bin, in cm
    type (string): type of file input
    title (string): title of the plot


    Returns
    -------
    "title".png (image):2D histogram Plot
    '''
    print "Constructing " + title + " Plot"
    font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 32,
        }

    num_bins = array.shape[0]
    delta = bin_um/10000.0
    dmax = (delta * num_bins)/2.0

    x = np.zeros((num_bins+1, num_bins+1))
    y = np.zeros((num_bins+1, num_bins+1))
    for i in range(num_bins+1):
        xx = -dmax + i*delta
        for j in range(num_bins+1):
            yy = -dmax + j*delta
            x[i,j] = xx
            y[i,j] = yy

    # Intiating plot
    fig = plt.figure()
    fig.set_figwidth(26)
    fig.set_figheight(12.0)

    # Making plot
    ax = fig.add_subplot(1,1,1)
    p = ax.pcolormesh(x,y, array, cmap=cm.afmhot, vmin= array.min(), vmax=array.max())
    ax.set_xlabel("X (cm)", fontdict=font)
    ax.set_ylabel("Y (cm)", fontdict=font)

    xmin = round(x.min(),1)
    xmax = round(x.max(),1)
    ymin = round(y.min(),1)
    ymax = round(y.max(),1)

    ax.set_xlim(int(xmin)- 0.5, int(xmax)+0.5)
    ax.set_ylim(int(ymin)- 0.5, int(ymax)+0.5)

    plt.colorbar(p)
    if type == 'carlo':
        x = "Carlo"
    elif type == 'flash4':
        x = "Flash"
    elif type == 'mitcsv':
        x = 'MITCSV'
        
    ax.set_title(x + ": " + title,fontdict=font)
    fig.savefig(title+".png", format='png')

def err2D_plot(array, bin_um, type, title):
    '''
    Genereates the 2D histogram plot based on 2D array

    Parameters
    ----------
    array (2D array): The array that will be plotted on a 2D histogram
    bin_um (float): Length of the side of a bin, in cm
    type (string): type of file input
    title (string): title of the plot


    Returns
    -------
    "title".png (image):2D histogram Plot
    '''
    print "Constructing Counts/Bin Plot"
    font = {'family': 'serif',
        'color':  'black',
        'weight': 'normal',
        'size': 32,
        }
    x,y = ru.position(array, bin_um)
    fig = plt.figure()
    fig.set_figwidth(26)
    fig.set_figheight(12.0)

    # Counts/Bin
    vm = max(abs(array.min()), abs(array.max()))
    ax = fig.add_subplot(1,1,1)
    p = ax.pcolormesh(x,y, array, cmap = cm.RdYlGn, vmin = -vm,
                      vmax = vm)
    ax.set_xlabel("X (cm)",fontdict=font)
    ax.set_ylabel("Y (cm)",fontdict=font)

    xmin = round(x.min(),1)
    xmax = round(x.max(),1)
    ymin = round(y.min(),1)
    ymax = round(y.max(),1)

    ax.set_xlim(int(xmin)- 0.5, int(xmax)+0.5)
    ax.set_ylim(int(ymin)- 0.5, int(ymax)+0.5)

    plt.colorbar(p)
    if type == 'carlo':
        x = "Carlo"
    elif type == 'flash4':
        x = "Flash"
    elif type == 'mitcsv':
        x = 'MITCSV'
    ax.set_title(x + ": " + title + " Noise",fontdict=font)
    fig.savefig(title+" Noise.png", format='png')

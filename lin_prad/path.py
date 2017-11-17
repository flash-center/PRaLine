'''
This script is intended to calculate the actual Perpendicular Magntetic Field
of the Proton Radiography simulation
'''
import sys
import math
from re import match

import rad_ut as ru
from constants import M_PROTON_G, ESU, C, V_PER_E

import pandas as pd
import numpy as np


def mag_parse(fname, bin_um):
    '''
    Parses input file and Returns the 2D array relevant to the actual magnetic
    field for verfication purposes

    Parameters
    ----------
    fn(string): full filename (including path) of the proton detector file; e.g. "/home/myouts/blob.out", where "blob.out" is the basename
    bin_um(float): size of the square edge lengths with which to divide the detector for binning

    Returns
    -------
    Bperp(2D array of (x,y) tuple): Magnetic Perpendicular Integral
    '''
    # Data file
    fd = open(fname, 'r')

    line = fd.readline()
    while not match('^# Columns:', line):

        if match('^# rs:', line):
            s2d_cm = float(line.split()[2])

        if match('^# ri:', line):
            s2r_cm = float(line.split()[2])

        if match('^# raperture:', line):
            rap = float(line.split()[2])

        line = fd.readline()

    while match('^#', line):
        line = fd.readline()

    radius = rap * s2d_cm / s2r_cm  # radius of undeflected image of aperture at screen
    ru.dmax = 0.98 * radius / math.sqrt(2.0)  # half the width of the detector
    # number of bins per dimensions
    nbins = int(ru.dmax * 2 / (bin_um / 10000.0))
    ru.delta = 2.0 * ru.dmax / nbins  # width of a bin

    flux = np.zeros((nbins, nbins))  # num. of protons per bin
    Bperp = np.zeros((nbins, nbins, 2))  # B Integral
    J = np.zeros((nbins, nbins))

    nprot = 0
    while line:
        nprot += 1
        xx = float(line.split()[3])
        yy = float(line.split()[4])
        jj = float(line.split()[8])
        b0 = float(line.split()[9])
        b1 = float(line.split()[10])
        idx = ru.vec2idx((xx, yy))
        i = idx[0]
        j = idx[1]

        if (xx + ru.dmax) / ru.delta >= 0 and i < nbins and (yy + ru.dmax) / ru.delta >= 0 and j < nbins:
            flux[i, j] += 1
            Bperp[i, j, 0] += b0
            Bperp[i, j, 1] += b1
            J[i, j] += jj

        line = fd.readline()

    fd.close()

    avg_fluence = nprot / (math.pi * radius**2)
    im_fluence = flux.sum() / (4 * ru.dmax**2)

    for i in range(nbins):
        for j in range(nbins):
            try:
                Bperp[i, j, :] /= flux[i, j]
                J[i, j] /= flux[i, j]
            except ZeroDivisionError:
                print "Zero pixel, will screw everything up."
                raise ValueError

    return Bperp, J, avg_fluence, im_fluence


def parse_im(fname):
    # Data file
    fd = open(fname, 'r')

    line = fd.readline()
    while match('#', line):

        if match('# s2r_cm', line):
            s2r_cm = float(line.split()[2])

        if match('# s2d_cm', line):
            s2d_cm = float(line.split()[2])

        if match('# Ep_MeV', line):
            Ep_MeV = float(line.split()[2])

        if match('# bin_um', line):
            bin_um = float(line.split()[2])

        if match('# x-mask', line):
            x = float(line.split()[2])

        if match('# y-mask', line):
            y = float(line.split()[2])

        line = fd.readline()

    data = pd.read_csv(fname, header=None, comment='#',
                       delimiter=",").as_matrix()

    index = int(data.shape[1] / 3)
    num_bins = int(math.sqrt(index))

    flux2D = data[:, :index]
    flux2D.shape = (num_bins, num_bins)

    flux2D_ref = data[:, index:(index * 2)]
    flux2D_ref.shape = (num_bins, num_bins)

    mask = data[:, (index * 2):]
    mask.shape = (num_bins, num_bins)

    return flux2D, flux2D_ref, mask, s2r_cm, s2d_cm, Ep_MeV, bin_um

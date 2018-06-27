

'''
Acts as a wrapper that runs the anlysis of the input and plots the data in this project
'''
import sys
import math

from pradreader import reader
import algorithm as algo
import Bplot2 as plot
import image
import path
import rad_ut as ru

import numpy as np
import pandas as pd
import argparse as ap

def get_input_data():
        '''
        Command line options and variables
        '''
        parser = ap.ArgumentParser(
            description="This script is used to reconstruct the magnetic field of "
            "Proton Radiography experiment")
        parser.add_argument("-v", "--verbose", help="increase output verbosity",
                            action="store_true")
        parser.add_argument("input_file", type=str,
                            help="The filename including the path")
        #TODO: Implement masking tool
        # parser.add_argument("--x1", default=0, type=int,
        #                     help="the first percentage of the x interval e.g 10 percent DEFAULT:0")
        # parser.add_argument("--y1", default=0, type=int,
        #                     help="the first percentage of the y interval e.g 10 percent. DEFAULT:0")
        # parser.add_argument("--x2", default=0, type=int,
        #                     help="the latter percentage of the x interval e.g 30 percent. DEFAULT:0")
        # parser.add_argument("--y2", default=0, type=int,
        #                     help="the latter percentage of the y interval e.g 30 percent. DEFAULT:0")

        args = parser.parse_args()

        return args


def prad_wrap():
    '''
    Wrapper Function for Command line tool that returns various details and plots
    to analyze the experiment

    Parameters
    ----------
    filename(required): including path
    rtype(required): carlo, mitcsv, flash4
    bin_um(required): length of the bin in microns

    Returns
    -------
    files (string): flux and fluence contrast plots and other various plots if
    path integrated data is available.
    '''
    print "STARTING ANALYSIS AND PLOTTING..."
    # First Parameter: Path name of the file
    # Second Parameter: The type of experimental output
    # Input variables and options
    args = get_input_data()
    fn = args.input_file
    pr = reader.loadPRRp(fn)
    rtype = pr.rtype
    flux = pr.flux2D
    flux_ref = pr.flux2D_ref
    sr2_cm = pr.s2r_cm
    s2d_cm = pr.s2d_cm
    Ep_MeV = pr.Ep_MeV
    bin_um = pr.bin_um
    flux = flux.T
    flux_ref = flux_ref.T
    print("\n")

    flux_min =10.0
    # Protons per bin 2D Histogram
    image.hist2D_plot(flux, bin_um, rtype,"Flux")
    print "Mean counts per bin: %12.5E ; Std. Dev. Counts per bin: %12.5E" % (flux.mean(), flux.std())
    print "Max counts per bin: %d ; Min counts per bin: %d" % (flux.max(), flux.min())
    print "Number of bins with zero protons: %d" % (flux.size - flux[ flux>0 ].size)
    print "Number of bins with %d or fewer protons: %d\n" % (flux_min, flux.size - flux[ flux>flux_min ].size)

    # Fluence Distrubtion of protons at the screen 2D Histogram
    Src, fluc = algo.steady_state(flux, flux_ref)
    image.hist2D_plot(fluc, bin_um, rtype, "Fluence")
    Flpos = fluc[flux >= flux_min]
    print "Mean Fluct.: %12.5E ; Std. Dev. Fluct.: %12.5E" % (Flpos.mean(), Flpos.std())
    print "Max Fluct: %12.5E ; Min Fluct: %12.5E" % (Flpos.max(), Flpos.min())

if __name__=="__main__":
    prad_wrap()

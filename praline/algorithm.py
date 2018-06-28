
'''
Many of the variables  and equations are in reference to the paper,
Inferring Morphology and Strength of Magnetic Fields From Proton Radiographs,
This script will reconstruct magnetic fields given the data
from a Proton Radiography experiment.
'''
import sys
import math

import rad_ut as ru
from constants import M_PROTON_G, ESU, C, V_PER_E

from re import match
import numpy as np

def b_field(s2r_cm, s2d_cm, Ep_MeV):
    '''
    Calculates the Uniform Magnetic field.

    Parameters
    ----------
    s2r_cm (float): Distance from the proton source to the detector, in cm
    s2d_cm (float): Distance from the proton source to the interaction region, in cm
    Ep_MeV (float): Kinetic energy

    Returns
    -------
    Bconst (float): Calculates the B, uniform magnetic field strength
    '''
    v = math.sqrt(2 * (Ep_MeV * V_PER_E) /
                  M_PROTON_G)  # Velocity of Proton

    Bconst = M_PROTON_G * C * v / \
        (ESU * (s2r_cm - s2d_cm))  # Uniform B Field Strength

    return Bconst

def steady_state(flux, flux_ref):
    '''
    Obtain the steady-state diffusion equation using Taylor Expansion
    This version does not assume flux and flux_ref are equal size in X and Y
    Uses NumPy broadcasting
    Re-written by Scott Feister 2018-06-27
    
    Parameters
    ----------
    flux (2D array): Number of protons per bin
    flux_ref (2D array): Number protons per bin without an interaction region

    Returns
    -------
    Lam (2D array): fluence contrast
    Src (2D array): Source term from multiplying the fluence contrast and exp(fluence contrast)
    '''
    Lam = np.zeros(flux.shape)

    ct = (flux_ref > 0) & (flux > 0) # Condition ensuring regions of zero flux are not computed
    
    # Taylor expansion
    Lam[ct] = 2.0 * ( 1.0 - np.sqrt(flux_ref[ct]/flux[ct]))
    
    # Source Term
    # RHS of the Steady-State Diffusion Equation
    Src = Lam * np.exp(Lam)

    return (Src, Lam)

def B_recon(flux, flux_ref, s2r_cm, s2d_cm, bin_um, Ep_MeV, tol_iter, max_iter):
    '''
    Produces a reconstructed magnetic field
    Assumes flux and flux_ref are equal size in X and Y
    
    Parameters
    ----------
    flux (2D array): Number of protons per bin
    flux_ref (2D array): Number protons per bin without an interaction region
    s2r_cm (float): Distance from the proton source to the detector, in cm
    s2d_cm (float): Distance from the proton source to the interaction region, in cm
    bin_um (float): Length of the side of a bin, in cm
    Ep_MeV (float): Kinetic Energy

    Returns
    -------
    BperpR (3D array of (x,y,2)): Reconstructed Magnetic Field components
    '''
    ru.delta = bin_um / 10000.0

    num_bins = flux_ref.shape[0]  # num_bins x num_bins
    # RHS of the Steady-State Diffusion Equation and Fluence Contrast
    Src, Lam = steady_state(flux, flux_ref)
    # The real component after Lam is transformed then convolved and then inversely transformed
    phi = ru.solve_poisson(Lam)
    # Uniform B Field Strength
    Bconst = b_field(s2r_cm, s2d_cm, Ep_MeV)
    # Iterate to solution
    print ("Gauss-Seidel Iteration...")
    GS = ru.Gauss_Seidel(phi, np.exp(Lam), ru.D, ru.O, Src,
                         talk=20, tol=tol_iter, maxiter=max_iter)
    # Multiplying by the area of the bin
    phi *= (ru.delta**2)
    # Reconstructed perpendicular B Fields
    BperpR = np.zeros((num_bins, num_bins, 2))
    # Reconstructed Lateral motion of proton
    deltaXR = np.zeros((num_bins, num_bins, 2))

    for i in range(num_bins):
        for j in range(num_bins):
            # Reconstructed Data
            deltaXR[i, j] = -ru.gradient(phi, (i, j))
            BperpR[i, j, 0] = Bconst * deltaXR[i, j, 1]
            BperpR[i, j, 1] = -Bconst * deltaXR[i, j, 0]

    print ("#L2 norm of residual = %12.5E ;  Number of Gauss-Seidel iterations = %d\n" % GS)
    return BperpR

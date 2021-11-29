# -*- coding: utf-8 -*-
"""
Created on Mon Nov 29 12:23:49 2021

@author: gabri
"""
import kwant
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

def make_low_energy_model(mu,t, u, Delta0, Delta1, Delta11):
    """
    Minimal low energy level for TRIPTOPS based on [Haim] equation 9.

    Parameters
    ----------
    mu : float
        Chemical potential.
    t : float
        Hopping parameter.
    u : float
        Spin-orbit coupling.
    Delta0 : float
        On-site paring potential.
    Delta1 : float
        Nearest neighbor singlet pairing potential.
    Delta11 : float
        Nearest neighbor triplet-component pairing potential.
        
    Returns
    -------
    syst : kwant.builder.Builder
        The representation for the minimal model.
        
    """
    sym = kwant.TranslationalSymmetry((-1,))
    syst = kwant.Builder(sym)
    lat = kwant.lattice.chain(norbs=2)      #two orbitals per site, one for each spin
    def onsite(site, mu):
        return -mu*np.eye(2)
    syst[lat(0)] = onsite
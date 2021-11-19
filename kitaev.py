# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 13:24:37 2021

@author: gabri
"""
import kwant
import numpy as np

def make_kitaev_chain_infinite(t, mu, Delta):
    """
    Create a an infinite translationally invariant Kitaev chain.
    """
    # Start with an empty tight-binding system. On each site, there
    # are now electron and hole orbitals, so we must specify the
    # number of orbitals per site 'norbs'.
    kitaev_chain_infinite = kwant.Builder(kwant.TranslationalSymmetry((-1,)))
    lat = kwant.lattice.chain(norbs=2)  
    # The superconducting order parameter couples electron and hole orbitals
    # on each site, and hence enters as an onsite potential
    kitaev_chain_infinite[lat(0)] = -mu * np.array([[1, 0], [0, -1]])
    # Hoppings
    kitaev_chain_infinite[kwant.HoppingKind((1,), lat)] = -t * np.array([[1, 0],[0, -1]]) - Delta * np.eye(2)
    
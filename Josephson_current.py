# -*- coding: utf-8 -*-
"""
Created on Mon Dec 13 10:55:20 2021

@author: gabri
"""
import kwant
import numpy as np
import matplotlib.pyplot as plt

#Pauli matrices
tau_z = np.array([[1, 0],[0, -1]])
tau_y = np.array([[0, -1j],[1j, 0]])

def onsite(onsite, mu):
    return -mu * tau_z

def hopping(site1, site2, t, Delta):
    return -t * tau_z - 1j * Delta * tau_y

def hopping_contact(site1, site2, t, phi):
    return -t * tau_z * np.exp(1j*np.pi*phi)

def make_Josephson_junction(t=1, mu=0, Delta=1, L=25, phi=1):
    """
    Create a 1D tight-binding model for
    Josephson's junction with magnetic flux.

    Parameters
    ----------
    t : float, optional
        Hopping parameter. The default is 1.
    mu : float optional
        Chemical potential. The default is 1.
    Delta : float, optional
        Superconducting gap. The default is 1.
    L : int, optional
        Number of sites in the chain. The default is 25.
    phi : float, optional
        Magnetic flux in units of quantum flux h/(2e).

    Returns
    -------
    kwant.builder.Builder
        The representation for the tight-binding model.

    """
    kitaev_chain_finite = kwant.Builder()
    lat = kwant.lattice.chain(norbs=2)  
    # The superconducting order parameter couples electron and hole orbitals
    # on each site, and hence enters as an onsite potential
    # There are L sites in each superconductor and a contact in zeroth position,
    kitaev_chain_finite[(lat(x) for x in range(-L, L+1))] = onsite
    # Hoppings
    kitaev_chain_finite[lat.neighbors()] = hopping
    kitaev_chain_finite[lat(0), lat(1)] = hopping_contact
    #kwant.plot(kitaev_chain_finite)
    return kitaev_chain_finite

phi = np.linspace(-np.pi, np.pi, 1000)
def Josephson_current(phi):
    fundamental_energy = []
    for phi in phi:
        params = dict(t=1, mu=0, Delta=1, L=25, phi=phi)
        kitaev = make_Josephson_junction(t=1, mu=0, Delta=1, L=25, phi=phi)
        kitaev = kitaev.finalized()
        H = kitaev.hamiltonian_submatrix(params=params)
        eigenvalues, eigenvectors = np.linalg.eig(H)
        fundamental_energy.append(np.sum(eigenvalues) / 2)
    current = np.diff(fundamental_energy)
    return current

current = Josephson_current(phi)
fig, ax = plt.subplots()
ax.plot(phi[:-1], current)
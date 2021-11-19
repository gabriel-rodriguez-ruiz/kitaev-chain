# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 13:24:37 2021

@author: gabri
"""
import kwant
import numpy as np
import matplotlib.pyplot as plt

#Pauli matrices
tau_z = np.array([[1, 0],[0, -1]])
tau_y = np.array([[0, -1j],[1j, 0]])
tau_x = np.array([[0, 1],[1, 0]])

def make_kitaev_chain_finite(t=1, mu=1, Delta=1, L=25):
    """
    Create a finite translationally invariant Kitaev chain.
    """
    # Start with an empty tight-binding system. On each site, there
    # are now electron and hole orbitals, so we must specify the
    # number of orbitals per site 'norbs'.
    kitaev_chain_finite = kwant.Builder()
    lat = kwant.lattice.chain(norbs=2)  
    # The superconducting order parameter couples electron and hole orbitals
    # on each site, and hence enters as an onsite potential
    def onsite(onsite, mu):
        return -mu * tau_z
    kitaev_chain_finite[(lat(x) for x in range(L))] = onsite
    # Hoppings
    def hop(site1, site2, t, Delta):
        return -t * tau_z - 1j * Delta * tau_y
    kitaev_chain_finite[lat.neighbors()] = hop   
    return kitaev_chain_finite    

#kitaev = make_kitaev_chain_finite().finalized()
#energies = np.linalg.eigvals(kitaev.hamiltonian_submatrix())

def plot_spectrum(syst, mu):
    fig = kwant.plotter.spectrum(syst, ("mu", mu), params=dict(t=1, Delta=1))
    fig.canvas.manager.set_window_title("Bandas de energía")    
    plt.xlabel(r"$\frac{\mu}{t}$")
    plt.ylabel(r"Energía")

def main():
    syst = make_kitaev_chain_finite()
    # Check that the system looks as intended.
    #kwant.plot(syst)
    
    # Finalize the system.
    syst = syst.finalized()
        
    mu = np.linspace(-4, 4)
    # We should see the energy bands.
    plot_spectrum(syst, mu)
    
if __name__ == '__main__':
    main()

    
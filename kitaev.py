# -*- coding: utf-8 -*-
"""
Created on Fri Nov 19 13:24:37 2021

@author: gabri
"""
import kwant
import numpy as np
import matplotlib.pyplot as plt

def make_kitaev_chain_finite(t=1, mu=1, Delta=1, L=25):
    """
    Create a finite translationally invariant Kitaev chain.

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

    Returns
    -------
    kwant.builder.Builder
        The representation for the Kitaev chain.

    """
    #Pauli matrices
    tau_z = np.array([[1, 0],[0, -1]])
    tau_y = np.array([[0, -1j],[1j, 0]])
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

def plot_spectrum(syst, mu):
    """
    Plot the spectrum by changing the parameter 'mu'.

    Parameters
    ----------
    syst : kwant.builder.FiniteSystem
        Finite Kitaev chain.
    mu : np.array
        Values for the chemical potential.

    Returns
    -------
    fig : plt.figure

    """
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
    
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
    Create a finite Kitaev chain.

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

def make_kitaev_chain_infinite(t=1, mu=1, Delta=1):
    """
    Create an infinite translationally invariant Kitaev chain.

    Parameters
    ----------
    t : float, optional
        Hopping parameter. The default is 1.
    mu : float optional
        Chemical potential. The default is 1.
    Delta : float, optional
        Superconducting gap. The default is 1.

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
    lat = kwant.lattice.chain(norbs=2)
    # The symetry
    sym = kwant.TranslationalSymmetry((-1,))
    kitaev_chain_infinite = kwant.Builder(sym)

    # The superconducting order parameter couples electron and hole orbitals
    # on each site, and hence enters as an onsite potential
    def onsite(site, mu):
        return -mu * tau_z
    kitaev_chain_infinite[lat(0)] = onsite
    # Hoppings
    def hop(site1, site2, t, Delta):
        return -t * tau_z - 1j * Delta * tau_y
    kitaev_chain_infinite[kwant.HoppingKind((1,), lat)] = hop   
    return kitaev_chain_infinite    

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
    fig.canvas.manager.set_window_title("Bandas de energía de una cadena finita")    
    plt.xlabel(r"$\frac{\mu}{t}$")
    plt.ylabel(r"Energía")

def main():
    kitaev_finite = make_kitaev_chain_finite()
    # Check that the system looks as intended.
    #kwant.plot(syst)
    
    # Finalize the system.
    kitaev_finite = kitaev_finite.finalized()
        
    mu = np.linspace(-4, 4)
    # We should see the energy bands.
    plot_spectrum(kitaev_finite, mu)
    
    # In the case of an infinite kitaev chain
    kitaev_infinite = make_kitaev_chain_infinite()
    kitaev_infinite = kitaev_infinite.finalized()
    fig = kwant.plotter.bands(kitaev_infinite, params=dict(t=1, mu=1, Delta=1))
    fig.canvas.manager.set_window_title("Bandas de energía de una cadena finita")    
    
if __name__ == '__main__':
    main()
    
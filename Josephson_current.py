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

def hopping(site1, site2, t, Delta, phi):
    if (site1.pos == [0.0] and site2.pos == [-1.0]) or (site1.pos == [-1.0] and site2.pos == [0.0]):
        return -t * tau_z * np.exp(1j*phi)
    else:
        return -t * tau_z - 1j * Delta * tau_y
    
def hopping_all_gauge(site1, site2, t, Delta, phi):
    if (site1.pos == [0.0] and site2.pos == [-1.0]) or (site1.pos == [-1.0] and site2.pos == [0.0]):
        return -t * tau_z
    else:
        if site2.pos[0]<-1:
            return -t * tau_z + - 1j * Delta * tau_y * np.exp(-1j*phi)
        else:
            return -t * tau_z + - 1j * Delta * tau_y * np.exp(1j*phi)

def make_Josephson_junction(t=1, mu=0, Delta=1, L=25, phi=0):
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
    # There are L sites in each superconductor
    kitaev_chain_finite[(lat(x) for x in range(-L, L))] = onsite
    # Hoppings
    kitaev_chain_finite[lat.neighbors()] = hopping_all_gauge
    # The hopping in the contact has a phase dependence.
    #kitaev_chain_finite[lat(-1), lat(0)] = hopping_contact
    #kitaev_chain_finite[lat(0), lat(-1)] = hopping_contact
    #kwant.plot(kitaev_chain_finite)
    return kitaev_chain_finite

def Josephson_current(syst, params):
    fundamental_energy = []
    for phi in params["phi"]:
        params["phi"] = phi
        H = syst.hamiltonian_submatrix(params=params)
        eigenvalues, eigenvectors = np.linalg.eig(H)
        fundamental_energy.append(np.sum(eigenvalues, where=eigenvalues>0) / 2)
    #print(fundamental_energy)
    current = np.diff(fundamental_energy)
    return current

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
    fig = kwant.plotter.spectrum(syst, ("mu", mu), params=dict(t=1, Delta=1, phi=2))
    fig.canvas.manager.set_window_title("Bandas de energía de una cadena finita")    
    plt.xlabel(r"$\frac{\mu}{t}$")
    plt.ylabel(r"Energía")

def site_color(site):
    if site.pos==[0.0] or site.pos==[-1.0]:
        return 'red'
    else:
        return 'blue'
    
def hop_color(site1, site2):
    #print(site1.pos, site2.pos)
    if site1.pos==[0.0] and site2.pos==[-1.0]:
        return 'red'
    else:
        return 'blue'
    
mu = np.linspace(-4, 4)
kitaev = make_Josephson_junction()
kwant.plot(kitaev, site_color=site_color, hop_color=hop_color)
kitaev = kitaev.finalized()
params = dict(t=1, mu=0, Delta=1, L=25, phi=0)
plot_spectrum(kitaev, mu)
# We should see the energy bands

phi = np.linspace(-np.pi, np.pi, 1000)
params["phi"] = phi
current = Josephson_current(kitaev, params)
#energy = Josephson_current(kitaev, params)
plt.figure("Corriente Josephson", clear=True)
plt.title("Corriente Josephson")
plt.xlabel(r"$\frac{\Phi}{\Phi_0}$")
plt.ylabel("J")
plt.xticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi],
           [r"$-\pi$", r"$-\frac{\pi}{2}$", "0", r"$\frac{\pi}{2}$", "$\pi$"])
plt.grid()    
plt.plot(phi[:-1], current)

#fig, ax = plt.subplots()
#ax.plot(phi[:-1], current)
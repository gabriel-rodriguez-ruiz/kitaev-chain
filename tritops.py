# -*- coding: utf-8 -*-
"""
Created on Mon Nov 29 12:23:49 2021

@author: gabri
"""
import kwant
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

sigma_x = np.array([[0, 1], [1, 0]])
sigma_y = np.array([[0, -1j], [1j, 0]])
sigma_z = np.array([[1, 0], [0, -1]])

tau_x = np.array([[0, 1], [1, 0]])
tau_y = np.array([[0, -1j], [1j, 0]])
tau_z = np.array([[1, 0], [0, -1]])

def make_low_energy_model_finite(mu=0, t=10, u=-5j, Delta_0=0, Delta_1=4, Delta_11=0, L=25):
    """
    Minimal low energy level for TRITOPS based on [Haim] equation 9.

    Parameters
    ----------
    mu : float
        Chemical potential.
    t : float
        Hopping parameter.
    u : float
        Spin-orbit coupling.
    Delta_0 : float
        On-site paring potential.
    Delta_1 : float
        Nearest neighbor singlet pairing potential.
    Delta_11 : float
        Nearest neighbor triplet-component pairing potential.
    L : int, optional
        Number of sites in the chain. The default is 25.
 
    Returns
    -------
    syst : kwant.builder.Builder
        The representation for the minimal model.
        
    """
    syst = kwant.Builder()
    lat = kwant.lattice.chain(norbs=4)      #four orbitals per site, psi=(c+, c-, c*+, c*-)
    def onsite(site, mu, Delta_0):
        return -mu/2 * np.kron(sigma_z, np.eye(2)) + 1/2j * Delta_0 * np.kron(sigma_x, sigma_y)
    syst[(lat(x) for x in range(L))] = onsite
    def hopping(site1, site2, t, u, Delta_1, Delta_11):     #np.kron() is the tensor product
        return (-t * np.kron(sigma_z, np.eye(2)) + 1j * u * np.kron(sigma_z, sigma_z)
                + 1/2 * (Delta_1 + Delta_11) * np.kron(sigma_x, 1j*sigma_y))
    syst[lat.neighbors()] = hopping
    return syst

def make_low_energy_model_infinite(mu=0, t=10, u=-5j, Delta_0=0, Delta_1=4, Delta_11=0):
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
    Delta_0 : float
        On-site paring potential.
    Delta_1 : float
        Nearest neighbor singlet pairing potential.
    Delta_11 : float
        Nearest neighbor triplet-component pairing potential.
        
    Returns
    -------
    syst : kwant.builder.Builder
        The representation for the minimal model.
        
    """
    sym = kwant.TranslationalSymmetry((-1,))
    syst = kwant.Builder(sym)
    lat = kwant.lattice.chain(norbs=4)      #four orbitals per site, psi=(c+, c-, c*+, c*-)
    def onsite(site, mu, Delta_0):
        return -mu/2 * np.kron(sigma_z, np.eye(2)) + 1/2j * Delta_0 * np.kron(sigma_x, sigma_y)
    syst[lat(0)] = onsite
    def hopping(site1, site2, t, u, Delta_1, Delta_11):     #np.kron() is the tensor product
        return (-t * np.kron(sigma_z, np.eye(2)) + 1j * u * np.kron(sigma_z, sigma_z)
                + 1/2 * (Delta_1 + Delta_11) * np.kron(sigma_x, 1j*sigma_y))
    syst[kwant.HoppingKind((1,), lat)] = hopping
    return syst

def plot_spectrum(syst, mu, params):
    """
    Plot the spectrum by changing the parameter 'mu'.

    Parameters
    ----------
    syst : kwant.builder.FiniteSystem
        Finite Kitaev chain.
    mu : np.array
        Values for the chemical potential.
    params : dict
        The rest of the parameters, which will be kept constant.
    
    Returns
    -------
    fig : plt.figure

    """
    fig = kwant.plotter.spectrum(syst, ("mu", mu), params=params)
    fig.canvas.manager.set_window_title("Bandas de energía de un triptop finito")    
    plt.xlabel(r"$\mu$")
    plt.ylabel(r"Energía")
    plt.ylim((-6,6))

# The function to be called anytime a slider's value changes
def update(value):
    global fig, tritops_infinite
    u = -1j*value       #value is Rashba parameter
    new_bands = kwant.physics.Bands(tritops_infinite, params=dict(mu=0, t=10, u=u, Delta_0=0, Delta_1=4, Delta_11=0))
    i = 0
    for line in fig.axes[0].lines:
        line.set_ydata([new_bands(k)[i] for k in line.get_xdata()])
        i += 1
    fig.canvas.draw()
    fig.canvas.flush_events()

def main():
    global tritops_infinite, tritops_finite, Rashba_slider, fig
    
    tritops_finite = make_low_energy_model_finite(L=50)
    # Check that the system looks as intended.
    #kwant.plot(syst)
    
    # Finalize the system.
    tritops_finite = tritops_finite.finalized()
        
    mu = np.linspace(-20, 20, 1000)
    # We should see the energy bands.
    plot_spectrum(tritops_finite, mu, params=dict(t=10, u=-5j, Delta_0=0, Delta_1=4, Delta_11=0))
    
    tritops_infinite = make_low_energy_model_infinite()
    tritops_infinite = tritops_infinite.finalized()
    fig = kwant.plotter.bands(tritops_infinite, momenta=1000, params=dict(mu=0, t=10, u=-5j, Delta_0=0, Delta_1=4, Delta_11=0))
    fig.canvas.manager.set_window_title("Bandas de energía de un tritops")
    # Make a horizontal slider to control the frequency.
    ax_Rashba = plt.axes([0.25, 0.1, 0.65, 0.03])
    plt.subplots_adjust(left=0.25, bottom=0.25)
    Rashba_slider = Slider(ax=ax_Rashba, label=r'$\lambda_R$', valmin=-10, valmax=10, valinit=5)

if __name__ == '__main__':
    main()
    Rashba_slider.on_changed(update)

# -*- coding: utf-8 -*-
"""
Created on Mon Nov 29 12:23:49 2021

@author: gabri
"""
import kwant
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
import os

my_path = os.path.abspath(__file__) # Figures out the absolute path for you in case your working directory moves around.


sigma_x = np.array([[0, 1], [1, 0]])
sigma_y = np.array([[0, -1j], [1j, 0]])
sigma_z = np.array([[1, 0], [0, -1]])

tau_x = np.array([[0, 1], [1, 0]])
tau_y = np.array([[0, -1j], [1j, 0]])
tau_z = np.array([[1, 0], [0, -1]])

def onsite(site, mu, Delta_0):
    return -mu * np.kron(sigma_z, np.eye(2)) - Delta_0 * np.kron(sigma_y, sigma_y)

def hopping(site1, site2, t, lambda_R, Delta_1):     #np.kron() is the tensor product
    return (-t * np.kron(sigma_z, np.eye(2)) - lambda_R * np.kron(sigma_z, 1j*sigma_y)
            - Delta_1 * np.kron(sigma_y, sigma_y))

def hopping_y(site1, site2, t, lambda_R, Delta_1):     #np.kron() is the tensor product
    return (-t * np.kron(sigma_z, np.eye(2)) + lambda_R * np.kron(1j*np.eye(2), sigma_x)
            - Delta_1 * np.kron(sigma_y, sigma_y))

def make_low_energy_model_finite(mu=0, t=10, lambda_R=5, Delta_0=0, Delta_1=2, L=25):
    """
    Minimal low energy level for TRITOPS based on [Zhang] equation 1
    but in 1D. A factor 1/2 is left out. H=1/2 psi* H_BdG psi

    Parameters
    ----------
    mu : float
        Chemical potential.
    t : float
        Hopping parameter.
    lambda_R : float
        Rashba spin-orbit coupling.
    Delta_0 : float
        On-site paring potential.
    Delta_1 : float
        Nearest neighbor pairing potential.
    L : int, optional
        Number of sites in the chain. The default is 25.
 
    Returns
    -------
    syst : kwant.builder.Builder
        The representation for the minimal model.
        
    """
    syst = kwant.Builder()
    lat = kwant.lattice.chain(norbs=4)      #four orbitals per site, psi=(c+, c-, c*+, c*-)
    syst[(lat(x) for x in range(L))] = onsite
    syst[lat.neighbors()] = hopping
    return syst

def make_low_energy_model_infinite(mu=0, t=10, lambda_R=5, Delta_0=0, Delta_1=2):
    """
    Minimal low energy level for TRITOPS based on [Zhang] equation 1.
    A factor 1/2 is left out. H=1/2 psi* H_BdG psi

    Parameters
    ----------
    mu : float
        Chemical potential.
    t : float
        Hopping parameter.
    lambda_R : float
        Rashba spin-orbit coupling.
    Delta_0 : float
        On-site paring potential.
    Delta_1 : float
        Nearest neighbor pairing potential.
      
    Returns
    -------
    syst : kwant.builder.Builder
        The representation for the minimal model.
        
    """
    sym = kwant.TranslationalSymmetry((-1,))
    syst = kwant.Builder(sym)
    lat = kwant.lattice.chain(norbs=4)      #four orbitals per site, psi=(c+, c-, c*+, c*-)
    syst[lat(0)] = onsite
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
    fig.canvas.manager.set_window_title("Bandas de energía de un tritops finito")    
    plt.xlabel(r"$\mu$")
    plt.ylabel(r"Energía")
    lambda_R = params["lambda_R"]
    plt.xticks([-4*lambda_R, -2*lambda_R, 0, 2*lambda_R, 4*lambda_R], [r"$-4\lambda_R$", r"$-2\lambda_R$",
                                       "0", r"$2\lambda_R$", r"$4\lambda_R$"])
    Delta_1 = params["Delta_1"]
    plt.yticks([-2*Delta_1, -Delta_1, 0, Delta_1, 2*Delta_1], [r"$-2\Delta_1$", r"$-\Delta_1$",
                                   0 ,r"$\Delta_1$", r"$2\Delta_1$"])
    plt.ylim((-3*Delta_1, 3*Delta_1))
    return fig
    
# The function to be called anytime a slider's value changes
def update(value):
    global fig, tritops_infinite
    # retrieve the values from the sliders
    x = [Rashba_slider.val, mu_slider.val]
    new_bands = kwant.physics.Bands(tritops_infinite, params=dict(mu=x[1], t=10, lambda_R=x[0], Delta_0=0, Delta_1=2))
    i = 0
    for line in fig.axes[0].lines:
        line.set_ydata([new_bands(k)[i] for k in line.get_xdata()])
        i += 1
    fig.canvas.draw()
    fig.canvas.flush_events()

def make_low_energy_model_ribbon(mu=0, t=10, lambda_R=5, Delta_0=0, Delta_1=2, L=25):
    """
    Minimal low energy level for TRITOPS based on [Zhang] equation 1 for ribbon of width L.
    A factor 1/2 is left out. H=1/2 psi* H_BdG psi

    Parameters
    ----------
    mu : float
        Chemical potential.
    t : float
        Hopping parameter.
    lambda_R : float
        Rashba spin-orbit coupling.
    Delta_0 : float
        On-site paring potential.
    Delta_1 : float
        Nearest neighbor pairing potential.
    L : float
        Width of the ribbon.
    
    Returns
    -------
    ribbon : kwant.builder.Builder
        The representation for the ribbon model.
        
    """
    #Create a 2D template
    sym = kwant.TranslationalSymmetry((1,0), (0,1))
    syst = kwant.Builder(sym)
    lat = kwant.lattice.square(1, norbs=4)
    syst[lat(0,0)] = onsite
    syst[kwant.HoppingKind((1, 0), lat)] = hopping
    syst[kwant.HoppingKind((0, 1), lat)] = hopping_y
    #Fill the target ribbon inside the template syst
    ribbon = kwant.Builder(kwant.TranslationalSymmetry([1,0]))
    ribbon.fill(syst, shape=(lambda site: 0 <= site.pos[1] < L), start=[0, 0])
    #syst[kwant.HoppingKind((1, 0), lat)] = hopping
    #syst[lat.neighbors()] = hopping 
    return ribbon

def main():
    global tritops_infinite, tritops_finite, Rashba_slider, fig, mu_slider
    
    tritops_finite = make_low_energy_model_finite(L=300)
    # Check that the system looks as intended.
    #kwant.plot(syst)
    
    # Finalize the system.
    tritops_finite = tritops_finite.finalized()
        
    mu = np.linspace(-30, 30)
    # We should see the energy bands.
    params = dict(t=10, lambda_R=5, Delta_0=0, Delta_1=2)
    fig = plot_spectrum(tritops_finite, mu, params=params)
    fig.savefig("../Images/Bandas de un tritops 2D finito.png")
    
    tritops_infinite = make_low_energy_model_infinite()
    tritops_infinite = tritops_infinite.finalized()
    fig = kwant.plotter.bands(tritops_infinite, momenta=1000, params=dict(mu=0, t=10, lambda_R=5, Delta_0=0, Delta_1=2))
    plt.xticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi],
               [r"$-\pi$", r"$-\frac{\pi}{2}$", "0", r"$\frac{\pi}{2}$", "$\pi$"])
    fig.canvas.manager.set_window_title("Bandas de energía de un tritops")
    # Make a horizontal slider to control the frequency.
    plt.subplots_adjust(left=0.25, bottom=0.25)
    ax_Rashba = fig.add_axes([0.25, 0.1, 0.65, 0.03])
    Rashba_slider = Slider(ax=ax_Rashba, label=r'$\lambda_R$', valmin=-30, valmax=30, valinit=5)
    #plt.subplots_adjust(left=0.25, bottom=0.25)
    ax_mu = fig.add_axes([0.25, 0.15, 0.65, 0.03])
    mu_slider = Slider(ax=ax_mu, label=r'$\mu$', valmin=-10, valmax=10, valinit=0)
    fig.savefig("../Images/Bandas de un tritops 2D infinito.png")

def main_2D():
    
    tritops_ribbon = make_low_energy_model_ribbon(L=100)
    # Check that the system looks as intended.
    #kwant.plot(syst)
    
    # Finalize the system.
    tritops_ribbon = tritops_ribbon.finalized()
        
    # We should see the energy bands.
    params = dict(mu=-10, t=10, lambda_R=5, Delta_0=-2, Delta_1=2)
    #mu=-10 si Delta_0=-2
    #mu=10 si Delta_0=2
    fig = kwant.plotter.bands(tritops_ribbon, momenta=1000, params=params)
    #fig.savefig("../Images/Bandas de una cinta tritops Delta0=-2.png")
    plt.title(f"Cinta con mu={params['mu']}, t={params['t']}, lambda_R={params['lambda_R']}, Delta_0={params['Delta_0']}, Delta_1={params['Delta_1']}")
    plt.xticks([-np.pi, -np.pi/2, 0, np.pi/2, np.pi],
               [r"$-\pi$", r"$-\frac{\pi}{2}$", "0", r"$\frac{\pi}{2}$", "$\pi$"])
    Delta_1 = params["Delta_1"]
    t = params["t"]
    plt.yticks([-t, -2*Delta_1, -Delta_1, 0, Delta_1, 2*Delta_1, t], ["-t", r"$-2\Delta_1$", r"$-\Delta_1$",
                                   0 ,r"$\Delta_1$", r"$2\Delta_1$", "t"])
    plt.ylim((-t, t))
    
if __name__ == '__main__':
    #main()
    main_2D()
    #Rashba_slider.on_changed(update)
    #mu_slider.on_changed(update)

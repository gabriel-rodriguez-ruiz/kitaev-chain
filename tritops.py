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

def make_low_energy_model(mu=0, t=10, u=-5j, Delta_0=0, Delta_1=2, Delta_11=0):
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
        return -mu * np.eye(4) + 1/2j * Delta_0 * np.kron(sigma_y, sigma_x)
    syst[lat(0)] = onsite
    def hopping(site1, site2, t, u, Delta_1, Delta_11):     #np.kron() is the tensor product
        return (-t * np.kron(np.eye(2), sigma_z) + 1j * u * np.kron(sigma_z, sigma_z)
                + 1/2 * (Delta_1 + Delta_11) * np.kron(1j * sigma_y, sigma_x))
    syst[kwant.HoppingKind((1,), lat)] = hopping
    return syst

# The function to be called anytime a slider's value changes
def update(value):
    global fig, tritops
    u = -1j*value       #value is Rashba parameter
    new_bands = kwant.physics.Bands(tritops, params=dict(mu=0, t=10, u=u, Delta_0=0, Delta_1=2, Delta_11=0))
    i = 0
    for line in fig.axes[0].lines:
        line.set_ydata([new_bands(k)[i] for k in line.get_xdata()])
        i += 1
    fig.canvas.draw()
    fig.canvas.flush_events()

def main():
    global tritops, Rashba_slider, fig
    tritops = make_low_energy_model()
    tritops = tritops.finalized()
    fig = kwant.plotter.bands(tritops, params=dict(mu=0, t=10, u=-5j, Delta_0=0, Delta_1=2, Delta_11=0))
    fig.canvas.manager.set_window_title("Bandas de energía de un triptops")
    # Make a horizontal slider to control the frequency.
    ax_Rashba = plt.axes([0.25, 0.1, 0.65, 0.03])
    plt.subplots_adjust(left=0.25, bottom=0.25)
    Rashba_slider = Slider(ax=ax_Rashba, label=r'$\lambda_R$', valmin=-10, valmax=10, valinit=5)

if __name__ == '__main__':
    main()
    Rashba_slider.on_changed(update)

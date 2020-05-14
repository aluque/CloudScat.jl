import sys
import scipy.constants as co

import numpy as np
from matplotlib import pyplot as plt
import h5py

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage[varg]{txfonts}')
plt.rc('axes', titlesize=54)
plt.rc('font', family='serif', size=12)

FNAME = "first_hitting.h5"
FOUT = "first_hitting.pdf"

def main():
    plt.figure(figsize=(5, 4))
    plot_mc()
    plot_model()
    
    plt.xlim([0, 10])
    plt.ylim([0, 500])
    plt.xlabel("time (ms)")
    plt.ylabel("normalized photon flux (s$^{-1}$)")
    plt.grid()
    plt.legend(loc="upper right")

    plt.savefig(FOUT)
    print(FOUT, "saved")
    
def plot_mc():
    fp = h5py.File(FNAME, "r")
    nobs = 1
    
    mu, w = np.polynomial.legendre.leggauss(5)
    tl = None
    t = None

    S = np.zeros_like(w)
    dt = 0.0
    
    for (i, wi) in enumerate(w):
        tl0 = np.array(fp[f"obs{i+1:05d}/timeline"])
        
        if tl is None:
            tl = wi * tl0
        else:
            tl += wi * tl0

        if t is None:
            t = np.array(fp[f"obs{i+1:05d}/t"])
            dt = t[1] - t[0]

        S[i] = sum(tl0) * dt
        
    delay = (400 * co.kilo) / co.c

    tshift = t - delay

    d = 400 * co.kilo
    tl *= np.pi * d**2

    dt = tshift[1] - tshift[0]
    
    plt.plot(tshift / co.milli, tl, c='k', lw=1.0, label="Monte Carlo")
    
    fp.close()
    
    
def plot_model():
    # Parameters from the Mie solver
    g = 0.8744118234466889
    omega0 = 0.9999971807871738
    Qext = 2.0399006732764797 

    Nd = 1e8
    R = 1e-5
    L = 5 * co.kilo
    
    nu = co.c * Nd * Qext * np.pi * R**2
    D = co.c**2 / (3 * nu * (1 - g * omega0))

    tauD = L**2 / (4 * D)
    tauA = 1 / (nu * (1 - omega0))

    d = (400 - 15) * co.kilo

    # Geometric factor
    S = (4 * np.pi * d**2)

    # Lambert factor
    lamb = 4
    
    print(f"omega0 = {omega0}")
    print(f"1 - omega0 = {1-omega0}")
    print(f"Qext = {Qext}")
    print(f"2 - Qext = {2-Qext}")
    
    print(f"D = {D}")
    print(f"tauD = {tauD}")
    print(f"tauA = {tauA}")
    
    t = np.linspace(1e-7, 10e-3, 500)
    
    F = (np.exp(-t / tauA - tauD / t) * (t / tauD)**(-1.5)
         / (np.sqrt(np.pi) * tauD))

    plt.plot(t / co.milli, F, c='r', label="Diffusion model")
    
    
if __name__ == '__main__':
    main()

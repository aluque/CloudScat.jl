import sys
import scipy.constants as co

import numpy as np
from matplotlib import pyplot as plt
import h5py

FNAME = "first_hitting_lowerbnd.h5"
FOUT = "first_hitting_lowerbnd.pdf"

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage[varg]{txfonts}')
plt.rc('axes', titlesize=54)
plt.rc('font', family='serif', size=12)

def main():
    plt.figure(figsize=(5, 4))

    plot_mc()
    plot_model()
    
    plt.xlim([0, 5])
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
        print(f"obs{i+1:05d}/timeline")
        
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

    plt.plot(tshift / co.milli, tl, c='k', lw=1.0, label="Monte Carlo")
    print("Total MC emissions:", sum(tl) * dt)
    
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
    
    w = 8 * co.kilo
    L = 5 * co.kilo

    print(f"tauD = {tauD}")
    print(f"tauA = {tauA}")
    
    t = np.linspace(1e-7, 10e-3, 500)
    
    F = (np.exp(-t / tauA - tauD / t) * (t / tauD)**(-1.5)
         / (np.sqrt(np.pi) * tauD))

    S = [np.zeros_like(t)]
    for k in range(1, 3):
        S.append(S[k - 1] +
                 ((2 * np.pi * D / w**2) * k * np.sin(np.pi * L * k / w)
                  * np.exp(-np.pi**2 * D * k**2 * t / w**2)))

    tauS = w**2 / (D * np.pi**2)
    tauT = 1 / (1 / tauS + 1 / tauA)

    plt.plot(t / co.milli, F, c='r', label="No lower boundary")
    plt.plot(t / co.milli, S[1], c='b', label="One-term")
    plt.plot(t / co.milli, S[2], c='b', ls="--", label="Two-term")
    
    
if __name__ == '__main__':
    main()

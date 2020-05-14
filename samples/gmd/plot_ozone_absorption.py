import sys
import os
import scipy.constants as co

import numpy as np
from matplotlib import pyplot as plt
import h5py
from scipy.interpolate import interp1d

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage[varg]{txfonts}')
plt.rc('axes', titlesize=54)
plt.rc('font', family='serif', size=12)

FOUT = "ozone_absorption.pdf"
PROFILES = [
    ("Tropical", "Hale", "Tropical", "#ff7700"),
    ("Tropical", "WarrenBrandt", "Tropical (ice)", "#ff7700"),
    ("MidLat_Summer", "Hale", "Midlatitude Summer", "r"),
    ("MidLat_Winter", "Hale", "Midlatitude Winter", "b"),
    ("US_Std_1976", "Hale", "US Standard Atmosphere", "g")]
PROFILES_PATH = "profiles"


def main():
    plt.figure(figsize=(8, 6))
    plt.subplots_adjust(left=0.075, right=0.95, wspace=0.35)
    
    plt.subplot(1, 2, 1)
    o3margins = plot_all_profiles()
    
    plt.subplot(1, 2, 2)
    plot_photo(o3margins)

    plt.savefig(FOUT)
    #plt.show()
    print(FOUT, "saved")


def plot_all_profiles():
    o3margins = []
    for prof, refindex, label, color in PROFILES:
        o3margin = plot_profile(os.path.join(PROFILES_PATH, f"{prof}.dat"),
                                label, color)

        o3margins.append(o3margin)
        
    plt.xlabel("O$_3$ concentration (10$^{18}$ m$^{-3}$)")
    plt.ylabel("altitude (km)")
    plt.ylim([0, 40])
    plt.xlim([0, 7])

    # The cloud location
    plt.fill_between([0, 10], 7, 12, alpha=0.2, facecolor='k')

    plt.axhline(10, color='k', ls='--', lw=1)
    
    plt.grid()
    plt.legend(loc="upper right")
    
    return o3margins


def plot_profile(fname, label, color):
    d = np.loadtxt(fname, skiprows=4)
    z = d[:, 1] * co.kilo
    T = d[:, 3]
    o3 = d[:, 5] * 1e-5 * co.atm / co.k / T

    o3margin = (interp1d(z, o3)(1.0e4),
                interp1d(z, o3)(1.5e4))
    if "ice" not in label:
        plt.plot(o3 / 1e18, z / co.kilo, label=label, color=color)

    return o3margin


def plot_photo(o3margins):
    for i, (prof, refindex, label, color) in enumerate(PROFILES):
        plot_photo_for_profile(f"ozone_absorption_{prof}_{refindex}.h5",
                               refindex, label, color,
                               o3margins[i])
        
    plt.xlim([0, 1.5])

    plt.xlabel("time (ms)")
    plt.ylabel("normalized photometer signal (m$^{-2}$s$^{-1}$)")
    plt.semilogy()
    plt.ylim([1e-15, 2e-10])
    
    plt.grid()
    plt.legend(loc="upper right")

    
def plot_photo_for_profile(fname, refindex, label, color, o3margin):

    fp = h5py.File(fname, "r")
    obs = 1
    
    tl = np.array(fp[f"obs{obs:05d}/timeline"])        
    t = np.array(fp[f"obs{obs:05d}/t"])
    obs = fp[f"obs{obs:05d}"]

    tshift = t - obs.attrs["delay"]
    dt = t[1] - t[0]
    
    sigma = 3.145E-19 * co.centi**2
    alpha = [sigma * no3 * co.c for no3 in o3margin]
    
    ls = '-' if refindex != "WarrenBrandt" else '--'
    plt.plot(tshift / co.milli, tl, ls=ls, lw=1.5, label=label, c=color)

    i0 = np.argmax(tl) + 3
    f = 1.
    if refindex != "WarrenBrandt":
        plt.fill_between(tshift / co.milli,
                         f * tl[i0] * np.exp(-alpha[1] * (tshift - tshift[i0])),
                         f * tl[i0] * np.exp(-alpha[0] * (tshift - tshift[i0])),
                         facecolor=color, edgecolor=color, lw=0.5, alpha=0.15)
    
    fp.close()
        
    
if __name__ == '__main__':
    main()

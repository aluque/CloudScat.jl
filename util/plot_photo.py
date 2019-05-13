import sys
import scipy.constants as co

import numpy as np
from matplotlib import pyplot as plt
from scipy import stats
from scipy.optimize import curve_fit
from scipy.integrate import trapz

import h5py

def get_parser():
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument("input",
                        help="HDF5 input file1")
    parser.add_argument("--observer",
                        type=int,
                        default=1,
                        help="Observer number (1-based)")
    return parser


def main():
    parser = get_parser()                        
    args = parser.parse_args()

    plotone(args.input, args.observer)
        
    
def plotone(fname, obs):
    fp = h5py.File(fname, "r")
    
    tl = np.array(fp[f"obs{obs:05d}/timeline"])
    t = np.array(fp[f"obs{obs:05d}/t"])
    g = fp[f"obs{obs:05d}"]

    decay = decay_time(t - g.attrs["delay"], tl)

    tshift = t - g.attrs["delay"]

    plt.plot(tshift / co.milli,
             tl, label="shift=%f" % g.attrs["shift"])

    sigma, l0, z = lognorm_fit(tshift, tl)

    # def lnorm(t, sigma, l0, A):
    #     return A * np.exp(-(np.log(tshift) - l0)**2 / (2 * sigma**2))

    # popt, pcov = curve_fit(lnorm, tshift, tl, p0=[sigma, l0, np.exp(z)])
    # sigma, l0, A = popt

    A = np.exp(z)
    t0 = np.exp(l0)
    itotal = A * np.sqrt(2 * np.pi) * np.exp(sigma**2 / 2) * t0 * sigma

    dt = t[1] - t[0]
    itotal0 = sum(tl) * dt

    print(f"    sigma = {sigma}")
    print(f"    t0 = {np.exp(l0) / co.milli} ms")
    print(f"    A  = {A}")
    print(f"    itotal  = {itotal} ({itotal0})")

    plt.plot(tshift / co.milli,
             A * np.exp(-(np.log(tshift) - l0)**2 / (2 * sigma**2)),
             c='k', alpha=0.5, lw=0.75)

    print("shift = %g km ==>  decay = %g ms" %
          (g.attrs["shift"] / co.kilo,
           decay / co.milli))
        
    plt.xlabel("Time (ms)")
    plt.ylabel("photons / m$^2$ / s / source photon")
    #plt.xlim([1e-2, 15])
    #plt.loglog()
    plt.show()

    
def decay_time(t, s):
    imax = np.argmax(s)

    t = t[imax:] - t[imax]
    s = s[imax:] / s[imax]

    flt = s > 0.1
    a, b, _, _, _ = stats.linregress(t[flt], np.log(s[flt]))

    return -1 / a


def lognorm_fit(t, s):
    imax = np.argmax(s)
    
    flt = np.logical_and(s > 0.05 * s[imax], t > 0)

    v = np.vander(np.log(t[flt]), 3)
    a, b, c = np.linalg.lstsq(v, np.log(s[flt]))[0]

    sigma = 1 / np.sqrt(-2 * a)
    l0 = -0.5 * b / a
    z = c - 0.25 * b**2 / a

    return sigma, l0, z
    
    
if __name__ == '__main__':
    main()

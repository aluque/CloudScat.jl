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

    parser.add_argument("--nofit",
                        action="store_true",
                        help="Remove the fits",
                        default=False)

    return parser


def main():
    parser = get_parser()                        
    args = parser.parse_args()

    plotone(args.input, args.observer, fits=not args.nofit)
        
    
def plotone(fname, obs, fits=True):
    fp = h5py.File(fname, "r")
    
    tl = np.array(fp[f"obs{obs:05d}/timeline"])
    t = np.array(fp[f"obs{obs:05d}/t"])
    obs = fp[f"obs{obs:05d}"]
    params = fp["parameters"]

    n, Qext, r = (params.attrs[k] for k in ("nscat", "Qext", "radius"))
    g, cloud_top, z, om0 = (params.attrs[k] for k in ("g", "cloud_top",
                                                      "source_altitude",
                                                      "ω₀"))
    
    L = cloud_top - z
    lmbd = 1 / (n * Qext * np.pi * r**2)
    D = lmbd * co.c / (3 * (1 - g * om0))
    tau = L**2 / (4 * D)
    alpha = tau * (1 - om0) * co.c / lmbd
    
    print("From parameters:")
    print(f"  g     = {g}")
    print(f"  tau   = {tau}")
    print(f"  alpha = {alpha}")
    
    tshift = t - obs.attrs["delay"]

    #norm, tstart, tauf, alphaf = exittime_fit(tshift, tl)

    norm, tstart, tauf, alphaf = exittime_fit2(tshift, tl)
    
    print("From fit:")
    print(f"  tau   = {tauf}")
    print(f"  alpha = {alphaf}")

    plt.plot(tshift / co.milli,
             tl, label="Model" % obs.attrs["shift"])

    if fits:
        plt.plot(tshift / co.milli,
                 norm * exittime(tshift, tstart, tauf, alphaf),
                 label="Fit" % obs.attrs["shift"], color="k")

        plt.plot(tshift / co.milli,
                 norm * exittime(tshift, 0.0, tau, alpha),
                 label="Analytical" % obs.attrs["shift"], color="r")

        
    plt.xlabel("Time (ms)")
    plt.ylabel("photons / m$^2$ / s / source photon")
    plt.legend()
    
    #plt.xlim([1e-2, 15])
    #plt.loglog()
    plt.show()


def exittime(t, tstart, tau, alpha):
    t1 = t - tstart
    w = (np.exp(2 * np.sqrt(alpha)) *
         np.sqrt(tau / t1**3) * np.exp(-tau / t1)
         * np.exp(-alpha * t1 / tau) / np.sqrt(np.pi))
    return np.where(w > 0, w, 0)

    
def exittime_fit(t, s):
    """ Fit a signal `s` sampled at times `t` to a exit-time times absorption
    curve. Returns:
    norm: the multiplying factor.
    t0: the starting time of the signal
    tau: the characteristic exit time
    alpha: the absroption rate in terms of tau.
    """
    # Assuming uniform sampling
    dt = t[1] - t[0]
    
    norm = sum(s)

    tavg = sum(t * s) / norm
    tdev2 = sum((t - tavg)**2 * s) / norm
    tskew = sum((t - tavg)**3 * s) / norm / tdev2**1.5
    
    alpha = (3 / (np.sqrt(2) * tskew))**4
    tau = np.sqrt(2 * alpha**1.5 * tdev2)

    t0 = tavg - tau / np.sqrt(alpha)
    
    return dt * norm, t0, tau, alpha
    

def exittime_fit2(t, s):
    dt = t[1] - t[0]
    norm = sum(s) * dt
    
    popt, pcov = curve_fit(exittime, t, s / norm, p0=[0.0, 5e-4, 0.05]) 

    return (norm, *popt)


if __name__ == '__main__':
    main()

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
    parser.add_argument("input", nargs='+',
                        help="HDF5 input file(s)")

    parser.add_argument("--observer",
                        type=int,
                        default=1,
                        help="Observer number (1-based)")

    parser.add_argument("--L",
                        action="store",
                        type=float,
                        help="Use this distance and plot the model",
                        default=None)
    
    parser.add_argument("--xlim", "-x",
                        help="Limits of the x-axis in ms (x0:x1)", 
                        action='store', default=None)

    parser.add_argument("--ylim", "-y",
                        help="Limits of the z-axis (y0:y1)", 
                        action='store', default=None)

    parser.add_argument("--output", "-o",
                        action="store",
                        help="Output file",
                        default=None)
    return parser


def main():
    parser = get_parser()                        
    args = parser.parse_args()

    for fname in args.input:
        plotone(fname, args.observer)

        if args.L is not None:
            plotmodel(fname, args.observer, args.L)

    
    xlim = [0, 5]
    if args.xlim is not None:
        xlim = [float(v) for v in args.xlim.split(':')]

    plt.xlim(xlim)
    
    if args.ylim is not None:
        ylim = [float(v) for v in args.ylim.split(':')]
        plt.ylim(ylim)
    
    plt.xlabel("Time (ms)")
    plt.ylabel("photons / m$^2$ / s / source photon")
    plt.legend()
        
    if args.output is not None:
        plt.savefig(args.output)
    else:
        plt.show()
        

def plotmodel(fname, obs, L):
    fp = h5py.File(fname, "r")
    
    tl = np.array(fp[f"obs{obs:05d}/timeline"])
    t = np.array(fp[f"obs{obs:05d}/t"])
    obs = fp[f"obs{obs:05d}"]
    params = fp["parameters"]

    g, n, Qext, r, om0 = (params.attrs[k] for k in ("g",
                                                    "nscat",
                                                    "Qext",
                                                    "radius",
                                                    "ω₀"))
    
    lmbd = 1 / (n * Qext * np.pi * r**2)
    D = lmbd * co.c / (3 * (1 - g * om0))
    tau = L**2 / (4 * D)
    alpha = tau * (1 - om0) * co.c / lmbd
    nu = (1 - om0) * co.c / lmbd
    
    tshift = t - obs.attrs["delay"]
    R = obs.attrs["delay"] * co.c

    # The normalization here is not 1 / 4 \pi R^2 because of Lambert's cosine
    # law.  
    norm = 1 / (np.pi * R**2)
    signal = exittime(tshift, 0.0, tau, alpha)

    dt = tshift[1] - tshift[0]
    print("Norm of the analytical expression:", dt * sum(signal))
    
    plt.plot(tshift / co.milli, norm * signal,
             label="Point source" % obs.attrs["shift"], color="r")

    fp.close()

    
def plotone(fname, obs):
    fp = h5py.File(fname, "r")
    
    tl = np.array(fp[f"obs{obs:05d}/timeline"])
    t = np.array(fp[f"obs{obs:05d}/t"])
    obs = fp[f"obs{obs:05d}"]
    params = fp["parameters"]

    tshift = t - obs.attrs["delay"]

    dt = tshift[1] - tshift[0]
    R = obs.attrs["delay"] * co.c

    print("Norm of the MC results:", dt * sum(tl) * (4 * np.pi * R**2))

    plt.plot(tshift / co.milli,
             tl, label=fname)

    fp.close()
 

def exittime(t, tstart, tau, alpha):
    t1 = t - tstart
    w = (np.exp(2 * np.sqrt(alpha)) *
         np.sqrt(tau / t1**3) * np.exp(-tau / t1)
         * np.exp(-alpha * t1 / tau) / np.sqrt(np.pi))
    return np.where(w > 0, w, 0)




if __name__ == '__main__':
    main()

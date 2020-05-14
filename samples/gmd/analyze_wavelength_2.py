import sys
import string

from itertools import product
import scipy.constants as co

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from scipy import stats

import h5py

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage[varg]{txfonts}')
plt.rc('axes', titlesize=54)
plt.rc('font', family='serif', size=14)

FOUT = "wavelength.pdf"

def main():
    obs = 1
    plt.figure(figsize=(7, 9))
    plt.subplots_adjust(top=0.98, bottom=0.1, hspace=0.075)
    lst = list(product([10, 12], [10, 20])) + [[0, 20]]
    
    for i, (h, R) in enumerate(lst):
        for lmbd in (337, 777):
            analyze(h, R, lmbd)
        
    
def analyze(h, R, lmbd):        
    print(f"h = {h}, R = {R}, lmbd = {lmbd}")

    if h != 0:
        fname = f"wavelength_2_{lmbd}nm_{h}km_{R}um.h5"
    else:
        fname = f"wavelength_2_extended_{lmbd}nm_{R}um.h5"

    fp = h5py.File(fname, "r")
    obs = 1
    
    # Note that the image is transposed wrt the julia array.
    img = np.array(fp[f"obs{obs:05d}/image"])

    width, height = img.shape
    x, y = np.arange(width), np.arange(height)
    v = img[:, height // 2]

    def m(x):
        return np.sum(v * x) / np.sum(v)

    
    width = np.sqrt(m(x**2) - m(x)**2)

    print(f"width = {width}")
    print()
    
    
def plot_map(ax, h, R, lmbd):
    if h != 0:
        fname = f"wavelength_2_{lmbd}nm_{h}km_{R}um.h5"
    else:
        fname = f"wavelength_2_extended_{lmbd}nm_{R}um.h5"

    fp = h5py.File(fname, "r")
    obs = 1
    
    # Note that the image is transposed wrt the julia array.
    img = np.array(fp[f"obs{obs:05d}/image"])

    width, height = img.shape
    ax.pcolormesh(img[492:532, 492:532], cmap="gnuplot2", rasterized=True)
    noxticks(ax)
    noyticks(ax)
    ax.tick_params('both', length=2, width=0.5, which='major')
    ax.axhline(512 - 492, lw=0.75, c="#777777")
    
    ax.text(0.03, 0.05, f"{lmbd} nm", color="w",
            transform=ax.transAxes)   

def move_down(ax):
    [left, bottom, width, height] = ax.get_position().bounds
    ax.set_position([left, bottom - 0.04, width, height])
    
def noxticks(ax):
    """ Remove xticks from the plot. """
    loc = ax.get_xticks()
    ax.set_xticklabels(['' for l in loc])

def noyticks(ax):
    """ Remove xticks from the plot. """
    loc = ax.get_yticks()
    ax.set_yticklabels(['' for l in loc])
    
    
if __name__ == '__main__':
    main()

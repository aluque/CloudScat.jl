import sys
import os
from warnings import warn
import re
import numpy as np
from matplotlib import pyplot as plt
import scipy.constants as co

from spacepy import pycdf

import h5py

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage[varg]{txfonts}')
plt.rc('axes', titlesize=54)
plt.rc('font', family='serif', size=14)

MMIA_FNAME = "mmia_triggered_level1_observation_time_2019-11-22_08-43-05_observation_id_27386_group_id_212080.cdf"
MC_FNAME = "asim_27386.h5"
FRAME = 1
FOUT = "asim_27386.pdf"

def main():
    plt.figure(figsize=(10, 5.))
    plt.subplots_adjust(left=0.05, right=0.95, top=0.95, bottom=0.12, wspace=0.2)

    plt.subplot(1, 2, 1)
    x, y = plot_mmia()

    plt.subplot(1, 2, 2)
    plot_mc(x, y)
    
    #plt.show()
    plt.savefig(FOUT)

    
def plot_mc(x, y):
    fp = h5py.File(MC_FNAME, "r")
    obs = 1
    
    # Note that the image is transposed wrt the julia array.
    img = np.array(fp[f"obs{obs:05d}/image"])

    width, height = img.shape

    plt.pcolormesh(x, y, img[y[0]:y[-1] + 1, x[0]:x[-1] + 1],
                   cmap="gnuplot2", rasterized=True)
    
    cbar = plt.colorbar()
    cbar.set_label("Normalized integrated radiance (m$^{-2}$sr$^{-1}$)")
    plt.xlabel("px")
    plt.ylabel("py")

    plt.axis('scaled')

    
def plot_mmia():
    cdf = pycdf.CDF(MMIA_FNAME)

    if (not cdf['CHU1Data_exists'][FRAME] or
        not cdf['CHU2Data_exists'][FRAME]):
        print("Error: data is not available", file=sys.stderr)
        sys.exit(-1)

    
    row0, row1 = [cdf[s][FRAME] for s in
                  ('chu_minimum_row', 'chu_maximum_row')]
    col0, col1 = [cdf[s][FRAME] for s in
                  ('chu_minimum_column', 'chu_maximum_column')]

    nrows = row1 - row0 + 1
    ncols = col1 - col0 + 1

    dark_pixels = 1024 + 32 - col1
    chus = np.array([cdf["CHU%d_photon_flux" % i][FRAME, :]\
                     .reshape((nrows, ncols))[:, :-dark_pixels]
                     for i in (1, 2)])
        
    
    # We substract 16 to remove the initial dark reference pixels
    x = np.r_[col0:col1 + 1 - dark_pixels] - 16
    y = np.r_[row0:row1 + 1]

    # Item #0 is the 337 data
    plt.pcolormesh(x, y, chus[0], cmap="gnuplot2", rasterized=True)
    cbar = plt.colorbar()
    cbar.set_label("Radiance (µW m$^{-2}$sr$^{-1}$)")
    plt.xlabel("px")
    plt.ylabel("py")
    plt.axis('scaled')

    return x, y


if __name__ == '__main__':
    main()

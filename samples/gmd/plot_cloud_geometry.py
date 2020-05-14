import sys
import scipy.constants as co

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from scipy import stats

import h5py

plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage[varg]{txfonts}')
plt.rc('axes', titlesize=54)
plt.rc('font', family='serif', size=16)

FNAME = "cloud_geometry.h5"
FOUT = "cloud_geometry.pdf"

def main():
    obs = 1
    
    fp = h5py.File(FNAME, "r")
    
    # Note that the image is transposed wrt the julia array.
    img = np.array(fp[f"obs{obs:05d}/image"])

    width, height = img.shape
    x, y = np.arange(width), np.arange(height)

    kwargs = {}

    xlim = [894, 1024]
    xfilter = np.logical_and(xlim[0] < x, x <= xlim[1])
    img = img[:, xfilter]
    x = x[xfilter]
        
    ylim = [894, 1024]
    yfilter = np.logical_and(ylim[0] < y, y <= ylim[1])
    img = img[yfilter, :]
    y = y[yfilter]

    clim = [1e-13, 1e-8]
    kwargs['vmin'], kwargs['vmax'] = clim        
    
    kwargs['norm'] = LogNorm()
    plt.subplots_adjust(bottom=0.125, top=0.925)
    
    plt.pcolormesh(x, y, img, **kwargs, cmap="gnuplot2", rasterized=True)
    cbar = plt.colorbar(extend='min')
    # cbar.set_label("normalized radiance (s$^{\mathdefault{-1}}$ sr$^{\mathdefault{-1}}$)")
    cbar.set_label("normalized signal (sr$^{-1}$)")

    for (x, y, lbl) in [(970, 965, "1"),
                        (940, 930, "2"),
                        (990, 1000, "3")]:
        plt.text(x, y, f"\Huge{{{lbl}}}", color="#ffffff")
        
    plt.axis('scaled')
    plt.xlabel("px")
    plt.ylabel("py")

    plt.savefig(FOUT)
    #plt.show()
    
    
    
if __name__ == '__main__':
    main()

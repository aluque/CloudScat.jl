import sys
import scipy.constants as co

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from scipy import stats

import h5py

def get_parser():
    import argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument("input",
                        help="HDF5 input file1")

    parser.add_argument("--observer",
                        type=int,
                        default=1,
                        help="Observer number")

    parser.add_argument("--xlim", "-x",
                        help="Limits of the x-axis (x0:x1)", 
                        action='store', default=None)

    parser.add_argument("--ylim", "-y",
                        help="Limits of the z-axis (y0:y1)", 
                        action='store', default=None)

    parser.add_argument("--log", action='store_true',
                        help="Logarithmic scale?")

    parser.add_argument("--clim", "-c",
                        help="Limits of the color axis (c0:c1)", 
                        action='store', default=None)

    parser.add_argument("--output", "-o",
                        action="store",
                        help="Output file",
                        default=None)

    return parser


def main():
    parser = get_parser()                        
    args = parser.parse_args()
    obs = args.observer
    
    fp = h5py.File(args.input, "r")
    
    # Note that the image is transposed wrt the julia array.
    img = np.array(fp[f"obs{obs:05d}/image"])

    plt.figure(f"obs{obs:05d}")

    width, height = img.shape
    x, y = np.arange(width), np.arange(height)
    
    kwargs = {}

    if args.xlim is not None:
        xlim = [float(v) for v in args.xlim.split(':')]
        xfilter = np.logical_and(xlim[0] < x, x < xlim[1])
        img = img[:, xfilter]
        x = x[xfilter]
        
    if args.ylim is not None:
        ylim = [float(v) for v in args.ylim.split(':')]
        yfilter = np.logical_and(ylim[0] < y, y < ylim[1])
        img = img[yfilter, :]
        y = y[yfilter]

    if args.clim is not None:
        clim = [float(v) for v in args.clim.split(':')]
        kwargs['vmin'], kwargs['vmax'] = clim

        
    if args.log:
        kwargs['norm'] = LogNorm()

    plt.pcolormesh(x, y, img, **kwargs)
    cbar = plt.colorbar()

    plt.title(args.input)

    plt.axis('scaled')
    plt.xlabel("px")
    plt.ylabel("py")

    if args.output is not None:
        plt.savefig(args.output)
    else:
        plt.show()

    
    
if __name__ == '__main__':
    main()

#! /usr/bin/env python

import numpy as np
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import healpy as hp
from matplotlib.ticker import MultipleLocator, StrMethodFormatter
import sys

if __name__=="__main__":
    if len(sys.argv) < 3:
        sys.exit("Provide the freq and fits filename")
        
    filename = sys.argv[2]
    freq = sys.argv[1]

    T = hp.read_map(filename)

#    hp.mollview(T, title='Temperature Map at {0} GHz'.format(freq), unit='K', coord='G',
#            norm='hist')
#    hp.graticule() # add meridians and parallels
#
#    plt.savefig('{0} GHz Temperature Map.pdf'.format(freq))

    # Compute the power spectrum
    LMAX = 2048
    cl = hp.anafast(T, lmax=LMAX)
    l = np.arange(len(cl))

    Dl = l * (l + 1)* cl / (2*np.pi)
    Dl *= 1e6 # Let's get the units right. Want it in MicroKelvin

    xlabelFormatter = StrMethodFormatter("{x:d}")
    xticks = [r"{0:d}".format(i) for i in [2, 10, 50, 100, 500, 1000, 1500, 2000]]
    yticks = [r"{0:d}".format(i) for i in [0, 1000, 2000, 3000, 4000, 5000, 6000]]
    ylabelFormatter = StrMethodFormatter("{x:3.0f}")
    fig, ax = plt.subplots(figsize=(12,12))
    ax.semilogx(l, Dl)
    ax.set_xlabel(r'l')
    ax.set_ylabel(r'$l (l+1) C_l / 2 \pi\ [\mu K^2]$')
    # ax.xaxis.set_major_formatter(xlabelFormatter)
    ax.set_xticklabels(xticks)
    ax.set_yticklabels(yticks)
    # ax.yaxis.set_major_formatter(ylabelFormatter)

    plt.savefig('TT Power Spectrum at {0} GHz.pdf'.format(freq))

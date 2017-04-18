#! /usr/bin/env python

import numpy as np
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import healpy as hp
import sys

if __name__=="__main__":
    filename = sys.argv[2]
    freq = sys.argv[1]

    T = hp.read_map(filename)

    hp.mollview(T, title='Temperature Map at {0} GHz'.format(freq), unit='K', coord='G',
            norm='hist')
    hp.graticule() # add meridians and parallels

    plt.savefig('{0} GHz Temperature Map.pdf'.format(freq))

    # Compute the power spectrum
    LMAX = 2048
    cl = hp.anafast(T, lmax=LMAX)
    l = np.arange(len(cl))

    Dl = l * (l + 1)* cl / (2*np.pi)

    fig, ax = plt.subplots(figsize=(12,12))
    ax.semilogx(l, Dl)
    ax.set_xlabel(r'l')
    ax.set_ylabel(r'$l (l+1) C_l / 2 \pi [\mu K^2]$')

    plt.savefig('TT Power Spectrum at {0} GHz.pdf'.format(freq))

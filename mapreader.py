#! /usr/bin/env python

import numpy as np
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import healpy as hp

filename = "HFI_SkyMap_100_2048_R2.02_full.fits"
T = hp.read_map(filename)

hp.mollview(T, title='Temperature Map at 100 GHz', unit='mK', coord='G',
        norm='hist')
hp.graticule() # add meridians and parallels

plt.savefig('100 GHz Temperature Map.pdf')

# Compute the power spectrum
LMAX = 2048
cl = hp.anafast(T, lmax=LMAX)
l = np.arange(LMAX)

cl = l * (l + 1)* cl / (2*np.pi)

fig, ax = plt.subplots(figsize=(12,12))
ax.plot(l, cl)
ax.set_xlabel(r'l')
ax.set_ylabel(r'$l (l+1) C_l / 2 \pi [\mu K^2]$')

plt.savefig('TT Power Spectrum.pdf')

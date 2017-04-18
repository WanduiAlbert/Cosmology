#! /usr/bin/env python

import numpy as np
import healpy as hp
import matplotlib
matplotlib.use('Agg')
from matplotlib import cm, colors

import matplotlib.pyplot as plt
import h5py

cmb_data = h5py.File("planck_100_longlat_0p5deg.mat", 'r')

print cmb_data.keys()

phi = np.array(cmb_data['phi']).reshape(-1)
theta = np.array(cmb_data['theta']).reshape(-1)
T = np.array(cmb_data['p100/T']).transpose()


# fig, ax = plt.subplots(figsize=(12,12), subplot_kw=dict(projection='mollweide'))
# ax.pcolor(phi, theta, T, cmap=cm.jet)
# ax.xaxis.set_visible(False)
# plt.savefig('Temperature_100GHz.pdf')

phig, thetag = np.meshgrid(phi, theta)
T[np.abs(thetag) < 0.3] = 0

fig, ax = plt.subplots(figsize=(12,12), subplot_kw=dict(projection='mollweide'))
ax.pcolormesh(phig, thetag, T, cmap=cm.jet)
# ax.xaxis.set_visible(False)
plt.savefig('Temperature_100GHz_masked.pdf')

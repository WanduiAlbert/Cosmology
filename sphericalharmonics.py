#! /usr/bin/env python


import numpy as np
import matplotlib
matplotlib.use('TkAgg')
from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, StrMethodFormatter
import scipy.special as sp
from scipy.misc import factorial

rad2deg = 180/np.pi

# Number of phi values and the angular increment
nphi = 90
da = 2*np.pi/nphi
# Calculate theta and phi values
phi = np.linspace(0, 2*np.pi, 100)
theta = np.linspace(0, np.pi, 100)
phi, theta = np.meshgrid(phi, theta)

def sph_harmonic(l, m, theta, phi):
    P_lm = sp.lpmv(m,l,np.cos(theta))
    prefactor = np.sqrt((1.0 * (2*l + 1) * factorial(l-m))/(4*np.pi * factorial(l+m)))

    Y_lm = prefactor * P_lm * np.exp(1j * m * phi)

    return Y_lm

Y_lm = sph_harmonic(2, 0, theta, phi).real
Y_sup = np.sqrt(3./8)*sph_harmonic(2,2,theta,phi) + np.sqrt(3./8)*sph_harmonic(2,-2,theta,phi) -\
0.5*sph_harmonic(2,0,theta,phi)


##################################################################################################
r = np.abs(Y_sup)
# The Cartesian coordinates of the spheres
x = r * np.sin(theta) * np.cos(phi)
y = r * np.sin(theta) * np.sin(phi)
z = r * np.cos(theta)

N = r/r.max() # Color scaling for mapping on to the harmonics
fig, ax = plt.subplots(subplot_kw=dict(projection='3d'), figsize=(12,10))
im = ax.plot_surface(x, y, z, rstride=1, cstride=1, facecolors=cm.jet(N))
m = cm.ScalarMappable(cmap=cm.jet)
m.set_array(r)    # Assign the unnormalized data array to the mappable
                  #so that the scale corresponds to the values of R
fig.colorbar(m, shrink=0.8);
ax.set_axis_off()

plt.show()
##################################################################################################
# Calculate the spherical harmonic Y(l,m) and normalize to [0,1]
fcolors = Y_sup.real
fmax, fmin = fcolors.max(), fcolors.min()
fcolors = (fcolors - fmin)/(fmax - fmin)

x = np.sin(theta) * np.cos(phi)
y = np.sin(theta) * np.sin(phi)
z = np.cos(theta)

# # Set the aspect ratio to 1 so our sphere looks spherical
# fig2 = plt.figure(figsize=plt.figaspect(1.))
fig2 = plt.figure(figsize=(12,12))
ax2 = fig2.add_subplot(111, projection='3d')
ax2.plot_surface(x, y, z,  rstride=1, cstride=1, facecolors=cm.seismic(fcolors))
m2 = cm.ScalarMappable(cmap=cm.seismic)
m2.set_array(fcolors)    # Assign the unnormalized data array to the mappable
                  #so that the scale corresponds to the values of R
fig2.colorbar(m2, shrink=0.8);
# # Turn off the axis planes
ax2.set_axis_off()

plt.show()

# ylabels = map(str, range(-90, 90, 15))
# yFormatter = StrMethodFormatter(r"{x:2.1f}")
# fig, ax = plt.subplots(figsize=(12,12), subplot_kw=dict(projection='mollweide'))
# ax.pcolor(phi - np.pi, np.pi/2 - theta, np.abs(Y_sup)**2);
# ax.xaxis.set_visible(False)
# ax.yaxis.set_ticklabels(ylabels);
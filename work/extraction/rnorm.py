#!/usr/bin/env python

"""
Make plots about R-matrix normalization
"""

import os, sys, glob
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table

from desispec.io import findfile, read_frame

os.environ['SPECPROD']='fuji'
night, expid, camera = 20210101, 70350, 'z1'
frame = read_frame(findfile('frame', night, expid, camera))

#- Convolve with an [OII]-like double Gaussian
wave = frame.wave
rnorm = np.sum(frame.R[0].data, axis=0)

#- Lorentzian kernal
n = 31
xx = np.arange(n)
gg = np.zeros(n)
sigma = 2.0/0.8
x0 = n/2 - 2.5/0.8
x1 = n/2 + 2.5/0.8
gg += np.exp(-(xx-x0)**2 / (2*sigma**2))
gg += np.exp(-(xx-x1)**2 / (2*sigma**2))
gg /= np.sum(gg)

gg_rnorm = np.convolve(rnorm, gg, mode='same')

#- Convolve with an [OIII]-like Lorentzian

xx = np.arange(31)
fwhm = 2.7/0.8    #- typical [OIII] FWHM in units of wavelength bins
x0 = len(xx)/2
lz = (fwhm/2)**2 / ((fwhm/2)**2 + (xx-x0)**2)
lz /= np.sum(lz)

lz_rnorm = np.convolve(rnorm, lz, mode='same')

ii = (9200 <= frame.wave) & (frame.wave < 9400)

fig, (ax1, ax2) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [1,3]},figsize=(5,3))

ax1.plot(frame.wave[ii], frame.flux[0,ii], 'k-')
ax1.set_ylim(0, 1.1*np.max(frame.flux[0,ii]))
ax1.set_ylabel(r'counts / $\AA$')
ax1.axes.xaxis.set_ticklabels([])

ax2.plot(frame.wave[ii], rnorm[ii], label='$\delta$ function')
ax2.plot(frame.wave[ii], lz_rnorm[ii], label='[OIII]-like Lorentzian')
ax2.plot(frame.wave[ii], gg_rnorm[ii], label='[OII]-like Double Gaussian')
ax2.legend(loc='upper left',fontsize="small")

ax2.set_ylim(0.8, 1.2)
ax2.set_ylabel('R normalization bias')
ax2.set_xlabel(r'Wavelength [$\AA$]')
ax2.grid(axis='y')

fig.tight_layout()
outfile = 'rnorm-bias.pdf'
plt.savefig(outfile)
print(f'Wrote {outfile}')

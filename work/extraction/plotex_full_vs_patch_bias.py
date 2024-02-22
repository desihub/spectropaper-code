#!/usr/bin/env python
# coding: utf-8

#- Plot extraction bias from divide-and-conquer vs. full extractions
#- Writes extraction_full_vs_patch_bias.pdf

#- non-DESI imports
import os, sys, glob, time

import matplotlib.pyplot as plt
import numpy as np
import fitsio

from astropy.visualization import ZScaleInterval

plt.ion();

#- DESI imports
import desispec.io
from desispec.io import findfile, read_image, read_frame
from desispec.interpolation import resample_flux
from desispec.resolution import Resolution
from specter.psf import load_psf
from specter.extract import ex2d

#- Pick an example image
specprod = 'fuji'
os.environ['SPECPROD'] = specprod
reduxdir = os.path.expandvars(f'$DESI_ROOT/spectro/redux/{specprod}/preproc/')
night = 20210101
expid = 70350
camera = 'r1'

dw = 0.8
ww = 7300 + np.arange(20)*dw
specrange = (25,30)

print(f'Using {night=} {expid=} {camera=} {specrange=} wavelengths {ww[0]} to {ww[0]}')

#- Read inputs
print('Reading inputs')
img = read_image(findfile('preproc', night, expid, camera))
psf = load_psf(findfile('psf', night, expid, camera))

# Compare results from full extraction with those from divide-and-conquer

dw = 0.8
nwave = 200
wave = 7300 + np.arange(nwave)*dw
nspec = 5

print('Extracting spectra')
flux, ivar, rdat = ex2d(img.pix, img.ivar*(img.mask==0), psf, 25, nspec, wave, wavesize=nwave)
# plt.plot(wave, flux[0], 'k-')

subwave = list()
subflux = list()
subivar = list()
npatch = 4
for ww in np.array_split(wave, npatch):
    f1, iv1, rdat1 = ex2d(img.pix, img.ivar*(img.mask==0), psf, 25, nspec, ww)
    subwave.append(ww)
    subflux.append(f1)
    subivar.append(iv1)
    # plt.plot(ww, f1[0], '-')


print('Making plots')
#plt.figure(figsize=(6,3))
plt.figure(figsize=(5,2.5))

#- Flux fiber 25
ax = plt.subplot(221)
for w1, f1, iv1 in zip(subwave, subflux, subivar):
    plt.plot(w1, f1[0], lw=1)

plt.plot(wave, flux[0], 'k:', lw=1)
plt.ylim(200, 1200)
plt.ylabel(r'flux [e/$\mathrm{\AA}$]')
plt.title('fiber 25')
ax.axes.xaxis.set_ticklabels([])

#- Flux fiber 29
ax = plt.subplot(222)
for w1, f1, iv1 in zip(subwave, subflux, subivar):
    plt.plot(w1, f1[4], lw=1)
plt.ylim(200, 1200)

plt.plot(wave, flux[4], 'k:', lw=1)
plt.title('fiber 29')
ax.axes.xaxis.set_ticklabels([])
ax.axes.yaxis.set_ticklabels([])

#- Flux differences
ax = plt.subplot(223)
for w1, f1, iv1 in zip(subwave, subflux, subivar):
    ii = np.isin(wave, w1)
    plt.plot(w1, (flux[0][ii] - f1[0]) * np.sqrt(iv1[0]), lw=1)

plt.ylim(-0.05, 0.05)
plt.xlabel(r'Wavelength [$\mathrm{\AA}$]')
plt.ylabel('$\Delta$ flux / $\sigma$')

ax = plt.subplot(224)
for w1, f1, iv1 in zip(subwave, subflux, subivar):
    ii = np.isin(wave, w1)
    plt.plot(w1, (flux[4][ii] - f1[4]) * np.sqrt(iv1[4]), lw=1)

plt.ylim(-0.05, 0.05)
plt.xlabel(r'Wavelength [$\mathrm{\AA}$]')
ax.axes.yaxis.set_ticklabels([])

plt.tight_layout()
plt.savefig('extraction_full_vs_patch_bias.pdf')

print('Wrote extraction_full_vs_patch_bias.pdf')

#!/usr/bin/env python
# coding: utf-8

# Compare the choice of extraction resolution to the pixel scale and instrument resolution
# Writes extraction_wavelengths.pdf

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

print(f'Using {night=} {expid=}')

print('Reading inputs')
camera_psf = dict()
wavelengths = dict()
for camera in ['b1', 'r1', 'z1']:
    camera_psf[camera] = load_psf(findfile('psf', night, expid, camera))
    wavelengths[camera] = fitsio.read(findfile('frame', night, expid, camera), 'WAVELENGTH')


i=100
colors = dict(b1='C0', r1='C3', z1='C7')
plt.figure(figsize=(6,4))
for camera, psf in camera_psf.items():
    wmin = wavelengths[camera][0]
    wmax = wavelengths[camera][-1]
    print(f'{camera} {wmin:.1f} {wmax:.1f}')
    ww = psf.wavelength(i, np.arange(0, psf.npix_y))
    plt.plot(ww, np.gradient(ww), color=colors[camera], lw=3)
    plt.plot(ww, psf.ysigma(i, ww), '--', color=colors[camera], lw=3)

plt.text(3800, 1.2, 'PSF Gaussian sigma')

plt.plot([3600, 9824], [0.8, 0.8], 'k-', lw=1)
plt.text(3800, 0.85, 'Extraction wavelength step')

plt.ylabel(r'Pixel & PSF scale [$\mathrm{\AA}$]')
plt.xlabel(r'Wavelength [$\mathrm{\AA}$]')
plt.ylim(0.0, 1.4)
plt.text(3800, 0.5, 'Pixel scale')
plt.savefig('extraction_wavelengths.pdf')

print('Wrote extraction_wavelengths.pdf')


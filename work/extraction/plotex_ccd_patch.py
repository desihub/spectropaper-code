#!/usr/bin/env python
# coding: utf-8

# Make plot of CCD extraction patch
# Writes extraction_ccd_patch.pdf

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

pixpad_frac=0.8
wavepad_frac=0.2
ny, nx = spotsize = psf.pix(specrange[0], ww[0]).shape
extra_ypix = int(round(pixpad_frac * ny))

x1, x2, y1, y2 = psf.xyrange(specrange, ww)
y1a = y1 - extra_ypix
y2a = y2 + extra_ypix

x2a = psf.xyrange(specrange[1], ww)[1]

wmin = psf.wavelength(specrange[0], y1a-int(ny*wavepad_frac))
wmax = psf.wavelength(specrange[0], y2a+int(ny*wavepad_frac))

waveall = 7300-50*dw + np.arange(200)*dw
ii = (wmin <= waveall) & (waveall <= wmax)
waveall = waveall[ii]
jj = (waveall < ww[0]) | (ww[-1] < waveall)
waveextra = waveall[jj]

x0, x3 = x1-15, x2+15
y0, y3 = y1a-10, y2a+10
patch0 = img.pix[y1:y2, x1:x2] * (img.ivar[y1:y2, x1:x2] > 0) * (img.mask[y1:y2, x1:x2]==0)
patch1 = img.pix[y0:y3, x0:x3] * (img.ivar[y0:y3, x0:x3] > 0) * (img.mask[y0:y3, x0:x3]==0)
vmin, vmax = ZScaleInterval().get_limits(patch1)

print('Making plots')
plt.figure(figsize=(4,4))
plt.imshow(patch1, extent=(x0-0.5,x3-0.5,y0-0.5,y3-0.5), vmin=vmin, vmax=vmax, cmap='Greys')
# plt.imshow(patch0, extent=(x1-0.5,x2-0.5,y1-0.5,y2-0.5), vmin=vmin, vmax=vmax, cmap='Blues')

plt.plot([x1,x2a,x2a,x1,x1], [y1a,y1a,y2a,y2a,y1a], 'r-')
plt.plot([x1,x2,x2,x1,x1], [y1,y1,y2,y2,y1], 'b-')

for i in range(*specrange):
    x, y = psf.xy(i, waveextra)
    plt.plot(x, y, 'rx')
    x, y = psf.xy(i, ww)
    plt.plot(x, y, 'b.')

x, y = psf.xy(specrange[1], waveall)
plt.plot(x, y, 'rx')

plt.xlim(x0,x3); plt.ylim(y0,y3)
plt.xlabel('CCD x pixel'); plt.ylabel('CCD y pixel')
plt.savefig('extraction_ccd_patch.pdf')

print('Wrote extraction_ccd_patch.pdf')


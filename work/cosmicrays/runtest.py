#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import fitsio
from desispec.io import read_image
from desispec.maskbits import ccdmask
from hacked_cosmics import reject_cosmic_rays_ala_sdss
print(ccdmask)

def tune_axes(a) :
    a.axes.xaxis.set_ticklabels([])
    a.axes.yaxis.set_ticklabels([])
    a.tick_params(bottom = False)
    a.tick_params(left = False)

img = read_image("preproc-r0-00124034.fits")
img.mask[img.mask&ccdmask.COSMIC>0] -= ccdmask.COSMIC  # remove COSMIC mask

if 0 :
    rejected_0 = reject_cosmic_rays_ala_sdss(img,nsig=6.,cfudge=3.,c2fudge=0.5,niter=100,step=0)
    rejected_1 = reject_cosmic_rays_ala_sdss(img,nsig=6.,cfudge=3.,c2fudge=0.5,niter=100,step=1)
    rejected_2 = reject_cosmic_rays_ala_sdss(img,nsig=6.,cfudge=3.,c2fudge=0.5,niter=100,step=2)
    rejected_3 = reject_cosmic_rays_ala_sdss(img,nsig=6.,cfudge=3.,c2fudge=0.5,niter=100,step=3)
    fitsio.write("rej0.fits",rejected_0.astype(int),clobber=True)
    fitsio.write("rej1.fits",rejected_1.astype(int),clobber=True)
    fitsio.write("rej2.fits",rejected_2.astype(int),clobber=True)
    fitsio.write("rej3.fits",rejected_3.astype(int),clobber=True)
else :
    rejected_0 = fitsio.read("rej0.fits")
    rejected_1 = fitsio.read("rej1.fits")
    rejected_2 = fitsio.read("rej2.fits")
    rejected_3 = fitsio.read("rej3.fits")

b0=2150
e0=b0+100
b1=1390
e1=b1+100

b0=2327
e0=b0+200
#b1=1850
b1=916
e1=b1+200

b0=2060
e0=b0+200
#b1=1850
b1=1470
e1=b1+200

# having a hard time to find examples with step 4 helps ...
b0=3190
e0=b0+100
b1=420
e1=b1+100

plt.tick_params(bottom = False)
plt.tick_params(left = False)

cmap="gray_r"
plt.figure("cosmics-example",figsize=(8,3.5))
ny=2;nx=5
p=1
textsize=18
a=plt.subplot(ny,nx,p); p+=1
plt.imshow(img.pix[b0:e0,b1:e1],vmin=-10,vmax=100,origin="lower",cmap=cmap)
tune_axes(a)
a=plt.subplot(ny,nx,p); p+=1
a.set_title("step 1", size=textsize)
plt.imshow(img.pix[b0:e0,b1:e1]*(rejected_0[b0:e0,b1:e1]==0),vmin=-10,vmax=100,origin="lower",cmap=cmap)
tune_axes(a)
a=plt.subplot(ny,nx,p); p+=1
a.set_title("step 2", size=textsize)
plt.imshow(img.pix[b0:e0,b1:e1]*(rejected_1[b0:e0,b1:e1]==0),vmin=-10,vmax=100,origin="lower",cmap=cmap)
tune_axes(a)
a=plt.subplot(ny,nx,p); p+=1
a.set_title("step 3", size=textsize)
plt.imshow(img.pix[b0:e0,b1:e1]*(rejected_2[b0:e0,b1:e1]==0),vmin=-10,vmax=100,origin="lower",cmap=cmap)
tune_axes(a)
a=plt.subplot(ny,nx,p); p+=1
a.set_title("step 4", size=textsize)
plt.imshow(img.pix[b0:e0,b1:e1]*(rejected_3[b0:e0,b1:e1]==0),vmin=-10,vmax=100,origin="lower",cmap=cmap)
tune_axes(a)

p=7
a=plt.subplot(ny,nx,p); p+=1
plt.imshow(rejected_0[b0:e0,b1:e1],origin='lower',cmap=cmap)
tune_axes(a)
a=plt.subplot(ny,nx,p); p+=1
plt.imshow(rejected_1[b0:e0,b1:e1],origin='lower',cmap=cmap)
tune_axes(a)
a=plt.subplot(ny,nx,p); p+=1
plt.imshow(rejected_2[b0:e0,b1:e1],origin='lower',cmap=cmap)
tune_axes(a)
a=plt.subplot(ny,nx,p); p+=1
plt.imshow(rejected_3[b0:e0,b1:e1],origin='lower',cmap=cmap)
tune_axes(a)
plt.tight_layout()
plt.show()

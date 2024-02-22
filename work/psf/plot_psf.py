#!/usr/bin/env python

import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
import specter.psf
import sys
import argparse
import string
import os.path
from   scipy.signal import fftconvolve
from   desispec.io import read_image

import matplotlib.pylab as pylab

fontsize=12
params = {'axes.labelsize': fontsize,
          'axes.titlesize': fontsize,
          'xtick.labelsize': fontsize,
          'ytick.labelsize': fontsize}
pylab.rcParams.update(params)

def readpsf(filename) :
    try :
        psftype=pyfits.open(filename)[0].header["PSFTYPE"]
    except KeyError :
        psftype=""
    print("PSF Type=",psftype)
    if psftype=="GAUSS-HERMITE" :
        return specter.psf.GaussHermitePSF(filename)
    elif psftype=="SPOTGRID" :
        return specter.psf.SpotGridPSF(filename)

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--psf', type = str, default = None, required = True,
                    help = 'path of psf file')
parser.add_argument('--preproc', type = str, default = None, required = True,
                    help = 'path of preproc arc lamp file')
parser.add_argument('--fiber', type = int, default = None, required = True,
                    help = 'fiber for psf1')
parser.add_argument('--wavelength', type = float, default = 6000., required = False,
                    help = 'wavelength')
parser.add_argument('-o','--output', type = str, default = None, required = False,
                    help = 'path to output image (png) file')
#parser.add_argument('--no-pixel-convolution', action = "store_true",
#                    help = 'do not convolve PSFs with pixel size')

args        = parser.parse_args()


psf=readpsf(args.psf)
xy=psf.xy(args.fiber,args.wavelength)

preproc=read_image(args.preproc)

print("for psf, xy=",xy)


xx,yy,xypix=psf.xypix(args.fiber,args.wavelength)

oversample=10
x1d=np.linspace(xx.start-0.5,xx.stop-0.5,(xx.stop-xx.start)*oversample)
y1d=np.linspace(yy.start-0.5,yy.stop-0.5,(yy.stop-yy.start)*oversample)
x=np.tile(x1d,(y1d.size,1))
y=np.tile(y1d,(x1d.size,1)).T
#y=x.T
fpix=psf._value(x,y,args.fiber,args.wavelength)
fpix /= np.sum(fpix)

fig=plt.figure("psf-r1",figsize=(8,6))


mx=5
my=5

#a0=plt.subplot(111)
plt.text(0.5,0.8,"fiber #%d lambda=%dA"%(args.fiber,args.wavelength),fontsize=10,color="k",horizontalalignment='center')#,transform=a0.transAxes)


a=plt.subplot(231)
a.set_title("PSF model",size=fontsize)
a.imshow(fpix,origin='lower',interpolation="nearest")
for i in range(1,xx.stop-xx.start) :
    a.axvline(i*oversample,color="white",alpha=0.3)
for i in range(1,yy.stop-yy.start) :
    a.axhline(i*oversample,color="white",alpha=0.3)
a.axis("off")



#plt.text(-hw+0.3,-hw+0.8,"fiber #%d lambda=%dA"%(args.fiber,args.wavelength),fontsize=10,color="white")
#plt.text(-hw+0.3,-hw+0.1,"(x,y)=(%4.1f,%4.1f)"%(xy[0],xy[1]),fontsize=10,color="white")
#plt.xlabel("x ccd")
#plt.ylabel("y ccd")

psfcore=(xypix>np.max(xypix)/3.)*xypix
#psfcore/=np.sum(psfcore)
flux=np.sum(preproc.pix[yy,xx]*psfcore)/np.sum(psfcore**2)
vmax=np.max(preproc.pix[yy,xx])
a=plt.subplot(232)
a.set_title("integrated PSF",size=fontsize)
a.imshow(flux*xypix,origin='lower',interpolation="nearest",vmin=0,vmax=vmax)
a.axis("off")


a=plt.subplot(233)
a.set_title("data",size=fontsize)
a.imshow(preproc.pix[yy,xx],origin='lower',interpolation="nearest",vmin=0,vmax=vmax)
a.axis("off")


a=plt.subplot(212)
dx=x1d-xy[0]
x=x1d
y=np.repeat(xy[1],x.shape)
fpix=psf._value(x,y,args.fiber,args.wavelength)
norm=np.max(fpix)
plt.plot(dx,fpix/norm)
if args.fiber>0 :
    fpix2=psf._value(x,y,args.fiber-1,args.wavelength)
    norm=np.max(fpix2)
    plt.plot(dx,fpix2/norm,":")
if args.fiber<psf.nspec-1 :
    fpix2=psf._value(x,y,args.fiber+1,args.wavelength)
    norm=np.max(fpix2)
    plt.plot(dx,fpix2/norm,":")

plt.xlabel("pixels",fontsize=fontsize)
#plt.ylabel("normalized PSF profile along a CCD row")
#plt.yscale("log")
plt.grid()
plt.tight_layout()
plt.show()

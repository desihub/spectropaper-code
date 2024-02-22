#!/usr/bin/env python

import numpy as np
import astropy.io.fits as pyfits
import pylab
import specter.psf
import sys
import argparse
import string
from   scipy.signal import fftconvolve
import matplotlib.pyplot as plt
import argparse

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


parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
description='''Measure the resolution of each fiber as a function of wavelength and save it in an ASCII file.
''')
parser.add_argument('-i','--input', type = str, default = None, required = True, help = 'psf fits file')
parser.add_argument('-o','--output', type = str, default = None, required = False, help = 'output ASCII file')
parser.add_argument('--plot', action='store_true')


args = parser.parse_args()

hw=5.
n1d=51
x1d=np.linspace(-hw,hw,51)
x=np.tile(x1d,(n1d,1))
y=x.T
kernel=(np.abs(x1d)<0.5).astype(float)
kernel/=np.sum(kernel)

filename=args.input
psf=readpsf(filename)
wave=np.linspace(psf.wmin,psf.wmax,15)

fibers=np.arange(psf.nspec)
if len(fibers)>41 :
    fibers=np.linspace(0,len(fibers)-1,30).astype(int)

nfibers=len(fibers)
FWHM=np.zeros((nfibers,wave.size))

for ifib,fiber in enumerate(fibers) :
    for u,w in enumerate(wave) :
        xy=psf.xy(fiber,w)

        #fpix=psf._value(x+xy[0],y+xy[1],fiber,w+1.) # just to check orientation
        fpix=psf._value(x+xy[0],y+xy[1],fiber,w)

        if False :
            for i in range(n1d) :
                fpix[i]=fftconvolve(fpix[i],kernel, mode='same')
            for j in range(n1d) :
                fpix[:,j]=fftconvolve(fpix[:,j],kernel, mode='same')

        proj=np.sum(fpix,axis=1)
        proj/=np.max(proj)
        #plt.plot(x1d,proj)
        #plt.show()
        x0=np.interp(0.5,proj[:n1d//2-1],x1d[:n1d//2-1])
        x1=np.interp(0.5,proj[n1d//2+1:][::-1],x1d[n1d//2+1:][::-1])
        FWHM[ifib,u]=x1-x0

if args.output :
    tmp=np.zeros((wave.size,nfibers+1))
    tmp[:,0]=wave
    for i in range(len(fibers)):
        tmp[:,i+1]=FWHM[i]
    np.savetxt(args.output,tmp)
    print("wrote",args.output)

if args.plot :
    for i in range(nfibers):
        plt.plot(wave,FWHM[i],alpha=0.5)
    plt.grid()
    plt.show()

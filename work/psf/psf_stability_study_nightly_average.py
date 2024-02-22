#!/usr/bin/env python

import numpy as np
import fitsio
import astropy.io.fits as pyfits
from  astropy.table import Table
import matplotlib.pyplot as plt
import specter.psf
import sys
import argparse
import os.path
from desispec.util import parse_fibers
from desiutil.log  import get_logger

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
parser.add_argument('--psf', type = str, nargs = "*", default = None, required = True,
                    help = 'path of psf files')
parser.add_argument('--fibers', type = str, default = None, required = False,
                    help = 'defines from_to which fiber to work on. (ex: --fibers=50:60,4 means that only fibers 4, and fibers from 50 to 60 (excluded) will be plotted)')
parser.add_argument('-o','--output', type = str, default = None, required = False,
                    help = 'path to output ascii file')


args        = parser.parse_args()
log = get_logger()
psfs=[]
nights=[]
for filename in args.psf :
    log.info("reading %s"%filename)
    psfs.append(readpsf(filename))

    # the following few lines are pretty bad ... because no info on exposure is saved in PSF file for andes run
    vals=os.path.basename(filename).replace(".fits","").split("-")
    night=int(vals[-1])
    cam=vals[-2]
    nights.append(night)

wmin=psfs[0]._wmin_all
wmax=psfs[0]._wmax_all
nw=20
waves=np.linspace(wmin+(wmax-wmin)/nw/2.,wmax-(wmax-wmin)/nw/2.,nw)
fibers=parse_fibers(args.fibers)
if fibers is None :
        fibers = np.arange(psfs[0].nspec)

cumul_night = []
cumul_fiber = []
cumul_wave = []
cumul_xc = []
cumul_yc = []
cumul_xsig = []
cumul_ysig = []
cumul_emission_line_coeff = []
cumul_continuum_coeff     = []


for fiber in fibers :
    for wave in waves :
        images = []
        i0 = []
        i1 = []

        for psf,night in zip(psfs,nights) :


            xx, yy, ccdpix = psf.xypix(fiber,wave)
            this_xc=psf.x(fiber,wave)
            this_yc=psf.y(fiber,wave)


            xprof = np.mean(ccdpix,axis=0)
            dx  = np.arange(xx.start,xx.stop)-this_xc
            yprof = np.mean(ccdpix,axis=1)
            dy  = np.arange(yy.start,yy.stop)-this_yc
            this_xsig = np.sqrt(np.sum(xprof*dx**2)/np.sum(xprof))
            this_ysig = np.sqrt(np.sum(yprof*dy**2)/np.sum(yprof))

            cumul_night.append(night)
            cumul_fiber.append(fiber)
            cumul_wave.append(wave)
            cumul_xc.append(this_xc)
            cumul_yc.append(this_yc)
            cumul_xsig.append(this_xsig)
            cumul_ysig.append(this_ysig)

            images.append(ccdpix)
            i1.append(xx.start)
            i0.append(yy.start)

        mi0 = int(np.min(i0))
        mi1 = int(np.min(i1))
        n=len(images)
        # add a margin and apply offset if necessary
        nimages=np.zeros((n,images[0].shape[0]+3,images[0].shape[1]+3))
        for j in range(n) :
            nimages[j,i0[j]-mi0:i0[j]-mi0+images[j].shape[0],i1[j]-mi1:i1[j]-mi1+images[j].shape[1]] = images[j]

        images=nimages
        n=images.shape[0]
        mimage=np.mean(images,axis=0)
        x=np.tile(np.arange(mimage.shape[0]),(mimage.shape[1],1))    # check with visual inspection of ccd image
        y=np.tile(np.arange(mimage.shape[1]),(mimage.shape[0],1)).T  # y is wavelength axis

        for j in range(n) :
            this_delta_ratio_emission_line = np.sum(images[j]*mimage)/np.sum(images[j]**2)-1
            pmimage=np.sum(mimage,axis=0) # projection to get 1D PSF along cross-dispersion for continuum fit normalization
            pimage=np.sum(images[j],axis=0)
            this_delta_ratio_continuum = np.sum(pimage*pmimage)/np.sum(pimage**2)-1

            cumul_emission_line_coeff.append(this_delta_ratio_emission_line)
            cumul_continuum_coeff.append(this_delta_ratio_continuum)


if args.output :
    table=Table()
    table["night"]=cumul_night
    table["fiber"]=cumul_fiber
    table["wavelength"]=cumul_wave
    table["emission_line_coeff"]=cumul_emission_line_coeff
    table["continuum_coeff"]=cumul_continuum_coeff
    table["x"]=cumul_xc
    table["y"]=cumul_yc
    table["xsig"]=cumul_xsig
    table["ysig"]=cumul_ysig
    table.write(args.output,overwrite=True)
    print("wrote",args.output)

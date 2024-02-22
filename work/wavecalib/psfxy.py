#!/usr/bin/env python

import numpy as np
import sys,os
import argparse
import matplotlib.pyplot as plt
import fitsio
import astropy.io.fits as pyfits
from astropy.table import Table

from desispec.io import read_xytraceset

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,
description='''Measure the resolution of each fiber as a function of wavelength and save it in an ASCII file.''')
parser.add_argument('-i','--infile', type = str, default = None, required = True, nargs="*", help = 'psf fits files')
parser.add_argument('-o','--outfile', type = str, default = None, required = True, help = 'output csv table')
parser.add_argument('--plot', action='store_true')



args = parser.parse_args()

wave  = None
fiber = None

keys=["MJD-OBS","HUMIDITY","PRESSURE","EXPID"]
okeys=["RTEMP","RHUMID"]

cumul_wave=[]
cumul_fiber=[]
cumul_x=[]
cumul_y=[]
cumul_head={k:[] for k in keys}
for k in okeys :
    cumul_head[k]=[]



for filename in args.infile :
    night=filename.split("/")[-3]
    expid=filename.split("/")[-2]
    unit=int(os.path.basename(filename).split("-")[-2][-1])
    print(filename,night,expid,unit)
    preproc=filename.replace("exposures","preproc").replace("fit-psf","preproc")
    raw="/global/cfs/cdirs/desi/spectro/data/{}/{}/desi-{}.fits.fz".format(night,expid,expid)
    print(raw)
    if not os.path.isfile(raw) :
        print(raw,"doesnt exist???")
        sys.exit(12)
    try:
        #head = fitsio.read_header(raw,1)
        head = pyfits.open(raw)[1].header
    except :
        print("failed to read header in hdu 1 of",raw)
        print(sys.exc_info())
        continue
    try :
        t=fitsio.read(raw,"SPECTCONS")
    except :
        print("failed to read SPECTCONS hdu in",raw)
        print(sys.exc_info())
        continue
    try :
        tset = read_xytraceset(filename)
    except :
        print("failed to read",filename)
        print(sys.exc_info())
        continue
    if wave is None :
        wave = (tset.wavemin+tset.wavemax)/2.
        fiber = tset.nspec//2
        print(wave,fiber)
    x=tset.x_vs_wave(fiber,wave)
    y=tset.y_vs_wave(fiber,wave)
    print(x,y)

    cumul_wave.append(wave)
    cumul_fiber.append(fiber)
    cumul_x.append(x)
    cumul_y.append(y)
    for k in keys :
        if k in head :
            cumul_head[k].append(head[k])
        else :
            cumul_head[k].append(0.)
    for k in okeys :
        if k in t.dtype.names :
            cumul_head[k].append(t[k][t["unit"]==unit][0])
        else :
            cumul_head[k].append(0.)
t=Table()
t["wave"]=cumul_wave
t["fiber"]=cumul_fiber
t["x"]=cumul_x
t["y"]=cumul_y
for k in keys :
    t[k]=cumul_head[k]
for k in okeys :
    t[k]=cumul_head[k]
t.write(args.outfile,overwrite=True)
print("wrote",args.outfile)


if args.plot :

    plt.figure("x-vs-mjd")
    plt.plot(t["mjd"],t["x"],"-o")
    plt.figure("y-vs-mjd")
    plt.plot(t["mjd"],t["y"],"-o")
    plt.show()

#!/usr/bin/env python

import sys,os
import argparse
import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt



plt.figure("fiber-cross-talk-vs-wave")
am1=plt.subplot(121,title="fiber-1")
ap1=plt.subplot(122,title="fiber+1")
#am2=plt.subplot(223,title="fiber-2")
#ap2=plt.subplot(224,title="fiber+2")

plt.figure("fiber-cross-talk-vs-fiber")
afm1=plt.subplot(121,title="fiber-1")
afp1=plt.subplot(122,title="fiber+1")

if len(sys.argv)<2 :
    print(sys.argv[0],"fiber-crosstalk-*.fits")
    sys.exit(12)

for filename in sys.argv[1:] :
    print(filename)
    c=filename.split(".")[0].split("-")[-1]
    #filename="fiber-crosstalk-{}.fits".format(c)
    if not os.path.isfile(filename) : continue
    table=Table.read(filename)
    print(table.dtype.names)

    wave=table["WAVELENGTH"]

    biw=np.where((wave>=9600)&(wave<=9700))[0]
    #biw=np.where((wave>=4000)&(wave<=5000))[0]
    #biw=np.where((wave>=5500)&(wave<=7000))[0]

    bfiber=(np.arange(20)+0.5)*25
    bval=np.zeros(20)
    berr=np.zeros(20)
    sw=np.zeros(wave.size)
    swx=np.zeros(wave.size)
    for bundle in range(20) :
        vals=table["CROSSTALK-B{:02d}-F-1".format(bundle)]
        ivar=table["CROSSTALKIVAR-B{:02d}-F-1".format(bundle)]
        sw += ivar
        swx += ivar*vals
        swb=np.sum(ivar[biw])
        if swb>0 :
            bval[bundle]=np.sum(ivar[biw]*vals[biw])/swb
            berr[bundle]=np.sqrt(1/swb)
    vals=swx/(sw+(sw==0))

    am1.plot(wave,vals,label=c+" mean of all bundles")
    ok=(bval!=0)
    afm1.errorbar(bfiber[ok],bval[ok],berr[ok],fmt="o")

    wave=table["WAVELENGTH"]
    bval=np.zeros(20)
    berr=np.zeros(20)
    sw=np.zeros(wave.size)
    swx=np.zeros(wave.size)
    for bundle in range(20) :
        vals=table["CROSSTALK-B{:02d}-F+1".format(bundle)]
        ivar=table["CROSSTALKIVAR-B{:02d}-F+1".format(bundle)]
        sw += ivar
        swx += ivar*vals
        swb=np.sum(ivar[biw])
        if swb>0 :
            bval[bundle]=np.sum(ivar[biw]*vals[biw])/swb
            berr[bundle]=np.sqrt(1/swb)
    vals=swx/(sw+(sw==0))

    ap1.plot(wave,vals,label=c+" mean of all bundles")
    ok=(bval!=0)
    afp1.errorbar(bfiber[ok],bval[ok],berr[ok],fmt="o")





    #vals=table["CROSSTALK-B{:02d}-F+1".format(bundle)]
    #ok=(vals!=0)
    #if np.sum(ok)>0 : ap1.plot(wave[ok],table["CROSSTALK-B{:02d}-F+1".format(bundle)][ok],label=c)

    #ap1.plot(table["WAVELENGTH"],table["CROSSTALK+1"],label=c)
    #am2.plot(table["WAVELENGTH"],table["CROSSTALK-2"],label=c)
    #ap2.plot(table["WAVELENGTH"],table["CROSSTALK+2"],label=c)

am1.grid()
ap1.grid()
#am2.grid()
#ap2.grid()
afm1.grid()
afp1.grid()
afm1.set_ylim([0,0.012])
afp1.set_ylim([0,0.012])

am1.legend(loc="upper left")

plt.show()

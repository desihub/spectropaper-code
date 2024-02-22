#!/usr/bin/env python

"""

"""

import os,sys,glob
import argparse
import numpy as np
import astropy.io.fits as pyfits
import matplotlib.pyplot as plt
from desispec.util import parse_fibers
from desispec.maskbits import specmask
from desispec.fluxcalibration import isStdStar

def reshape(itmp) :
    tmp=np.array(itmp)
    tmp=tmp.reshape((tmp.shape[0]*tmp.shape[1],tmp.shape[2]))
    return tmp


parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i','--infile', type = str, default = None, required = True, nargs="*",
                    help = 'path to frame')
parser.add_argument('-t','--title', type = str, default = "skyres", required = False)
parser.add_argument('--show-sky',action='store_true')
parser.add_argument('--show-all-fibers',action='store_true')

args = parser.parse_args()


flux=[]
ivar=[]
sky=[]

for filename in args.infile :
    print(filename)
    h=pyfits.open(filename)
    wave=h["WAVELENGTH"].data
    fibermap=h["FIBERMAP"].data
    fibers=np.where(fibermap["OBJTYPE"]=="SKY")[0] # sky fibers

    print("fibers=",fibers)
    if fibers.size == 0 :
        print("no selected fibers among sky fibers of",filename)
        continue
    camera=h[0].header["camera"].strip()
    print(camera)
    # find corresponding sky
    sky_filename = filename.replace("sframe-","sky-").replace("-modmask","")

    if not os.path.isfile(sky_filename) :
        print("no",sky_filename)
        continue
    print("USING SKY: ",sky_filename)
    h_sky=pyfits.open(sky_filename)

    maskout = specmask.mask('SOMEBADPIX|ALLBADPIX|COSMIC|LOWFLAT|BADFIBERFLAT|BAD2DFIT|NODATA')

    for fiber in fibers :

        fflux=h[0].data[fiber]
        fivar=h["ivar"].data[fiber]*((h["mask"].data[fiber]&maskout)==0)
        fsky=h_sky["SKY"].data[fiber]
        sky_ivar=h_sky["ivar"].data[fiber]*(h_sky["mask"].data[fiber]==0)

        if np.sum(fivar>0) < 10 :
            print("skip fiber",fiber)
            continue
        flux.append(h[0].data[fiber])
        ivar.append(fivar)
        sky.append(fsky)

    h.close()
    h_sky.close()


flux=np.array(flux)
ivar=np.array(ivar)
sky=np.array(sky)
res = flux.copy()

meansky=np.zeros(res.shape[1])
meanres=np.zeros(res.shape[1])
rms=np.zeros(res.shape[1])
rms_err=np.zeros(res.shape[1])
erms=np.zeros(res.shape[1])
erms03=np.zeros(res.shape[1])
erms1=np.zeros(res.shape[1])
erms2=np.zeros(res.shape[1])
erms3=np.zeros(res.shape[1])
#erms10=np.zeros(res.shape[1])
add_read_noise=0.

for j in range(res.shape[1]) :
    ok=np.where(ivar[:,j]>1/1000.**2)[0]
    if ok.size==0 :
        continue
    meansky[j]=np.median(sky[ok,j])
    meanres[j]=np.median(res[ok,j])
    rms[j]=np.sqrt(np.mean(res[ok,j]**2))

    if 1 : # use mean for expected variance
        erms[j]=np.sqrt(np.mean(1./ivar[ok,j])+add_read_noise**2)
        rms_err[j]=erms[j]/np.sqrt(2*ok.size)
        erms03[j]=np.sqrt(np.mean(1./ivar[ok,j]+(0.003*sky[ok,j])**2)+add_read_noise**2)
        erms1[j]=np.sqrt(np.mean(1./ivar[ok,j]+(0.01*sky[ok,j])**2)+add_read_noise**2)
        erms2[j]=np.sqrt(np.mean(1./ivar[ok,j]+(0.02*sky[ok,j])**2)+add_read_noise**2)
        erms3[j]=np.sqrt(np.mean(1./ivar[ok,j]+(0.03*sky[ok,j])**2)+add_read_noise**2)
        #erms10[j]=np.sqrt(np.mean(1./ivar[ok,j]+(0.1*sky[ok,j])**2)+add_read_noise**2)
    else : # use median for expected variance to avoid fibers with cosmics boosting the variance
        erms[j]=np.sqrt(np.mean(1./ivar[ok,j])+add_read_noise**2)
        rms_err[j]=erms[j]/np.sqrt(2*ok.size)
        erms03[j]=np.sqrt(np.mean(1./ivar[ok,j]+(0.003*sky[ok,j])**2)+add_read_noise**2)
        erms1[j]=np.sqrt(np.mean(1./ivar[ok,j]+(0.01*sky[ok,j])**2)+add_read_noise**2)
        erms2[j]=np.sqrt(np.mean(1./ivar[ok,j]+(0.02*sky[ok,j])**2)+add_read_noise**2)
        erms3[j]=np.sqrt(np.mean(1./ivar[ok,j]+(0.03*sky[ok,j])**2)+add_read_noise**2)
        #erms10[j]=np.sqrt(np.mean(1./ivar[ok,j]+(0.1*sky[ok,j])**2)+add_read_noise**2)



owave=wave.copy()
ok=np.where((rms>0))[0]
wave=wave[ok]
meanres=meanres[ok]
meansky=meansky[ok]
rms=rms[ok]
erms=erms[ok]
erms03=erms03[ok]
erms1=erms1[ok]
erms2=erms2[ok]
erms3=erms3[ok]
#erms10=erms10[ok]

figname="residuals"
if args.title is not None :
    figname=args.title
plt.figure(figname)

a0=plt.subplot(211)
a1=plt.subplot(223)
a2=plt.subplot(224)

nsig=3

select0=np.repeat(True,wave.shape)
select1=(wave>7500)&(wave<7600)
select2=(wave>9300)&(wave<9400)

for a,s in zip([a0,a1,a2],[select0,select1,select2]) :

    #for i in range(res.shape[0]) :
    #    a.plot(wave[s],res[i,s]*(ivar[i,s]>0),c="purple",alpha=0.3)

    a.fill_between(wave[s],erms[s]-nsig*rms_err[s],erms[s]+nsig*rms_err[s],color="gray",alpha=0.7,label="expected noise ($\pm {} \sigma$)".format(nsig))
    a.plot(wave[s],erms2[s]+nsig*rms_err[s],":",c="gray",alpha=0.8,label=None)
    a.plot(wave[s],erms1[s]+nsig*rms_err[s],":",c="gray",alpha=0.8,label=None)

    a.plot(wave[s],rms[s],c="C0",label="r.m.s. < expected + 1% sky")

    ii=np.where(rms[s]>erms3[s]+nsig*rms_err[s])[0]
    step=(wave[1]-wave[0])/2.
    if ii.size>0 :
        label="r.m.s. > expected + 3% sky".format(nsig)
        for i in ii :
            x=[wave[s][i]-step,wave[s][i],wave[s][i]+step]
            y=np.interp(x,wave[s],rms[s])
            a.plot(x,y,"-",c="red",label=label)
            label=None

    ii=np.where(rms[s]>erms2[s]+nsig*rms_err[s])[0]
    step=(wave[1]-wave[0])/2.
    if ii.size>0 :
        label="r.m.s. > expected + 2% sky".format(nsig)
        for i in ii :
            x=[wave[s][i]-step,wave[s][i],wave[s][i]+step]
            y=np.interp(x,wave[s],rms[s])
            #a.plot(x,y,"-",c="red",label=label)
            label=None

    ii=np.where((rms[s]>erms1[s]+nsig*rms_err[s])&(rms[s]<erms2[s]+nsig*rms_err[s]))[0]
    if ii.size>0 :
        label="r.m.s. > expected + 1% sky".format(nsig)
        for i in ii :
            x=[wave[s][i]-step,wave[s][i],wave[s][i]+step]
            y=np.interp(x,wave[s],rms[s])
            a.plot(x,y,"-",c="C1",label=label)
            label=None

a0.legend(loc="upper left",fontsize="medium",title="sky residuals r.m.s.")
a0.set_xlabel("wavelength (A)")
a0.set_ylabel("electrons/A")
a0.grid()
a1.set_xlabel("wavelength (A)")
a1.set_ylabel("electrons/A")
a1.grid()
a2.set_xlabel("wavelength (A)")
a2.set_ylabel("electrons/A")
a2.grid()
plt.tight_layout()
plt.show()

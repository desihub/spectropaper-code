#!/usr/bin/env python

import sys,os
import astropy.io.fits as pyfits
from astropy.table import Table
import fitsio
import numpy as np
import glob

import matplotlib
import matplotlib.pyplot as plt

from scipy.ndimage.filters import median_filter
from scipy.signal import fftconvolve,correlate

from desiutil.log import get_logger
from desispec.io import findfile,read_fiberflat
from desispec.io.fiberflat_vs_humidity import read_fiberflat_vs_humidity

log=get_logger()

filename=sys.argv[1]
fiberflat,humidity,wave,header = read_fiberflat_vs_humidity(filename)

cam=header["CAMERA"]
fiber=250
fig=plt.figure("fiberflat-vs-humidity",figsize=(5,3))
a=plt.subplot(122)

colors = plt.cm.jet(np.linspace(0,1,humidity.size))

for index in range(humidity.size) :
    flat=fiberflat[index]
    flat /= np.mean(flat,axis=0)[None,:]
    a.plot(wave,flat[fiber],"-",c=colors[index])

norm=matplotlib.colors.Normalize(vmin=np.min(humidity),vmax=np.max(humidity))
cb=plt.colorbar(plt.cm.ScalarMappable(norm=norm,cmap=plt.cm.jet))
cb.set_label("Humidity (%)")
a.set_xlabel("Wavelength (A)")
a.set_ylabel("flatfield variation")
a.grid()



#fig=plt.figure("variation-from-cross-correlation")

a=plt.subplot(121)

#fig.suptitle(f"{cam} fiber #{fiber}")
# average of cross-correlated flat
iref=humidity.size//2

width=300
crossdflat=np.zeros((fiberflat.shape[0],fiberflat.shape[2]))
crossdflat[iref]=fiberflat[iref,fiber]-median_filter(fiberflat[iref,fiber],width)
dx=np.zeros(humidity.size)
u=np.arange(fiberflat.shape[2])
for i in range(humidity.size) :
    if i!=iref :
        dflat=fiberflat[i,fiber]-median_filter(fiberflat[i,fiber],width)
        corr=correlate(dflat,crossdflat[iref],"same")
        j=np.argmax(corr)
        # find barycenter of n pts
        hw=5
        x=np.arange(j-hw,j+hw+1)
        y=corr[j-hw:j+hw+1]
        c=np.polyfit(x,y,2)
        val=-c[1]/2./c[0]
        #val=np.sum(corr[j-hw:j+hw+1]*np.arange(j-hw,j+hw+1))/np.sum(corr[j-hw:j+hw+1])
        dx[i]=val-corr.size//2
        print(i,"dx=",dx[i])
        #plt.clf()
        #plt.plot(corr[j-10:j+11])
        #plt.axvline(val-(j-10))
        #plt.show()
        crossdflat[i]=np.interp(u,u-dx[i],dflat)
    #plt.plot(crossdflat[i],alpha=0.4)

dwdx=np.mean(np.gradient(wave))

dflat=np.median(crossdflat,axis=0)
margin=70
dflat[:margin]=0.
dflat[-margin:]=0.
#plt.subplot(1,2,1)
#plt.plot(wave,dflat,color="k")
#plt.xlabel("Wavelength (A)")
#plt.ylabel("Flatfield variation")
#plt.subplot(1,2,2)
#c=np.polyfit(humidity,dx,1)
#pol=np.poly1d(c)
a.plot(humidity,dx*dwdx,"o",color="k")
hh=np.linspace(5.,45.,100)
power=0.1
y=(humidity**power-humidity[iref]**power)
s=np.sum(y*dx)/np.sum(y**2)
y=s*(hh**power-humidity[iref]**power)
a.plot(hh,y*dwdx,"-",color="gray")
a.grid()
a.set_xlabel("Humidity (%)")
a.set_ylabel("Wavelength shift (A)")
plt.tight_layout()
plt.show()

sys.exit(0)
plt.figure()
noabsflat=np.zeros((fiberflat.shape[0],fiberflat.shape[2]))
for i in range(humidity.size) :
   noabsflat[i]=fiberflat[i,fiber]-np.interp(u,u+dx[i],dflat)
   plt.plot(noabsflat[i],alpha=0.4)
noabsflat=np.median(noabsflat,axis=0)
plt.plot(noabsflat,color="k")

betterflat=np.zeros((fiberflat.shape[0],fiberflat.shape[2]))

s=2.
hw=int(3*s+1)
x=np.linspace(-hw,hw,2*hw+1)
kern=np.exp(-x**2/2/s**2)
kern/=np.sum(kern)

for i in range(humidity.size) :
    betterflat[i]=noabsflat+fftconvolve(fiberflat[i,fiber]-noabsflat,kern,mode="same")


plt.figure()

i1=0
i2=i1+2
x1=humidity[i1]
x2=humidity[i2]
flat1=fiberflat[i1,fiber]
flat2=fiberflat[i2,fiber]
i3=i1+1
x=humidity[i3]

interflat=((x2-x)*fiberflat[i1,fiber]+(x-x1)*fiberflat[i2,fiber])/(x2-x1)

dx1=dx[i1]
dx2=dx[i2]
flat1=fiberflat[i1,fiber]-np.interp(u,u+dx1,dflat)
flat2=fiberflat[i2,fiber]-np.interp(u,u+dx2,dflat)
dxx=((x2-x)*dx1+(x-x1)*dx2)/(x2-x1)
interflat2=((x2-x)*flat1+(x-x1)*flat2)/(x2-x1) + np.interp(u,u+dxx,dflat)

#flat1=fiberflat[i1,fiber]/(1+np.interp(u,u+dx1,dflat))
#flat2=fiberflat[i2,fiber]/(1+np.interp(u,u+dx2,dflat))
#dxx=((x2-x)*dx1+(x-x1)*dx2)/(x2-x1)
#interflat3=np.exp(((x2-x)*np.log(flat1)+(x-x1)*np.log(flat2))/(x2-x1) + np.log(1+np.interp(u,u+dxx,dflat)))




#interflat3=((x2-x)*betterflat[i1]+(x-x1)*betterflat[i2])/(x2-x1)


#flat1=fftconvolve(fiberflat[i1,fiber]/noabsflat,kern,mode="same")
#flat2=fftconvolve(fiberflat[i2,fiber]/noabsflat,kern,mode="same")
#interflat4=np.exp(((x2-x)*np.log(flat1)+(x-x1)*np.log(flat2))/(x2-x1))*noabsflat


plt.plot(wave,fiberflat[i3,fiber],label="data")
plt.plot(wave,interflat,"--",label="interflat1")
plt.plot(wave,interflat2,"--",label="interflat2")
#plt.plot(wave,interflat3,"--",label="interflat3")
#plt.plot(wave,interflat4,"--",label="interflat4")
plt.legend()
plt.show()

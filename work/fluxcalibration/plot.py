#!/usr/bin/env python

#import astropy.units as units

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table


def plotprof(x,y,color="k",fmt="o",plot=True) :
    bins=np.exp(np.linspace(0.,np.log(np.max(x)*1.001),8))
    h1,_=np.histogram(x,bins=bins)
    hx,_=np.histogram(x,bins=bins,weights=x)
    hy,_=np.histogram(x,bins=bins,weights=y)
    hy2,_=np.histogram(x,bins=bins,weights=y**2)
    kk=(h1>1)
    mx=hx[kk]/h1[kk]
    my=hy[kk]/h1[kk]
    rmsy=np.sqrt((hy2[kk]/h1[kk]-my**2)*(h1[kk]/(h1[kk]-1)))
    if plot :
        plt.errorbar(mx,my,rmsy,fmt=fmt,color=color)
    for i in range(mx.size) :
        print("{} mx={:.3f} my={:.3f} rmsy={:.3f}".format(i,mx[i],my[i],rmsy[i]))



t=Table.read("fluxes.csv")
print("read fluxes.csv")



ii=(t["IMAGING_FLUX_R"]>10)&(t["IMAGING_FIBERFLUX_R"]>1)&(t["SPECTRO_FLUX_R"]>10)

# rm 2 outliers
bad = ((t["IMAGING_FLUX_R"]>200)&(t["SPECTRO_FLUX_R"]<20))
bad |= ((t["IMAGING_FLUX_R"]<20)&(t["SPECTRO_FLUX_R"]>200))
print("rm",np.sum(bad),"outliers")
ii &= ~bad



plt.figure()
plt.plot(t["IMAGING_FLUX_R"][ii],t["IMAGING_FIBERFLUX_R"][ii]/t["IMAGING_FLUX_R"][ii],".")

jj=(ii)&(t["IMAGING_FIBERFLUX_R"]>0.77*t["IMAGING_FLUX_R"])


t["MAG_R"]=-2.5*np.log10(t["IMAGING_FLUX_R"])+22.5
t["FIBERMAG_R"]=-2.5*np.log10(t["IMAGING_FIBERFLUX_R"])+22.5


fig=plt.figure("flux-ratio",figsize=(5,3))
#fig.suptitle("Flux ratio in r-band (spectroscopy/imaging)")
if 0 :
    a=plt.subplot(221)
    a.plot(t["IMAGING_FLUX_R"][ii],t["SPECTRO_FLUX_R"][ii],".")
    a.plot(t["IMAGING_FLUX_R"][jj],t["SPECTRO_FLUX_R"][jj],".")
    a.grid()
    a=plt.subplot(222)
    a.plot(t["IMAGING_FIBERFLUX_R"][ii],t["SPECTRO_FIBERFLUX_R"][ii],".")
    a.plot(t["IMAGING_FIBERFLUX_R"][jj],t["SPECTRO_FIBERFLUX_R"][jj],".")
    a.grid()

a=plt.subplot(221)
a.plot(t["FIBERMAG_R"][ii],t["SPECTRO_FLUX_R"][ii]/t["IMAGING_FLUX_R"][ii],".")
a.plot(t["FIBERMAG_R"][jj],t["SPECTRO_FLUX_R"][jj]/t["IMAGING_FLUX_R"][jj],".")
#x=t["FIBERMAG_R"][jj]
#y=t["SPECTRO_FLUX_R"][jj]/t["IMAGING_FLUX_R"][jj]
#plotprof(x,y,color="k",fmt="o")


a.set_ylim([0.5,1.5])
a.set_xlabel("FIBERMAG_R")
a.set_ylabel("FLUX_R ratio")
a.grid()

a=plt.subplot(222)
a.plot(t["IMAGING_FIBERFLUX_R"][ii]/t["IMAGING_FLUX_R"][ii],t["SPECTRO_FLUX_R"][ii]/t["IMAGING_FLUX_R"][ii],".")
a.plot(t["IMAGING_FIBERFLUX_R"][jj]/t["IMAGING_FLUX_R"][jj],t["SPECTRO_FLUX_R"][jj]/t["IMAGING_FLUX_R"][jj],".")
a.set_ylim([0.,2])
a.set_xlabel("Fraction of flux in fiber")
a.set_ylabel("FLUX_R ratio")
a.grid()

a=plt.subplot(223)
a.plot(t["FIBERMAG_R"][ii],t["SPECTRO_FIBERFLUX_R"][ii]/t["IMAGING_FIBERFLUX_R"][ii],".")
a.plot(t["FIBERMAG_R"][jj],t["SPECTRO_FIBERFLUX_R"][jj]/t["IMAGING_FIBERFLUX_R"][jj],".")
x=t["FIBERMAG_R"][ii]
y=t["SPECTRO_FIBERFLUX_R"][ii]/t["IMAGING_FIBERFLUX_R"][ii]
plotprof(x,y,color="k",fmt="o",plot=False)
a.set_ylim([0.5,1.5])
a.set_xlabel("FIBERMAG_R")
a.set_ylabel("FIBERFLUX_R ratio")

a.grid()

a=plt.subplot(224)
a.plot(t["IMAGING_FIBERFLUX_R"][ii]/t["IMAGING_FLUX_R"][ii],t["SPECTRO_FIBERFLUX_R"][ii]/t["IMAGING_FIBERFLUX_R"][ii],".")
a.plot(t["IMAGING_FIBERFLUX_R"][jj]/t["IMAGING_FLUX_R"][jj],t["SPECTRO_FIBERFLUX_R"][jj]/t["IMAGING_FIBERFLUX_R"][jj],".")
a.set_ylim([0.,2])
a.set_xlabel("Fraction of flux in fiber")
a.set_ylabel("FIBERFLUX_R ratio")
a.grid()

#plt.xlabel("FIBERFLUX_R (imaging)")
#plt.ylabel("FIBERFLUX_R (spectroscopy)")

plt.tight_layout()
plt.show()

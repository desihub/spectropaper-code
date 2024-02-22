#!/usr/bin/env python

import os,sys
import glob
import fitsio
import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
from desitarget.targetmask import desi_mask

ofilename="lrg-tsnr2-deltachi2-guadalupe.csv"

print("reading",ofilename)
t=Table.read(ofilename)
x=t["MEDIAN_TSNR2_LRG"]
y=t["MEDIAN_DELTACHI2"]

dx=x-np.mean(x)
dy=y-np.mean(y)
corr=np.mean(dx*dy)/np.sqrt(np.mean(dx**2)*np.mean(dy**2))
print("correlation coefficent=",corr)

mx=np.mean(x)
ii=(x>mx)
dx=x[ii]-np.mean(x[ii])
dy=y[ii]-np.mean(y[ii])
corr2=np.mean(dx*dy)/np.sqrt(np.mean(dx**2)*np.mean(dy**2))
print("correlation coefficent (for half with TSNR>={}) = {}".format(mx,corr2))


snr2time = fitsio.read_header(os.environ["DESIMODEL"]+"/data/tsnr/tsnr-ensemble-lrg.fits")["SNR2TIME"]
print("snr2time=",snr2time)

ref_snr2 = 1000./snr2time

plt.figure("deltachi2-vs-tsnr2",figsize=(5,4))
a=plt.subplot(111,title="Median per exposure for LRG targets")
plt.plot(t["MEDIAN_TSNR2_LRG"],t["MEDIAN_DELTACHI2"],".")
plt.text(0.02,0.94,"Correlation coefficient = {:.2f}".format(corr),transform=a.transAxes)
plt.xlabel("LRG TSNR$^2$")
plt.ylabel("redrock $\Delta \chi^2$")
plt.axvline(ref_snr2,color="k",linestyle="--")
plt.grid()

x=t["MEDIAN_TSNR2_ELG"]
dx=x-np.mean(x)
dy=y-np.mean(y)
corr_bis=np.mean(dx*dy)/np.sqrt(np.mean(dx**2)*np.mean(dy**2))
print("correlation coefficent (wrong TSNR)=",corr_bis)
mx=np.mean(x)
ii=(x>mx)
dx=x[ii]-np.mean(x[ii])
dy=y[ii]-np.mean(y[ii])
corr2=np.mean(dx*dy)/np.sqrt(np.mean(dx**2)*np.mean(dy**2))
print("correlation coefficent (wrong TSNR) (for half with TSNR>={}) = {}".format(mx,corr2))

plt.figure("deltachi2-vs-wrong-tsnr2")
a=plt.subplot(111,title="Median per exposure for LRG targets")
plt.plot(t["MEDIAN_TSNR2_ELG"],t["MEDIAN_DELTACHI2"],".")
plt.text(0.02,0.94,"Correlation coefficient = {:.2f}".format(corr_bis),transform=a.transAxes)
plt.xlabel("ELG TSNR$^2$ (for LRG targets)")
plt.ylabel("redrock $\Delta \chi^2$")
plt.grid()

plt.tight_layout()

plt.show()

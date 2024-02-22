#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import astropy.io.fits as pyfits

plt.figure("sky-throughput-correction",figsize=(5,3))
labels={"b":"blue","r":"red","z":"NIR"}
for cam,color in zip(["b","r","z"],["b","r","k"]) :
    ifile="throughput-correction-{}.csv".format(cam)
    t=Table.read(ifile)
    corr=t["CORR"]
    err=t["ERR"]
    bins=np.linspace(0.91,1.076,50)
    ok=(err>0)&(err<0.03)
    histo,_ = np.histogram(corr[ok],bins=bins)
    center=bins[:-1]+(bins[1]-bins[0])/2.
    #plt.plot(center,histo,"-",label="{} cameras\nrms={:4.3f}, $\sigma$={:4.3f}".format(labels[cam],np.std(corr[ok]),np.sqrt(np.mean(err[ok]**2))),color=color)
    plt.step(center,histo,where="mid",label="{} cameras\nrms={:4.3f}, $\sigma$={:4.3f}".format(labels[cam],np.std(corr[ok]),np.sqrt(np.mean(err[ok]**2))),color=color)
    #plt.hist(corr[ok],bins=bins,alpha=0.4,label="{} cameras\nrms={:4.3f}, $\sigma$={:4.3f}".format(labels[cam],np.std(corr[ok]),np.sqrt(np.mean(err[ok]**2))),edgecolor=color,color=None)
plt.legend(loc="upper left")
plt.grid()
plt.locator_params(axis='x', nbins=3)
plt.xlabel("Throughput Correction")

plt.tight_layout()
plt.show()

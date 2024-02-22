#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import astropy.io.fits as pyfits

#import glob
#ifiles=np.sort(glob.glob("epsilon_wave_*.csv"))

ifiles=sys.argv[1:]
print(ifiles)
x=[]
y=[]

for ifile in ifiles :
    t=Table.read(ifile)
    ii=(t["WAVE"]<8600)|(t["WAVE"]>8645) # avoid broad blend of lines
    x.append(t["WAVE"][ii])
    y.append(t["SIGMA"][ii])
x=np.hstack(x)
y=np.hstack(y)

plt.figure("sky-epsilon-lambda",figsize=(5,3))
#plt.plot(x,y,".",alpha=0.7,color="gray")

# median sigma per line to avoid effect of cosmic rays
xx=[]
yy=[]
ee=[]
for tx in np.unique(x) :
    ii=(x==tx)
    nmeas=np.sum(ii)
    if nmeas>3 :
        xx.append(tx)
        med=np.median(y[ii])
        #rms=1.48*np.median(np.abs(y-med))
        ok=ii&(np.abs(y-med)<0.2)
        yy.append(np.mean(y[ok]))
        ee.append(np.std(y[ok])/np.sqrt(np.sum(ok)))

#plt.plot(x,y,".")
#plt.plot(xx,yy,"o")
plt.errorbar(xx,yy,ee,fmt="o")


plt.xlabel("Wavelength ($\AA$)")
plt.ylabel("$\epsilon_{\lambda}$ ($\AA$)")
plt.grid()
plt.tight_layout()
plt.show()

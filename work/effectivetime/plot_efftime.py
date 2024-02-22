#!/usr/bin/env python

import os,sys
import glob
import fitsio
import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
#from desitarget.targetmask import desi_mask

from desispec.efftime import compute_efftime

t=Table.read("/global/cfs/cdirs/desi/spectro/redux/guadalupe/exposures-guadalupe.csv")
ii=(t["FAFLAVOR"]=="maindark")&(t["EFFTIME_GFA"]>0)

t=t[ii]
snr2time = fitsio.read_header(os.environ["DESIMODEL"]+"/data/tsnr/tsnr-ensemble-lrg.fits")["SNR2TIME"]
print("snr2time=",snr2time)

plt.figure("efftime",figsize=(5,4))
plt.plot(t["EFFTIME_GFA"],t["LRG_EFFTIME_DARK"],".")
#ok=(t["LRG_EFFTIME_DARK"]>800)&(t["EFFTIME_GFA"]>0)&(t["EFFTIME_GFA"]<2000)
#plt.plot(t["EFFTIME_GFA"][ok],t["LRG_EFFTIME_DARK"][ok]/t["EFFTIME_GFA"][ok],".")
#print(np.mean(t["LRG_EFFTIME_DARK"][ok]/t["EFFTIME_GFA"][ok]))
a=np.sum(t["LRG_EFFTIME_DARK"]*t["EFFTIME_GFA"])/np.sum(t["EFFTIME_GFA"]**2)
print(a)
u=[0,1400]
plt.plot(u,u,color="gray")

x=t["EFFTIME_GFA"]
y=t["LRG_EFFTIME_DARK"]
dx=x-np.mean(x)
dy=y-np.mean(y)
corr=np.mean(dx*dy)/np.sqrt(np.mean(dx**2)*np.mean(dy**2))
print("corr=",corr)

plt.xlabel(r"$T_{r}$ (seconds)")
plt.ylabel(r"$T_{spec}$ (seconds)")
plt.grid()
plt.tight_layout()

#plt.figure()
#plt.plot(t["EFFTIME_ETC"],t["LRG_EFFTIME_DARK"],".")

#plt.plot(t["EFFTIME_ETC"][ok],t["LRG_EFFTIME_DARK"][ok]/t["EFFTIME_ETC"][ok],".")
#print(np.mean(t["LRG_EFFTIME_DARK"][ok]/t["EFFTIME_ETC"][ok]))

#efftime_dark , efftime_bright , efftime_backup = compute_efftime(t)

#plt.figure()
#plt.plot(t["EFFTIME_GFA"],efftime_dark,".")

plt.show()

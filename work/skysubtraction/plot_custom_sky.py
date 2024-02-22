#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
import astropy.io.fits as pyfits

z=Table.read("sky-z0.csv")
covar=pyfits.open("sky-covar-z0.fits")[0].data
ccovar=pyfits.open("csky-covar-z0.fits")[0].data

err=np.diag(covar).copy()
cerr=np.diag(ccovar).copy()
wave=z["WAVE"]
plt.figure("sky-convolved-deconvolved",figsize=(5,3))
plt.errorbar(wave,z["SKY"]/1e4,err/1e4,fmt="o",label=r"deconvolved sky model $S$")
plt.errorbar(wave,z["CSKY"]/1e4,cerr/1e4,fmt="o",label=r"reconvolved sky model $\tilde{S}$")

plt.plot(wave,np.array(z["SKY"])/1e4,"-",c="gray",alpha=0.8)
plt.plot(wave,np.array(z["CSKY"])/1e4,"-",c="gray",alpha=0.8)
plt.xlabel(r"Wavelength ($\AA$)")
plt.ylabel(r"Counts ($10^4 \, e^- \, \AA^{-1}$)")
plt.legend(loc="upper left")
plt.xlim([8881,8929])
plt.grid()
plt.tight_layout()

corr_covar = covar.copy()
for i in range(covar.shape[0]) :
    corr_covar[i,:] /= np.sqrt(covar[i,i])
    corr_covar[:,i] /= np.sqrt(covar[i,i])
corr_ccovar = ccovar.copy()
for i in range(ccovar.shape[0]) :
    corr_ccovar[i,:] /= np.sqrt(ccovar[i,i])
    corr_ccovar[:,i] /= np.sqrt(ccovar[i,i])

if 0 :
    plt.figure("correlation-matrix")
    print(np.min(corr_covar),np.max(corr_covar))
    plt.subplot(1,2,1,title="deconvolved sky model")
    plt.imshow(corr_covar[1000:1010,1000:1010],origin="lower",vmin=-1,vmax=1)
    plt.subplot(1,2,2,title="reconvolved sky model")
    plt.imshow(corr_ccovar[1000:1010,1000:1010],origin="lower",vmin=-1,vmax=1)
    plt.colorbar()
    plt.tight_layout()

#fig=plt.figure("sky-spectrum-correlation",figsize=(5,5))
fig, ax = plt.subplots(2, 1, sharex=True,figsize=(5,5),num="sky-spectrum-correlation")
nn=covar.shape[0]
i=np.arange(nn-4)
a1=ax[0]
a2=ax[1]
a1.set_title("deconvolved spectrum",fontsize="medium")
#a1=plt.subplot(211,title="deconvolved spectrum")
a1.plot(wave[i],corr_covar[i,i+1],label=r"$\Delta=1$")
a1.plot(wave[i],corr_covar[i,i+2],label=r"$\Delta=2$")
a1.plot(wave[i],corr_covar[i,i+3],label=r"$\Delta=3$")
a1.set_ylabel("$corr[i,i+\Delta]$")
a1.legend(loc="upper left")
a1.set_ylim([-1.09,1.59])
a1.grid()
a2.set_title("reconvolved spectrum",fontsize="medium")
#a2=plt.subplot(212,title="reconvolved spectrum")
a2.plot(wave[i],corr_ccovar[i,i+1],label=r"$\Delta=1$")
a2.plot(wave[i],corr_ccovar[i,i+2],label=r"$\Delta=2$")
a2.plot(wave[i],corr_ccovar[i,i+3],label=r"$\Delta=3$")
#a2.legend(loc="upper left")
a2.set_ylabel("$corr[i,i+\Delta]$")
a2.set_xlabel(r"Wavelength ($\AA$)")
a2.set_ylim([-0.099,0.099])
a2.grid()
plt.tight_layout()
plt.show()

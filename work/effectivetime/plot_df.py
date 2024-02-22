#!/usr/bin/env python

import os
import glob
import astropy.io.fits as pyfits
from astropy.convolution import convolve, Box1DKernel
import matplotlib.pyplot as plt

scaled_df=["BGS","Lya"]
#scaled_df=["BGS"]

plt.figure("tsnr-deltaf",figsize=(5,3))
smooth=1
for i,target in enumerate(["LRG","ELG","QSO","Lya","BGS"]) :
    filename="{}/data/tsnr/tsnr-ensemble-{}.fits".format(os.environ["DESIMODEL"],target.lower())
    print(filename)
    color="C{}".format(i)
    h=pyfits.open(filename)
    label=target
    if target in scaled_df :
        label += r" $(\times 0.1)$"

    for b in "BRZ" :
        flux = convolve(h["DFLUX_"+b].data[0], Box1DKernel(smooth), boundary='extend')
        if target in scaled_df :
            flux /= 10.
        plt.plot(h["WAVE_"+b].data,flux,color=color,label=label)
        label=None
plt.legend(loc="upper left")#,fontsize="x-small")
plt.xlabel("Wavelength ($\AA$)")
plt.ylabel(r"$< (\delta F)^2 > ^{1/2} \ \ (10^{-17} \, \rm{erg} \, \rm{s}^{-1} \, \rm{cm}^{-2} \, \AA^{-1})$")
plt.grid()
plt.tight_layout()
plt.show()

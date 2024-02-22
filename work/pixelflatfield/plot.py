#!/usr/bin/env python


import numpy as np
import fitsio
import matplotlib.pyplot as plt


def tune_axes(a) :
    a.axes.xaxis.set_ticklabels([])
    a.axes.yaxis.set_ticklabels([])
    a.tick_params(bottom = False)
    a.tick_params(left = False)

img1=fitsio.read("mean-r8.fits")[1950:2250,3240:3540]
img2=fitsio.read("pixflat-sm2-r8-20210816.fits.gz")[1950:2250,3240:3540]

mm=np.median(img1)
plt.figure("pixelflat")
a=plt.subplot(121)
a.imshow(img1,vmin=mm*0.8,vmax=mm*1.2,origin="lower",cmap="gray")
tune_axes(a)
a=plt.subplot(122)
a.imshow(img2,vmin=0.95,vmax=1.05,origin="lower",cmap="gray")
tune_axes(a)
plt.tight_layout()
plt.show()

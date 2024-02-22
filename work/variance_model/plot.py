#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import fitsio

img1=fitsio.read("preproc-r0-20220307-e-without-variance-model.fits","IVAR")[1200:1350,3300:3450]
img2=fitsio.read("preproc-r0-20220307-e-with-variance-model.fits","IVAR")[1200:1350,3300:3450]

def tune_axes(a) :
    a.axes.xaxis.set_ticklabels([])
    a.axes.yaxis.set_ticklabels([])
    a.tick_params(bottom = False)
    a.tick_params(left = False)


plt.figure("variance-model",figsize=(8,4))

vmin=0
vmax=0.2
a=plt.subplot(121)
a.imshow(img1,vmin=vmin,vmax=vmax,origin="lower",cmap="gray")
tune_axes(a)
a=plt.subplot(122)
a.imshow(img2,vmin=vmin,vmax=vmax,origin="lower",cmap="gray")
tune_axes(a)
plt.tight_layout()
plt.show()

#!/usr/bin/env python


import numpy as np
import matplotlib.pyplot as plt
import fitsio

img1=fitsio.read("preproc-r3-00125691-nocorr.fits")[1170:1370,3220:3420]
img2=fitsio.read("preproc-r3-00125691-withcorr.fits")[1170:1370,3220:3420]




def tune_axes(a) :
    a.axes.xaxis.set_ticklabels([])
    a.axes.yaxis.set_ticklabels([])
    a.tick_params(bottom = False)
    a.tick_params(left = False)


plt.figure("negative-trails-sm6-r3",figsize=(8,4))
vmin=-10
vmax=10
textsize=18
a=plt.subplot(121)
a.set_title("without correction", size=textsize)
a.imshow(img1,vmin=vmin,vmax=vmax,origin="lower",cmap="gray_r")
tune_axes(a)
a=plt.subplot(122)
a.set_title("with correction", size=textsize)
a.imshow(img2,vmin=vmin,vmax=vmax,origin="lower",cmap="gray_r")
tune_axes(a)
plt.tight_layout()
plt.show()

#!/usr/bin/env python

import numpy as np
import os,sys
import matplotlib.pyplot as plt
from astropy.table import Table


mmin=58840
#mmin=0

if 0 :
    print("read telemetry...")
    #telem=Table.read("telemetry.csv")
    #telem.write("telemetry.fits")
    telem=Table.read("telemetry.fits")
    print("done")

camera="blue" ; s=0.4
#camera="red" ; s=0.2
#camera="nir" ; s=0.05

t={}

#for cam in ["b8","r8","z8"] :
for cam in ["b4","r4","z4"] :
#for cam in ["r0","r1","r2","r3","r4","r5","r6","r7","r8","r9"] :
#for cam in ["z0","z1","z2","z3","z4","z5","z6","z7","z8","z9"] :

    filename="psfxy-{}.csv".format(cam)
    if os.path.isfile(filename) :
        t[cam]=Table.read(filename)
        selection=(t[cam]["MJD-OBS"]>mmin)&(t[cam]["x"]!=0)
        t[cam]=t[cam][selection]


#############################################################################
whats=["x","y"]
legend_label=None

fig,aa=plt.subplots(2,3,sharex=True,num="dxdy",gridspec_kw={'hspace': 0,'wspace': 0},figsize=(5,4))
colors=["blue","red","brown"]
camera=["Blue","Red","NIR"]
for i,cam in enumerate(t) :
    #unit=int(cam[-1])
    #color="C{}".format(unit)
    #color="k"

    days=t[cam]["MJD-OBS"]-mmin
    ii=(days<60)
    days=days[ii]
    x=t[cam]["x"][ii]
    y=t[cam]["y"][ii]
    dx=x-np.median(x)
    dy=y-np.median(y)


    a0=aa[0,i]
    a1=aa[1,i]
    a0.plot(days,dx,".",color=colors[i],alpha=0.8)
    a1.plot(days,dy,".",color=colors[i],alpha=0.8)

    if i==0 :
        a0.set_ylabel(r"$\delta X$ (pixels)")
        a1.set_ylabel(r"$\delta Y$ (pixels)")
    else :
        a0.yaxis.set_ticklabels([])
        a1.yaxis.set_ticklabels([])

    a0.set_title(camera[i],color=colors[i])
    a0.grid()
    a1.grid()
    a0.set_ylim([-0.7,0.7])
    a1.set_ylim([-0.199,0.199])

    print("cam {} rms(dx)={:4.3f} max(|dx|)={:4.3f} rms(dy)={:4.3f} max(|dy|)={:4.3f}".format(cam,np.std(dx),np.max(np.abs(dx)),np.std(dy),np.max(np.abs(dy))))

aa[1,1].set_xlabel("days")
#aa[2,0].set_xlabel("days")
#aa[2,1].set_xlabel("days")
plt.tight_layout()
plt.show()

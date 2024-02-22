#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import sys

fiberdiam=107.
pixsize=15.
fiber_spacing=230 #um
coll_fnumber = 3.57 # table 5.1
cam_fnumber = 1.7 # table 5.5 p 87

fiberdiam_pixels = fiberdiam * cam_fnumber / coll_fnumber / pixsize

print("demag=",cam_fnumber / coll_fnumber)
print("fiberdiam= {} um = {} pixels".format(fiberdiam * cam_fnumber / coll_fnumber,fiberdiam_pixels))

plt.figure("fwhm",figsize=(5,3))

tmp=np.loadtxt("fwhm-b1-20200315.txt")
bwave=tmp[:,0]
bfwhm=tmp[:,1:]
mean_bfwhm=np.mean(bfwhm,axis=1)
min_bfwhm=np.min(bfwhm,axis=1)
max_bfwhm=np.max(bfwhm,axis=1)
ok=(bwave>3600)
bwave=bwave[ok]
mean_bfwhm=mean_bfwhm[ok]
min_bfwhm=min_bfwhm[ok]
max_bfwhm=max_bfwhm[ok]

plt.plot(bwave,min_bfwhm,":",c="b")
plt.plot(bwave,max_bfwhm,":",c="b")
plt.plot(bwave,mean_bfwhm,"-",c="b")
plt.fill_between(bwave,min_bfwhm,max_bfwhm,color="b",alpha=0.3)


tmp=np.loadtxt("fwhm-r1-20200315.txt")
rwave=tmp[:,0]
rfwhm=tmp[:,1:]
mean_rfwhm=np.mean(rfwhm,axis=1)
min_rfwhm=np.min(rfwhm,axis=1)
max_rfwhm=np.max(rfwhm,axis=1)
ok=(rwave>5900)
rwave=rwave[ok]
mean_rfwhm=mean_rfwhm[ok]
min_rfwhm=min_rfwhm[ok]
max_rfwhm=max_rfwhm[ok]

plt.plot(rwave,min_rfwhm,":",c="r")
plt.plot(rwave,max_rfwhm,":",c="r")
plt.plot(rwave,mean_rfwhm,c="r")
plt.fill_between(rwave,min_rfwhm,max_rfwhm,color="r",alpha=0.3)


tmp=np.loadtxt("fwhm-z1-20200315.txt")
zwave=tmp[:,0]
zfwhm=tmp[:,1:]
mean_zfwhm=np.mean(zfwhm,axis=1)
min_zfwhm=np.min(zfwhm,axis=1)
max_zfwhm=np.max(zfwhm,axis=1)
ok=(zwave>7600)&(zwave<9800)
zwave=zwave[ok]
mean_zfwhm=mean_zfwhm[ok]
min_zfwhm=min_zfwhm[ok]
max_zfwhm=max_zfwhm[ok]

plt.plot(zwave,min_zfwhm,":",c="brown")
plt.plot(zwave,max_zfwhm,":",c="brown")
plt.plot(zwave,mean_zfwhm,c="brown")
plt.fill_between(zwave,min_zfwhm,max_zfwhm,color="brown",alpha=0.3)


#plt.axhline(fiberdiam_pixels,linestyle="--",color="k")


plt.grid()
plt.xlabel("Wavelength ($\AA$)")
plt.ylabel("FWHM (pixels)")
plt.tight_layout()
plt.show()

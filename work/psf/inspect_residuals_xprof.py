#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits

from desispec.io import read_image
from desispec.io import read_xytraceset



if len(sys.argv)<4 :
    print(sys.argv[0],"res-z1-00055705-10.fits preproc-z1-00055705.fits psf-z1-00055705-10.fits")
    sys.exit(1)

res_filename=sys.argv[1]
preproc_filename=sys.argv[2]
psf_filename=sys.argv[3]

#cam="z1"
#res_filename="res-{}-00055705-10.fits".format(cam)
#preproc_filename="preproc-{}-00055705.fits".format(cam)
#psf_filename="psf-{}-00055705-10.fits".format(cam)

h=pyfits.open(res_filename)
model=h["MODEL"].data
preproc=read_image(preproc_filename)
tset = read_xytraceset(psf_filename)



ii=(model.ravel()>1)
d = preproc.pix.ravel()[ii]
m = model.ravel()[ii]
ivar = preproc.ivar.ravel()[ii]*(preproc.mask.ravel()[ii]!=0)
pull = (d-m)*np.sqrt(ivar)
pull = h["PULL"].data.ravel()[ii]


image=preproc.pix
ivar=preproc.ivar*(preproc.mask==0)
n0=image.shape[0]
n1=image.shape[1]
py=np.mean(model[:,n1//2-10:n1//2+10],axis=1)
j=np.argmax(py[1500:2500])+1500
x=np.arange(n1)
ii=np.where(model[j]>0)[0]
b=ii[0]-3
e=ii[-1]+4

hw=2
x=x[b:e]
pd=np.mean(image[j-hw:j+hw+1,b:e],axis=0)
pm=np.mean(model[j-hw:j+hw+1,b:e],axis=0)
#plt.figure()
fig,aa=plt.subplots(2,1,sharex=True, gridspec_kw={'height_ratios': [1, 1], 'hspace': 0},figsize=(5.5,3.7))
fig.suptitle("Mean of CCD rows {} to {} (electrons per pixel)".format(j-hw,j+hw),size=11)
#aa[0].set_title("Mean of CCD rows {} to {} (electrons per pixel)".format(j-hw,j+hw))
a=aa[0]
a.plot(x,pd,label="data")
a.plot(x,pm,"--",label="model")
a.grid()
#a.set_ylabel("Electrons per pixel")
#a.legend(title="Mean of CCD rows {} to {}".format(j-hw,j+hw))
a=aa[1]#plt.subplot(212)
a.plot(x,pd,label="data")
a.plot(x,pm,"--",label="model")
a.set_yscale("log")
a.set_ylim([0.1,1.1*np.max(pm)])
a.set_xlabel("CCD column (pixels)")
#a.set_ylabel("mean electrons/pixel")
#a.legend(loc="lower center",title="Mean of CCD rows {} to {}".format(j-hw,j+hw),fontsize=8)

a.grid()
plt.tight_layout
plt.show()

#plt.plot(py)
plt.show()
sys.exit(12)
px=np.mean(model,axis=0)
b=np.where(px>0)[0][0]
e=np.where(px>0)[0][-1]
print(b,e)
data_py  = np.sum((ivar[:,b:e]>0)*image[:,b:e],axis=1)/np.sum(ivar[:,b:e]>0,axis=1)
model_py = np.sum((ivar[:,b:e]>0)*model[:,b:e],axis=1)/np.sum(ivar[:,b:e]>0,axis=1)


y=np.arange(n0)
wave=tset.wave_vs_y(fiber=250,y=y)

fig=plt.figure(figsize=(8,3))
ax1 = fig.add_subplot(111)

#fig,ax1 = plt.subplots()#constrained_layout=True)

#def y2w(yy):
#    return np.interp(yy,y,wave)
#def w2y(ww) :
#    return np.interp(ww,wave,y)
#ax2 = ax1.secondary_xaxis('top', functions=(y2w,w2y))

#def tick_function(yy):
#    V=np.interp(yy,y,wave)
#    return ["%3.0f" % v for v in V]
def tick_function(ww):
    V=np.interp(ww,wave,y)
    return ["%3.0f" % v for v in V]

ax1.plot(wave,data_py,"-",label="data",alpha=1,)
ax1.plot(wave,model_py,"--",label="model",alpha=1)
ax1.grid()
ax1.set_xlabel(r"Wavelength ($\AA$)")
ax1.set_ylabel("mean in pixel row (electrons/pix)")
ax1.legend(loc="upper right")
ax1.set_ylim([0,2*np.max(model_py)])
if False :
    new_tick_locations = np.array([0,1000,2000,3000,4000])
    ax2 = ax1.twiny()
    ax2.set_xlim(ax1.get_xlim())
    ax2.set_xticks(new_tick_locations)
    ax2.set_xticklabels(tick_function(new_tick_locations))
    #ax2.set_xlabel(r"Wavelength ($\AA$)")
    ax2.set_xlabel("CCD row (pixels)")


zoomwave=[float(v) for v in sys.argv[4:]]
print(zoomwave)
#for [3651.,4360.,5947.]
for i,zw in enumerate(zoomwave) :
    ax3=fig.add_subplot(2,len(zoomwave)+1,i+1)
    ii=np.abs(wave-zw)<5.
    ax3.plot(wave[ii],data_py[ii],"-",label="data",alpha=1)
    ax3.plot(wave[ii],model_py[ii],"--",label="model",alpha=1)
    ax3.set_yticks([])
    ax3.set_xticks([zw,])
plt.tight_layout()


plt.show()

#!/usr/bin/env python

import sys,os
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
from astropy.table import Table

from desispec.io import read_image
from desispec.io import read_xytraceset



if len(sys.argv)<4 :
    print(sys.argv[0],"res-z1-00055705-10.fits preproc-z1-00055705.fits psf-z1-00055705-10.fits 7604 8780 9803")
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


if False :
    ii=np.where(model.ravel()>1)[0]
    d = preproc.pix.ravel()[ii]
    m = model.ravel()[ii]
    ivar = preproc.ivar.ravel()[ii]*(preproc.mask.ravel()[ii]!=0)
    pull = (d-m)*np.sqrt(ivar)
    pull = h["PULL"].data.ravel()[ii]

    bins=np.linspace(0,10000,100)
    h1,bins=np.histogram(m,bins=bins)
    hx,bins=np.histogram(m,bins=bins,weights=m)
    hy,bins=np.histogram(m,bins=bins,weights=pull)
    hy2,bins=np.histogram(m,bins=bins,weights=pull**2)
    ok=(h1>10)
    h1=h1[ok]
    hx=hx[ok]
    hy=hy[ok]
    hy2=hy2[ok]
    mx=(hx/h1)
    my=(hy/h1)
    rmsy=np.sqrt((hy2/h1)-my**2)


    plt.plot(m,pull,".",alpha=0.1,color="gray")

    plt.plot(mx,rmsy,"-",color="k")
    plt.plot(mx,-rmsy,"-",color="k")
    plt.axhline(1,linestyle=":",color="blue")
    plt.axhline(-1,linestyle=":",color="blue")
    plt.grid()

table_filename=os.path.basename(res_filename).split(".")[0]+".csv"
print(f"table filename='{table_filename}'")

if not os.path.isfile(table_filename) :
    print("computing the spectra ...")
    image=preproc.pix
    ivar=preproc.ivar*(preproc.mask==0)
    n0=image.shape[0]
    n1=image.shape[1]
    px=np.mean(model,axis=0)
    b=np.where(px>0)[0][0]
    e=np.where(px>0)[0][-1]
    print(b,e)
    data_py  = np.sum((ivar[:,b:e]>0)*image[:,b:e],axis=1)/np.sum(ivar[:,b:e]>0,axis=1)
    model_py = np.sum((ivar[:,b:e]>0)*model[:,b:e],axis=1)/np.sum(ivar[:,b:e]>0,axis=1)
    y=np.arange(n0)
    wave=tset.wave_vs_y(fiber=250,y=y)
    t=Table()
    t["wave"]=wave
    t["data"]=data_py
    t["model"]=model_py
    t.write(table_filename,overwrite=False)
else :
    print("reading precomputed spectra")
    t=Table.read(table_filename)
    wave=t["wave"]
    data_py=t["data"]
    model_py=t["model"]

    zenodo=False ; zenodo_filename="../../zenodo/figure-14c.fits"

    if zenodo:
        t2=Table()
        t2["WAVELENGTH"]=wave
        t2["MEASURED_FLUX"]=data_py
        t2["MODEL_FLUX"]=model_py
        t2.write(zenodo_filename)
fig=plt.figure(figsize=(6,2.5))
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

# some lines
lines=[(3651.198,0.,"Hg"),]
lines.append((4359.560,0.,"Hg"))
lines.append((4047.708,0.,"Hg"))
lines.append((4801.254,0,"Cd"))
lines.append((5087.239,0,"Cd"))
lines.append((5462.268,0.,"Hg"))
lines.append((5771.210,-30,"Hg"))
#lines.append((5792.276,0,"Hg"))
lines.append((5854.110,+30.,"Ne"))
#lines.append((5946.481,0.,"Ne"))
lines.append((6144.763,0.,"Ne"))
lines.append((6404.018,0.,"Ne"))
lines.append((6508.325,0.,"Ne"))
lines.append((6680.120,0.,"Ne"))
lines.append((6931.379,-20,"Ne"))
lines.append((6967.352,+20,"Ar"))

lines.append((7034.352,0.,"Ne"))
lines.append((7247.163,0.,"Ne"))
lines.append((7440.947,0.,"Ne"))
lines.append((7603.638,-30.,"Kr "))
lines.append((7637.208,30.," Ar"))
lines.append((7950.362,0,"Ar"))
#lines.append((8233.896,0.,"Xe"))
lines.append((8426.963,0.,"Ar"))
lines.append((9125.471,0.,"Ar"))
lines.append((9227.030,0,"Ar"))
lines.append((9660.435,0.,"Ar"))
lines.append((8779.161,-30.,"Kr"))
lines.append((8821.832,30.,"Xe"))
lines.append((9802.384,0.,"Xe"))

for line in lines :
    lw=line[0]
    if lw>wave[0] and lw<wave[-1] :
        print("adding label for line",line)
        li=np.argmin(np.abs(wave-lw))
        lf=np.max(data_py[li-10:li+11])
        plt.text(lw+line[1],lf*1.05,line[2],horizontalalignment="center",verticalalignment="bottom")
ax1.grid()
ax1.set_xlabel(r"Wavelength ($\AA$)")
ax1.set_ylabel("Electrons per pixel")
ax1.legend(loc="upper right")
ax1.set_ylim([0,2.7*np.max(model_py)])
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

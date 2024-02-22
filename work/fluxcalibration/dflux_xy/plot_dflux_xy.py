#!/usr/bin/env python

import argparse
import os,sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from desispec.fiberfluxcorr import flat_to_psf_flux_correction
#flat_to_psf_flux_correction(fibermap)


def tune_axes(a) :
    a.axes.xaxis.set_ticklabels([])
    a.axes.yaxis.set_ticklabels([])
    a.axes.yaxis.set_visible(False)
    a.tick_params(bottom = False)
    a.tick_params(left = False)

parser = argparse.ArgumentParser()
parser.add_argument('-i','--infiles', type=str, required=True, nargs="*",
                        help = 'Input calibstars files')
parser.add_argument('--what', type=str, required=False, default="broadband",
                        help = 'for figures')
parser.add_argument('--no-colorbar', action='store_true',
                        help = '')

args = parser.parse_args()

xx=[]
yy=[]
ff=[]
zz=[]


et=Table.read("/global/homes/j/jguy/redux/daily/tsnr-exposures.fits")
et=Table.read("./extra-exposures3.fits")
e2i={e:i for i,e in enumerate(et["EXPID"])}

pp1=[]
pp2=[]

for filename in args.infiles :
    #print(filename)

    expid=int(filename.split("-")[-1].split(".")[0])

    index=e2i[expid]
    seeing=et["SEEING_GFA"][index]

    #adc1phi=et["ADC1PHI"][index]
    #adc2phi=et["ADC2PHI"][index]
    #if adc2phi > adc1phi + 180 : adc2phi -= 360.
    #if adc1phi > adc2phi + 180 : adc1phi -= 360.
    #pp1.append(adc1phi)
    #pp2.append(adc2phi)
    #phi=(adc1phi+adc2phi)/2.
    #cp=np.cos(phi/180*np.pi)
    #sp=np.sin(phi/180*np.pi)
    #M=np.array([[cp,-sp],[sp,cp]])


    print(filename,expid,seeing)
    if seeing>2. : continue

    t=Table.read(filename)
    #ok=t["VALID"]>0
    ok=np.abs(t["RCALIBFRAC"]-1)<0.2

    ff.append(t["FIBER"][ok])
    zz.append(t["RCALIBFRAC"][ok])

    x=t["X"][ok]
    y=t["Y"][ok]
    #xy=np.array([x,y])
    #xy=M.dot(xy)
    #x=xy[0]
    #y=xy[1]
    xx.append(x)
    yy.append(y)

    #adc1phi=et["ADC1PHI"][index]
    #plot_dflux_seeing.py:    adc2phi=et["ADC2PHI"][index]


x=np.hstack(xx)
y=np.hstack(yy)
f=np.hstack(ff)
z=np.hstack(zz)

fibermap=Table()
fibermap["FIBER_X"]=x
fibermap["FIBER_Y"]=y
corr=flat_to_psf_flux_correction(fibermap)
z *= corr

if 0 :
    plt.figure("flux-ratio-vs-fiber-"+args.what)
    plt.subplot(111,title=args.what+" flux ratio (spectro. / photom.) for std. stars (after flatfield)")
    plt.plot(f,z,".")
    plt.grid()
    for s in range(1,10) :
        plt.axvline(500*s,color="gray")
    plt.xlabel("fiber")
    plt.ylabel("normalized flux ratio")

if 0 :
    plt.figure("flux-ratio-vs-rfocal-"+args.what)
    r=np.sqrt(x**2+y**2)
    plt.subplot(111,title=args.what+" flux ratio (spectro. / photom.) for std. stars (after flatfield)")
    plt.plot(r,z,".")
    #plt.plot(r,corr,".")
    plt.grid()
    plt.xlabel("r_focal (mm)")
    plt.ylabel("normalized flux ratio")


vmin=0.8+0.0001
vmax=1.2-0.0001

if 0 :
    plt.figure("flux-ratio-vs-xy-"+args.what)
    plt.subplot(111,title=args.what+" flux ratio (spectro. / photom.)\n for std. stars (after flatfield)")
    plt.scatter(x,y,c=z,vmin=vmin,vmax=vmax)
    plt.colorbar()
    plt.xlabel("x_focal (mm)")
    plt.ylabel("y_focal (mm)")


figsize=(4,4)
if not args.no_colorbar:
    figsize=(4.8,4)
plt.figure("flux-ratio-vs-binned-xy-"+args.what,figsize=figsize)
bins1d=np.linspace(-410,410,50)

h1,bx,by = np.histogram2d(x,y,bins=(bins1d,bins1d))
hx,bx,by = np.histogram2d(x,y,weights=x,bins=(bins1d,bins1d))
hy,bx,by = np.histogram2d(x,y,weights=y,bins=(bins1d,bins1d))
hz,bx,by = np.histogram2d(x,y,weights=z,bins=(bins1d,bins1d))

mx=hx/h1
my=hy/h1
mz=hz/h1

a=plt.subplot(111)#,title=args.what) #+" flux ratio (spectro. / photom.)\n for std. stars (after flatfield)")
toto=a.imshow(mz.T,origin="lower",extent=(bins1d[0],bins1d[-1],bins1d[0],bins1d[-1]),vmin=vmin,vmax=vmax,cmap="gray")
a.axis("off")
if not args.no_colorbar:
    #plt.figure("flux-ratio-vs-binned-xy-colorbar",figsize=(2,4))
    cb=plt.colorbar(toto) #ax=a)
    from matplotlib import ticker
    tick_locator = ticker.MaxNLocator(nbins=4)
    cb.locator = tick_locator
    cb.update_ticks()
    cb.ax.tick_params(labelsize=17)
#a.set_xlabel("x_focal (mm)")
#a.set_ylabel("y_focal (mm)")
tune_axes(a)
plt.tight_layout()

t=Table()
mxb=mx
myb=my
mx=(bx[1:]+bx[:-1])/2.
my=(by[1:]+by[:-1])/2.
t["X"]=np.tile(mx,(mx.size,1)).T.ravel() # center of bin
t["Y"]=np.tile(my,(my.size,1)).ravel()
#t["XB"]=mxb.ravel() # barycenter of bin, for debugging
#t["YB"]=myb.ravel()
mz[np.isnan(mz)]=0
t["R"]=mz.ravel()
ofilename="flux-ratio-vs-xy-"+args.what+".csv"
t.write(ofilename,overwrite=True)
print("wrote",ofilename)


plt.show()

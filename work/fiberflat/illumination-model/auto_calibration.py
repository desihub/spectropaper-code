#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import os
import matplotlib as mpl
import sys

def fit_intensity(x,y,flux) :
    h=np.zeros((4,flux.size))
    h[0]+=1.
    h[1]+=x
    h[2]+=y
    h[3]+=(x**2+y**2)
    A=h.dot(h.T)
    B=np.sum(flux*h,axis=1)
    Ai=np.linalg.inv(A)
    c=Ai.dot(B)
    return c[0] # value at center of fov

def myplot(what,title) :
    fig=plt.figure(title)
    vmin=np.min(what)
    vmax=np.max(what)
    axcb  = fig.add_axes([0.85, 0.05, 0.05, 0.9])
    aximg = fig.add_axes([0.05, 0.05, 0.75, 0.9])
    aximg.scatter(fiber_x,fiber_y,c=(what-vmin)/(vmax-vmin),edgecolors="face",cmap=cmap,s=5)
    norm=mpl.colors.Normalize(vmin=vmin,vmax=vmax)
    cb=mpl.colorbar.ColorbarBase(axcb, cmap=cmap,norm=norm,orientation='vertical')


def myplot4(what,title) :
    fig=plt.figure(title)
    for l in range(4) :
        vmin=np.min(what[l])
        vmax=np.max(what[l])
        offset=0.5*np.array([l%2,(3-l)//2,0,0])
        axcb  = fig.add_axes(offset+0.45*np.array([0.85, 0.05, 0.05, 0.9]))
        aximg = fig.add_axes(offset+0.45*np.array([0.05, 0.05, 0.75, 0.9]))
        aximg.set_title("lamp %d"%l)
        aximg.axis("off")
        aximg.scatter(fiber_x,fiber_y,c=(what[l]-vmin)/(vmax-vmin),edgecolors="face",cmap=cmap,s=5)
        norm=mpl.colors.Normalize(vmin=vmin,vmax=vmax)
        cb=mpl.colorbar.ColorbarBase(axcb, cmap=cmap,norm=norm,orientation='vertical')

plot=False
fiber_flux={}

for source in ["sky","lamp_0","lamp_1","lamp_2","lamp_3"] :
    filename = "fiber_illumination_fraction_%s.txt"%source
    print("reading %s"%filename)
    
    tmp=np.loadtxt(filename).T
    #tmp=tmp[:,:100] # for debugging
    fiber_id=tmp[0].astype(int)
    fiber_x=tmp[1]
    fiber_y=tmp[2]
    fiber_theta=tmp[3]
    fiber_phi=tmp[4]
    fiber_flux[source]=tmp[5]
    fiber_flux[source] /= np.max(fiber_flux[source])

r=np.sqrt(fiber_x**2+fiber_y**2)
   
nfibers=fiber_id.size
fiber_spectro=fiber_id//500
spectros=np.unique(fiber_spectro)

nlamps=4

flux=np.zeros((nlamps,nfibers))
for l in range(nlamps) :
    flux[l]=fiber_flux["lamp_%d"%l]


true_fiber_transmission = np.random.normal(size=nfibers)*0.3+1
true_fiber_transmission[true_fiber_transmission<0.5]=0.2
true_fiber_transmission[true_fiber_transmission>2]=2.
# fiber transmission has an average of one per spectro by definition
for spec in spectros :
    true_fiber_transmission[fiber_spectro==spec] /= np.mean(true_fiber_transmission[fiber_spectro==spec])

true_spectro_transmission = np.random.normal(size=spectros.size)*0.2+1
true_spectro_transmission[true_spectro_transmission<0.5]=0.5
true_spectro_transmission[true_spectro_transmission>2]=2.
print("true_spectro_transmission=",true_spectro_transmission)

true_lamp_intensity = np.exp(np.random.normal(size=nlamps)*0.9)
true_lamp_intensity /= np.sum(true_lamp_intensity)
print("true_lamp_intensity=",true_lamp_intensity)

# measured flux affected by spectro transmission , fiber transmission , fluctuation of lamp intensity
meas_flux=flux.copy()
for lamp in range(nlamps) :
    meas_flux[lamp] *= true_lamp_intensity[lamp]*true_fiber_transmission
for spec in spectros :
    meas_flux[:,fiber_spectro==spec] *= true_spectro_transmission[spec]

meas_flux /= np.mean(meas_flux)

if plot :
    cmap=plt.cm.RdYlBu
    myplot4(flux,"lamps")
    myplot4(meas_flux,"data")


    
# fit intensities using polynomial fit
estimated_lamp_intensity=np.zeros(nlamps)
for lamp in range(nlamps) :
    estimated_lamp_intensity[lamp] = fit_intensity(fiber_x,fiber_y,meas_flux[lamp])

estimated_lamp_intensity /= np.sum(estimated_lamp_intensity)
rms_intensity_error=np.std(estimated_lamp_intensity-true_lamp_intensity)
max_intensity_error=np.max(np.abs(estimated_lamp_intensity-true_lamp_intensity))
print("meas lamp intensity= %s max err=%f"%(estimated_lamp_intensity,max_intensity_error))


sum_flux = np.sum((1./estimated_lamp_intensity)*meas_flux.T,axis=1)
meas_fiber_transmission=sum_flux/fiber_flux["sky"]
for spec in spectros :
    meas_fiber_transmission[fiber_spectro==spec] /= np.mean(meas_fiber_transmission[fiber_spectro==spec])
residuals=meas_fiber_transmission-true_fiber_transmission

if plot :
    myplot(meas_fiber_transmission,title="meas. fiber transmission")
    myplot(residuals,title="residuals")

for spectro in spectros :
    ok=np.where(fiber_spectro==spectro)[0]
    skyfibers = (np.random.uniform(size=ok.size)*ok.size).astype(int)
    c=np.polyfit(r[ok][skyfibers],residuals[ok][skyfibers],2)
    pol=np.poly1d(c)
    residuals[ok] -= pol(r[ok])
print("RESULT fiberflat rms,max= %f %f intensity max,rms= %f %f"%(np.std(residuals),np.max(np.abs(residuals)),rms_intensity_error,max_intensity_error))

if plot :
    myplot(residuals,title="residuals-bis")
    plt.show()

    

#!/usr/bin/env python

import matplotlib.pyplot as plt
from matplotlib import ticker
import numpy as np
import os
import matplotlib as mpl

lamps=[0]
ny=2
nx=2
index=1



cmap=plt.cm.gray

fiber_flux={}

for source in ["lamp_0","lamp_1","lamp_2","lamp_3","sky"] :
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


#fiber_flux["lamps0"]=(fiber_flux["lamp_0"])/fiber_flux["sky"]
#fiber_flux["(lamp 1)/sky"]=(fiber_flux["lamp_1"])/fiber_flux["sky"]
#fiber_flux["(lamp 2)/sky"]=(fiber_flux["lamp_2"])/fiber_flux["sky"]
#fiber_flux["(lamps 0+1)/sky"]=(fiber_flux["lamp_0"]+fiber_flux["lamp_1"])/fiber_flux["sky"]
#fiber_flux["(lamps 0+1+2)/sky"]=(fiber_flux["lamp_0"]+fiber_flux["lamp_1"]+fiber_flux["lamp_2"])/fiber_flux["sky"]
fiber_flux["lamps0123"]=(fiber_flux["lamp_0"]+fiber_flux["lamp_1"]+fiber_flux["lamp_2"]+fiber_flux["lamp_3"])/fiber_flux["sky"]
fiber_flux["lamps0123"] /= np.max(fiber_flux["lamps0123"])

for source in ["lamp_0","lamps0123"] :
#for source in ["(lamp 0)/sky","(lamp 1)/sky","(lamp 2)/sky"] :

    print("plotting %s"%source)



    fiber_flux[source] /= np.max(fiber_flux[source])

    min_flux=np.min(fiber_flux[source])
    max_flux=np.max(fiber_flux[source])
    print(min_flux,max_flux)

    if False : # old version
        fig=plt.figure(source,figsize=(10,8))
        axcb  = fig.add_axes([0.85, 0.05, 0.05, 0.9])
        aximg = fig.add_axes([0.05, 0.05, 0.75, 0.9])
        frac = (fiber_flux[source]-min_flux)/(max_flux-min_flux)
        aximg.scatter(fiber_x,fiber_y,c=frac,edgecolors="face",cmap=cmap)

        norm=mpl.colors.Normalize(vmin=min_flux,vmax=max_flux)
        cb=mpl.colorbar.ColorbarBase(axcb, cmap=cmap,norm=norm,orientation='vertical')

    if True : # new version
        fig=plt.figure(source,figsize=(4.9,4))
        axcb  = fig.add_axes([0.8, 0.05, 0.05, 0.9])
        aximg = fig.add_axes([0.05, 0.05, 0.75, 0.9])
        #frac = (fiber_flux[source]-min_flux)/(max_flux-min_flux)

        if source=="lamp_0" :
            max_flux=1
            min_flux=0.75
        else :
            max_flux=1
            min_flux=0.992
        aximg.scatter(fiber_x,fiber_y,c=fiber_flux[source],s=7,edgecolors="face",cmap=cmap,vmin=min_flux,vmax=max_flux)
        aximg.axis("off")
        #aximg.xaxis.set_ticklabels([])
        #aximg.yaxis.set_ticklabels([])
        #aximg.xaxis.set_ticks([])
        #aximg.yaxis.set_ticks([])

        norm=mpl.colors.Normalize(vmin=min_flux,vmax=max_flux)
        cb=mpl.colorbar.ColorbarBase(axcb,norm=norm,orientation='vertical',cmap=cmap)
        tick_locator = ticker.MaxNLocator(nbins=5)
        cb.locator = tick_locator
        cb.update_ticks()
        cb.ax.tick_params(labelsize=14)
        #plt.tight_layout()
        ofilename=source+".pdf"
        fig.savefig(ofilename,format="pdf")

    if source=="lamp0" : continue

    fig=plt.figure("%s-theta"%source)
    plt.plot(fiber_theta,fiber_flux[source],"o")
    c=np.polyfit(fiber_theta,fiber_flux[source],2)
    print("coefficients=",c)
    pol=np.poly1d(c)
    ii=np.argsort(fiber_theta)
    plt.plot(fiber_theta[ii],pol(fiber_theta[ii]),"-",color="r")

    plt.xlabel("theta (deg)")
    plt.ylabel("inhomogeneity (w.r.t. sky)")
    plt.grid()

    fig=plt.figure("%s-phi"%source)
    thetaref=1.
    ok=np.where(np.abs(fiber_theta-thetaref)<3)[0]
    #plt.plot(fiber_phi[ok],fiber_flux[source][ok]/pol(fiber_theta[ok]),"o")
    h0,bins=np.histogram(fiber_phi[ok],bins=40)
    hphi,junk=np.histogram(fiber_phi[ok],bins=bins,weights=fiber_phi[ok])
    hy,junk=np.histogram(fiber_phi[ok],bins=bins,weights=fiber_flux[source][ok]/pol(fiber_theta[ok]))
    hphi /= h0
    hy /= h0
    plt.plot(hphi,hy,"o")
    plt.ylim([0.99,1.01])
    plt.xlabel("phi (for theta ~ %2.1f deg)"%thetaref)
    plt.ylabel("inhomogeneity (w.r.t. sky)")
    plt.grid()

plt.show()

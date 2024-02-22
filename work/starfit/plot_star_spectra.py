#!/usr/bin/env python

import sys,os
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
from desispec.io import read_frame
from desispec.interpolation import resample_flux
from astropy.table import Table

# /project/projectdirs/desi/spectro/redux/andes/exposures/20200315/00055643/stdstars-0-00055643.fits

filename=sys.argv[1]

stdstars=pyfits.open(filename)
stdstars.info()

fibers=stdstars["FIBERS"].data
swave=stdstars["WAVELENGTH"].data
sflux=stdstars["FLUX"].data
sdata=stdstars["METADATA"].data

frames=[]
for band in "brz" :
    frame_filename = filename.replace("stdstars-","cframe-"+band)
    if os.path.isfile(frame_filename) :
        frame=read_frame(frame_filename)
        frames.append(frame)

colors=["C0","red","brown"]


zenodo=False

if zenodo :
    zenodo_data_dict = dict()
    zenodo_data_dict["WAVELENGTH"]=[]
    zenodo_data_dict["MEASURED_FLUX"]=[]
    zenodo_data_dict["CAMERA"]=[]
    zenodo_model_dict = dict()
    zenodo_model_dict["WAVELENGTH"]=[]
    zenodo_model_dict["MODEL_FLUX"]=[]
    zenodo_model_dict["CAMERA"]=[]

for starindex,fiber in enumerate(fibers) :

    fig=plt.figure("star-spectrum",figsize=(5,4))

    aa=[]
    aa.append(plt.axes([0.15,0.42,0.8,0.5]))
    for i in range(4) :
        aa.append(plt.axes([0.15+(0.8/4)*i,0.05,0.8/4,0.25]))



    #a00=plt.subplot(111)
    #aa=[ plt.subplot(211), plt.subplot(245), plt.subplot(246), plt.subplot(247), plt.subplot(248) ]

    #gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1])gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1])
    selections=[[1000.,10000.],[1000,3800],[4050,4400],[6550,6600],[8450,8750]]

    first=True
    for a,selection in zip (aa,selections) :
        mx=[]
        my=[]
        for frame,color in zip(frames,colors) :
            ii=np.where((frame.ivar[fiber]>0)&(frame.mask[fiber]==0)&(frame.wave>=selection[0])&(frame.wave<=selection[1]))[0]
            if ii.size < 10 : continue
            err=1./np.sqrt(frame.ivar[fiber,ii])
            a.plot(frame.wave[ii],frame.flux[fiber,ii],"-",color=color,alpha=0.8)
            #plt.errorbar(frame.wave[ii],frame.flux[fiber,ii],err,fmt=".-",color=color)
            if zenodo :
                zenodo_data_dict["WAVELENGTH"].append(frame.wave[ii])
                zenodo_data_dict["MEASURED_FLUX"].append(frame.flux[fiber,ii])
                zenodo_data_dict["CAMERA"].append(np.repeat(frame.meta["CAMERA"][0].upper(),ii.size))

            tmp = resample_flux(frame.wave,swave,sflux[starindex])
            tmp = frame.R[fiber].dot(tmp)
            norm = np.sum(tmp[ii]*frame.flux[fiber,ii])/np.sum(tmp[ii]**2)

            if first:
                blabla="fiber #{}\n".format(fiber)
                blabla+="$T_{eff} = %4.0f\,$K\n"%(sdata["TEFF"][starindex])
                blabla+="$\log g = %2.1f$\n"%(sdata["LOGG"][starindex])
                blabla+="[Fe/H] = %2.1f\n"%(sdata["FEH"][starindex])
                props = dict(boxstyle='round', facecolor='white', alpha=0.5)
                plt.text(0.98,0.96,blabla,fontsize=10, bbox=props,verticalalignment='top', horizontalalignment='right', transform=a.transAxes)
                a.set_xlabel("Wavelength ($\AA$)")
                a.set_ylabel("Flux ($10^{-17}$ erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)")
            else:
                #blabla=r"model$\,\times {:3.2f}$".format(norm)
                #props = dict(boxstyle='round', facecolor='white', alpha=0.5)
                #plt.text(0.98,0.0,blabla,fontsize=10, verticalalignment='bottom', horizontalalignment='right', transform=a.transAxes)
                a.get_yaxis().set_visible(False)
                #a.get_xaxis().set_visible(False)
                a.get_xaxis().set_ticks([])
                mmwave=int(np.mean(frame.wave[ii])/10)*10
                a.set_xlabel(f"${mmwave}$ $\AA$")

                #a.xaxis.set_major_locator(plt.MaxNLocator(1))

            ii=(frame.wave>=selection[0])&(frame.wave<=selection[1])
            mx.append(frame.wave[ii])
            my.append(tmp[ii]*norm)
            if zenodo :
                zenodo_model_dict["WAVELENGTH"].append(frame.wave[ii])
                zenodo_model_dict["MODEL_FLUX"].append(tmp[ii]*norm)
                zenodo_model_dict["CAMERA"].append(np.repeat(frame.meta["CAMERA"][0].upper(),np.sum(ii)))


        if zenodo :
            t=Table()
            for k in zenodo_data_dict :
                t[k]=np.hstack(zenodo_data_dict[k])
            t.write("../../zenodo/figure-36-data.fits")
            t=Table()
            for k in zenodo_model_dict :
                t[k]=np.hstack(zenodo_model_dict[k])
            t.write("../../zenodo/figure-36-model.fits")
            print("exit for zenodo")
            sys.exit(0)

        #mx=np.hstack(mx)
        #my=np.hstack(my)
        for mmx,mmy in zip(mx,my) :
            a.plot(mmx,mmy,"-",c="k",alpha=0.5)

        first=False


    #plt.box()
    #plt.xlabel("Wavelength ($\AA$)")
    #plt.ylabel("Flux ($10^{-17}$ ergs s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)")
    #plt.tight_layout()
    plt.show()

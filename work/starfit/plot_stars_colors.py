#!/usr/bin/env python

import os
import sys
import numpy as np

import fitsio
import matplotlib.pyplot as plt

filenames=sys.argv[1:]
fileid=[]
snr=[]
chi2=[]
data_gr=[]
model_gr=[]
logg=[]
teff=[]
feh=[]
xmin=12.
xmax=-12.
for fid,filename in enumerate(filenames) :

    t=fitsio.read(filename,"METADATA")
    print(t.dtype.names)
    #plt.plot(t["DATA_G-R"],t["MODEL_G-R"],"o",label='{}'.format(os.path.basename(filename)))
    data_gr.append(t["DATA_G-R"])
    model_gr.append(t["MODEL_G-R"])
    fileid.append(fid*np.ones(len(t)))
    logg.append(t["LOGG"])
    teff.append(t["TEFF"])
    feh.append(t["FEH"])
    if "BLUE_SNR" in t.dtype.names :
        snr.append(t["BLUE_SNR"])
    chi2.append(t["CHI2DOF"])

data_gr=np.hstack(data_gr)
model_gr=np.hstack(model_gr)
fileid=np.hstack(fileid)
logg=np.hstack(logg)
teff=np.hstack(teff)
feh=np.hstack(feh)

xmin=np.min(data_gr)
xmax=np.max(data_gr)
snr=np.hstack(snr)

print("Std of color difference                   = {:4.3f}".format(np.sqrt(np.mean((data_gr[model_gr!=0]-model_gr[model_gr!=0])**2))))


selection=(data_gr>-100)&(data_gr<0.35)&(model_gr!=0)&(snr>0)

plt.figure("stdstar-color-color",(5,3))
a=plt.subplot(111)
a.plot([xmin,xmax],[xmin,xmax],c="gray")

cmap = plt.get_cmap("tab10")
colors = cmap(np.linspace(0,1,len(filenames)))
color="C0"
for fid,filename in enumerate(filenames) :
    label=os.path.basename(filename)
    ms=6
    #plt.plot(data_gr[(fileid==fid)],model_gr[(fileid==fid)],"o",c=colors[fid],markersize=ms,alpha=0.6,fillstyle='none',label=None)
    a.plot(data_gr[selection&(fileid==fid)],model_gr[selection&(fileid==fid)],"o",c=color,markersize=ms,alpha=0.6,label=label)

print("mean model - data (G-R) = ",np.mean(model_gr[selection]-data_gr[selection]))
print("rms  model - data (G-R) = ",np.std(model_gr[selection]-data_gr[selection]))
print("number of stars = ",np.sum(selection))

a.set_xlabel("Data $g-r$")
a.set_ylabel("Model $g-r$")
blabla = "exposure #00055643 (900 sec)\n"
blabla += "r.m.s. $\Delta (g-r) = %4.3f$"%(np.std(model_gr[selection]-data_gr[selection]))
props = dict(boxstyle='round', facecolor='white', alpha=0.5)
plt.text(0.05,0.96,blabla,fontsize=10, bbox=props,verticalalignment='top', horizontalalignment='left', transform=a.transAxes)
a.grid()
plt.tight_layout()
plt.show()
#a.legend()

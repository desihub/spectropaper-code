#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

from desispec.fibercrosstalk import read_crosstalk_parameters,eval_crosstalk



params = read_crosstalk_parameters()


wave=np.linspace(8000,9800,500)
fibers=np.array([250])
dfiber=1

waves={}
waves["b"]=np.linspace(3600.0,5800.0,100)
waves["r"]=np.linspace(5760.0,7620.0,100)
waves["z"]=np.linspace(7520.0,9824.0,100)


plt.figure("fibercrosstalk",figsize=(5,3))

for p,dfiber in enumerate([1,2]) :
    a=plt.subplot(2,1,p+1)
    for s in range(10) :
        color=f"C{s}"
        for b in ["b","r","z"] :
            camera=f"{b}{s}"
            wave=waves[b]
            xtalk = eval_crosstalk(camera,wave,fibers,dfiber,params,apply_scale=True,nfiber_per_bundle=25)
            a.plot(wave,xtalk[0],"-",color=color)
            #if b=="z" : plt.plot(wave,xtalk[1],"--",color=color)
    plt.grid()
    plt.text(0.1,0.85,f"cross-talk between fibers N and N+{dfiber}",transform=a.transAxes,bbox=dict(facecolor='white', alpha=1))
    #a.set_ylabel(f"Cross-talk at +{dfiber} fiber")

a.set_xlabel("Wavelength ($\AA$)")

plt.tight_layout()
plt.show()

#!/usr/bin/env python


import numpy as np
import matplotlib.pyplot as plt
import fitsio
from specter.psf.gausshermite import GaussHermitePSF

hw=21

#cam="b2" ; title="SM3 "+cam[0] ; expids=np.arange(10381,10390,dtype=int)
#cam="r2" ; title="SM3 "+cam[0] ; expids=np.arange(10381,10390,dtype=int)
#cam="z2" ; title="SM3 "+cam[0] ; expids=np.arange(10381,10390,dtype=int)
cam="z1" ; title="SM1 "+cam[0] ; expids=np.arange(7859,7868,dtype=int)
#cam="z1" ; title="SM1 "+cam[0] ; expids=np.arange(7859,7860,dtype=int)





psf_filename="/global/cfs/cdirs/desi/users/jguy/fiber_crosstalk/psf-{}.fits".format(cam)
psf=GaussHermitePSF(psf_filename)

psf._polyparams['HSIZEX']=hw
for k in ['TAILXSCA','TAILYSCA','TAILAMP','TAILCORE','TAILINDE'] :
    print("{} = {}".format(k,psf.coeff[k].eval(0,9000.)))

waves=None
if cam[0]=="b" :
    waves = np.array([4047.5,4801.3,5087.2,5462.2])
if cam[0]=="r" :
    waves = np.array([5854,6600.7,6680.,6931.3,7441.])
if cam[0]=="z" :
    waves = np.array([7603.6384,8192.3082,8300.3907,9125.4719,9660.435,9660.4,9754.435])

xx=np.linspace(-hw,hw,hw*2*2+1)
wprof={} # data
wmprof={} # model

images=[]

for e in expids :
    print(e)
    images.append(fitsio.read("/global/cfs/cdirs/desi/users/jguy/fiber_crosstalk/preproc-{}-{:08d}.fits".format(cam,e)))

for wave in waves :
    profs=[]
    mprofs=[]
    for fiber in range(psf.nspec) :
        x,y = psf.xy(fiber,wave)
        print(wave,fiber,x,y)
        b = int(x-hw)
        e = int(x+hw+2)
        ii = np.arange(b,e)
        j  = int(y+0.5)
        for img in images :
            prof = np.interp(xx,ii-x,img[j,b:e])
            profs.append(prof)

        tx,ty,pix = psf.xypix(fiber,wave)
        mprof = np.interp(xx,np.arange(tx.start,tx.stop)-x,pix[j-ty.start,:])
        mprof *= ( np.sum(prof)/np.sum(mprof) )
        mprofs.append(mprof)

    profs=np.array(profs)
    prof=np.sum(profs,axis=0)
    prof/=np.max(prof)
    wprof[int(wave)]=prof
    mprofs=np.array(mprofs)
    mprof=np.sum(mprofs,axis=0)
    mprof/=np.max(mprof)
    wmprof[int(wave)]=mprof


plt.figure("psf-cross-profile-"+title.lower().replace(" ","-"))
plt.subplot(111,title=title)

colors = plt.cm.jet(np.linspace(0,1,len(waves)))
y=np.zeros(waves.shape)
my=np.zeros(waves.shape)
for i,wave in enumerate(waves) :
    prof=wprof[int(wave)]
    plt.plot(xx,prof,color=colors[i],label="{:4.1f}A".format(wave))
    mprof=wmprof[int(wave)]
    #plt.plot(xx,mprof,"--",color=colors[i],alpha=0.8)#,color=colors[i],label="{:4.1f}A".format(wave))
    y[i]=np.interp(7.5,xx,prof)
    my[i]=np.interp(7.5,xx,mprof)

plt.ylim([5e-6,1.5])
plt.grid()
plt.legend(loc="upper right")
plt.yscale('log')
plt.xlabel("pixels")
plt.tight_layout()
if 0 :
    plt.figure("psf-tails-vs-wave-"+title.lower().replace(" ","-"))
    plt.subplot(111,title=title)
    plt.plot(waves,y,"-o")
    plt.plot(waves,my,"--")
    plt.xlabel("wavelength (A)")
    plt.ylabel("flux frac")
    plt.grid()

plt.show()

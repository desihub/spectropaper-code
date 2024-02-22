#!/usr/bin/env python

import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt
import argparse
import os.path
from desispec.util import parse_fibers

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i', '--infile', type = str, default = None, required = True, nargs="*",
                    help = 'path to fiber flat file')
parser.add_argument('--spectra', action = 'store_true')
parser.add_argument('--fibers', type=str, default = None, required = False,
                    help = 'defines from_to which fiber to work on. (ex: --fibers=50:60,4 means that only fibers 4, and fibers from 50 to 60 (excluded) will be plotted)')

args        = parser.parse_args()

x=[]
y=[]
z=[]

if args.fibers is not None :
    fibers = parse_fibers(args.fibers)
else :
    fibers = None

ntot=0
n1=0
for filename in args.infile :

    h=fits.open(filename)
    h.info()
    wave=h["WAVELENGTH"].data

    ntot += h[0].data.size
    n1   += np.sum(h[0].data==1.)

    if args.spectra :
        figname=os.path.basename(filename).split(".")[0]
        fig = plt.figure(figname,figsize=(5,2.5))
        if fibers is None :
            fibers = np.arange(h[0].data.shape[0],dtype=int)
        for fiber in fibers :
            ok=np.where(h["MASK"].data[fiber]==0)[0]
            tx=wave[ok]
            ty=h[0].data[fiber,ok]
            #plt.plot(tx,ty)
            # plot one point out of N to make the figure lighter
            i=np.arange(tx.size)
            ii=(i%2==0)
            plt.plot(tx[ii],ty[ii])

        plt.grid()
        plt.ylim([0.6,1.3])
        plt.xlabel(r"Wavelength ($\AA$)")
        plt.ylabel("Relative transmission")
        plt.tight_layout()

    if "FIBERMAP" in h :
        h[0].data *= (h["MASK"].data == 0 )
        z.append( np.median(h[0].data,axis=1) )
        x.append( h["FIBERMAP"].data["FIBERASSIGN_X"] )
        y.append( h["FIBERMAP"].data["FIBERASSIGN_Y"] )

print("fraction of pixels with flat==1 :",float(n1)/ntot)

if False and len(x)>0 :
    x=np.hstack(x)
    y=np.hstack(y)
    z=np.hstack(z)
    vmin=0.7
    vmax=1.3
    plt.figure("focalplane")
    plt.scatter(x,y,c=z,vmin=vmin,vmax=vmax)
    plt.colorbar()
    plt.figure()
    plt.hist(z[(z>=vmin)&(z<=vmax)],bins=50)
    plt.tight_layout()
    print("median          =",np.median(z))
    print("median in range =",np.median(z[(z>=vmin)&(z<=vmax)]))
    print("mean            =",np.mean(z))
    print("mean in range   =",np.mean(z[(z>=vmin)&(z<=vmax)]))


plt.tight_layout(rect=[0.0,-0.06,1.02,1.03])
ofilename=figname+".pdf"
fig.savefig(ofilename,format="pdf") ; print("wrote",ofilename)
plt.show()

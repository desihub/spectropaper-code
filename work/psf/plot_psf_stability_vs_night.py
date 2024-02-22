#!/usr/bin/env python

import numpy as np
from  astropy.table import Table
import matplotlib.pyplot as plt
import sys
import argparse
import os.path

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i','--infile', type = str, default = None, required = True,
                    help = 'path to csv table')
args        = parser.parse_args()

table = Table.read(args.infile)

nights=table["night"]
waves=np.unique(table["wavelength"])
waves=waves[2:-2]
print(waves)
fibers=np.unique(table["fiber"])
print(fibers)

fig1=plt.figure(num="psf-stability",figsize=(13,4))

ny=2
nx=2
a_x = plt.subplot(ny,nx,1)
a_y = plt.subplot(ny,nx,2)
a_xsig = plt.subplot(ny,nx,3)
a_ysig = plt.subplot(ny,nx,4)

fig2=plt.figure(num="psf-stability-emm",figsize=(6,3))
a_emm = plt.subplot(111)
lstyle="."
color="k"
alpha=0.5

for wave in waves :
#for wave in [ waves[waves.size//2] ] :
    selection = (table["wavelength"]==wave)
    for fiber in fibers :
        selection &= (table["fiber"]==fiber)
        if np.sum(selection)==0 : continue
        mx = np.mean(table["x"][selection])
        my = np.mean(table["y"][selection])
        a_x.plot(nights[selection],table["x"][selection]-mx,lstyle,color=color,alpha=alpha)
        a_y.plot(nights[selection],table["y"][selection]-my,lstyle,color=color,alpha=alpha)
        mxsig=np.mean(table["xsig"][selection])
        mysig=np.mean(table["ysig"][selection])
        a_xsig.plot(nights[selection],table["xsig"][selection]/mxsig-1.,lstyle,color=color,alpha=alpha)
        a_ysig.plot(nights[selection],table["ysig"][selection]/mysig-1.,lstyle,color=color,alpha=alpha)
        memm=np.mean(table["emission_line_coeff"][selection])
        a_emm.plot(nights[selection],table["emission_line_coeff"][selection]-memm,lstyle,color=color,alpha=alpha)

a_x.set_ylabel("$\delta X$ (pixels)")
a_y.set_ylabel("$\delta Y$ (pixels)")
#a_xsig.set_ylabel("$\delta \, \ln (\sigma_X) $")
#a_ysig.set_ylabel("$\delta \, \ln (\sigma_Y) $")
a_xsig.set_ylabel("$\delta \, \sigma_X / \sigma_X $")
a_ysig.set_ylabel("$\delta \, \sigma_Y / \sigma_Y $")
a_x.grid()
a_y.grid()
a_xsig.grid()
a_ysig.grid()
a_xsig.set_xlabel("nights")
a_ysig.set_xlabel("nights")


a_emm.grid()
a_emm.set_xlabel("nights")
a_emm.set_ylabel(r"$\delta \, F / F$")


plt.tight_layout()

filename=os.path.basename(args.infile).split(".")[0]+".pdf"
fig1.savefig(filename)
print("wrote",filename)
filename=os.path.basename(args.infile).split(".")[0]+"-emm.pdf"
fig2.savefig(filename)
print("wrote",filename)

plt.show()

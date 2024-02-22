#!/usr/bin/env python

import sys,os
import argparse
import matplotlib.pyplot as plt
import numpy as np
import fitsio

#from pkg_resources import resource_exists, resource_filename

#from desispec.util import parse_fibers
from desispec.qproc.io import read_qframe,write_qframe
from desimeter.io import load_metrology


#from desispec.io import read_fibermap,read_frame
#from desispec.interpolation import resample_flux
#from desispec.fluxcalibration import isStdStar

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i','--infile', type = str, default = None, required = True, nargs="*",
                    help = 'path to one or several frame fits files')
#parser.add_argument('-o','--outfile', type = str, default = None, required = True,help = 'outfile')

metro = load_metrology()

moving_posids = ["M02780", "M03168", "M03353", "M03848", "M04015", "M04062",
                 "M04390", "M04527", "M04653", "M04662", "M04663", "M04666",
                 "M04674", "M04684", "M04760", "M04787", "M05225", "M05252",
                 "M05528", "M05692", "M05701", "M05704", "M05725", "M05726",
                 "M05781", "M05811", "M05814", "M05816", "M05843", "M05845",
                 "M05862", "M05889", "M05897", "M05930", "M05977", "M05980",
                 "M05983", "M06080", "M06157", "M06160", "M06163", "M06167",
                 "M06171", "M06609", "M06666", "M06711", "M06713", "M06734",
                 "M06763", "M06765", "M06767", "M06769", "M06777", "M07017",
                 "M07056", "M07269", "M07337", "M07343", "M07345", "M07477",
                 "M07509", "M07808", "M07840", "M07892", "M07907", "M07916",
                 "M07918", "M07935", "M07942", "M07943", "M07946", "M07957",
                 "M07959", "M07962", "M07963", "M07967", "M07971", "M08017",
                 "M08026", "M08054", "M08056", "M08060", "M08071", "M08096",
                 "M08135", "M08143", "M08147", "M08150", "M08151", "M08154",
                 "M08158", "M08162", "M08180", "M08182", "M08189", "M08190",
                 "M08192", "M08193"]


disabled_posids = ["M02956", "M03090", "M03848", "M04062", "M04343", "M04515",
                   "M04529", "M05218", "M05703", "M05842", "M05862", "M05871",
                   "M05888", "M05925", "M05928", "M05977", "M05980", "M06167",
                   "M06612", "M07364"]

moving_from_metrology=np.in1d(metro["DEVICE_ID"],moving_posids)
disabled_from_metrology=np.in1d(metro["DEVICE_ID"],disabled_posids)


if False :
    plt.figure("metrology")
    plt.plot(metro["X_FP"],metro["Y_FP"],".",color="gray")
    plt.plot(metro["X_FP"][moving_from_metrology],metro["Y_FP"][moving_from_metrology],"o",color="green")
    plt.plot(metro["X_FP"][disabled_from_metrology],metro["Y_FP"][disabled_from_metrology],"X",color="red")

args   = parser.parse_args()

amps = []

#plt.figure("mflux")

for filename in args.infile :
    if not os.path.isfile(filename) :
        print("missing",filename)
        continue
    print("adding",filename)
    frame = read_qframe(filename)

    print(frame.fibermap.dtype.names)

    #print(list(frame.fibermap["LOCATION"]))
    #print(list(metro["LOCATION"][moving]))

    moving=np.in1d(list(frame.fibermap["LOCATION"]),list(metro["LOCATION"][moving_from_metrology]))
    disabled=np.in1d(list(frame.fibermap["LOCATION"]),list(metro["LOCATION"][disabled_from_metrology]))

    if False :
        plt.figure("frame")
        x=frame.fibermap["FIBERASSIGN_X"]
        y=frame.fibermap["FIBERASSIGN_Y"]
        plt.plot(x,y,".",color="gray")
        plt.plot(x[moving],y[moving],"o",color="green")
        plt.plot(x[disabled],y[disabled],"x",color="red")
        plt.show()

    mwave = np.mean(frame.wave,axis=0)
    rflux = np.zeros(frame.flux.shape)
    for i in range(frame.flux.shape[0]) :
        rflux[i]=np.interp(mwave,frame.wave[i],frame.flux[i],left=0,right=0)
    mflux = np.median(rflux,axis=0)

    peaks = np.zeros(mflux.shape)
    peaks[1:-1] = (mflux[1:-1]>2000)&(mflux[1:-1]>mflux[:-2])&(mflux[1:-1]>mflux[2:])
    peakindices=np.where(peaks>0)[0]

    amp=[]
    for i in range(rflux.shape[0]) :
        ampi=[]
        for j in peakindices :
            a = np.sum(mflux[j-4:j+5]**2)
            b = np.sum(mflux[j-4:j+5]*rflux[i,j-4:j+5])
            ampi.append(b/a)
        amp.append(np.median(ampi))

    #plt.plot(mwave,rflux[i],alpha=0.5)
    #plt.plot(mwave,mflux)
    #plt.plot(mwave[peakindices],mflux[peakindices],"o")

    #plt.show()
    #sys.exit(12)

    #a = np.sum(mflux**2)
    #b = np.sum(rflux*mflux[None,:],axis=1)
    #amp = b/a
    amps.append(amp)

amps=np.array(amps)

mamp=np.median(amps,axis=0)
goodfibers=np.where(mamp>0.5)[0] # ignore broken fibers
amps=amps[:,goodfibers]
mamp=mamp[goodfibers]
amps/=mamp[None,:]

movingfibers=np.where(moving[goodfibers])[0]
disabledfibers=np.where(disabled[goodfibers])[0]

plt.figure("amp")
for amp in amps :
    plt.plot(amp[movingfibers],".",color="green")
    plt.plot(amp[disabledfibers],".",color="gray")


plt.figure("rms")
rmsamp=np.std(amps,axis=0)
maxamp=np.max(amps,axis=0)
minamp=np.min(amps,axis=0)
plt.plot(rmsamp[movingfibers])
plt.plot(maxamp[movingfibers]-1)
plt.plot(minamp[movingfibers]-1)

print("mean rms(moving)=",np.mean(rmsamp[movingfibers]))
print("mean rms(disabled)=",np.mean(rmsamp[disabledfibers]))
print("rms from moves = ",np.sqrt(np.mean(rmsamp[movingfibers]**2)-np.mean(rmsamp[disabledfibers]**2)))


plt.grid()

plt.show()

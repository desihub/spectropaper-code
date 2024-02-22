#!/usr/bin/env python

import numpy as np
import os,sys
import matplotlib.pyplot as plt
from astropy.table import Table


mmin=58840
#mmin=0


print("read telemetry...")
#telem=Table.read("telemetry.csv")
#telem.write("telemetry.fits")
telem=Table.read("telemetry.fits")
print("done")

camera="blue" ; s=0.4
#camera="red" ; s=0.2
#camera="nir" ; s=0.05

t={}

for cam in ["b0","b1","b2","b3","b4","b5","b6","b7","b8","b9"] :
#for cam in ["r0","r1","r2","r3","r4","r5","r6","r7","r8","r9"] :
#for cam in ["z0","z1","z2","z3","z4","z5","z6","z7","z8","z9"] :

    filename="psfxy-{}.csv".format(cam)
    if os.path.isfile(filename) :
        t[cam]=Table.read(filename)
        selection=(t[cam]["MJD-OBS"]>mmin)&(t[cam]["x"]!=0)
        t[cam]=t[cam][selection]


print("averaging")
for k in ["room_pressure","space_temp1","space_temp2","space_temp3"] :
    y=telem[k]
    n=y.size
    h=15
    tmp=np.zeros(n)
    nn=np.zeros(n)
    for d in range(-h,h+1) :
        tmp[h:n-h] += y[h+d:n-h+d]
        nn[h:n-h] += (y[h+d:n-h+d]!=0.)
    tmp[nn>0] /= nn[nn>0]
    telem["mean_"+k]=tmp

#plt.figure("pressure-vs-humidity")
#ii=telem["mean_room_pressure"]>0
#plt.plot(telem["mjd"][ii]-mmin,telem["mean_room_pressure"][ii]*5000-20)
#plt.plot(telem["mjd"][ii]-mmin,telem["space_humidity"][ii])

#tab=list(t.values())[0]
#telem["total_pressure"] = telem["mean_room_pressure"]

if 0 :
    plt.figure("pressure")
    tt=Table.read("psfxy-r3.csv")
    y1=tt["PRESSURE"]
    plt.plot(tt["MJD-OBS"]-mmin,y1,"o")

    ii=telem["mean_room_pressure"]>0
    y2=telem["mean_room_pressure"][ii]
    y2=telem["pressure"]
    plt.plot(telem["mjd"]-mmin,y2)

if 0 :
    plt.figure("humidity")
    plt.plot(telem["mjd"]-mmin,telem["humidity"])
    plt.plot(telem["mjd"]-mmin,telem["space_humidity"])

if 1 : # averaging
    for unit in range(10) :
        for k0 in ["blue_camera_temp","blue_camera_humidity","red_camera_temp","red_camera_humidity","nir_camera_temp","nir_camera_humidity"] :
            k=k0+"_{}".format(unit)
            y=telem[k]
            n=y.size
            h=15
            tmp=np.zeros(n)
            nn=np.zeros(n)
            for d in range(-h,h+1) :
                tmp[h:n-h] += y[h+d:n-h+d]
                nn[h:n-h] += (y[h+d:n-h+d]!=0.)
            tmp[nn>0] /= nn[nn>0]
            # fill blanks
            ii=np.arange(n)
            tmp[nn==0] = np.interp(ii[nn==0],ii[nn>0],tmp[nn>0])
            telem["mean_"+k]=tmp

if 0 :
    plt.figure("temperatures")
    for unit in range(10) :
        plt.plot(telem["mjd"]-mmin,telem["mean_red_camera_temp_{}".format(unit)])


#############################################################################
plt.figure(camera+"-dx-vs-humidy")
legend_label=None
for k in ["space_humidity"] :
    y=telem[k]
    ii=y!=0
    m=np.mean(y[ii])
    r=np.std(y[ii])
    ry=s*(y-m)/r
    legend_label="({}-{:3.2f})/{:3.2f}".format(k,m,r/s)
    plt.plot(telem["mjd"][ii]-mmin,ry[ii],alpha=0.5)

for cam in t :
    unit=int(cam[-1])
    color="C{}".format(unit)
    x=t[cam]["x"]
    ii=(t[cam]["MJD-OBS"]-mmin<46.5)
    dx=x-np.median(x[ii])
    #label="dX "+cam
    plt.plot(t[cam]["MJD-OBS"]-mmin,dx,"o",color=color,alpha=0.8)

plt.xlabel("days")
plt.ylabel("delta X (pixels)")
plt.legend([legend_label],loc="upper left",fontsize="small")
plt.grid()
#############################################################################

#############################################################################
plt.figure(camera+"-dy-vs-humidy")
for k in ["space_humidity"] :
    y=-telem[k]
    m=np.mean(y)
    r=np.std(y)
    ry=s*(y-m)/r
    legend_label="({}-{:3.2f})/{:3.2f}".format(k,m,r/s)
    plt.plot(telem["mjd"]-mmin,ry,alpha=0.5)

for cam in t :
    unit=int(cam[-1])
    color="C{}".format(unit)
    y=t[cam]["y"]
    ii=(t[cam]["MJD-OBS"]-mmin<46.5)
    dy=y-np.median(y[ii])
    #label="dY "+cam
    plt.plot(t[cam]["MJD-OBS"]-mmin,dy,"o",color=color,alpha=0.8)
plt.xlabel("days")
plt.ylabel("dY (pixels)")
plt.legend([legend_label],loc="upper left",fontsize="small")
plt.grid()
#############################################################################

#############################################################################
fig,aa = plt.subplots(nrows=10,ncols=1,sharex=True,num=camera+"-dx-vs-temp")

for cam in t :
    unit=int(cam[-1])
    a=aa[unit]
    color="C{}".format(unit)

    x=t[cam]["x"]
    ii=(t[cam]["MJD-OBS"]-mmin<46.5)
    dx=x-np.median(x[ii])
    a.plot(t[cam]["MJD-OBS"]-mmin,dx,"o",color=color,alpha=0.8)

    temp=telem["mean_red_camera_temp_{}".format(unit)]
    y=temp-np.median(temp)
    a.plot(telem["mjd"]-mmin,y,alpha=0.8,color=color)
    a.grid()

plt.xlabel("days")
plt.ylabel("dX (pixels)")
plt.legend([legend_label],loc="upper left",fontsize="small")

#############################################################################

#############################################################################
fig,aa = plt.subplots(nrows=10,ncols=1,sharex=True,num=camera+"-dy-vs-temp")

for cam in t :
    unit=int(cam[-1])
    a=aa[unit]
    color="C{}".format(unit)

    y=t[cam]["y"]
    ii=(t[cam]["MJD-OBS"]-mmin<46.5)
    dy=y-np.median(y[ii])
    a.plot(t[cam]["MJD-OBS"]-mmin,dy,"o",color=color,alpha=0.8)

    temp=telem["mean_red_camera_temp_{}".format(unit)]
    y=temp-np.median(temp)
    a.plot(telem["mjd"]-mmin,y,alpha=0.8,color=color)
    a.grid()

plt.xlabel("days")
plt.ylabel("dY (pixels)")
plt.legend([legend_label],loc="upper left",fontsize="small")

#############################################################################

plt.show()

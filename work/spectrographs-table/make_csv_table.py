#!/usr/bin/env python


import numpy as np
import matplotlib.pyplot as plt
import sys,os
from astropy.table import Table
import fitsio

filenames=sys.argv[1:]

cameras=[]
for s in range(10) :
    for b in 'brz' :
        cameras.append(f"{b}{s}")

amps=["A","B","C","D"]
res={}
for k in ["NIGHT","EXPID","CAMERA","AMP","RDNOISE","GAIN"] :
    res[k]=[]

for filename in filenames :
    t=Table.read(filename)
    night=int(os.path.basename(filename).replace("exposure_table_","").replace(".csv",""))
    ii=np.where((t["OBSTYPE"]=="dark")&(np.abs(t["EXPTIME"]-300)<2.))[0]
    if ii.size == 0 : continue
    i=ii[0]
    #night=t["NIGHT"]
    expid=t["EXPID"][i]

    for cam in cameras :
        efilename="/global/cfs/cdirs/desi/spectro/redux/daily/preproc/{}/{:08d}/preproc-{}-{:08d}.fits.gz".format(night,expid,cam,expid)
        if not os.path.isfile(efilename) : continue

        print(night,expid,cam)

        head=fitsio.read_header(efilename)
        for amp in amps :
            k="OBSRDN"+amp
            if not k in head : continue
            rdnoise=float(head[k])
            gain=float(head["GAIN"+amp])
            res["NIGHT"].append(night)
            res["EXPID"].append(expid)
            res["CAMERA"].append(cam)
            res["AMP"].append(amp)
            res["RDNOISE"].append(rdnoise)
            res["GAIN"].append(gain)

t=Table()
for k in res.keys() :
    t[k]=np.array(res[k])

ofilename="rdnoise-desi-202206.csv"
t.write(ofilename,overwrite=True)
print("wrote",ofilename)
print(t)

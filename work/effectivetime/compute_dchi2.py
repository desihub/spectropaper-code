#!/usr/bin/env python

import os,sys
import glob
import fitsio
import numpy as np
from astropy.table import Table
import matplotlib.pyplot as plt
from desitarget.targetmask import desi_mask

ofilename="lrg-tsnr2-deltachi2-guadalupe.csv"

if True :

    print("computing",ofilename)

    filenames=[]
    t=Table.read("/global/cfs/cdirs/desi/spectro/redux/guadalupe/tiles-guadalupe.csv")
    ii=(t["FAFLAVOR"]=="maindark")
    t=t[ii]
    tiles=sorted(t["TILEID"])

    res={}
    for k in ["TILEID","EXPID","MEDIAN_TSNR2_LRG","MEDIAN_TSNR2_ELG","MEDIAN_DELTACHI2","N_LRG"] :
        res[k]=[]

    for tile in tiles :
        expids=[os.path.basename(edir) for edir in sorted(glob.glob(f"/global/cfs/cdirs/desi/spectro/redux/guadalupe/tiles/perexp/{tile}/*"))]
        print(tile,expids)
        for expid in expids :
            xx=[]
            yy=[]
            zz=[]
            for s in range(10) :
                filename=f"/global/cfs/cdirs/desi/spectro/redux/guadalupe/tiles/perexp/{tile}/{expid}/coadd-{s}-{tile}-exp{expid}.fits"
                if not os.path.isfile(filename) :
                    continue
                scores=fitsio.read(filename,"SCORES")
                #print(scores.dtype.names)
                for target in ["LRG"] :
                    scores[f"TSNR2_{target}"]=scores[f"TSNR2_{target}_B"]+scores[f"TSNR2_{target}_R"]+scores[f"TSNR2_{target}_Z"]
                filename2=filename.replace("coadd-","redrock-")
                redshifts=fitsio.read(filename2,"REDSHIFTS")
                fibermap=fitsio.read(filename2,"FIBERMAP")
                ii=((fibermap["DESI_TARGET"]&desi_mask.LRG)>0)&(redshifts["DELTACHI2"]<6000)&(fibermap["FLUX_R"]>0.1)
                ii &= (fibermap["COADD_FIBERSTATUS"]==0)
                xx += list(scores["TSNR2_LRG"][ii])
                yy += list(redshifts["DELTACHI2"][ii])
                zz += list(scores["TSNR2_ELG"][ii])
            if len(xx)<100 : continue

            res["TILEID"].append(tile)
            res["EXPID"].append(expid)
            res["MEDIAN_TSNR2_LRG"].append(np.median(xx))
            res["MEDIAN_DELTACHI2"].append(np.median(yy))
            res["MEDIAN_TSNR2_ELG"].append(np.median(zz))
            res["N_LRG"].append(len(xx))
            nentries=len(res["EXPID"])
            print(nentries)
        #if nentries>50 : break
        if nentries%50==0 :
            t=Table()
            for k in res.keys() :
                t[k]=np.array(res[k])
            t.write(ofilename,overwrite=True)
            print("wrote",ofilename)

    t=Table()
    for k in res.keys() :
        t[k]=np.array(res[k])
    t.write(ofilename,overwrite=True)
    print("wrote",ofilename)

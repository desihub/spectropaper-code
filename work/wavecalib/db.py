#!/usr/bin/env python

import numpy as np
import psycopg2
from astropy.table import Table
from astropy.time import Time

# see also https://replicator.desi.lbl.gov/TV3/app/T/index
comm=psycopg2.connect(host="db.replicator.dev-cattle.stable.spin.nersc.org",port=60042,database='desi_dev',user='desi_reader',password='reader')
cx=comm.cursor()

def dbquery(bla) :
    cx.execute(bla)
    keys=[col.name for col in cx.description]
    print(keys)
    print(" fetching data...")
    tmp=cx.fetchall()
    print(" fill table...")
    res=Table()
    for i,k in enumerate(keys) :
        res[k]=[entry[i] for entry in tmp]
    print(" done")
    return res


print("telemetry from shack")
res=dbquery("select room_pressure,space_temp1,space_temp2,space_temp3,space_temp4,space_humidity,time_recorded from shack_wec where time_recorded between date '2019-10-01' and date '2020-03-17'")
res["mjd"]=Time(res['time_recorded']).mjd

if 1 :
    print("adding telemetry from tower")
    tmp=dbquery("select pressure,temperature,humidity,time_recorded from environmentmonitor_tower where time_recorded between date '2019-10-01' and date '2020-03-17'")

    print(" converting time ...")
    mjd=Time(tmp['time_recorded']).mjd
    print(" done")

    for k in tmp.dtype.names :
        if k=='time_recorded' : continue
        print("adding {} (with time interpolation)".format(k))
        res[k]=np.interp(res["mjd"],mjd,tmp[k])

if 1 :
    print("adding spectrographs_sensors telemetry")
    tmp=dbquery("select blue_camera_temp,blue_camera_humidity,red_camera_temp,red_camera_humidity,nir_camera_temp,nir_camera_humidity,time_recorded,unit from spectrographs_sensors where time_recorded between date '2019-10-01' and date '2020-03-17'")

    print(" converting time ...")
    mjd=Time(tmp['time_recorded']).mjd
    print(" done")

    units=np.unique(tmp["unit"])
    print("units=",units)
    for unit in np.unique(units) :
        for k in tmp.dtype.names :
            if k=='time_recorded' : continue
            if k=='unit' : continue
            k2=k+"_{}".format(unit)
            print("adding {} (with time interpolation)".format(k2))
            ii=(tmp["unit"]==unit)
            res[k2]=np.interp(res["mjd"],mjd[ii],tmp[k][ii])

res.write("telemetry.csv",overwrite=True)

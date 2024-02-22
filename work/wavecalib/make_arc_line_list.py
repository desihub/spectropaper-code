#!/usr/bin/env python

import os
import numpy as np

filename=os.environ["SPECEXDATA"]+"/specex_linelist_desi.txt"

res={}
res["lamp"]=[]
res["wavelength"]=[]

with open(filename) as file :

    for line in file.readlines() :
        if line[0]=="#" : continue
        vals=line.strip().split(" ")
        try :
            lamp=str(vals[0][:-1])
            wave=float(vals[1])
        except :
            continue
        #print(lamp,wave)
        res["lamp"].append(lamp)
        res["wavelength"].append(wave)
for k in res.keys() :
    res[k]=np.hstack(res[k])
ii=np.argsort(res["wavelength"])
for k in res.keys() :
    res[k]=res[k][ii]

if 1 : # a real table
    nlines=len(res["wavelength"])
    ncols=7
    nrows=nlines//ncols+1

    #text="""
    #\\begin{table}
    #\\centering
    #\\footnotesize
    #\\begin{tabular}{llll}
    #\\multicolumn{4}{c}{Wavelength(Lamp) ($\AA$)}\\\\
    #"""
    #print(text)
    for r in range(nrows) :
        line=""
        for c in range(ncols) :
            i=c+r*ncols
            #i=r+c*nrows
            if i<nlines :
                if c>0 :
                    line += " & "
                line+="{:.3f}({})".format(res["wavelength"][i],res["lamp"][i])
        line += " \\\\"
        print(line)
else : # a simple list
    line=""
    for i in range(len(res["wavelength"])) :
        line+="{:.3f}({}) ".format(res["wavelength"][i],res["lamp"][i])
    print(line)

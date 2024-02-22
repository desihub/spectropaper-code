#!/usr/bin/env python


import numpy as np
import matplotlib.pyplot as plt
import sys,os
from astropy.table import Table
import fitsio
from desispec.calibfinder import sm2sp

t=Table.read("rdnoise-desi-202206.csv")

sms=["SM%d"%i for i in range(1,11)]

text="""\\begin{table}[h]
\\centering
\\small
\\begin{tabular}{ccccc}
Spectrograph & \multicolumn{4}{c}{Read noise (electrons)} \\\\
- Camera & A & B & C & D \\\\
\\hline
"""
print(text)

for sm in sms :
    sp=sm2sp(sm.lower()).upper()
    petal=int(sp.replace("SP",""))
    #print("\\\\")
    #print("\multicolumn{5}{c}{Spectrograph %s (Petal %d)}\\\\"%(sm,petal))
    #print("\\hline")
    line=""
    cameras=[f"{b}{petal}" for b in ["b","r","z"]]
    for cam in cameras :
        line=f"{sm} - {cam}"
        for amp in ["A","B","C","D"] :
            ii=(t["CAMERA"]==cam)&(t["AMP"]==amp)
            #print(t[ii])
            rdnoise=np.median(t["RDNOISE"][ii])
            line+=" & {:1.2f}".format(rdnoise)
        line += "\\\\"
        print(line)
text="""\\end{tabular}
\\caption{CCD read noise for each camera and amplifier (average values from June 2022).}
\\label{table:spectrographs-ccd-rdnoise-table}
\\end{table}
"""
print(text)

if (1) :
    text="""\\begin{table}
\\centering
\\small
\\begin{tabular}{ccccc}
Spectrograph & \multicolumn{4}{c}{Gain (elec./ADU)} \\\\
- Camera & A & B & C & D \\\\
\\hline"""
    print(text)

    for sm in sms :
        sp=sm2sp(sm.lower()).upper()
        petal=int(sp.replace("SP",""))
        #print("\multicolumn{5}{c}{Spectrograph %s (Petal %d)}\\\\"%(sm,petal))
        #print("\\hline")
        line=""
        cameras=[f"{b}{petal}" for b in ["b","r","z"]]
        for cam in cameras :
            line=f"{sm} - {cam}"
            for amp in ["A","B","C","D"] :
                ii=(t["CAMERA"]==cam)&(t["AMP"]==amp)
                gain=np.median(t["GAIN"][ii])
                line+=" & {:1.2f}".format(gain)
            line += "\\\\"
            print(line)
    text="""\\end{tabular}
    \\caption{CCD electronic gains for each camera and amplifier.}
    \\label{table:spectrographs-ccd-gain-table}
    \\end{table}
    """
    print(text)

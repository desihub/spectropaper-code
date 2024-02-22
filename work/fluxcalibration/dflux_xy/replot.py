#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table

fig=plt.figure("flux-ratio-vs-xy")

for i,filename in enumerate(["flux-ratio-vs-xy-4500A-5500A.csv","flux-ratio-vs-xy-6000A-7300A.csv","flux-ratio-vs-xy-8500A-9800A.csv"]) :
    a=plt.subplot(1,3,i+1)
    t=Table.read(filename)
    ok=(t["R"]>0)
    a.scatter(t["X"][ok],t["Y"][ok],c=t["R"][ok],vmin=0.8,vmax=1.2,cmap="gray")
    a.axis("off")

plt.show()

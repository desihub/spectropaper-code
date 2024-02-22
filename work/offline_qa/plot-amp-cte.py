#!/usr/bin/env python

"""
Plot median of columns above/below amp boundaries
"""

import os, sys
import argparse
import numpy as np
import fitsio

import matplotlib.pyplot as plt

p = argparse.ArgumentParser()
p.add_argument('-n', '--night', type=int,help='YEARMMDD night')
p.add_argument('-e', '--expid', type=int, help='Exposure ID')
p.add_argument('-c', '--camera', type=str, help='Camera')
p.add_argument('-i', '--input', help='input preproc file')
p.add_argument('--nrow', type=int, default=21,
        help='number of rows to include in median')
p.add_argument('--xminmax', nargs=2, type=int, default=(0, 0),
        help='x (column) range to plot')
p.add_argument('--debias',action='store_true')

args = p.parse_args()

n = args.nrow
xmin, xmax = args.xminmax

if args.input is None:
    from desispec.io import findfile
    args.input = findfile('preproc', night=args.night, expid=args.expid,
            camera=args.camera)

img = fitsio.read(args.input, 'IMAGE')

ny, nx = img.shape

if xmax == 0 :
    xmin = nx//2 -300
    xmax = nx//2 + 300

above = np.median(img[ny//2:ny//2+n, xmin:xmax], axis=0)
below = np.median(img[ny//2-n:ny//2, xmin:xmax], axis=0)

if args.debias :
    margin=100
    bias1 = np.median(img[ny//2:ny//2+n,0:margin])
    bias2 = np.median(img[ny//2:ny//2+n,nx-margin:nx])
    print("above bias = {} {}".format(bias1,bias2))
    above[:nx//2-xmin] -= bias1
    above[nx//2-xmin:] -= bias2
    bias1 = np.median(img[ny//2-n:ny//2,0:margin])
    bias2 = np.median(img[ny//2-n:ny//2,nx-margin:nx])
    print("below bias = {} {}".format(bias1,bias2))
    below[:nx//2-xmin] -= bias1
    below[nx//2-xmin:] -= bias2




xx = np.arange(xmin, xmax)

title="qa-ctedet"
fig=plt.figure(title,figsize=(9,4))
ax=plt.subplot(111)
ax.plot(xx, below, label='below',color="C1")
ax.plot(xx, above, label='above',color="C0")
ax.axvline(nx//2,linestyle="--",color="k")
handles, labels = ax.get_legend_handles_labels()
ax.legend(reversed(handles), reversed(labels),loc="upper right")
plt.title(f'median of {n} rows above/below CCD amp boundary')
plt.ylim(-5,min(80,max(np.max(above),np.max(below))))
plt.xlim(xmin, xmax)
plt.xlabel('CCD column')
plt.grid()
plt.tight_layout()
fig.savefig(title+".pdf",format="pdf")
plt.show()

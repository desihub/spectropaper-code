#!/usr/bin/env python

import sys,os
import argparse
import matplotlib.pyplot as plt
import numpy as np
import fitsio

#from pkg_resources import resource_exists, resource_filename

#from desispec.util import parse_fibers
from desispec.qproc.io import read_qframe,write_qframe

#from desispec.io import read_fibermap,read_frame
#from desispec.interpolation import resample_flux
#from desispec.fluxcalibration import isStdStar

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('-i','--infile', type = str, default = None, required = True, nargs="*",
                    help = 'path to one or several frame fits files')
parser.add_argument('-o','--outfile', type = str, default = None, required = True,
                    help = 'outfile')


args   = parser.parse_args()

frame = None
flux = []

for filename in args.infile :
    if not os.path.isfile(filename) :
        print("missing",filename)
        continue
    print("adding",filename)
    tmp_frame = read_qframe(filename)
    if frame is None : frame = tmp_frame
    flux.append(frame.flux)
flux = np.array(flux)
frame.flux = np.median(flux,axis=0)

write_qframe(args.outfile,frame)
print("wrote",args.outfile)

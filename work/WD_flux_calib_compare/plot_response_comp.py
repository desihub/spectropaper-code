import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from astropy.table import Table


#####################################################
################## Plotting section #################
#####################################################
matplotlib.rc('font',**{'size' : 28, 'family':'serif',
'serif':['serif']})
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['ps.usedistiller'] = 'xpdf'
#matplotlib.rcParams['xtick.major.pad'] = 22
matplotlib.rcParams['ytick.major.pad'] = 22
#####################################################

#fig = plt.figure(figsize = (12,10))
fig = plt.figure(figsize = (9,8))
axesb = fig.add_axes([0,0.72,1,0.28])
axesr = fig.add_axes([0,0.36,1,0.28])
axesz = fig.add_axes([0,0   ,1,0.28])

##################################################
xmajorLocator   = MultipleLocator(250)
xmajorFormatter = FormatStrFormatter('%d')
xminorLocator   = MultipleLocator(50)
axesb.xaxis.set_major_locator(xmajorLocator)
axesb.xaxis.set_major_formatter(xmajorFormatter)
axesb.xaxis.set_minor_locator(xminorLocator)

ymajorLocator   = MultipleLocator(0.04)
ymajorFormatter = FormatStrFormatter('%.2f')
yminorLocator   = MultipleLocator(0.01)
#axesb.yaxis.set_major_locator(ymajorLocator)
#axesb.yaxis.set_major_formatter(ymajorFormatter)
#axesb.yaxis.set_minor_locator(yminorLocator)
##################################################

##################################################
xmajorLocator   = MultipleLocator(250)
xmajorFormatter = FormatStrFormatter('%d')
xminorLocator   = MultipleLocator(50)
axesr.xaxis.set_major_locator(xmajorLocator)
axesr.xaxis.set_major_formatter(xmajorFormatter)
axesr.xaxis.set_minor_locator(xminorLocator)

ymajorLocator   = MultipleLocator(0.04)
ymajorFormatter = FormatStrFormatter('%.2f')
yminorLocator   = MultipleLocator(0.01)
#axesr.yaxis.set_major_locator(ymajorLocator)
#axesr.yaxis.set_major_formatter(ymajorFormatter)
#axesr.yaxis.set_minor_locator(yminorLocator)
##################################################

##################################################
xmajorLocator   = MultipleLocator(250)
xmajorFormatter = FormatStrFormatter('%d')
xminorLocator   = MultipleLocator(50)
#axesz.xaxis.set_major_locator(xmajorLocator)
#axesz.xaxis.set_major_formatter(xmajorFormatter)
#axesz.xaxis.set_minor_locator(xminorLocator)

ymajorLocator   = MultipleLocator(0.04)
ymajorFormatter = FormatStrFormatter('%.2f')
yminorLocator   = MultipleLocator(0.01)
axesz.yaxis.set_major_locator(ymajorLocator)
axesz.yaxis.set_major_formatter(ymajorFormatter)
axesz.yaxis.set_minor_locator(yminorLocator)
##################################################


axesb.set_xlim(3590, 5810)
axesr.set_xlim(5750, 7630)
axesz.set_xlim(7510, 9835)
axesb.set_ylim(0.925,1.075)
axesr.set_ylim(0.925,1.075)
axesz.set_ylim(0.925,1.075)
axesz.set_xlabel('Wavelength[\AA]')
#axesb.set_ylabel('B-arm WD residual')
#axesr.set_ylabel('R-arm WD residual')
#axesz.set_ylabel('Z-arm WD residual')
axesb.set_ylabel('b')
axesr.set_ylabel('r')
axesz.set_ylabel('z')


axesb.axhline(1, lw = 0.5, ls = '--', color = '0.5', zorder = 0)
axesr.axhline(1, lw = 0.5, ls = '--', color = '0.5', zorder = 0)
axesz.axhline(1, lw = 0.5, ls = '--', color = '0.5', zorder = 0)

axesb.axvspan(4500,4700, color = '0.7', zorder = 0)

tmp = [axesb, axesr, axesz]
for axes in tmp:

  lines= [3656.65 ,3657.25,3658.04 ,3658.65 ,3659.41 ,3660.32 ,3661.27,3662.22 ,3663.41 ,3664.65 ,3666.08 ,3667.73,3669.45 ,3671.32 ,3673.81 ,3676.376,3679.370,3682.823,3686.831,3691.551,3697.157,3703.859,3711.978,3721.946,3734.369,3750.151,3770.633,3797.874,3835.397,3889.024,3970.075,4101.734,4340.472,4861.35,6562.79]
  for line in lines:
    axes.axvline(line, 0.1, 0.17, lw = 2, color = 'C1', alpha = 1, zorder = 1)

  lines = [8545.38,8598.39,8665.02,8750.46,8862.8,9015.3 ,9229.7 ,9546.2]
  for line in lines:
    axes.axvline(line, 0.1, 0.17, lw = 2, color = 'C1', alpha = 1, zorder = 1)

wb, fb = np.genfromtxt('Fuji_response_compare_b_arm.dat', unpack = True)
wr, fr = np.genfromtxt('Fuji_response_compare_r_arm.dat', unpack = True)
wz, fz = np.genfromtxt('Fuji_response_compare_z_arm.dat', unpack = True)

axesb.plot(wb, fb, color = '0.2', lw = 1.5, alpha = 1, zorder = 100)
axesr.plot(wr, fr, color = '0.2', lw = 1.5, alpha = 1, zorder = 100)
axesz.plot(wz, fz, color = '0.2', lw = 1.5, alpha = 1, zorder = 100)

plt.savefig('Fuji_calib_compare_paper.pdf', bbox_inches = 'tight')

zenodo = False
if zenodo :
  t = Table()
  t["WAVELENGTH"] = np.hstack([wb,wr,wz])
  t["RATIO"] = np.hstack([fb,fr,fz])
  t["CAMERA"]=np.hstack([np.repeat("B",fb.size),np.repeat("R",fr.size),np.repeat("Z",fz.size)])
  t.write("../../zenodo/figure-41.fits")

plt.show()
plt.close()

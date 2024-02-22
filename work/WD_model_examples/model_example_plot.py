import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import scipy.interpolate
import glob
from astropy.table import Table

def factor(wb, fb, mw, mf):
  check_f_spec=fb[(wb>4500.) & (wb<4700.)]
  mf[np.isnan(mf)] = 0.0
  check_f_model = mf[(mw > 4500) & (mw < 4700)]
  adjust = np.average(check_f_model)/np.average(check_f_spec)
  return adjust

spec1_b = "J095529.86+025913.47_39627859702057428-20210407-00083857-2152-b.dat"
spec1_r = "J095529.86+025913.47_39627859702057428-20210407-00083857-2152-r.dat"
spec1_z = "J095529.86+025913.47_39627859702057428-20210407-00083857-2152-z.dat"
spec2_b = "J095725.07+004151.71_39627805352266243-20210407-00083857-0520-b.dat"
spec2_r = "J095725.07+004151.71_39627805352266243-20210407-00083857-0520-r.dat"
spec2_z = "J095725.07+004151.71_39627805352266243-20210407-00083857-0520-z.dat"
spec3_b = "J095833.13+013049.20_39627823475854146-20210407-00083857-0625-b.dat"
spec3_r = "J095833.13+013049.20_39627823475854146-20210407-00083857-0625-r.dat"
spec3_z = "J095833.13+013049.20_39627823475854146-20210407-00083857-0625-z.dat"
spec4_b = "J095902.77+013848.21_39627829519843567-20210407-00083857-0917-b.dat"
spec4_r = "J095902.77+013848.21_39627829519843567-20210407-00083857-0917-r.dat"
spec4_z = "J095902.77+013848.21_39627829519843567-20210407-00083857-0917-z.dat"

model_1 = "17276.9_7.897_82.88.model"
model_2 = "30117.1_7.879_135.6.model"
model_3 = "11952.1_8.220_39.85.model"
model_4 = "12363.9_8.291_36.76.model"

w1b, f1b, i1b = np.genfromtxt(spec1_b, unpack = True)
w1r, f1r, i1r = np.genfromtxt(spec1_r, unpack = True)
w1z, f1z, i1z = np.genfromtxt(spec1_z, unpack = True)
w2b, f2b, i2b = np.genfromtxt(spec2_b, unpack = True)
w2r, f2r, i2r = np.genfromtxt(spec2_r, unpack = True)
w2z, f2z, i2z = np.genfromtxt(spec2_z, unpack = True)
w3b, f3b, i3b = np.genfromtxt(spec3_b, unpack = True)
w3r, f3r, i3r = np.genfromtxt(spec3_r, unpack = True)
w3z, f3z, i3z = np.genfromtxt(spec3_z, unpack = True)
w4b, f4b, i4b = np.genfromtxt(spec4_b, unpack = True)
w4r, f4r, i4r = np.genfromtxt(spec4_r, unpack = True)
w4z, f4z, i4z = np.genfromtxt(spec4_z, unpack = True)

mw1, mf1, = np.genfromtxt(model_1, unpack = True)
mw2, mf2, = np.genfromtxt(model_2, unpack = True)
mw3, mf3, = np.genfromtxt(model_3, unpack = True)
mw4, mf4, = np.genfromtxt(model_4, unpack = True)


adjust1 = factor(w1b, f1b, mw1, mf1)
adjust2 = factor(w2b, f2b, mw2, mf2)
adjust3 = factor(w3b, f3b, mw3, mf3)
adjust4 = factor(w4b, f4b, mw4, mf4)

#####################################################
################## Plotting section #################
#####################################################
matplotlib.rc('font',**{'size' : 20, 'family':'serif',
'serif':['serif']})
matplotlib.rcParams['text.usetex'] = True
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['ps.usedistiller'] = 'xpdf'
matplotlib.rcParams['xtick.major.pad'] = 16
matplotlib.rcParams['ytick.major.pad'] = 16
#####################################################
#fig = plt.figure(figsize= (7,8))
fig = plt.figure(figsize= (6.5,7))


##### AXES1 #####

axes1 = fig.add_axes([0,0,1,0.25])
axes1.set_xlabel('Wavelength (\AA)')

##################################################
xmajorLocator   = MultipleLocator(1000)
xmajorFormatter = FormatStrFormatter('%d')
xminorLocator   = MultipleLocator(200)
axes1.xaxis.set_major_locator(xmajorLocator)
axes1.xaxis.set_major_formatter(xmajorFormatter)
axes1.xaxis.set_minor_locator(xminorLocator)

ymajorLocator   = MultipleLocator(25)
ymajorFormatter = FormatStrFormatter('%d')
yminorLocator   = MultipleLocator(5)
#axes1.yaxis.set_major_locator(ymajorLocator)
#axes1.yaxis.set_major_formatter(ymajorFormatter)
#axes1.yaxis.set_minor_locator(yminorLocator)
##################################################

axes1.set_ylim(0.1, 120)
axes1.set_xlim(3500, 9999)
axes1.axvspan(4500,4700, color = '0.8')

axes1.text(-0.14, +1.25, r'Flux (10$^{-17}$\,erg\,s$^{-1}$\,cm$^{-2}$\,\AA$^{-1}$)', transform=axes1.transAxes, rotation = 90)
axes1.plot(w1b, f1b, color = 'C0', linewidth = 0.6)
axes1.plot(w1r, f1r, color = 'C1', linewidth = 0.6)
axes1.plot(w1z, f1z, color = 'C3', linewidth = 0.6)
axes1.plot(mw1, mf1/adjust1, color = '0.3', linewidth = 1)

zenodo=False

if zenodo :
  t=Table()
  t["WAVELENGTH"]=np.hstack([w1b,w1r,w1z])
  t["MEASURED_FLUX"]=np.hstack([f1b,f1r,f1z])
  t["CAMERA"]=np.hstack([np.repeat("B",f1b.size),np.repeat("R",f1r.size),np.repeat("Z",f1z.size)])
  t.write("../../zenodo/figure-40a-data.fits")
  #t.write("../../zenodo/figure-40a-data.csv")
  t=Table()
  t["WAVELENGTH"]=mw1
  t["MODEL_FLUX"]=mf1/adjust1
  t.write("../../zenodo/figure-40a-model.fits")
  #t.write("../../zenodo/figure-40a-model.csv")




##### AXES2 #####

axes2 = fig.add_axes([0,0.25,1,0.25])

##################################################
xmajorLocator   = MultipleLocator(1000)
xmajorFormatter = FormatStrFormatter('%d')
xminorLocator   = MultipleLocator(200)
axes2.xaxis.set_major_locator(xmajorLocator)
axes2.xaxis.set_major_formatter(xmajorFormatter)
axes2.xaxis.set_minor_locator(xminorLocator)

ymajorLocator   = MultipleLocator(25)
ymajorFormatter = FormatStrFormatter('%d')
yminorLocator   = MultipleLocator(5)
axes2.yaxis.set_major_locator(ymajorLocator)
axes2.yaxis.set_major_formatter(ymajorFormatter)
axes2.yaxis.set_minor_locator(yminorLocator)
##################################################

axes2.set_ylim(0.1, 135)
axes2.set_xlim(3500, 10000)
axes2.axvspan(4500,4700, color = '0.8')

axes2.xaxis.set_ticklabels([])
axes2.plot(w2b, f2b, color = 'C0', linewidth = 0.6)
axes2.plot(w2r, f2r, color = 'C1', linewidth = 0.6)
axes2.plot(w2z, f2z, color = 'C3', linewidth = 0.6)
axes2.plot(mw2, mf2/adjust2, color = '0.3', linewidth = 1)

if zenodo :
  t=Table()
  t["WAVELENGTH"]=np.hstack([w2b,w2r,w2z])
  t["MEASURED_FLUX"]=np.hstack([f2b,f2r,f2z])
  t["CAMERA"]=np.hstack([np.repeat("B",f2b.size),np.repeat("R",f2r.size),np.repeat("Z",f2z.size)])
  t.write("../../zenodo/figure-40b-data.fits")
  #t.write("../../zenodo/figure-40b-data.csv")
  t=Table()
  t["WAVELENGTH"]=mw2
  t["MODEL_FLUX"]=mf2/adjust2
  t.write("../../zenodo/figure-40b-model.fits")
  #t.write("../../zenodo/figure-40b-model.csv")


##### AXES3 #####

axes3 = fig.add_axes([0,0.5,1,0.25])

##################################################
xmajorLocator   = MultipleLocator(1000)
xmajorFormatter = FormatStrFormatter('%d')
xminorLocator   = MultipleLocator(200)
axes3.xaxis.set_major_locator(xmajorLocator)
axes3.xaxis.set_major_formatter(xmajorFormatter)
axes3.xaxis.set_minor_locator(xminorLocator)

ymajorLocator   = MultipleLocator(25)
ymajorFormatter = FormatStrFormatter('%d')
yminorLocator   = MultipleLocator(5)
axes3.yaxis.set_major_locator(ymajorLocator)
axes3.yaxis.set_major_formatter(ymajorFormatter)
axes3.yaxis.set_minor_locator(yminorLocator)
##################################################

axes3.set_ylim(0.1, 148)
axes3.set_xlim(3500, 10000)
axes3.axvspan(4500,4700, color = '0.8')

axes3.xaxis.set_ticklabels([])
axes3.plot(w3b, f3b, color = 'C0', linewidth = 0.6)
axes3.plot(w3r, f3r, color = 'C1', linewidth = 0.6)
axes3.plot(w3z, f3z, color = 'C3', linewidth = 0.6)
axes3.plot(mw3, mf3/adjust3, color = '0.3', linewidth = 1)

if zenodo :
  t=Table()
  t["WAVELENGTH"]=np.hstack([w3b,w3r,w3z])
  t["MEASURED_FLUX"]=np.hstack([f3b,f3r,f3z])
  t["CAMERA"]=np.hstack([np.repeat("B",f3b.size),np.repeat("R",f3r.size),np.repeat("Z",f3z.size)])
  t.write("../../zenodo/figure-40c-data.fits")
  #t.write("../../zenodo/figure-40c-data.csv")
  t=Table()
  t["WAVELENGTH"]=mw3
  t["MODEL_FLUX"]=mf3/adjust3
  t.write("../../zenodo/figure-40c-model.fits")
  #t.write("../../zenodo/figure-40c-model.csv")

##### AXES4 #####

axes4 = fig.add_axes([0,0.75,1,0.25])

##################################################
xmajorLocator   = MultipleLocator(1000)
xmajorFormatter = FormatStrFormatter('%d')
xminorLocator   = MultipleLocator(200)
axes4.xaxis.set_major_locator(xmajorLocator)
axes4.xaxis.set_major_formatter(xmajorFormatter)
axes4.xaxis.set_minor_locator(xminorLocator)

ymajorLocator   = MultipleLocator(15)
ymajorFormatter = FormatStrFormatter('%d')
yminorLocator   = MultipleLocator(3)
axes4.yaxis.set_major_locator(ymajorLocator)
axes4.yaxis.set_major_formatter(ymajorFormatter)
axes4.yaxis.set_minor_locator(yminorLocator)
##################################################

axes4.set_ylim(0.1, 60)
axes4.set_xlim(3500, 10000)
axes4.axvspan(4500,4700, color = '0.8')
axes4.xaxis.set_ticklabels([])
axes4.plot(w4b, f4b, color = 'C0', linewidth = 0.6)
axes4.plot(w4r, f4r, color = 'C1', linewidth = 0.6)
axes4.plot(w4z, f4z, color = 'C3', linewidth = 0.6)
axes4.plot(mw4, mf4/adjust4, color = '0.3', linewidth = 1)

if zenodo :
  t=Table()
  t["WAVELENGTH"]=np.hstack([w4b,w4r,w4z])
  t["MEASURED_FLUX"]=np.hstack([f4b,f4r,f4z])
  t["CAMERA"]=np.hstack([np.repeat("B",f4b.size),np.repeat("R",f4r.size),np.repeat("Z",f4z.size)])
  t.write("../../zenodo/figure-40d-data.fits")
  #t.write("../../zenodo/figure-40d-data.csv")
  t=Table()
  t["WAVELENGTH"]=mw4
  t["MODEL_FLUX"]=mf4/adjust4
  t.write("../../zenodo/figure-40d-model.fits")
  #t.write("../../zenodo/figure-40d-model.csv")

plt.savefig('WD_model_fit_examples.pdf', bbox_inches = 'tight')
plt.show()
plt.close()

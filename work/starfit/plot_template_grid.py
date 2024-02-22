#!/usr/bin/env python

import os,sys
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pyfits
from desispec.io.filters import load_legacy_survey_filter
import astropy.units as units

filename=os.path.join(os.getenv('DESI_BASIS_TEMPLATES'), 'stdstar_templates_v2.2.fits')
h=pyfits.open(filename)
h.info()
t=h["METADATA"].data
print(t.dtype.names)
ok=np.where((t["TEFF"]>4900)&(t["TEFF"]<9100))[0]
t=t[ok]
wave=h[2].data
flux=h[0].data[ok]
fluxunits = 1e-17 * units.erg / units.s / units.cm**2 / units.Angstrom
g_filter=load_legacy_survey_filter(band="G",photsys="S")
r_filter=load_legacy_survey_filter(band="R",photsys="S")
z_filter=load_legacy_survey_filter(band="Z",photsys="S")
gmag=g_filter.get_ab_magnitude(flux*fluxunits,wave)
rmag=r_filter.get_ab_magnitude(flux*fluxunits,wave)
zmag=z_filter.get_ab_magnitude(flux*fluxunits,wave)
gr=gmag-rmag
rz=rmag-zmag

print("wave[0]=",wave[0])
print("wave[-1]=",wave[-1])
print("c dwave/wave= {} km/s".format(np.mean(np.gradient(np.log(wave)))*2.997e5))


stdstars=(gr>0)&(gr<0.35)&(rz<0.2)

#color0="gray"
#color1="k"
color0="C0"
color1="C1"

alpha=1


fig=plt.figure("template-grid",figsize=(5,5))
a1=plt.subplot(221)
a1.plot(t["TEFF"],t["FEH"],"o",color=color0)
a1.plot(t["TEFF"][stdstars],t["FEH"][stdstars],"o",color=color1)
a1.grid()
a1.set_xlabel("$T_{eff}$ (K)")
a1.set_ylabel("[Fe/H]")

a2=plt.subplot(222)
a2.plot(t["LOGG"],t["FEH"],"o",alpha=alpha,color=color0)
a2.plot(t["LOGG"][stdstars],t["FEH"][stdstars],"o",alpha=alpha,color=color1)
a2.grid()
a2.set_ylabel("[Fe/H]")
a2.set_xlabel("$\log_{10} g$")

a3=plt.subplot(223)
a3.plot(t["TEFF"],t["LOGG"],"o",alpha=alpha,color=color0)
a3.plot(t["TEFF"][stdstars],t["LOGG"][stdstars],"o",alpha=alpha,color=color1)
a3.grid()
a3.set_xlabel("$T_{eff}$ (K)")
a3.set_ylabel("$\log_{10} g$")

a4=plt.subplot(224)
a4.plot(t["TEFF"],gr,"o",alpha=alpha,color=color0)
a4.plot(t["TEFF"][stdstars],gr[stdstars],"o",alpha=alpha,color=color1)
c=np.polyfit(t["TEFF"],gr,1)
print(c)
pol=np.poly1d(c)
x=np.array([4000,9000])
plt.plot(x,pol(x),":",color="gray")
a4.grid()
a4.set_xlabel("$T_{eff}$ (K)")
a4.set_ylabel("g-r (AB mags)")
a4.axhline(0,color="gray",linestyle="--")
a4.axhline(0.35,color="gray",linestyle="--")

#a4=plt.subplot(224)
#plt.colorbar(toto,ax=a3)

plt.tight_layout()

plt.show()

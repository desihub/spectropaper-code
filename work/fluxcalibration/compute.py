#!/usr/bin/env python

#import astropy.units as units

import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table

from desispec.io import read_spectra
from desispec.interpolation import resample_flux
from desispec.io.filters import load_legacy_survey_filter
from desispec.coaddition import coadd_cameras
from desispec.magnitude import compute_broadband_flux,ab_flux_in_ergs_s_cm2_A

from desitarget.targets import main_cmx_or_sv

# choose a main survey dark tile with a speed of ~1
# 3393,main,dark,dark,maindark,1,1591.0,309.664,1.226,1006.5,941.1,1098.0,1000.0,obsend,1014.2,941.1,995.7,849.0,dark,0.85,20210708

t=Table.read("/global/cfs/cdirs/desi/spectro/redux/guadalupe/tiles-guadalupe.csv")
selection=(t["SURVEY"]=="main")&(t["PROGRAM"]=="dark")&(t["EFFTIME_SPEC"]>900)&(t["EXPTIME"]<1100)
ii=np.where(selection)[0][:3]
filenames=[]
for i in ii :
    tile=t["TILEID"][i]
    night=t["LASTNIGHT"][i]
    for sp in range(10) :
        filenames.append(f"/global/cfs/cdirs/desi/spectro/redux/guadalupe/tiles/cumulative/{tile}/{night}/coadd-{sp}-{tile}-thru{night}.fits")
print(filenames)


decam_r=load_legacy_survey_filter("r","S")

# trim
ii=np.where(decam_r.response>0.05*np.max(decam_r.response))[0]
decam_r_wavelength = decam_r.wavelength[ii]
decam_r_response   = decam_r.response[ii]

#plt.plot(decam_r_wavelength,decam_r_response)
#plt.show()

res={}
for k in ["TILEID","NIGHT","FIBER","IMAGING_FLUX_R","IMAGING_FIBERFLUX_R","SPECTRO_FLUX_R","SPECTRO_FIBERFLUX_R"] :
    res[k]=[]

for filename in filenames :
    print(filename)
    spectra=read_spectra(filename)
    coadd_spectra=coadd_cameras(spectra)

    # keep only entries from first expid
    expid=spectra.exp_fibermap["EXPID"][0]
    fmap=spectra.exp_fibermap[spectra.exp_fibermap["EXPID"]==expid]
    # match target ids
    t2i={t:i for i,t in enumerate(fmap["TARGETID"])}
    ii=[t2i[t] for t in coadd_spectra.fibermap["TARGETID"]]
    fmap=fmap[ii]
    # add back "PSF_TO_FIBER_SPECFLUX"
    coadd_spectra.fibermap["PSF_TO_FIBER_SPECFLUX"] = fmap["PSF_TO_FIBER_SPECFLUX"]
    spectra=coadd_spectra

    # avoid QSOs that can be variable
    target_colnames, target_masks, survey = main_cmx_or_sv(spectra.fibermap)
    desi_target = spectra.fibermap[target_colnames[0]]  # (SV1_)DESI_TARGET
    desi_mask   = target_masks[0]

    #selection_mask = desi_mask.mask("LRG|ELG|ELG_LOP|ELG_HIP|STD_FAINT|STD_WD|STD_BRIGHT|BGS_ANY")
    #selection      = np.where((desi_target & selection_mask)>0)[0]
    selection      =  np.where(spectra.fibermap["OBJTYPE"]=="TGT")[0]
    wave=spectra.wave["brz"]
    flux=spectra.flux["brz"][selection]
    ivar=spectra.ivar["brz"][selection]
    mask=spectra.mask["brz"][selection]

    wmin=min(wave[0],decam_r.wavelength[0])
    wmax=max(wave[-1],decam_r.wavelength[-1])
    nn=int((wmax-wmin)/10.+1)
    twave=np.linspace(wmin,wmax,nn)
    denominator = compute_broadband_flux(twave,ab_flux_in_ergs_s_cm2_A(twave),decam_r_wavelength,decam_r_response)
    flux_scale  = 10**(0.4*22.5)/denominator


    #fluxunits = 1e-17 * units.erg / units.s / units.cm**2 / units.Angstrom



    spec_rflux  = np.zeros(flux.shape[0])
    for i in range(flux.shape[0]) :
        #print(i)
        tflux,_ = resample_flux(twave,wave,flux[i],ivar=ivar[i]*(mask[i]==0))
        bbflux = 1e-17 * compute_broadband_flux(twave,tflux,decam_r_wavelength,decam_r_response)

        if False and bbflux>0 : # cross-checking
            rmag  = -2.5*np.log10(bbflux/denominator)
            rmag2 =  decam_r.get_ab_magnitude(tflux * fluxunits, twave)
            print(i,rmag,rmag2)

        spec_rflux[i]   = bbflux * flux_scale

    nspec=selection.size
    tile=spectra.meta["TILEID"]
    res["TILEID"].append( np.repeat(tile,nspec))
    night=spectra.meta["NIGHT"]
    res["NIGHT"].append( np.repeat(night,nspec))
    res["FIBER"].append( list(spectra.fibermap["FIBER"][selection]) )
    res["IMAGING_FLUX_R"].append( list(spectra.fibermap["FLUX_R"][selection]) )
    res["IMAGING_FIBERFLUX_R"].append( list(spectra.fibermap["FIBERFLUX_R"][selection]) )
    res["SPECTRO_FLUX_R"].append( spec_rflux )
    res["SPECTRO_FIBERFLUX_R"].append( list(spec_rflux*spectra.fibermap["PSF_TO_FIBER_SPECFLUX"][selection]) )


    #plt.plot(spectra.fibermap["FLUX_R"],spec_rflux,".",color="C0")
    #plt.plot(spectra.fibermap["FIBERFLUX_R"][selection],spec_rflux*spectra.fibermap["PSF_TO_FIBER_SPECFLUX"][selection],".",color="C1")

t=Table()
for k in res.keys() :
    v = np.hstack(res[k])
    print(k,len(v))
    t[k] = v
t.write("fluxes.csv",overwrite=True)
print("wrote fluxes.csv")

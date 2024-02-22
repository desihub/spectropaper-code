#!/usr/bin/env python

from glob import glob
import os
import fitsio
import numpy as np
from astropy.table import Table
from desitarget.targetmask import desi_mask
from desispec.maskbits import fibermask
from desispec.night_qa import  get_surveys_night_expids, get_ctedet_night_expid
from desispec.io import findfile, read_spectra, read_table
from desispec.interpolation import resample_flux
import redrock.templates
from desiutil.log import get_logger
import matplotlib.pyplot as plt

plots = "ctedet,skyzfiber,elgspectrum"
# plots = "ctedet"

#outdir = "/global/cfs/cdirs/desi/users/raichoor/spectro_paper"
outdir = "."
prod = "/global/cfs/cdirs/desi/spectro/redux/daily"
fig_fmt = "pdf"


def get_infos():
    tileid, night, archivedate = 6611, 20220524, 20220525
    petal, camera = 3, "z"
    elg_tid, elg_fiber = 39628523475831150, 1727
    return tileid, night, archivedate, petal, camera, elg_tid, elg_fiber


if "ctedet" in plots.split(","):
    tileid, night, archivedate, petal, camera, elg_tid, elg_fiber = get_infos()
    # https://github.com/desihub/desispec/blob/063d71e471c8c39dbbfbba0d0bf48dc6a77a092a/py/desispec/night_qa.py#L466-L531
    #ctedet_expid = get_ctedet_night_expid(night, prod)
    ctedet_expid = 136562
    nrow = 21
    ylim = (-5, 10)
    clim = (-5, 5)
    petcam_xmin, petcam_xmax = None, None
    fig, ax = plt.subplots(figsize=(30, 5))
    #
    fn = findfile("preproc", night, ctedet_expid, camera+str(petal),
            specprod_dir=prod)
    ax.set_title(
        "NIGHT={}, EXPID={}, PETAL={}, CAMERA={} ; Median of {} rows above/below CCD amp boundary".format(
            night, ctedet_expid, petal, camera, nrow,
        )
    )
    # AR read
    with fitsio.FITS(fn) as fx:
        img = fx["IMAGE"].read()
    ny, nx = img.shape
    if petcam_xmin is None:
        petcam_xmin = 0
    if petcam_xmax is None:
        petcam_xmax = nx
    above = np.median(img[ny // 2: ny // 2 + nrow, petcam_xmin : petcam_xmax], axis=0)
    below = np.median(img[ny // 2 - nrow : ny // 2, petcam_xmin : petcam_xmax], axis=0)
    xx = np.arange(petcam_xmin, petcam_xmax)
    # AR plot 1d median
    ax.plot(xx, above, alpha=0.5, label="above (AMPC : x < {}; AMPD : x > {}".format(nx // 2 - 1, nx // 2 -1))
    ax.plot(xx, below, alpha=0.5, label="below (AMPA : x < {}; AMPB : x > {}".format(nx // 2 - 1, nx // 2 -1))
    ax.legend(loc=2)
    # AR amplifier x-boundary
    ax.axvline(nx // 2 - 1, color="k", ls="--")
    ax.set_xlabel("CCD column")
    ax.set_ylabel("Electrons")
    ax.set_xlim(petcam_xmin, petcam_xmax)
    ax.set_ylim(ylim)
    ax.grid()
    plt.savefig(os.path.join(outdir, "qa-ctedet.{}".format(fig_fmt)), bbox_inches="tight")
    plt.close()


if "skyzfiber" in plots.split(","):
    tileid, night, archivedate, petal, camera, elg_tid, elg_fiber = get_infos()
    # AR copying/adapting https://github.com/desihub/desispec/blob/6d23718c4e0a96f1b8e3790332f530e07521dc36/py/desispec/night_qa.py#L680-L779
    _, tileids, _ = get_surveys_night_expids(night)
    # tileids = [10474]
    #
    dchi2_threshold=9
    group="cumulative"
    log = get_logger()
    #
    # AR safe
    tileids = np.unique(tileids)
    # AR gather all infos from the redrock*fits files
    fibers, zs, dchi2s, faflavors = [], [], [], []
    nfn = 0
    for tileid in tileids:
        # AR main backup/bright/dark ?
        faflavor = None
        fns = sorted(
            glob(
                os.path.join(
                    os.getenv("DESI_ROOT"),
                    "spectro",
                    "data",
                    "{}".format(night),
                    "*",
                    "fiberassign-{:06d}.fits*".format(tileid),
                )
            )
        )
        if len(fns) > 0:
            hdr = fitsio.read_header(fns[0], 0)
            if "FAFLAVOR" in hdr:
                faflavor = hdr["FAFLAVOR"]
        log.info("identified FAFLAVOR for {}: {}".format(tileid, faflavor))
        # AR
        tmp = findfile("redrock", night=night, tile=tileid, groupname=group, spectrograph=0, specprod_dir=prod)
        tiledir = os.path.dirname(tmp)
        fns = sorted(glob(os.path.join(tiledir, f"redrock-?-{tileid}-*{night}.fits*")))
        nfn += len(fns)
        for fn in fns:
            fm = fitsio.read(fn, ext="FIBERMAP", columns=["OBJTYPE", "FIBER", "TARGETID"])
            rr = fitsio.read(fn, ext="REDSHIFTS", columns=["Z", "DELTACHI2"])
            sel = fm["OBJTYPE"] == "SKY"
            log.info("selecting {} / {} SKY fibers in {}".format(sel.sum(), len(rr), fn))
            fibers += fm["FIBER"][sel].tolist()
            zs += rr["Z"][sel].tolist()
            dchi2s += rr["DELTACHI2"][sel].tolist()
            faflavors += [faflavor for x in range(sel.sum())]
    fibers, zs, dchi2s, faflavors = np.array(fibers), np.array(zs), np.array(dchi2s), np.array(faflavors, dtype=str)
    # AR plot
    # plot_faflavors = ["all", "mainbackup", "mainbright", "maindark"]
    plot_faflavor = "maindark"
    ylim = (-1.1, 1.1)
    # ylim = (0, 2)
    yticks = np.array([0, 0.1, 0.25, 0.5, 1, 2, 3, 4, 5, 6])
    fig, ax = plt.subplots()
    if plot_faflavor == "all":
        faflavor_sel = np.ones(len(fibers), dtype=bool)
        title = "NIGHT = {}\nAll tiles ({} fibers)".format(night, len(fibers))
    else:
        faflavor_sel = faflavors == plot_faflavor
        # title = "NIGHT = {}\nFAFLAVOR={} ({} fibers)".format(night, plot_faflavor, faflavor_sel.sum())
        title = "NIGHT={}, FAFLAVOR={}, PETAL={}".format(night, plot_faflavor, petal)
    if faflavor_sel.sum() < 5000:
        alpha = 0.3
    else:
        # alpha = 0.1
        alpha = 0.5
    """
    for sel, selname, color in zip(
        [
            (faflavor_sel) & (dchi2s < dchi2_threshold),
            (faflavor_sel) & (dchi2s > dchi2_threshold),
        ],
        [
            "OBJTYPE=SKY and DELTACHI2<{}".format(dchi2_threshold),
            "OBJTYPE=SKY and DELTACHI2>{}".format(dchi2_threshold),
        ],
        ["orange", "b"]
    ):
    """
    for sel, selname, color in zip([faflavor_sel], ["OBJTYPE=SKY"], ["orange"]):
        # ax.scatter(fibers[sel], np.log10(0.1 + zs[sel]), c=color, s=1, alpha=alpha, label="{} ({} fibers)".format(selname, sel.sum()))
        ax.scatter(fibers[sel], np.log10(0.1 + zs[sel]), c=color, s=5, alpha=alpha, label=selname)
    ax.grid()
    ax.set_title(title)
    ax.set_xlabel("FIBER")
    # ax.set_xlim(-100, 5100)
    ax.set_xlim(1500, 2000)
    ax.set_ylabel("Z")
    ax.set_ylim(ylim)
    ax.set_yticks(np.log10(0.1 + yticks))
    ax.set_yticklabels(yticks.astype(str))
    # ax.legend(loc=2, markerscale=10)
    ax.legend(loc=2, markerscale=5)
    # AR highlihgt problematic region
    # AR display outer edge
    angs = np.linspace(2 * np.pi, 0, 1000)
    xcen, ycen = 1725, np.log10(0.1 + 1.2)
    xrad, yrad = 35, 0.20
    dxs = xcen + xrad * np.cos(angs)
    dys = ycen + yrad * np.sin(angs)
    ax.plot(dxs, dys, color="k", zorder=1)
    plt.savefig(os.path.join(outdir, "qa-skyzfiber.{}".format(fig_fmt)), bbox_inches="tight")
    plt.close()


# AR adapted from https://github.com/desihub/desispec/blob/master/bin/plot_spectra
if "elgspectrum" in plots.split(","):
    tileid, night, archivedate, petal, camera, elg_tid, elg_fiber = get_infos()
    rebin = 32
    spectra_fn = os.path.join(prod, "tiles", "archive", str(tileid), str(archivedate), "spectra-{}-{}-thru{}.fits".format(elg_fiber // 500, tileid, night))
    redrock_fn = spectra_fn.replace("spectra", "redrock")
    rest_frame = False
    ylim = (-1, 2)
    spec = read_spectra(spectra_fn)
    j = np.where(spec.fibermap["TARGETID"]==elg_tid)[0][0]
    wavescale=1.
    #
    templates = dict()
    for filename in redrock.templates.find_templates():
        tx = redrock.templates.Template(filename)
        templates[(tx.template_type, tx.sub_type)] = tx
    #
    lines = {
        'Ha'      : 6562.8,
        'Hb'       : 4862.68,
        'Hg'       : 4340.464,
        'Hd'       : 4101.734,
        'OIII-b'       :  5006.843,
        'OIII-a'       : 4958.911,
        'MgII'    : 2799.49,
        'OII'         : 3728,
        'CIII'  : 1909.,
        'CIV'    : 1549.06,
        'SiIV'  : 1393.76018,
        'LYA'         : 1215.67,
        'LYB'         : 1025.72
    }
    #
    redshifts=read_table(redrock_fn, "REDSHIFTS")
    model_flux=dict()
    line="TARGETID={}".format(elg_tid)
    if redshifts is not None :
        j=np.where(redshifts["TARGETID"]==elg_tid)[0][0]
        line += " Z={} SPECTYPE={} ZWARN={}".format(redshifts["Z"][j],redshifts["SPECTYPE"][j],redshifts["ZWARN"][j])
        zval=redshifts["Z"][j]
        tx = templates[(redshifts['SPECTYPE'][j], redshifts['SUBTYPE'][j])]
        for band in spec.bands:
            model_flux[band] = np.zeros(spec.wave[band].shape)
            coeff = redshifts['COEFF'][j][0:tx.nbasis]
            model = tx.flux.T.dot(coeff).T
            mx = resample_flux(spec.wave[band], tx.wave*(1+redshifts['Z'][j]), model)
            k=np.where(spec.fibermap["TARGETID"]==elg_tid)[0][0]
            model_flux[band] = spec.R[band][k].dot(mx)
    #
    fig, ax = plt.subplots()
    for b in spec._bands :
        i=np.where(spec.ivar[b][j]*(spec.mask[b][j]==0)>1./100.**2)[0]
        if rebin is not None and rebin>0:
            rwave=np.linspace(spec.wave[b][0],spec.wave[b][-1],spec.wave[b].size//rebin)
            rflux,rivar = resample_flux(rwave,spec.wave[b],spec.flux[b][j],ivar=spec.ivar[b][j]*(spec.mask[b][j]==0))
        else:
            rwave = spec.wave[b][i]
            rflux = spec.flux[b][j, i]
            rivar = spec.ivar[b][j, i]
        """
        if args.errors:
            plt.fill_between(wavescale*rwave, rflux-1./np.sqrt(rivar),
                             rflux+1./np.sqrt(rivar), alpha=0.5)
            line, = plt.plot(wavescale*rwave, rflux)
            plt.plot(wavescale*rwave, 1/np.sqrt(rivar)-2,
                     color=line.get_color(), linestyle='--')
        else:
            plt.plot(wavescale*rwave, rflux)
        """
        ax.plot(wavescale*rwave, rflux, label="{}-camera".format(b))
    # model
    for band in spec.bands:
        ax.plot(wavescale*spec.wave[band],model_flux[band],"-",alpha=0.6, color="k")
        for elem in lines :
            line=(1+zval)*lines[elem]
            if line>spec.wave[band][0] and line<spec.wave[band][-1] :
                ax.axvline(wavescale*line,color="k",linestyle="--",alpha=0.4)
                y=np.interp(wavescale*line,wavescale*spec.wave[band],model_flux[band])
                ax.text(wavescale*(line+60),ylim[0]+0.05 * (ylim[1] - ylim[0]),elem.split("-")[0],color="k")
    ax.set_title("ELG TARGETID={} (FIBER={}, Z={:.4f})".format(elg_tid, elg_fiber, zval))
    ax.set_xlabel(r"Wavelength [$\AA$]")
    ax.set_ylabel(r"Flux [$10^{-17}$ erg.s$^{-1}$.cm$^{-2}$.$\AA^{-1}$]")
    ax.set_xlim(3500, 10000)
    ax.set_ylim(ylim)
    ax.grid()
    ax.legend(loc=1)
    plt.savefig(os.path.join(outdir, "qa-elgspectrum.{}".format(fig_fmt)), bbox_inches="tight")
    plt.close()

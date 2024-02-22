#!/usr/bin/env python
# coding: utf-8

# ## Covariance from a simulated spectrum
# Writes extraction_covcorr.pdf and extraction_pull.pdf
#
# Re-extract N>>1 times with different noise realizations

#- non-DESI imports
import os, sys, glob, time

### get_ipython().run_line_magic('matplotlib', 'inline')
import matplotlib.pyplot as plt
import numpy as np
import fitsio

from astropy.visualization import ZScaleInterval

plt.ion();

#- DESI imports
import desispec.io
from desispec.io import findfile, read_image, read_frame
from desispec.interpolation import resample_flux
from desispec.resolution import Resolution
from specter.psf import load_psf
from specter.extract import ex2d


specprod = 'fuji'
os.environ['SPECPROD'] = specprod
reduxdir = os.path.expandvars(f'$DESI_ROOT/spectro/redux/{specprod}/preproc/')
night = 20210101
expid = 70350
camera = 'r1'

dw = 0.8
ww = 7300 + np.arange(50)*dw

print(f'Reading PSF for {night=} {expid=} {camera=}')
psf = load_psf(findfile('psf', night, expid, camera))

nspec = 1
nwave = len(ww)
# inphot = 10*np.ones((nspec, nwave))
inphot = np.tile(np.linspace(1,1000,nwave), nspec).reshape(nspec, nwave)
inphot[0,15] += 500
inphot[0,25] += 1000
inphot[0,35] += 10000

xyrange = xmin, xmax, ymin, ymax = psf.xyrange((0,nspec), ww)
simpix = psf.project(ww, inphot, xyrange=xyrange)


np.random.seed(0)
readnoise = 3.0
imgivar = 1/(simpix + readnoise**2)

np.random.seed(0)
multiflux = list()
multiivar = list()
multirdat = list()
ntest = 10000
#ntest = 10
t0 = time.time()
print(f'Simulating {ntest} extractions with different noise (this is slow)')
for i in range(ntest):
    if i%100 == 0 : print(f"{i+1}/{ntest}")
    simimg = np.random.poisson(simpix) + np.random.normal(scale=readnoise, size=simpix.shape)
    flux, ivar, rdat = ex2d(simimg, imgivar, psf, specmin=0, nspec=1, wavelengths=ww, xyrange=xyrange)
    multiflux.append(flux.ravel())
    multiivar.append(ivar.ravel())
    multirdat.append(rdat)

multiflux = np.array(multiflux)
multiivar = np.array(multiivar)

dt = time.time() - t0
print(f'{dt/60:.1f} minutes to process {ntest} extractions')


cov = np.cov(multiflux, rowvar=False)
corr = np.corrcoef(multiflux, rowvar=False)

#- detail: ivar only depends upon input pixel variance, not the noisy signal,
#- so every row of multiivar is actually the same; take the mean anyway in case
#- we update how to do the sims...
meanvar = np.mean(1/multiivar, axis=0)


print('Making plots')
fig=plt.figure(figsize=(5,3.5))

extent=(ww[0]-dw/2,ww[-1]-dw/2,ww[0]-dw/2,ww[-1]-dw/2)
ax = plt.subplot(111)
toto=ax.imshow(corr, vmin=-1, vmax=1, cmap='RdBu', extent=extent, origin="lower")
ax.set_title('Correlation')
ax.set_xlabel(r'Wavelength [$\AA$]')
ax.set_ylabel(r'Wavelength [$\AA$]')
#ax.axes.yaxis.set_ticklabels([])
plt.colorbar(toto)
plt.tight_layout()
fig.savefig('extraction_covcorr_top.pdf')

fig=plt.figure(figsize=(6,5))
ax = plt.subplot(311)
plt.plot(ww, np.diag(cov), label='observed variance')
plt.plot(ww, meanvar, '--', label='calculated variance')
plt.ylabel('Variance')
# plt.xlabel(r'Wavelength [$\AA$]')
plt.title(f'Extracted flux variance of {ntest} noise realizations')
ax.axes.xaxis.set_ticklabels([])
plt.legend()

alpha=0.8

ax = plt.subplot(312)
deltacov = (np.diag(cov) - meanvar)
plt.plot(ww, deltacov)
errvar = np.sqrt(2*meanvar**2 / (2*ntest))  # error on the variance, maybe
plt.plot(ww, errvar, color='k', alpha=alpha)
plt.plot(ww, -errvar, color='k', alpha=alpha)
plt.axhline(0, color='k', alpha=alpha)
plt.ylim(-100, 100)
plt.ylabel(r'$\Delta$ Variance')
ax.axes.xaxis.set_ticklabels([])

ax = plt.subplot(313)
deltacov = (np.diag(cov) - np.mean(1/multiivar, axis=0)) / np.diag(cov)
plt.plot(ww, deltacov)
plt.axhline(0, color='k', alpha=alpha)
plt.ylim(-0.05, 0.05)
plt.ylabel(r'$\Delta$(var) / var')
plt.xlabel(r'Wavelength [$\AA$]')

plt.tight_layout()

fig.savefig('extraction_covcorr_bottom.pdf')

print('Wrote extraction_covcorr_top.pdf and extraction_covcorr_bottom.pdf')

# ## Pull of extractions on an exactly aligned wavelength grid

# In[17]:


R = Resolution(rdat[0])
model = R.dot(inphot[0])/dw


# In[18]:


pull = (multiflux - model)*np.sqrt(multiivar)
np.mean(pull), np.std(pull), 1.4828*np.median(np.abs(pull - np.median(pull)))


# In[19]:


# plt.figure(figsize=(12,4))
#plt.figure(figsize=(8,4))
fig=plt.figure(figsize=(5,5))
ax = plt.subplot(311)
for i in range(min(100,ntest)):
    label = '2D extraction' if i==0 else None
    plt.plot(ww, multiflux[i], '-', color="C0", lw=1, alpha=alpha, label=label)

plt.plot(ww, model, '--', color="C1", label='Resolution convolved input')
# plt.plot(ww, rowflux, '-', color='C0', label='row-by-row')
plt.ylabel(r'e/$\AA$')
plt.legend(loc='upper left')
ax.axes.xaxis.set_ticklabels([])

#- log(flux)
ax = plt.subplot(312)
for i in range(min(100,ntest)):
    label = '2D extraction' if i==0 else None
    y = multiflux[i].copy()
    y[y<0] = np.NaN
    plt.plot(ww, np.log10(y), '-', color="C0", lw=1, alpha=alpha, label=label)

plt.plot(ww, np.log10(model), '--', color="C1", label='Resolution convolved input')
plt.ylabel(r'$\log_{10}$(e/$\AA$)')
ax.axes.xaxis.set_ticklabels([])

#- pull mean and std vs. wavelength
ax = plt.subplot(313)
plt.plot(ww, np.mean(pull, axis=0))
plt.plot(ww, np.std(pull, axis=0))
plt.axhline(0.0, linestyle=":", color='k', alpha=alpha)
plt.axhline(1.0, linestyle=":", color='k', alpha=alpha)
plt.text(ww[0], 0.1, 'mean(pull)')
plt.text(ww[0], 1.1, 'stddev(pull)')
plt.xlabel(r'Wavelength [$\AA$]')
plt.ylim(-0.2, 1.4)
plt.tight_layout()
fig.savefig('extraction_pull_top.pdf')


#- histogram pull
fig=plt.figure(figsize=(4,4))
ax = plt.subplot(111)
plt.hist(pull.ravel(), 50, (-5, 5))
plt.xlabel(r'pull = (flux - model)/$\sigma$')
xx = np.linspace(-5, 5)
yy = np.exp(-xx**2/2) * pull.size/np.sqrt(2*np.pi) / 5
plt.plot(xx, yy, '-',color="C1")
chimean = np.mean(pull)
chistd = np.std(pull)
plt.text(0.05, 0.90, f'mean = {chimean:.3f}', transform=ax.transAxes)
plt.text(0.05, 0.83, f'stddev = {chistd:.3f}', transform=ax.transAxes);
plt.tight_layout()
fig.savefig('extraction_pull_bottom.pdf')

print('Wrote extraction_pull_top.pdf and extraction_pull_bottom.pdf')

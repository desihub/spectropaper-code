# How to generate plots related to spectral extraction

Note that all these scripts/notebooks generate their plots in
spectropaper/work/extraction/ .  When you are happy with the plots,
then you copy them by hand into their final location in spectropaper/figures/ .
This ensures that re-running the script doesn't accidentally alter the
reference figure in git.

## Extraction plots

# Run the Extraction.ipynb notebook, then

Run the extraction plotting scripts then copy the files:

```
python plotex_ccd_patch.py
python plotex_wavelengths.py
python plotex_full_vs_patch_bias.py
python plotex_covcorr_pull.py

cp extraction*.pdf ../../figures
```

Note that `plotex_covcorr_pull.py` takes several minutes to run because it is
simulating N>>1 extractions to evaluate the covariance.

## R normalization bias plot

```
python rnorm.py
cp rnorm-bias.pdf ../../figures
```

## Interpolating the resolution matrix plots for appendix B

```
python plot_interpresolution.py
cp exbias*.pdf ../../figures
```


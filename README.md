assorted calculations for TESS mission extension

## `/src/`

* `visualize_survey_designs.py`: make maps (in Mollweide projection) of
  output from `tessmaps`, for possible pointing strategies of extended mission.

* `planet_yield_plots.py`: make plots (mostly bar charts) to assess the output
  of the [@mrtommyb](https://github.com/mrtommyb) detection simulations, at
  <https://github.com/mrtommyb/textended>.

* `period_barchart.py`: makes plot comparing "C3PO", "hemi", and "allsky"
  period distribution from the Huang et al (2018) detection simulation

* `ing_egr_snr.py`: run a toy ingress/egress calculation to understand what
  fraction of stars short cadence helps for. (it's a small fraction).

* `fraction_of_sky_observed.py`: calculate the fraction of the sky observed by
  different extended mission pointings.

## `/results/`

* `/visualize_survey_designs/`: maps of different sky-tiling strategies.

* `/planet_yield_plots/`: output from the @mrtommyb detection simulations.

## `/visualization_ext_fields/`

Old code from the Bouma et al (2017) extended mission white paper used to
generate 3d plots of different extended mission sky-tiling strategies.

# -*- coding: utf-8 -*-
"""
usage:
    python fraction_of_sky_observed.py >> ../results/fraction_of_sky_observed.txt
"""
from __future__ import division, print_function

import numpy as np, pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

from tessmaps.get_time_on_silicon import \
        given_cameras_get_stars_on_silicon as gcgss

from astropy import units as u, constants as const
from astropy.coordinates import SkyCoord

from numpy import array as nparr

import os
from glob import glob

if __name__=="__main__":

    fnames = [
        '../data/em2_v00_coords_observed_forproposal_S1_S26.csv',
        '../data/em2_v00_coords_observed_forproposal_S1_S56.csv',
        '../data/em2_v00_coords_observed_forproposal_S1_S97.csv',
    ]

    print(10*'-')
    for fname in fnames:

        df = pd.read_csv(fname, sep=';')

        obsd = (nparr(df['n_observations']) > 0)

        n_stars = len(df)
        n_obsd = len(df[obsd])

        print('{:s} --- {:.2f}% of sky observed'.
              format(os.path.basename(fname), 100*n_obsd/n_stars))

    pri = pd.read_csv('../data/em2_v00_coords_observed_forproposal_S1_S26.csv', sep=';')
    em1 = pd.read_csv('../data/em2_v00_coords_observed_forproposal_S27_S56.csv', sep=';')
    em2 = pd.read_csv('../data/em2_v00_coords_observed_forproposal_S57_S97.csv', sep=';')

    # One of our Primary Mission Objectives (PMOs) for EM1 was "Re-observe and
    # double the time coverage for 80% of the Prime Mission fields."  Did we do
    # this?  What fraction of the sky that was observed during the Prime
    # Mission was (or is scheduled to be) re-observed during EM1?

    sel = em1[em1.n_observations>0].index.isin(pri[pri.n_observations>0].index)

    n_pri = len(pri[pri.n_observations>0])
    n_reobs_em1 = len(em1[em1.n_observations>0][sel])

    print(10*'-')
    q1 = (
        'What fraction of the sky that was observed during the Prime '
        'Mission was (or is scheduled to be) re-observed during EM1?\n'
    )
    print('{:s}A: {:.2f}% of sky that was observed during Prime Mission '
          'has been or is scheduled to be re-observed during EM1'.
          format(q1, 100*n_reobs_em1/n_pri))

    # Another PMO was "Observe 70% of the portion of the sky that was not
    # observed during the Prime Mission, including the ecliptic."  Did we do
    # this?  What fraction of the sky that was not observed during the PM was
    # or will be observed during EM1?

    sel = em1[em1.n_observations>0].index.isin(pri[pri.n_observations==0].index)

    n_not_obs_pri = len(pri[pri.n_observations==0])
    n_obs_em1_but_not_pri = len(em1[em1.n_observations>0][sel])

    print(10*'-')
    q2 = (
        'What fraction of the sky that was not observed during the '
        'PM was or will be observed during EM1?\n'
    )
    print('{:s}A: {:.2f}% of sky that was not observed during Prime Mission '
          'has been or is scheduled to be re-observed during EM1'.
          format(q2, 100*n_obs_em1_but_not_pri/n_not_obs_pri))


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
        '../data/em2_v01_coords_observed_forproposal_S1_S26.csv',
        '../data/em2_v01_coords_observed_forproposal_S1_S55.csv',
        '../data/em2_v01_coords_observed_forproposal_S1_S97.csv',
        '../data/em2_v01_coords_observed_forproposal_S1_S123.csv',
    ]

    print(10*'-')
    for fname in fnames:

        df = pd.read_csv(fname, sep=';')

        obsd = (nparr(df['n_observations']) > 0)

        n_stars = len(df)
        n_obsd = len(df[obsd])

        print('{:s} --- {:.2f}% of sky observed'.
              format(os.path.basename(fname), 100*n_obsd/n_stars))

    pri = pd.read_csv('../data/em2_v01_coords_observed_forproposal_S1_S26.csv', sep=';')
    em1 = pd.read_csv('../data/em2_v01_coords_observed_forproposal_S27_S55.csv', sep=';')
    em2 = pd.read_csv('../data/em2_v01_coords_observed_forproposal_S56_S97.csv', sep=';')
    em3 = pd.read_csv('../data/em2_v01_coords_observed_forproposal_S98_S123.csv', sep=';')

    # One of our Primary Mission Objectives (PMOs) for EM1 was "Re-observe and
    # double the time coverage for 80% of the Prime Mission fields."  Did we do
    # this?  What fraction of the sky that was observed during the Prime
    # Mission was (or is scheduled to be) re-observed during EM1?

    sel = (
        (em1.n_observations > 0)
        &
        (pri.n_observations > 0)
    )

    n_pri = len(pri[pri.n_observations>0])
    n_reobs_em1 = len(em1[sel])

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

    sel = (
        (em1.n_observations > 0)
        &
        (pri.n_observations == 0)
    )

    n_not_obs_pri = len(pri[pri.n_observations==0])
    n_obs_em1_but_not_pri = len(em1[sel])

    print(10*'-')
    q2 = (
        'What fraction of the sky that was not observed during the '
        'PM was or will be observed during EM1?\n'
    )
    print('{:s}A: {:.2f}% of sky that was not observed during Prime Mission '
          'has been or is scheduled to be re-observed during EM1'.
          format(q2, 100*n_obs_em1_but_not_pri/n_not_obs_pri))
    print(10*'-')


    ##########################################
    q3 = (
        'What fraction of the sky has been observed twice with at least 12 '
        'months separation?'
    )
    n_12mo_pri = len(pri[pri.n_observations>=12])

    # in EM1: means reobserved, or n_obs>=12 in primary.
    sel = (
        ((em1.n_observations > 0) & (pri.n_observations > 0)) # reobserved
        |
        (pri.n_observations >= 12) # n_obs>=12 in primary
    )
    n_12mo_em1 = len(em1[sel])

    # in EM2: means reobserved from either EM1 or EM2, or n_obs>=12 in primary.
    sel = (
        ((em2.n_observations > 0) & (pri.n_observations > 0)) # reobserved pri->EM2
        |
        ((em2.n_observations > 0) & (em1.n_observations > 0)) # reobserved EM1->EM2
        |
        ((em1.n_observations > 0) & (pri.n_observations > 0)) # reobserved pri->EM1
        |
        (pri.n_observations >= 12) # n_obs>=12 in primary
    )
    n_12mo_em2 = len(em1[sel])

    # in EM3: means reobserved from either EM1 or EM2, or n_obs>=12 in primary.
    sel = (
        ((em3.n_observations > 0) & (pri.n_observations > 0)) # reobserved pri->EM3
        |
        ((em3.n_observations > 0) & (em1.n_observations > 0)) # reobserved EM1->EM3
        |
        ((em3.n_observations > 0) & (em2.n_observations > 0)) # reobserved EM2->EM3
        |
        ((em2.n_observations > 0) & (pri.n_observations > 0)) # reobserved pri->EM2
        |
        ((em2.n_observations > 0) & (em1.n_observations > 0)) # reobserved EM1->EM2
        |
        ((em1.n_observations > 0) & (pri.n_observations > 0)) # reobserved pri->EM1
        |
        (pri.n_observations >= 12) # n_obs>=12 in primary
    )
    n_12mo_em3 = len(em1[sel])


    print(q3)
    print(
        "A:\n"
        f"Primary mission S1-S26: {100*n_12mo_pri/n_stars:.2f}% \n"
        f"Primary+EM1 S1-S55: {100*n_12mo_em1/n_stars:.2f}% \n"
        f"Primary+EM1+EM2 S1-S97: {100*n_12mo_em2/n_stars:.2f}% \n"
        f"Primary+EM1+EM2+2year S1-S123: {100*n_12mo_em3/n_stars:.2f}% \n"
    )

    ##########################################
    q4 = (
        'What fraction of the sky has been observed at least once?'
    )

    sel = (
        (pri.n_observations > 0)
    )
    n_pri = len(pri[sel])

    sel = (
        (em1.n_observations > 0)
        |
        (pri.n_observations > 0)
    )
    n_em1 = len(em1[sel])

    sel = (
        (em2.n_observations > 0)
        |
        (em1.n_observations > 0)
        |
        (pri.n_observations > 0)
    )
    n_em2 = len(em2[sel])

    sel = (
        (em3.n_observations > 0)
        |
        (em2.n_observations > 0)
        |
        (em1.n_observations > 0)
        |
        (pri.n_observations > 0)
    )
    n_em3 = len(em3[sel])

    print(q4)
    print(
        "A:\n"
        f"Primary mission S1-S26: {100*n_pri/n_stars:.2f}% \n"
        f"Primary+EM1 S1-S55: {100*n_em1/n_stars:.2f}% \n"
        f"Primary+EM1+EM2 S1-S97: {100*n_em2/n_stars:.2f}% \n"
        f"Primary+EM1+EM2+2year S1-S123: {100*n_em3/n_stars:.2f}% \n"
    )



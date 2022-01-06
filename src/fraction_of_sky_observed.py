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

def main(versionstr=None):
    # e.g., "v01" for versionstr.
    assert versionstr[0] == 'v'
    assert len(versionstr) == 3

    # fnames = [
    #     f'../data/em2_{versionstr}_coords_observed_forproposal_S1_S26.csv',
    #     f'../data/em2_{versionstr}_coords_observed_forproposal_S1_S55.csv',
    #     f'../data/em2_{versionstr}_coords_observed_forproposal_S1_S97.csv',
    #     f'../data/em2_{versionstr}_coords_observed_forproposal_S1_S123.csv',
    # ]

    EM1end = 'S56' #if versionstr=='v00' else 'S55'
    EM2start = 'S57' #if versionstr=='v00' else 'S56'

    df_cum_em1 = pd.read_csv(f'../data/em2_{versionstr}_coords_observed_forproposal_S1_{EM1end}.csv')
    df_cum_em2 = pd.read_csv(f'../data/em2_{versionstr}_coords_observed_forproposal_S1_S97.csv')
    if versionstr != 'v00':
        df_cum_em3 = pd.read_csv(f'../data/em2_{versionstr}_coords_observed_forproposal_S1_S123.csv')

    n_stars = len(df_cum_em1)

    # print(10*'-')
    # for fname in fnames:

    #     df = pd.read_csv(fname, sep=';')

    #     obsd = (nparr(df['n_observations']) > 0)

    #     n_stars = len(df)
    #     n_obsd = len(df[obsd])

    #     namestr = (
    #         os.path.basename(fname).replace('_coords_observed_forproposal_','_')
    #     )
    #     print('{:s} --- {:.2f}% of sky observed for >=1 sector'.
    #           format(namestr, 100*n_obsd/n_stars))

    #     obsd = (nparr(df['n_observations']) >= 12)
    #     n_obsd = len(df[obsd])
    #     print('{:s} --- {:.2f}% of sky observed for >= 12 sectors'.
    #           format(namestr, 100*n_obsd/n_stars))

    #     obsd = (nparr(df['n_observations']) >= 24)
    #     n_obsd = len(df[obsd])
    #     print('{:s} --- {:.2f}% of sky observed for >= 24 sectors'.
    #           format(namestr, 100*n_obsd/n_stars))

    #     obsd = (nparr(df['n_observations']) >= 26)
    #     n_obsd = len(df[obsd])
    #     print('{:s} --- {:.2f}% of sky observed for >= 26 sectors'.
    #           format(namestr, 100*n_obsd/n_stars))

    #     obsd = (nparr(df['n_observations']) >= 36)
    #     n_obsd = len(df[obsd])
    #     print('{:s} --- {:.2f}% of sky observed for >= 36 sectors'.
    #           format(namestr, 100*n_obsd/n_stars))

    #     obsd = (nparr(df['n_observations']) >= 39)
    #     n_obsd = len(df[obsd])
    #     print('{:s} --- {:.2f}% of sky observed for >= 39 sectors'.
    #           format(namestr, 100*n_obsd/n_stars))

    #     print(10*'-')


    pri = pd.read_csv(f'../data/em2_{versionstr}_coords_observed_forproposal_S1_S26.csv', sep=';')
    em1 = pd.read_csv(f'../data/em2_{versionstr}_coords_observed_forproposal_S27_{EM1end}.csv', sep=';')
    em2 = pd.read_csv(f'../data/em2_{versionstr}_coords_observed_forproposal_{EM2start}_S97.csv', sep=';')
    try:
        em3 = pd.read_csv(f'../data/em2_{versionstr}_coords_observed_forproposal_S98_S123.csv', sep=';')
    except:
        em3 = None
    try:
        em1mid = pd.read_csv(f'../data/em2_{versionstr}_coords_observed_forproposal_S27_S47.csv', sep=';')
    except:
        em1mid = None

    # # One of our Primary Mission Objectives (PMOs) for EM1 was "Re-observe and
    # # double the time coverage for 80% of the Prime Mission fields."  Did we do
    # # this?  What fraction of the sky that was observed during the Prime
    # # Mission was (or is scheduled to be) re-observed during EM1?

    sel = (
        (em1.n_observations > 0)
        &
        (pri.n_observations > 0)
    )

    n_pri = len(pri[pri.n_observations>0])
    n_reobs_em1 = len(em1[sel])

    if em1mid is not None:
        sel = (
            (em1mid.n_observations > 0)
            &
            (pri.n_observations > 0)
        )
        n_reobs_em1mid = len(em1mid[sel])

    print(10*'-')
    q1 = (
        'What fraction of the sky that was observed during the Prime '
        'Mission was (or is scheduled to be) re-observed during EM1?\n'
    )
    print('{:s}A: {:.2f}% of sky that was observed during Prime Mission '
          'has been or is scheduled to be re-observed during EM1'.
          format(q1, 100*n_reobs_em1/n_pri))
    if em1mid is not None:
        print('A at S47: {:.2f}% of sky that was observed during Prime Mission '
              'has been reobserved.'.
              format(100*n_reobs_em1mid/n_pri))

    # # Another PMO was "Observe 70% of the portion of the sky that was not
    # # observed during the Prime Mission, including the ecliptic."  Did we do
    # # this?  What fraction of the sky that was not observed during the PM was
    # # or will be observed during EM1?

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
    q4 = (
        'Fraction of sky observed at least once'
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
        f"Primary+EM1 S1-{EM1end}: {100*n_em1/n_stars:.2f}% \n"
        f"Primary+EM1+EM2 S1-S97: {100*n_em2/n_stars:.2f}% \n"
        f"Primary+EM1+EM2+2year S1-S123: {100*n_em3/n_stars:.2f}% \n"
    )

    cum_sky_once_row = np.round([100*n_pri/n_stars, 100*n_em1/n_stars,
                        100*n_em2/n_stars, 100*n_em3/n_stars],2)

    ##########################################
    q3 = (
        'Fraction of the sky observed at least twice'
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
        f"Primary+EM1 S1-{EM1end}: {100*n_12mo_em1/n_stars:.2f}% \n"
        f"Primary+EM1+EM2 S1-S97: {100*n_12mo_em2/n_stars:.2f}% \n"
        f"Primary+EM1+EM2+2year S1-S123: {100*n_12mo_em3/n_stars:.2f}% \n"
    )

    cum_sky_twice_row = np.round([100*n_12mo_pri/n_stars, 100*n_12mo_em1/n_stars,
                         100*n_12mo_em2/n_stars, 100*n_12mo_em3/n_stars],2)

    ##########################################
    q5 = (
        'Fraction of the sky observed for >=19 sectors (>= 500 days)'
    )
    df_cum_em1 = pd.read_csv(f'../data/em2_{versionstr}_coords_observed_forproposal_S1_{EM1end}.csv', sep=';')
    df_cum_em2 = pd.read_csv(f'../data/em2_{versionstr}_coords_observed_forproposal_S1_S97.csv', sep=';')
    df_cum_em3 = pd.read_csv(f'../data/em2_{versionstr}_coords_observed_forproposal_S1_S123.csv', sep=';')

    NSECTOR = 19
    n_pri = len(pri[pri.n_observations>=NSECTOR])
    n_em1 = len(df_cum_em1[df_cum_em1.n_observations>=NSECTOR])
    n_em2 = len(df_cum_em2[df_cum_em2.n_observations>=NSECTOR])
    n_em3 = len(df_cum_em3[df_cum_em3.n_observations>=NSECTOR])

    print(q5)
    print(
        "A:\n"
        f"Primary mission S1-S26: {100*n_pri/n_stars:.2f}% \n"
        f"Primary+EM1 S1-{EM1end}: {100*n_em1/n_stars:.2f}% \n"
        f"Primary+EM1+EM2 S1-S97: {100*n_em2/n_stars:.2f}% \n"
        f"Primary+EM1+EM2+2year S1-S123: {100*n_em3/n_stars:.2f}% \n"
    )

    cum_19 = np.round([100*n_pri/n_stars, 100*n_em1/n_stars, 100*n_em2/n_stars,
              100*n_em3/n_stars],2)

    ##########################################
    q6 = (
        'Fraction of the sky observed for >=36 sectors (>= 1000 days)'
    )

    NSECTOR = 36
    n_pri = len(pri[pri.n_observations>=NSECTOR])
    n_em1 = len(df_cum_em1[df_cum_em1.n_observations>=NSECTOR])
    n_em2 = len(df_cum_em2[df_cum_em2.n_observations>=NSECTOR])
    n_em3 = len(df_cum_em3[df_cum_em3.n_observations>=NSECTOR])

    print(q6)
    print(
        "A:\n"
        f"Primary mission S1-S26: {100*n_pri/n_stars:.2f}% \n"
        f"Primary+EM1 S1-{EM1end}: {100*n_em1/n_stars:.2f}% \n"
        f"Primary+EM1+EM2 S1-S97: {100*n_em2/n_stars:.2f}% \n"
        f"Primary+EM1+EM2+2year S1-S123: {100*n_em3/n_stars:.2f}% \n"
    )

    cum_36 = np.round([100*n_pri/n_stars, 100*n_em1/n_stars, 100*n_em2/n_stars,
              100*n_em3/n_stars],2)

    data = np.array([cum_sky_once_row, cum_sky_twice_row, cum_19, cum_36])

    cols = 'Prime,EM1,EM2,EM3'.split(',')
    index = [q4,q3,q5,q6]

    df = pd.DataFrame(
        {k:v for k,v in zip(cols, data.T)}, index=index
    )
    print(df)
    outpath = f'../results/fraction_of_sky/{versionstr}.csv'
    df.to_csv(outpath)

    outpath = f'../results/fraction_of_sky/{versionstr}_EM2_EM3_only.csv'
    df[['EM2','EM3']].to_csv(outpath)




if __name__=="__main__":
    # main('v05')
    # main('v06')
    #main('v07')
    main('v09')
    #main('v00')

    # main('v01')
    # main('v02')
    # main('v03')
    # main('v04')

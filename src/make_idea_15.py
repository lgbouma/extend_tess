# -*- coding: utf-8 -*-
"""
NOTE:
    Sun Apr 19 09:32:51 2020
    this "shift" was never executed, because doing things with the "C3PO"-style
    true-northern pointings messes up the naive "latitude shift" implemented
    below.

----------

* shift everything "up" 1 degree in latitude towards the pole, except when
observing ecliptic

* shift everything forward 6 degrees in longitude, except when observing
ecliptic
"""
from __future__ import division, print_function

import numpy as np, pandas as pd
import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt

from numpy import array as nparr

def main():

    df = pd.read_csv("../data/idea_15_final_truenorth.csv", sep=";")
    raise NotImplementedError

    lon_shift = 6 # degrees
    lat_shift = 0 # degrees  # NOTE: a hack because otherwise breaks 24 deg separation

    scb_elat = nparr(df['spacecraft_boresight_elat'])
    scb_elon = nparr(df['spacecraft_boresight_elon'])
    sc_eroll = nparr(df['spacecraft_eroll'])

    c1_elat = nparr(df['cam1_elat'])
    c2_elat = nparr(df['cam2_elat'])
    c3_elat = nparr(df['cam3_elat'])
    c4_elat = nparr(df['cam4_elat'])

    c1_elon = nparr(df['cam1_elon'])
    c2_elon = nparr(df['cam2_elon'])
    c3_elon = nparr(df['cam3_elon'])
    c4_elon = nparr(df['cam4_elon'])

    # shift everything "up" 1 degree in latitude towards the pole, except when
    # observing ecliptic
    is_ecl = (sc_eroll != 0)
    is_north = (scb_elat > 0)
    is_south = (scb_elat < 0)

    scb_elat[~is_ecl & is_north] += lat_shift
    scb_elat[~is_ecl & is_south] -= lat_shift
    c1_elat[~is_ecl & is_north] += lat_shift
    c1_elat[~is_ecl & is_south] -= lat_shift
    c2_elat[~is_ecl & is_north] += lat_shift
    c2_elat[~is_ecl & is_south] -= lat_shift
    c3_elat[~is_ecl & is_north] += lat_shift
    c3_elat[~is_ecl & is_south] -= lat_shift

    # dumb way of dealing with the poles
    c4_elat[~is_ecl & is_north] -= lat_shift
    c4_elat[~is_ecl & is_south] += lat_shift

    # shift everything forward 6 degrees in longitude, except when observing
    # ecliptic
    scb_elon[~is_ecl] += lon_shift
    c1_elon[~is_ecl] += lon_shift
    c2_elon[~is_ecl] += lon_shift
    c3_elon[~is_ecl] += lon_shift
    c4_elon[~is_ecl] += lon_shift

    # because of the latitude shift
    c4_elon[~is_ecl] += 180

    outdf = pd.DataFrame({
        'start': nparr(df['start']),
        'end': nparr(df['end']),
        'orbit': nparr(df['sector']),
        'spacecraft_boresight_elat': scb_elat,
        'spacecraft_boresight_elon': np.mod(scb_elon,360),
        'spacecraft_eroll': sc_eroll,
        'cam4_elat': c4_elat,
        'cam4_elon': np.mod(c4_elon,360),
        'cam3_elat': c3_elat,
        'cam3_elon': np.mod(c3_elon,360),
        'cam2_elat': c2_elat,
        'cam2_elon': np.mod(c2_elon,360),
        'cam1_elat': c1_elat,
        'cam1_elon': np.mod(c1_elon,360),
        'comment': nparr(df['comment'])
    })

    outcol = ['start', 'end', 'orbit', 'spacecraft_boresight_elat',
              'spacecraft_boresight_elon', 'spacecraft_eroll', 'cam4_elat',
              'cam4_elon', 'cam3_elat', 'cam3_elon', 'cam2_elat', 'cam2_elon',
              'cam1_elat', 'cam1_elon', 'comment']

    outpath = "../data/idea_15_final_truenorth_shifted.csv"
    outdf[outcol].to_csv(outpath, index=False, sep=";")
    print('saved {}'.format(outpath))

if __name__=="__main__":
    main()

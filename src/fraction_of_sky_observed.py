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

    fnames = np.sort(glob('../data/*_observed_merged.csv'))

    for fname in fnames:

        df = pd.read_csv(fname, sep=';')

        obsd = (nparr(df['n_observations']) > 0)

        n_stars = len(df)
        n_obsd = len(df[obsd])

        print('{:s} --- {:.2f}% of sky observed'.
              format(os.path.basename(fname), 100*n_obsd/n_stars))

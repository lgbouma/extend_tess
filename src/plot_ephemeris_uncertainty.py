# -*- coding: utf-8 -*-
"""
make plot that shows why we're observing one region of the ecliptic, but not
others.
"""
from __future__ import division, print_function

import numpy as np, pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

from numpy import array as nparr
import os
from glob import glob
from parse import parse
from shutil import copyfile

from astropy.time import Time
from astropy import units as u

from scipy.stats import linregress
from scipy.optimize import curve_fit

from matplotlib.ticker import MaxNLocator
import matplotlib.transforms as transforms

from astrobase.timeutils import get_epochs_given_midtimes_and_period

def linear_model(xdata, m, b):
    return m*xdata + b


def plot_ephemeris_uncertainty(seed=42):

    primarystart = Time('2018-07-18')
    emstart = Time('2020-06-27')

    N_tra = 3 # in primary...
    sigma_tc_single_tra = (4*u.minute).to(u.hr).value
    period_yr = 10/365.25 # days
    ext_transit_epoch = 73 # 73 transits later... roughly 2 years
    _epoch = np.arange(0, 400, 1)

    # epochs [0,1,2,73] are observed
    primary_epochs = np.arange(0,N_tra)
    epochs = np.concatenate((np.atleast_1d(primary_epochs),
                             np.atleast_1d(np.array(ext_transit_epoch))))

    # true transit times of observed transits.
    true_transit_times = primarystart.byear + period_yr*epochs

    np.random.seed(seed)
    err_tc_hr = np.random.normal(0., sigma_tc_single_tra, size=N_tra+1)

    # add errors to true times. get observed transit times (in years)
    observed_transit_times = true_transit_times + err_tc_hr/(24*365.25)

    print((observed_transit_times-true_transit_times)*365.25*24*60)

    # fit t_mid_observed vs epoch. tmid = t_0 + P*E. The slope, "P" is in
    # years; so is std_err_period.
    (slope, intercept, r_value, p_value, std_err_period
    ) = linregress(epochs[:-1], observed_transit_times[:-1])

    # when refitting, you need to set the epoch zero-point to be the weighted
    # mid-time of the timeseries... note: this doesn't really change anything!
    # i did the fitting without this, and it changed nothing.
    epoch_shifted, _ = (
        get_epochs_given_midtimes_and_period(
            observed_transit_times, period_yr,
            err_t_mid=np.ones_like(epochs)*sigma_tc_single_tra/(24*365.25) ,
            verbose=True)
    )

    (slope_EM, intercept_EM, r_value_EM, p_value_EM, std_err_period_EM
    ) = linregress(epoch_shifted, observed_transit_times)

    ##########################################
    # EXTRA CHECK: manually checked that the period errors from these different
    # least squares routines agree to 1e-14 (years) in the "_EM" case, and
    # 1e-10 (years) in the prime-only case.
    popt, pcov = curve_fit(
        linear_model, epochs[:-1], observed_transit_times[:-1],
        p0=(period_yr, intercept), sigma=np.ones_like(epochs[:-1])*sigma_tc_single_tra/(24*365.25)
    )
    lsfit_period = popt[0]
    lsfit_period_err = pcov[0,0]**0.5
    lsfit_t0 = popt[1]
    lsfit_t0_err = pcov[1,1]**0.5

    popt, pcov = curve_fit(
        linear_model, epoch_shifted, observed_transit_times,
        p0=(period_yr, intercept), sigma=np.ones_like(epochs)*sigma_tc_single_tra/(24*365.25)
    )
    lsfit_period = popt[0]
    lsfit_period_err_EM = pcov[0,0]**0.5
    lsfit_t0 = popt[1]
    lsfit_t0_err_EM = pcov[1,1]**0.5
    # END EXTRA CHECK
    ##########################################

    f, ax = plt.subplots(figsize=(6,4))

    # true midtimes of transits at every epoch (regardless of observability)
    t_interp_byear = primarystart.byear + _epoch*period_yr

    # inferred midtimes from primary only. tmid = t_0 + period*epoch. note
    # "1_sig" in variable names is really 2 sigma.
    tmid_interp = intercept + slope*_epoch
    tmid_interp_plus_1_sig = intercept + (slope+1.96*std_err_period)*_epoch
    tmid_interp_minus_1_sig = intercept + (slope-1.96*std_err_period)*_epoch
    # inferred midtimes from extended + primary.
    tmid_interp_EM = intercept_EM + slope_EM*_epoch
    tmid_interp_plus_1_sig_EM = (
        intercept_EM + (slope_EM+1.96*std_err_period_EM)*_epoch
    )
    tmid_interp_minus_1_sig_EM = (
        intercept_EM + (slope_EM-1.96*std_err_period_EM)*_epoch
    )

    # width is the predicted tmid (+2sigma) minus the predicted tmid (-2sigma).
    prime = (
        tmid_interp_plus_1_sig - tmid_interp_minus_1_sig
    )*(365.25*24)

    extend = (
        tmid_interp_plus_1_sig_EM - tmid_interp_minus_1_sig_EM
    )*(365.25*24)

    ax.plot(t_interp_byear[prime>0], prime[prime>0], color='#1f77b4',
            zorder=-3, alpha=1, lw=3)
    sel = (extend>0) & (t_interp_byear > emstart.byear)
    ax.plot(t_interp_byear[sel], extend[sel], color='#ff7f0e',
            zorder=-2, alpha=1, lw=3)

    ax.set_xlabel('Year', fontsize='large')
    ax.set_ylabel('Transit mid-time uncertainty [hours]', fontsize='large')
    ax.set_yscale('log')

    ax.xaxis.set_major_locator(MaxNLocator(integer=True))

    ax.get_yaxis().set_tick_params(which='both', direction='in')
    ax.get_xaxis().set_tick_params(which='both', direction='in')

    ax.set_ylim((0.08,108))
    ax.set_xlim([2018,2029])

    if np.array_equal(
        ax.get_yticks(),
        np.array([1.e-03, 1.e-02, 1.e-01, 1.e+00, 1.e+01, 1.e+02, 1.e+03,
                  1.e+04])
    ):
        ax.set_yticklabels('0.001,0.01,0.1,1,10,100,1000,10000'.split(','))

    f.tight_layout()
    outpath = '../results/ephemeris_uncertainty.pdf'
    f.savefig(outpath, bbox_inches='tight')
    print('made {}'.format(outpath))


def get_ephemeris_uncertainty(seed=42):

    primarystart = Time('2018-07-18')
    emstart = Time('2020-06-27')

    N_tra = 3 # in primary...
    sigma_tc_single_tra = (4*u.minute).to(u.hr).value
    period_yr = 10/365.25 # days
    ext_transit_epoch = 73 # 73 transits later... roughly 2 years
    _epoch = np.arange(0, 400, 1)

    # epochs [0,1,2,73] are observed
    primary_epochs = np.arange(0,N_tra)
    epochs = np.concatenate((np.atleast_1d(primary_epochs),
                             np.atleast_1d(np.array(ext_transit_epoch))))

    # true transit times of observed transits.
    true_transit_times = primarystart.byear + period_yr*epochs

    np.random.seed(seed)
    err_tc_hr = np.random.normal(0., sigma_tc_single_tra, size=N_tra+1)

    # add errors to true times. get observed transit times (in years)
    observed_transit_times = true_transit_times + err_tc_hr/(24*365.25)

    print((observed_transit_times-true_transit_times)*365.25*24*60)

    # fit t_mid_observed vs epoch. tmid = t_0 + P*E. The slope, "P" is in
    # years; so is std_err_period.
    (slope, intercept, r_value, p_value, std_err_period
    ) = linregress(epochs[:-1], observed_transit_times[:-1])

    # when refitting, you need to set the epoch zero-point to be the weighted
    # mid-time of the timeseries... note: this doesn't really change anything!
    # i did the fitting without this, and it changed nothing.
    epoch_shifted, _ = (
        get_epochs_given_midtimes_and_period(
            observed_transit_times, period_yr,
            err_t_mid=np.ones_like(epochs)*sigma_tc_single_tra/(24*365.25) ,
            verbose=True)
    )

    (slope_EM, intercept_EM, r_value_EM, p_value_EM, std_err_period_EM
    ) = linregress(epoch_shifted, observed_transit_times)

    # true midtimes of transits at every epoch (regardless of observability)
    t_interp_byear = primarystart.byear + _epoch*period_yr

    # inferred midtimes from primary only. tmid = t_0 + period*epoch. note
    # "1_sig" in variable names is really 2 sigma.
    tmid_interp = intercept + slope*_epoch
    tmid_interp_plus_1_sig = intercept + (slope+1.96*std_err_period)*_epoch
    tmid_interp_minus_1_sig = intercept + (slope-1.96*std_err_period)*_epoch
    # inferred midtimes from extended + primary.
    tmid_interp_EM = intercept_EM + slope_EM*_epoch
    tmid_interp_plus_1_sig_EM = (
        intercept_EM + (slope_EM+1.96*std_err_period_EM)*_epoch
    )
    tmid_interp_minus_1_sig_EM = (
        intercept_EM + (slope_EM-1.96*std_err_period_EM)*_epoch
    )

    # width is the predicted tmid (+2sigma) minus the predicted tmid (-2sigma).
    prime = (
        tmid_interp_plus_1_sig - tmid_interp_minus_1_sig
    )*(365.25*24)

    extend = (
        tmid_interp_plus_1_sig_EM - tmid_interp_minus_1_sig_EM
    )*(365.25*24)

    sel = (extend>0) & (t_interp_byear > emstart.byear)

    return t_interp_byear[prime>0], prime[prime>0], t_interp_byear[sel], extend[sel]


def make_plot():

    f, ax = plt.subplots(figsize=(6,4))

    x1s, y1s, x2s, y2s = [], [], [], []
    for seed in np.arange(40,541):

        print(seed)
        x1,y1,x2,y2 = get_ephemeris_uncertainty(seed=seed)
        x1s.append(x1)
        y1s.append(y1)
        x2s.append(x2)
        y2s.append(y2)

        ax.plot(x1, y1, color='#1f77b4',
                zorder=-3, alpha=0.1, lw=0.5)
        ax.plot(x2, y2, color='#ff7f0e',
                zorder=-2, alpha=0.2, lw=0.5)

    # shape is nseeds x ntimes
    x1s = np.atleast_2d(x1s)
    y1s = np.atleast_2d(y1s)
    x2s = np.atleast_2d(x2s)
    y2s = np.atleast_2d(y2s)

    ax.plot(np.median(x1s,axis=0), np.median(y1s,axis=0), color='black',
            zorder=1, alpha=1, lw=4, ls='--')
    ax.plot(np.median(x2s,axis=0), np.median(y2s,axis=0), color='black',
            zorder=2, alpha=1, lw=4, ls='-')

    ax.set_xlabel('Year', fontsize='large')
    ax.set_ylabel('Transit mid-time uncertainty [hours]', fontsize='large')
    ax.set_yscale('log')

    ax.xaxis.set_major_locator(MaxNLocator(integer=True))

    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.get_yaxis().set_tick_params(which='both', direction='in')
    ax.get_xaxis().set_tick_params(which='both', direction='in')

    ax.set_ylim((0.08,108))
    ax.set_xlim([2018,2029])

    if np.array_equal(
        ax.get_yticks(),
        np.array([1.e-03, 1.e-02, 1.e-01, 1.e+00, 1.e+01, 1.e+02, 1.e+03,
                  1.e+04])
    ):
        ax.set_yticklabels('0.001,0.01,0.1,1,10,100,1000,10000'.split(','))

    f.tight_layout()
    outpath = '../results/ephemeris_uncertainty_randomseeds.pdf'
    f.savefig(outpath, bbox_inches='tight')
    print('made {}'.format(outpath))


def make_plot():

    f, ax = plt.subplots(figsize=(6,4))

    x1s, y1s, x2s, y2s = [], [], [], []
    for seed in np.arange(40,541):

        print(seed)
        x1,y1,x2,y2 = get_ephemeris_uncertainty(seed=seed)
        x1s.append(x1)
        y1s.append(y1)
        x2s.append(x2)
        y2s.append(y2)

        # ax.plot(x1, y1, color='#1f77b4',
        #         zorder=-3, alpha=0.1, lw=0.5)
        # ax.plot(x2, y2, color='#ff7f0e',
        #         zorder=-2, alpha=0.2, lw=0.5)

    # shape is nseeds x ntimes
    x1s = np.atleast_2d(x1s)
    y1s = np.atleast_2d(y1s)
    x2s = np.atleast_2d(x2s)
    y2s = np.atleast_2d(y2s)

    ax.plot(np.median(x1s,axis=0), np.median(y1s,axis=0), color='C0',
            zorder=1, alpha=1, lw=4, ls='-')
    ax.plot(np.median(x2s,axis=0), np.median(y2s,axis=0), color='C1',
            zorder=2, alpha=1, lw=4, ls='-')

    ax.set_xlabel('Year', fontsize='large')
    ax.set_ylabel('Transit mid-time uncertainty [hours]', fontsize='large')
    ax.set_yscale('log')

    ax.xaxis.set_major_locator(MaxNLocator(integer=True))

    ax.yaxis.set_ticks_position('both')
    ax.xaxis.set_ticks_position('both')
    ax.get_yaxis().set_tick_params(which='both', direction='in')
    ax.get_xaxis().set_tick_params(which='both', direction='in')

    ax.set_ylim((0.08,108))
    ax.set_xlim([2018,2029])

    if np.array_equal(
        ax.get_yticks(),
        np.array([1.e-03, 1.e-02, 1.e-01, 1.e+00, 1.e+01, 1.e+02, 1.e+03,
                  1.e+04])
    ):
        ax.set_yticklabels('0.001,0.01,0.1,1,10,100,1000,10000'.split(','))

    f.tight_layout()
    outpath = '../results/ephemeris_uncertainty_medians.pdf'
    f.savefig(outpath, bbox_inches='tight')
    print('made {}'.format(outpath))




if __name__=="__main__":

    seed = 60

    # plot_ephemeris_uncertainty(seed=seed)

    # make_all_showed_plot()

    make_plot()

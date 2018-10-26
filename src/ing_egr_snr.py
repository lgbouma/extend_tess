# -*- coding: utf-8 -*-
'''
given the TIC CTL, if you assume every star hosts a transiting, P=10 day
super-Earth, over a total of 5 years observing, what is the phase-folded SNR of
the ingress or egress phase?
'''
from __future__ import division, print_function

import numpy as np, pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from glob import glob
from astropy.coordinates import SkyCoord

from numpy import array as nparr

from tessmaps.get_data import make_prioritycut_ctl
from tessmaps import get_time_on_silicon as gts

def make_ctl_72():

    ctldir = '/Users/luke/local/TIC/CTL72/'

    subcols = ['RA', 'DEC', 'TESSMAG', 'TEFF', 'PRIORITY', 'RADIUS', 'MASS',
               'CONTRATIO', 'ECLONG', 'ECLAT', 'DIST', 'TICID', 'SPEC_LIST',
               'N_STAR', 'N_CONT', 'N_SKY', 'N_READ', 'N_SYS']

    make_prioritycut_ctl(datadir=ctldir,
                         prioritycut=0.0005,
                         subcols=subcols,
                         savpath='../data/TIC72_prioritycut.csv')

def get_nsectors_obsd():

    df = pd.read_csv('../data/TIC72_prioritycut.csv')

    coords = SkyCoord(nparr(df['RA']), nparr(df['DEC']), unit='deg')

    outdf = gts.get_time_on_silicon(coords)

    outname = '../data/TIC72_nsectors_obsd.csv'
    outdf.to_csv(outname, index=False)
    print('saved {:s}'.format(outname))


def merge_the_dfs():

    ticdf = pd.read_csv('../data/TIC72_prioritycut.csv')
    nsecdf = pd.read_csv('../data/TIC72_nsectors_obsd.csv')
    ndf = nsecdf[['ra', 'dec', 'elon', 'elat', 'total_sectors_obsd']]

    mdf = pd.merge(ticdf, ndf, left_on=['RA', 'DEC'], right_on=['ra', 'dec'])

    mdf[mdf['elat']<0].to_csv('../data/TIC72_obsd_southhemi.csv', index=False)


def calculate_ingress_snr_stats(rp=2, P_orb=10):
    '''
    assume every star hosts a transiting, P=10 day super-Earth, over a total of
    1 year observing, what is the phase-folded SNR of the ingress or egress
    phase?
    '''

    # from a single year of southern hemisphere observations
    df = pd.read_csv('../data/TIC72_obsd_southhemi.csv')

    from astropy import units as u, constants as c

    rp = rp*u.Rearth
    rstar = nparr(df['RADIUS'])*u.Rsun
    mstar = nparr(df['MASS'])*u.Msun

    # on average, the ingress has depth half the full depth
    signal = ((rp/rstar)**2).cgs.value / 2

    noise_1hr = np.sqrt(
        nparr(df['N_STAR'])**2 +
        nparr(df['N_CONT'])**2 +
        nparr(df['N_SKY'])**2 +
        nparr(df['N_READ'])**2 +
        nparr(df['N_SYS'])**2
    )

    # compute T_tot - T_full following Winn 2010 Eq 14 and 15.
    P_orb = P_orb*u.day
    k = (rp/rstar).cgs.value
    b = 0
    sini = 1

    a = ( P_orb**2 * c.G * mstar / (4*np.pi**2) )**(1/3)
    a = a.to(u.au)

    T_tot = P_orb/np.pi * np.arcsin(
        rstar/a * np.sqrt( (1+k)**2 - b**2  )/sini
    )
    T_full = P_orb/np.pi * np.arcsin(
        rstar/a * np.sqrt( (1-k)**2 - b**2  )/sini
    )
    T_ingr = (T_tot-T_full)/2
    T_ingr_hr = T_ingr.to(u.hr*u.rad).value

    # assume gaussian noise
    noise = noise_1hr / np.sqrt(T_ingr_hr)
    snr_ingr = signal/noise

    df['T_ingr_hr'] = T_ingr_hr
    df['snr_ingr'] = snr_ingr

    N_tra = nparr(df['total_sectors_obsd'])

    df['snr_ingr_pf'] = snr_ingr * np.sqrt(N_tra)

    savname = '../results/TIC72_ingress_snr.csv'
    df.to_csv(savname, index=False)
    print('saved {:s}'.format(savname))


def _TIC72_ingress_snr_distribution(x, yrstr='2yrs'):

    plt.style.use("dark_background")
    plt.close('all')
    f, ax = plt.subplots()

    bins = np.logspace(-1,2,7)
    n, _, patches = ax.hist(x, bins=bins)

    ax.set_yscale('log')
    ax.set_xscale('log')

    ax.set_xlabel('phased-folded SNR of ingress')
    ax.set_ylabel('number')
    ax.set_title(
        'Assume every star in TIC7.2 has P=10day, Rp=2Re transiting'
        '\nplanet. This is the distribution of ingress SNR ({:s})'.
        format(yrstr),
        fontsize='x-small'
    )

    f.tight_layout()
    f.savefig(
        '../results/TIC72_ingress_snr_distribution_{:s}.png'.format(yrstr),
        dpi=350
    )
    f.savefig(
        '../results/TIC72_ingress_snr_distribution_{:s}.pdf'.format(yrstr)
    )


def _TIC72_ingress_durn_distribution(x):

    plt.style.use("dark_background")
    plt.close('all')
    f, ax = plt.subplots()

    bins = np.logspace(0,2,9)
    n, _, patches = ax.hist(x, bins=bins)

    ax.set_yscale('log')
    ax.set_xscale('log')

    ax.set_xlabel('ingress duration [minute]')
    ax.set_ylabel('number')
    ax.set_title(
        'Assume every star in TIC7.2 has P=10day, Rp=2Re transiting'
        '\nplanet. This is the distribution of ingress duration.',
        fontsize='x-small'
    )

    f.tight_layout()
    f.savefig(
        '../results/TIC72_ingress_duration_distribution.png',
        dpi=350
    )
    f.savefig(
        '../results/TIC72_ingress_duration_distribution.pdf'
    )



def make_distribn_plots():

    df = pd.read_csv('../results/TIC72_ingress_snr.csv')

    # multiply by factor of two to get estimate for both hemispheres
    x = nparr(df['snr_ingr_pf'])
    x = np.repeat(x, 2)

    # anything that is not finite is probably because some of the stellar
    # parameters are missing.
    # NOTE this might be janky if the M dwarfs preferentially have missing
    # parameters...
    # FIXME FIXME FIXME
    x = x[(~np.isnan(x)) & (x>0)]

    _TIC72_ingress_snr_distribution(x, yrstr='2yrs')

    x = nparr(df['snr_ingr_pf']) * np.sqrt(2)
    x = x[(~np.isnan(x)) & (x>0)]
    x = np.repeat(x, 2)

    _TIC72_ingress_snr_distribution(x, yrstr='4yrs')

    T_ingr_hr = nparr(df['T_ingr_hr'])
    T_ingr_hr = T_ingr_hr[np.isfinite(T_ingr_hr)]

    _TIC72_ingress_durn_distribution(T_ingr_hr*60)


def tmag_v_teff_long_short_duration_ingress(T_ingr_split_min=8):

    df = pd.read_csv('../results/TIC72_ingress_snr.csv')

    x = nparr(df['TEFF'])
    y = nparr(df['TESSMAG'])
    T_ingr_hr = nparr(df['T_ingr_hr'])

    sel = np.isfinite(x) & np.isfinite(y) & np.isfinite(T_ingr_hr)

    x = x[sel]
    y = y[sel]
    T_ingr_hr = T_ingr_hr[sel]

    ########## 

    plt.close('all')
    plt.style.use("dark_background")
    f, ax = plt.subplots()

    ax.scatter(x[T_ingr_hr>=T_ingr_split_min/60],
               y[T_ingr_hr>=T_ingr_split_min/60],
               rasterized=True,
               label='T_ingr > {:d} minutes'.format(T_ingr_split_min),
               s=2, zorder=2, color='#E42E45')
    ax.scatter(x[T_ingr_hr<T_ingr_split_min/60],
               y[T_ingr_hr<T_ingr_split_min/60],
               rasterized=True,
               label='T_ingr < {:d} minutes'.format(T_ingr_split_min),
               s=2, zorder=1, color='#37B442')

    ax.legend(loc='best')

    ax.set_yscale('linear')
    ax.set_xscale('log')

    ax.set_xlim([1.2e4, 2e3])
    ax.set_ylim([17, 0])
    ax.tick_params(axis='x', which='minor')

    ax.set_xlabel('effective temperature [K]')
    ax.set_ylabel('T mag')
    ax.set_title(
        'Assume every star in TIC7.2 CTL has P=10day, Rp=2Re transiting'
        '\nplanet. Which have the longer ingress durations?',
        fontsize='x-small'
    )

    f.tight_layout()
    savname = (
        '../results/tmag_v_teff_long_short_duration_ingress_{:d}min.png'.
        format(T_ingr_split_min)
    )
    f.savefig(savname, dpi=350)
    print('saved {:s}'.format(savname))


def tmag_v_teff_high_low_ingress_snr(SNR_ingr_split=3):

    df = pd.read_csv('../results/TIC72_ingress_snr.csv')

    x = nparr(df['TEFF'])
    y = nparr(df['TESSMAG'])
    snr_ingr_pf = nparr(df['snr_ingr_pf'])

    sel = np.isfinite(x) & np.isfinite(y) & np.isfinite(snr_ingr_pf)

    x = x[sel]
    y = y[sel]
    snr_ingr_pf = snr_ingr_pf[sel]

    ########## 

    plt.close('all')
    plt.style.use("dark_background")
    f, ax = plt.subplots()

    ax.scatter(x[snr_ingr_pf>=SNR_ingr_split],
               y[snr_ingr_pf>=SNR_ingr_split],
               rasterized=True,
               label='SNR_ingr_pf > {:d}'.format(SNR_ingr_split),
               s=2, zorder=2, color='#E42E45')
    ax.scatter(x[snr_ingr_pf<SNR_ingr_split],
               y[snr_ingr_pf<SNR_ingr_split],
               rasterized=True,
               label='SNR_ingr_pf < {:d}'.format(SNR_ingr_split),
               s=2, zorder=1, color='#37B442')

    ax.legend(loc='best')

    ax.set_yscale('linear')
    ax.set_xscale('log')

    ax.set_xlim([1.2e4, 2e3])
    ax.set_ylim([17, 0])
    ax.tick_params(axis='x', which='minor')

    ax.set_xlabel('effective temperature [K]')
    ax.set_ylabel('T mag')
    ax.set_title(
        'If every star in TIC7.2 CTL has P=10day, Rp=2Re transiting'
        '\nplanet, which have the highest SNR during ingress/egress?',
        fontsize='small'
    )

    f.tight_layout()
    savname = (
        '../results/tmag_v_teff_high_low_ingress_snr_split{:d}.png'.
        format(SNR_ingr_split)
    )
    f.savefig(savname, dpi=350)
    print('saved {:s}'.format(savname))




if __name__ == "__main__":

    make_ctl = 0
    get_nsectors_obs = 0
    merge_dfs = 0
    calc_stats = 0
    make_distribn_plot = 1
    make_scatter_plot = 0

    if make_ctl:
        make_ctl_72()
    if get_nsectors_obs:
        get_nsectors_obsd()
    if merge_dfs:
        merge_the_dfs()
    if calc_stats:
        calculate_ingress_snr_stats()
    if make_distribn_plot:
        make_distribn_plots()
    if make_scatter_plot:
        tmag_v_teff_high_low_ingress_snr(SNR_ingr_split=2)
        tmag_v_teff_high_low_ingress_snr(SNR_ingr_split=3)
        tmag_v_teff_long_short_duration_ingress(T_ingr_split_min=6)
        tmag_v_teff_long_short_duration_ingress(T_ingr_split_min=8)
        tmag_v_teff_long_short_duration_ingress(T_ingr_split_min=10)

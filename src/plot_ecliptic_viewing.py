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

def make_data_dir():

    outdir='../data/vanderspek_ecliptic/'
    if not os.path.exists(outdir):
        os.mkdir(outdir)

    startdir = '../vanderspek_emi_eclpointing/'

    emangle_pngs = glob(os.path.join(startdir, 'emangles_S*_O*_d0_R90.png'))
    emi_pngs = glob(os.path.join(startdir, 'emi_S*_O*_d0_R90.png'))

    emangle_dats = glob(os.path.join(startdir, 'emangles_S*_O*_d0_R90.dat'))
    emi_dats = glob(os.path.join(startdir, 'emi_S*_O*_d0_R90.dat'))

    for fclass in [emangle_pngs, emi_pngs, emangle_dats, emi_dats]:

        for f in fclass:

            fmt = '{}_S{:d}_O{:d}_d0_R90.{}'
            res = parse(fmt,os.path.basename(f))

            dtype = res[0]
            sectornum = res[1]
            orbitnum = res[2]
            fext = res[3]

            outname = '{}_S{}_O{}_d0_R90.{}'.format(
                dtype, str(sectornum).zfill(3), str(orbitnum).zfill(3), fext
            )

            outpath = os.path.join(outdir, outname)

            copyfile(f, outpath)
            print('{}->{}'.format(f, outpath))


def calculate_efficiency_vs_time(df, proximity_cutoff_deg=40):
    """
    Given a critical separation angle, calculate "efficiency with which we
    would observe the ecliptic", where

        * efficiency = 1 if no body (earth/moon) is < the critical angle from
        any camera's boresight.

        * efficiency = 0.75 if one camera is knocked out, by either the moon,
        or the earth.

        * efficiency = 0.5 if two cameras are knocked out (ditto).

        * = 0.25 if three ...

    Call this the "em_efficiency" (earth-moon efficiency).

    Similarly, we can have the "earth efficiency", and the "moon efficiency"
    for each single object.
    """

    colnames = ['c1_e','c2_e','c3_e','c4_e', 'c1_m','c2_m','c3_m','c4_m']

    e_efficiency = (nparr(df['c1_e']>proximity_cutoff_deg).astype(int) +
                    nparr(df['c2_e']>proximity_cutoff_deg).astype(int) +
                    nparr(df['c3_e']>proximity_cutoff_deg).astype(int) +
                    nparr(df['c4_e']>proximity_cutoff_deg).astype(int)
                   )/4

    m_efficiency = (nparr(df['c1_m']>proximity_cutoff_deg).astype(int) +
                    nparr(df['c2_m']>proximity_cutoff_deg).astype(int) +
                    nparr(df['c3_m']>proximity_cutoff_deg).astype(int) +
                    nparr(df['c4_m']>proximity_cutoff_deg).astype(int)
                   )/4

    em_efficiency = (
            nparr( (df['c1_m']>proximity_cutoff_deg) & (df['c1_e']>proximity_cutoff_deg) ).astype(int) +
            nparr( (df['c2_m']>proximity_cutoff_deg) & (df['c2_e']>proximity_cutoff_deg) ).astype(int) +
            nparr( (df['c3_m']>proximity_cutoff_deg) & (df['c3_e']>proximity_cutoff_deg) ).astype(int) +
            nparr( (df['c4_m']>proximity_cutoff_deg) & (df['c4_e']>proximity_cutoff_deg) ).astype(int)
    )/4

    df['e_efficiency'] = e_efficiency
    df['m_efficiency'] = m_efficiency
    df['em_efficiency'] = em_efficiency


    return df


def make_timeseries_df(datadir, n_of_interest=85):
    """
    Make the DataFrame containing, as a function of time:
        * Earth to camera boresight angle
        * Moon to camera boresight angle
        * Sector number
        * Orbit number (assuming 1 orbit = 0.5 sectors).
        * "Efficiency".

    Kwargs:
        n_of_interest (int): number of files at which to cut, otherwise the
        timeseries Roland provided is a bit too long.
    """

    emi_dats = np.sort(glob(os.path.join(datadir, 'emi_S*_O*_d0_R90.dat'))
                      )[:n_of_interest]

    list_ = []
    n_entries = []

    sectornum = 1
    orbitnum= 9
    for emi_dat in emi_dats:
        df = pd.read_csv(emi_dat, skiprows=[0,1], names=[
            'sc_time', 'date', 'time', 'd_days', 'c1_e', 'c2_e', 'c3_e',
            'c4_e', 'c1_m', 'c2_m', 'c3_m', 'c4_m', 'e_dist', 'm_dist'],
            delim_whitespace=True)

        df['sectornum'] = sectornum

        orbitnums = np.ones_like(list(range(len(df))))*orbitnum
        orbitnums[int(np.floor(len(orbitnums)/2)):] += 1

        df['orbitnum'] = orbitnums

        list_.append(df)
        n_entries.append(len(df))

        sectornum += 1
        orbitnum += 2

    df = pd.concat(list_, axis=0, ignore_index=True)

    # ok, now make column of "efficiency" at every moment.
    df = calculate_efficiency_vs_time(df)

    return df


def plot_efficiency_vs_time_for_sectornum(df, sectornum=1):
    # look at efficiency vs time for one sector
    df = df[df['sectornum']==sectornum]

    plt.close('all')
    f,ax = plt.subplots(figsize=(12,3))

    ax.plot(df['sc_time'], df['em_efficiency'], label='EM efficiency',
            rasterized=True)
    ax.plot(df['sc_time'], df['e_efficiency'], label='E efficiency',
            rasterized=True)
    ax.plot(df['sc_time'], df['m_efficiency'], label='M efficiency',
            rasterized=True)

    ax.set_xlabel('spacecraft time [seconds] (two orbits shown)')
    ax.set_ylabel('efficiency')

    ax.legend(loc='best',fontsize='x-small')

    onum0 = 2*sectornum+7
    onum1 = onum0 + 1
    ax.set_title('orbits {} and {} (sector {}), ecliptic dec=0'.
                 format(onum0, onum1, sectornum))

    outdir = '../results/ecliptic_timing/efficiency_vs_time_by_sector/'
    outpath = os.path.join(
        outdir,'efficiency_vs_time_sector{}.png'.format(str(sectornum).zfill(3))
    )
    f.savefig(outpath,bbox_inches='tight')
    print('made {}'.format(outpath))



def plot_avgd_efficiency_vs_orbit_number_plot(df):

    ##########################################
    # calculate the average earth-moon efficiency in each orbit

    orbitnums = nparr(list(range(9,123)))

    em_efficiency = []
    for orbitnum in orbitnums:
        em_efficiency.append(
            np.mean(df[df['orbitnum']==orbitnum]['em_efficiency'])
        )
    em_efficiency = nparr(em_efficiency)

    ##########################################
    # get the anti-sun direction based on idea_2_SNSNS_hemi
    df_primary = pd.read_csv('../data/primary_mission.csv', sep=';')
    df_idea2 = pd.read_csv('../data/idea_2_SNSNS_hemi.csv', sep=';')
    df_dirns = pd.concat([df_primary, df_idea2], axis=0, ignore_index=True)
    df_dirns = df_dirns.sort_values(by='orbit')
    antisun_elon = nparr(df_dirns['spacecraft_boresight_elon'])

    tstart = Time(nparr(df_dirns['start']).astype(str), format='iso', scale='utc')
    half_orbit_days = (27.32/4) # days
    half_orbit_years = half_orbit_days / (365.25)
    tmid_orbit = tstart.byear + half_orbit_years

    ##########################################
    # make the plot

    f,ax = plt.subplots(figsize=(9,3))

    # cyclic options include: twilight, twilight_shifted, hsv
    cm = plt.cm.get_cmap('twilight_shifted')

    # make scatter points
    cs = ax.scatter(orbitnums, em_efficiency, c=antisun_elon, cmap=cm,
                    zorder=12)

    # make smoothed line...
    from scipy.signal import savgol_filter
    smooth_efficiency = savgol_filter(em_efficiency, 13, 3, mode='mirror')

    ax.plot(orbitnums, smooth_efficiency, lw=1, c='gray', zorder=2, alpha=0.5)

    cbar = f.colorbar(cs, pad=0.01, ticks=np.arange(0,360+60,60))
    cbar.ax.set_yticklabels(list(map(lambda x: str(x)+'$\!$$^\circ$',
                                     np.arange(0,360+60,60))))
    cbar.ax.tick_params(direction='in')

    ax.set_xlabel('Orbit number')
    ax.set_ylabel('Ecliptic observing efficiency')
    cbar.ax.set_ylabel('Anti-Sun ecliptic longitude')

    ax.set_ylim([0,1])

    # make the matching x-axis with time. and plot the corresponding lines
    ax2 = ax.twiny()
    ax2.set_xlabel('Year')
    ax2.scatter(tmid_orbit, em_efficiency, c=antisun_elon, cmap=cm, zorder=10,
                alpha=0)

    primarystart = Time('2018-07-18')
    emstart = Time('2020-06-27')
    emend = Time('2022-10-08')

    eclstart = Time('2021-10-04')
    eclend = Time('2022-02-18')
    # NOTE: we botched the start time by like an orbit in the original
    # selection window...
    eclstart += 27.3*u.day/2
    eclend += 27.3*u.day/2

    ax2.vlines(eclstart.byear, 0, 1, color='gray', zorder=-3, alpha=0.9, lw=1)
    ax2.vlines(eclend.byear, 0, 1, color='gray', zorder=-3, alpha=0.9, lw=1)
    ax2.fill_between(np.linspace(eclstart.byear,eclend.byear,100), 0, 1,
                     color='gray', zorder=-5, alpha=0.1)

    ax2.vlines(primarystart.byear, 0, 1, color='gray', zorder=-3, alpha=0.3,
               linestyles='--', lw=1)
    ax2.vlines(emstart.byear, 0, 1, color='gray', zorder=-3, alpha=0.3,
               linestyles='--', lw=1)
    ax2.vlines(emend.byear, 0, 1, color='gray', zorder=-3, alpha=0.3,
               linestyles='--', lw=1)

    # ax2.text(
    #     primarystart.byear + (emstart.byear - primarystart.byear)/2,
    #     0.08,
    #     'Prime',
    #     ha='center',va='top'
    # )
    # ax2.text(
    #     emstart.byear + (emend.byear - emstart.byear)/2,
    #     0.08,
    #     'Extended',
    #     ha='center',va='top'
    # )
    ax2.text(
        eclstart.byear + (eclend.byear - eclstart.byear)/2,
        0.3,
        'Ecliptic\nSurvey',
        ha='center',va='center',
        rotation=0,
        fontsize='small'
    )

    ax.get_yaxis().set_tick_params(which='both', direction='in')
    ax.get_xaxis().set_tick_params(which='both', direction='in')
    ax2.get_xaxis().set_tick_params(which='both', direction='in')

    outdir = '../results/ecliptic_timing/'
    outpath = os.path.join(
        outdir,'avgd_efficiency_vs_orbit_number.png'
    )
    f.tight_layout()
    f.savefig(outpath,bbox_inches='tight', dpi=400)
    print('made {}'.format(outpath))


if __name__ == "__main__":

    datadir = '../data/vanderspek_ecliptic/'

    if not os.path.exists(datadir):
        make_data_dir()

    df = make_timeseries_df(datadir)

    outdir = '../results/ecliptic_timing/efficiency_vs_time_by_sector/'
    for sectornum in range(1,58):
        outpath = os.path.join(
            outdir,'efficiency_vs_time_sector{}.png'.format(str(sectornum).zfill(3))
        )
        if not os.path.exists(outpath):
            plot_efficiency_vs_time_for_sectornum(df, sectornum=sectornum)

    plot_avgd_efficiency_vs_orbit_number_plot(df)

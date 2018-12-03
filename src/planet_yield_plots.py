# -*- coding: utf-8 -*-
"""
much of this adapted from
https://github.com/mrtommyb/textended/blob/master/code/make_detected_planet_catalog.ipynb
"""

from __future__ import division, print_function

import numpy as np, pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

from numpy import array as nparr
import os


def get_data(fpath):
    """
    available columns:

    ['Unnamed: 0', 'Unnamed: 0.1', 'TICID', 'RA', 'DEC', 'PLX', 'ECLONG',
    'ECLAT', 'V', 'J', 'Ks', 'TESSMAG', 'TEFF', 'RADIUS', 'MASS', 'CONTRATIO',
    'PRIORITY', 'isMdwarf', 'isGiant', 'isSubgiant', 'cosi', 'noise_level',
    'Nplanets', 'planetRadius', 'planetPeriod', 'starID', 'T0', 'ars', 'ecc',
    'omega', 'rprs', 'impact', 'duration', 'duration_correction',
    'transit_depth', 'transit_depth_diluted', 'has_transits', '1', '2', '3',
    '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16',
    '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28',
    '29', '30', '31', '32', '33', '34', '35', '36', '37', '38', '39', '40',
    '41', '42', '43', '44', '45', '46', '47', '48', '49', '50', '51', '52',
    '53', '54', '55', '56', '57', '58', '59', '60', '61', '62', '63', '64',
    '65', 'isObserved', 'Ntransits', 'Ntransits_primary', 'SNR', 'SNR_primary',
    'needed_for_detection', 'detected', 'needed_for_detection_primary',
    'detected_primary', 'insol', 'inOptimisticHZ']
    """

    df = pd.read_csv(fpath)

    return df

def radius_period_diagram(df, xlim=(0.3,500), ylim=(0.3,30), outdir=None):

    fig,ax = plt.subplots(figsize=(6,4))

    in_prime = nparr(df['detected_primary'])
    in_ext = nparr(df['detected']) & ~in_prime

    ax.scatter(df.loc[in_prime, 'planetPeriod'],
               df.loc[in_prime, 'planetRadius'],
               label='detected after primary', lw=0, alpha=0.8, s=10)

    ax.scatter(df.loc[in_ext, 'planetPeriod'],
               df.loc[in_ext, 'planetRadius'],
               label='detected after extended', lw=0, alpha=0.8, s=10)

    ax.set_ylabel('Radius [$R_\oplus$]')
    ax.set_xlabel('Period [days]')

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.legend(loc='lower right', fontsize='xx-small')

    ax.get_yaxis().set_tick_params(which='both', direction='in')
    ax.get_xaxis().set_tick_params(which='both', direction='in')

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    titlestr = os.path.basename(outdir).replace('_',' ')
    ax.set_title(titlestr, fontsize='x-small')

    savpath = os.path.join(outdir,'radius_period_diagram.png')
    fig.savefig(savpath, bbox_inches='tight', dpi=400)
    print('made {}'.format(savpath))


def radius_period_diagram_multitransiting(df, xlim=(0.3,500), ylim=(0.3,30),
                                          outdir=None):

    fig,ax = plt.subplots(figsize=(6,4))

    in_prime = nparr(df['detected_primary'])
    in_ext = nparr(df['detected']) & ~in_prime

    N_transiting_planets = len(df)
    N_systems = len(np.unique(df['TICID']))
    u, u_inds, inv_inds, counts = np.unique(df['TICID'],
                                            return_index=True,
                                            return_inverse=True,
                                            return_counts=True)
    multiplicity_function = counts     # length: number of systems
    df_multiplicity = counts[inv_inds] # length: number of transiting planets

    is_single = df_multiplicity == 1
    is_multi = df_multiplicity != 1

    ax.scatter(df.loc[in_prime & is_single, 'planetPeriod'],
               df.loc[in_prime & is_single, 'planetRadius'],
               label='detected after primary & single', alpha=0.3,
               facecolors='none', edgecolors='#1f77b4', lw=0.5 , zorder=1, s=5)

    ax.scatter(df.loc[in_ext, 'planetPeriod'],
               df.loc[in_ext, 'planetRadius'],
               label='detected after extended & single', lw=0.5, alpha=0.3,
               facecolors='none', edgecolors='#ff7f0e', zorder=2, s=5)

    ax.scatter(df.loc[in_prime & is_multi, 'planetPeriod'],
               df.loc[in_prime & is_multi, 'planetRadius'],
               label='detected after primary & multi', lw=0, alpha=0.8,
               facecolors='#1f77b4', s=10)

    ax.scatter(df.loc[in_ext & is_multi, 'planetPeriod'],
               df.loc[in_ext & is_multi, 'planetRadius'],
               label='detected after extended & multi', lw=0, alpha=0.8,
               facecolors='#ff7f0e', s=10)

    ax.set_ylabel('Radius [$R_\oplus$]')
    ax.set_xlabel('Period [days]')

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.legend(loc='lower right', fontsize='xx-small')

    ax.get_yaxis().set_tick_params(which='both', direction='in')
    ax.get_xaxis().set_tick_params(which='both', direction='in')

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    titlestr = os.path.basename(outdir).replace('_',' ')
    ax.set_title(titlestr, fontsize='x-small')

    savpath = os.path.join(outdir,'radius_period_diagram_multitransiting.png')
    fig.savefig(savpath, bbox_inches='tight', dpi=400)
    print('made {}'.format(savpath))


def radius_insolation_diagram(df, xlim=None, ylim=None, outdir=None):

    fig,ax = plt.subplots(figsize=(6,4))

    in_prime = nparr(df['detected_primary'])
    in_ext = nparr(df['detected']) & ~in_prime

    ax.scatter(df.loc[in_prime, 'insol'],
               df.loc[in_prime, 'planetRadius'],
               label='detected after primary', lw=0, alpha=0.8, s=10)

    ax.scatter(df.loc[in_ext, 'insol'],
               df.loc[in_ext, 'planetRadius'],
               label='detected after extended', lw=0, alpha=0.8, s=10)

    ax.set_ylabel('Radius [$R_\oplus$]')
    ax.set_xlabel('Insolation [$S_\oplus$]')

    ax.set_xscale('log')
    ax.set_yscale('log')

    ax.legend(loc='best', fontsize='xx-small')

    ax.get_yaxis().set_tick_params(which='both', direction='in')
    ax.get_xaxis().set_tick_params(which='both', direction='in')

    if xlim:
        ax.set_xlim(xlim)
    if ylim:
        ax.set_ylim(ylim)

    titlestr = os.path.basename(outdir).replace('_',' ')
    ax.set_title(titlestr, fontsize='x-small')

    savpath = os.path.join(outdir,'radius_insolation_diagram.png')
    fig.savefig(savpath, bbox_inches='tight', dpi=400)
    print('made {}'.format(savpath))


def autolabel(rects, ax):
    """
    Attach a text label above each bar displaying its height
    """
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2., height+(height*0.02),
                '%d' % int(height),
                ha='center', va='bottom')

def autolabel2(rects, ax):
    """
    Attach a text label above each bar displaying its height
    """
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2., height-(height*0.22),
                '%d' % int(height),
                ha='center', va='bottom')


def number_of_detections_vs_plradius_barchart(df, txt=None, ylim=None,
                                              justmultis=False, outstr='',
                                              yscale=None, outdir=None):

    in_pri = nparr(df['detected_primary'])
    after_ext = nparr(df['detected'])

    if outstr=='justmultis_':
        u, u_inds, inv_inds, counts = np.unique(df['TICID'],
                                                return_index=True,
                                                return_inverse=True,
                                                return_counts=True)
        df_multiplicity = counts[inv_inds]
        is_single = df_multiplicity == 1
        is_multi = df_multiplicity != 1

        after_ext &= is_multi
        in_pri &= is_multi

    elif outstr=='oneortwotraonly_':
        ntra_pri = df['Ntransits_primary']
        in_pri &= (  (ntra_pri==1) | (ntra_pri==2)  )
        after_ext &= (  (ntra_pri==1) | (ntra_pri==2)  )


    counts_ext = np.histogram(df.loc[after_ext, 'planetRadius'],
                              bins=[0,1.25,2,4,25])

    counts_pri = np.histogram(df.loc[in_pri, 'planetRadius'],
                              bins=[0,1.25,2,4,25])

    tl = ['$<1.25$', '$1.25-2.0$', '$2.0-4.0$', '$>4.0$']

    plt.close("all")
    fig, ax = plt.subplots(1,1,figsize=[6,6])

    if outstr not in ['oneortwotraonly_']:
        h_pri = ax.bar(np.arange(4), counts_pri[0], tick_label=tl,
                       label='after primary', zorder=2, color='#1f77b4')
    h_ext = ax.bar(np.arange(4), counts_ext[0], tick_label=tl,
                   label='after extended', zorder=1, color='#ff7f0e')

    ax.set_xlabel('Planet radii ($ R_\oplus$)')
    ax.set_ylabel('Number of planets')

    if yscale:
        ax.set_yscale(yscale)
    else:
        ax.set_yscale('log')

    if txt:
        ax.text(0.05, 0.8, txt, ha='left', va='center', fontsize='x-small',
                transform=ax.transAxes)

    autolabel(h_ext, ax)
    if outstr not in ['oneortwotraonly_']:
        autolabel2(h_pri, ax)

    ax.legend()

    if ylim:
        ax.set_ylim(ylim)

    ax.get_yaxis().set_tick_params(which='both', direction='in')
    ax.get_xaxis().set_tick_params(which='both', direction='in')

    titlestr = os.path.basename(outdir).replace('_',' ')
    ax.set_title(titlestr, fontsize='x-small')

    savpath = os.path.join(
        outdir,outstr+'number_of_detections_vs_plradius_barchart.png')
    fig.savefig(savpath, bbox_inches='tight', dpi=400)
    print('made {}'.format(savpath))


def number_of_detections_vs_period_barchart(df, txt=None, ylim=None,
                                            justmultis=False, outstr='',
                                            yscale=None, outdir=None):

    in_pri = nparr(df['detected_primary'])
    after_ext = nparr(df['detected'])

    if outstr=='justmultis_':
        u, u_inds, inv_inds, counts = np.unique(df['TICID'],
                                                return_index=True,
                                                return_inverse=True,
                                                return_counts=True)
        df_multiplicity = counts[inv_inds]
        is_single = df_multiplicity == 1
        is_multi = df_multiplicity != 1

        after_ext &= is_multi
        in_pri &= is_multi

    elif outstr=='oneortwotraonly_':
        ntra_pri = df['Ntransits_primary']
        in_pri &= (  (ntra_pri==1) | (ntra_pri==2)  )
        after_ext &= (  (ntra_pri==1) | (ntra_pri==2)  )


    counts_ext = np.histogram(df.loc[after_ext, 'planetPeriod'],
                              bins=[0,20,50,1000])

    counts_pri = np.histogram(df.loc[in_pri, 'planetPeriod'],
                              bins=[0,20,50,1000])

    tl = ['$<20$', '$20-50$', '$>50$']

    plt.close("all")
    fig, ax = plt.subplots(1,1,figsize=[6,6])

    if outstr not in ['oneortwotraonly_']:
        h_pri = ax.bar(np.arange(3), counts_pri[0], tick_label=tl,
                       label='after primary', zorder=2, color='#1f77b4')
    h_ext = ax.bar(np.arange(3), counts_ext[0], tick_label=tl,
                   label='after extended', zorder=1, color='#ff7f0e')

    ax.set_xlabel('Orbital period (days)')
    ax.set_ylabel('Number of planets')

    if yscale:
        ax.set_yscale(yscale)
    else:
        ax.set_yscale('log')

    if txt:
        ax.text(0.05, 0.8, txt, ha='left', va='center', fontsize='x-small',
                transform=ax.transAxes)

    autolabel(h_ext, ax)
    if outstr not in ['oneortwotraonly_']:
        autolabel2(h_pri, ax)

    ax.legend()

    if ylim:
        ax.set_ylim(ylim)

    ax.get_yaxis().set_tick_params(which='both', direction='in')
    ax.get_xaxis().set_tick_params(which='both', direction='in')

    titlestr = os.path.basename(outdir).replace('_',' ')
    ax.set_title(titlestr, fontsize='x-small')

    savpath = os.path.join(
        outdir,outstr+'number_of_detections_vs_period_barchart.png')
    fig.savefig(savpath, bbox_inches='tight', dpi=400)
    print('made {}'.format(savpath))



def number_of_detections_vs_plradius_barchart_Tmaglt(df, txt=None, Tmagcut=10,
                                                     ylim=None,
                                                     justmultis=False,
                                                     outstr='', outdir=None,
                                                     yscale='log'):

    in_pri = nparr(df['detected_primary'])
    after_ext = nparr(df['detected'])

    if outstr=='justmultis_':
        u, u_inds, inv_inds, counts = np.unique(df['TICID'],
                                                return_index=True,
                                                return_inverse=True,
                                                return_counts=True)
        df_multiplicity = counts[inv_inds]
        is_single = df_multiplicity == 1
        is_multi = df_multiplicity != 1

        after_ext &= is_multi
        in_pri &= is_multi

    sel = df['TESSMAG'] < Tmagcut

    counts_ext = np.histogram(df.loc[after_ext & sel, 'planetRadius'],
                              bins=[0,1.25,2,4,25])

    counts_pri = np.histogram(df.loc[in_pri & sel, 'planetRadius'],
                              bins=[0,1.25,2,4,25])

    tl = ['$<1.25$', '$1.25-2.0$', '$2.0-4.0$', '$>4.0$']

    plt.close("all")
    fig, ax = plt.subplots(1,1,figsize=[6,6])

    h_pri = ax.bar(np.arange(4), counts_pri[0], tick_label=tl,
                   label='after primary', zorder=2)
    h_ext = ax.bar(np.arange(4), counts_ext[0], tick_label=tl,
                   label='after extended', zorder=1)

    if txt:
        ax.text(0.05, 0.8, txt, ha='left', va='center', fontsize='x-small',
                transform=ax.transAxes)

    ax.set_xlabel('Planet radii ($ R_\oplus$)')
    ax.set_ylabel('Number of planets')

    ax.set_yscale(yscale)

    autolabel(h_ext, ax)
    autolabel2(h_pri, ax)

    ax.legend()

    if ylim:
        ax.set_ylim(ylim)

    ax.get_yaxis().set_tick_params(which='both', direction='in')
    ax.get_xaxis().set_tick_params(which='both', direction='in')

    titlestr = os.path.basename(outdir).replace('_',' ')
    ax.set_title(titlestr, fontsize='x-small')

    savpath = ( os.path.join(outdir,outstr+
        'number_of_detections_vs_plradius_barchart_Tmaglt{:d}.png'.
        format(Tmagcut))
    )
    fig.savefig(savpath, bbox_inches='tight', dpi=400)
    print('made {}'.format(savpath))


def number_of_detections_vs_ntra_barchart(df, txt=None, ylim=None, outstr='',
                                          outdir=None):

    in_pri = nparr(df['detected_primary'])
    after_ext = nparr(df['detected'])
    if outstr=='oneortwotraonly_':
        ntra_pri = df['Ntransits_primary']
        in_pri &= (  (ntra_pri==1) | (ntra_pri==2)  )
        after_ext &= (  (ntra_pri==1) | (ntra_pri==2)  )

    counts_ext = np.histogram(df.loc[after_ext, 'Ntransits'],
                              bins=[0.5,2.5,4.5, 8.5,15.5,27.5, 50, 10000])

    counts_pri = np.histogram(df.loc[in_pri, 'Ntransits_primary'],
                              bins=[0.5,2.5,4.5, 8.5,15.5,27.5, 50, 10000])

    tl = ['2', '3-4', '5-8', '9-15', '16-27', '27-50', '50+']

    plt.close("all")
    fig, ax = plt.subplots(1,1,figsize=[6,6])

    h_pri = ax.bar(np.arange(7), counts_pri[0], tick_label=tl,
                   label='after primary', zorder=2)
    h_ext = ax.bar(np.arange(7), counts_ext[0], tick_label=tl,
                   label='after extended', zorder=1)

    ax.set_xlabel('Number of transits')
    ax.set_ylabel('Number of planets')

    if txt:
        ax.text(0.05, 0.8, txt, ha='left', va='center', fontsize='x-small',
                transform=ax.transAxes)

    autolabel(h_ext, ax)
    autolabel2(h_pri, ax)

    ax.legend()

    if ylim:
        ax.set_ylim(ylim)

    titlestr = os.path.basename(outdir).replace('_',' ')
    ax.set_title(titlestr, fontsize='x-small')

    ax.get_yaxis().set_tick_params(which='both', direction='in')
    ax.get_xaxis().set_tick_params(which='both', direction='in')

    savpath = os.path.join(
        outdir,outstr+'number_of_detections_vs_ntra_barchart.png')
    fig.savefig(savpath, bbox_inches='tight', dpi=400)
    print('made {}'.format(savpath))


def number_of_detections_vs_insolation_barchart(df, txt=None, ylim=None,
                                                justmultis=False, outstr='',
                                                outdir=None):

    in_pri = nparr(df['detected_primary'])
    after_ext = nparr(df['detected'])

    if justmultis:
        u, u_inds, inv_inds, counts = np.unique(df['TICID'],
                                                return_index=True,
                                                return_inverse=True,
                                                return_counts=True)
        df_multiplicity = counts[inv_inds]
        is_single = df_multiplicity == 1
        is_multi = df_multiplicity != 1

        after_ext &= is_multi
        in_pri &= is_multi

        outstr = 'justmultis_'

    counts_ext = np.histogram(df.loc[after_ext, 'insol'],
                              bins=[0,0.5,2,1e6])

    counts_pri = np.histogram(df.loc[in_pri, 'insol'],
                              bins=[0,0.5,2,1e6])

    tl = ['$<0.5$', '$0.5-2$', '$>2$']

    plt.close("all")
    fig, ax = plt.subplots(1,1,figsize=[6,6])

    h_pri = ax.bar(np.arange(3), counts_pri[0], tick_label=tl,
                   label='after primary', zorder=2)
    h_ext = ax.bar(np.arange(3), counts_ext[0], tick_label=tl,
                   label='after extended', zorder=1)

    ax.set_xlabel('Planet insolation ($ S_\oplus$)')
    ax.set_ylabel('Number of planets')

    ax.set_yscale('log')

    if txt:
        ax.text(0.05, 0.8, txt, ha='left', va='center', fontsize='x-small',
                transform=ax.transAxes)

    autolabel(h_ext, ax)
    autolabel2(h_pri, ax)

    ax.legend()

    if ylim:
        ax.set_ylim(ylim)

    titlestr = os.path.basename(outdir).replace('_',' ')
    ax.set_title(titlestr, fontsize='x-small')

    ax.get_yaxis().set_tick_params(which='both', direction='in')
    ax.get_xaxis().set_tick_params(which='both', direction='in')

    savpath = ( os.path.join( outdir,
        '{:s}number_of_detections_vs_insolation_barchart.png'.format(outstr))
    )
    fig.savefig(savpath, bbox_inches='tight', dpi=400)
    print('made {}'.format(savpath))


def number_of_detections_vs_sptype_barchart(df, txt=None, ylim=None,
                                            justmultis=False, outstr='',
                                            outdir=None):


    in_pri = nparr(df['detected_primary'])
    after_ext = nparr(df['detected'])

    if justmultis:
        u, u_inds, inv_inds, counts = np.unique(df['TICID'],
                                                return_index=True,
                                                return_inverse=True,
                                                return_counts=True)
        df_multiplicity = counts[inv_inds]
        is_single = df_multiplicity == 1
        is_multi = df_multiplicity != 1

        after_ext &= is_multi
        in_pri &= is_multi

        outstr = 'justmultis_'

    counts_ext = np.histogram(df.loc[after_ext, 'TEFF'],
                              bins=[2285,3905,5310,5980,7330.,10050])

    counts_pri = np.histogram(df.loc[in_pri, 'TEFF'],
                              bins=[2285,3905,5310,5980,7330.,10050])

    tl = ['M', 'K', 'G', 'F', 'A']

    plt.close("all")
    fig, ax = plt.subplots(1,1,figsize=[6,6])

    h_pri = ax.bar(np.arange(5), counts_pri[0], tick_label=tl,
                   label='after primary', zorder=2)
    h_ext = ax.bar(np.arange(5), counts_ext[0], tick_label=tl,
                   label='after extended', zorder=1)

    ax.set_xlabel('Host spectral type')
    ax.set_ylabel('Number of planets')

    ax.set_yscale('linear')

    if txt:
        ax.text(0.05, 0.8, txt, ha='left', va='center', fontsize='x-small',
                transform=ax.transAxes)

    autolabel(h_ext, ax)
    autolabel2(h_pri, ax)

    ax.legend()

    if ylim:
        ax.set_ylim(ylim)

    titlestr = os.path.basename(outdir).replace('_',' ')
    ax.set_title(titlestr, fontsize='x-small')

    ax.get_yaxis().set_tick_params(which='both', direction='in')
    ax.get_xaxis().set_tick_params(which='both', direction='in')

    savpath = ( os.path.join( outdir,
        '{:s}number_of_detections_vs_sptype_barchart.png'.format(outstr))
    )
    fig.savefig(savpath, bbox_inches='tight', dpi=400)
    print('made {}'.format(savpath))


def new_planets_newsky_or_snrboost_vs_plradius(df, txt=None, ylim=None,
                                               outdir=None):
    """
    * Number of planets found in the EM
        * New Planets found due to increased SNR relative to Prime
        * New Planets found in regions of sky that Prime did not observe
    """

    in_pri = nparr(df['detected_primary'])
    after_ext = nparr(df['detected'])
    only_ext = after_ext & ~in_pri

    ntra_pri = nparr(df['Ntransits_primary'])
    ntra_ext = nparr(df['Ntransits'])

    from_newsky = (ntra_pri == 0)
    from_snrorntra = ~from_newsky

    bins = [0,1.25,2,4,25]
    counts_newsky = np.histogram(
        df.loc[only_ext & from_newsky, 'planetRadius'],
        bins=bins)

    counts_snrorntra = np.histogram(
        df.loc[only_ext & from_snrorntra, 'planetRadius'],
        bins=bins)

    tl = ['$<1.25$', '$1.25-2.0$', '$2.0-4.0$', '$>4.0$']

    plt.close("all")
    fig, ax = plt.subplots(1,1,figsize=[6,6])

    ind = np.arange(len(bins)-1) # x locations for groups
    width = 0.35                 # width of bars

    h_newsky = ax.bar(ind, counts_newsky[0], width, tick_label=tl,
                   label='from new sky', zorder=2)
    h_snrorntra = ax.bar(ind+width, counts_snrorntra[0], width, tick_label=tl,
                   label='from SNR or Ntra boost', zorder=1)

    if txt:
        ax.text(0.05, 0.8, txt, ha='left', va='center', fontsize='x-small',
                transform=ax.transAxes)

    ax.set_xticks(ind + width/2)

    ax.set_xlabel('Planet radii ($ R_\oplus$)')
    ax.set_ylabel('Number of planets')

    ax.set_yscale('linear')

    autolabel(h_newsky, ax)
    autolabel(h_snrorntra, ax)

    ax.legend()

    if ylim:
        ax.set_ylim(ylim)

    titlestr = os.path.basename(outdir).replace('_',' ')
    ax.set_title(titlestr, fontsize='x-small')

    ax.get_yaxis().set_tick_params(which='both', direction='in')
    ax.get_xaxis().set_tick_params(which='both', direction='in')

    savpath = ( os.path.join( outdir,
        'new_planets_newsky_or_snrboost_vs_plradius.png' )
    )
    fig.savefig(savpath, bbox_inches='tight', dpi=400)
    print('made {}'.format(savpath))


if __name__=="__main__":

    names = ['idea_1_SNE','idea_2_SNSNS','idea_3_SNNSN','idea_4_EC3PO']
    outdirs = [os.path.join('../results/planet_yield_plots/',dirname)
               for dirname in names]
    datanames = [n.split('_')[-1] for n in names]
    datapaths = [os.path.join(
        '../data/tommyb', 'detected_planet_catalog_{:s}-v3.csv.bz2'.
        format(dn)) for dn in datanames]

    one_thru_three = 0
    four_thru_eight = 0
    nine_thru_thirteen = 0
    fourteen_thru_X =1

    for datapath, outdir in zip(datapaths, outdirs):
        print(datapath)
        df = get_data(datapath)

        if not os.path.exists(outdir):
            os.mkdir(outdir)

        if one_thru_three:
            radius_period_diagram(df, outdir=outdir, xlim=(0.3,500),
                                  ylim=(0.3,30))
            radius_period_diagram_multitransiting(df, outdir=outdir,
                                                  xlim=(0.3,500),
                                                  ylim=(0.3,30))
            radius_insolation_diagram(df, outdir=outdir, xlim=(0.11, 2.5e5),
                                      ylim=(0.45,30))

        if four_thru_eight:

            number_of_detections_vs_plradius_barchart(
                df, txt='SNR>10 and $N_{\mathrm{tra,ext}} >= 3$',
                ylim=(1,3800), outdir=outdir)

            Tmagcut=10
            number_of_detections_vs_plradius_barchart_Tmaglt(
                df, Tmagcut=Tmagcut,
                txt='SNR>10 and $N_{\mathrm{tra,ext}} >= 3$ '+
                'and T<{:d}'.format(Tmagcut),
                ylim=(1,3800), outdir=outdir)

            number_of_detections_vs_ntra_barchart(
                df, txt='SNR>10 and $N_{\mathrm{tra,ext}} >= 3$',
                ylim=(0,1600), outdir=outdir)


            number_of_detections_vs_plradius_barchart(
                df, txt='SNR>10 and $N_{\mathrm{tra,ext}} >= 3$ and just multis',
                ylim=(1,3800), outstr='justmultis_', outdir=outdir)

            Tmagcut=10
            number_of_detections_vs_plradius_barchart_Tmaglt(
                df, Tmagcut=Tmagcut,
                txt='SNR>10 and $N_{\mathrm{tra,ext}} >= 3$ '+
                'and T<{:d} and just multis'.format(Tmagcut),
                ylim=(0,1000), outstr='justmultis_', outdir=outdir,
                yscale='linear')

        if nine_thru_thirteen:

            number_of_detections_vs_insolation_barchart(
                df, txt='SNR>10 and $N_{\mathrm{tra,ext}} >= 3$',
                outdir=outdir, ylim=(0.5, 6800))

            number_of_detections_vs_sptype_barchart(
                df, txt='SNR>10 and $N_{\mathrm{tra,ext}} >= 3$',
                outdir=outdir, ylim=(0,1800))

            new_planets_newsky_or_snrboost_vs_plradius(
                df,
                txt='new planets only\nSNR>10 and $N_{\mathrm{tra,ext}} >= 3$',
                outdir=outdir,
                ylim=(0.5,2100)
            )

            number_of_detections_vs_ntra_barchart(
                df,
                txt='SNR>10 and $N_{\mathrm{tra,pri}}$ is 1 or 2, and $N_{\mathrm{tra,ext}} >= 3$',
                ylim=(0,1600), outstr='oneortwotraonly_', outdir=outdir)

            number_of_detections_vs_plradius_barchart(
                df, txt='SNR>10 and $N_{\mathrm{tra,pri}}$ is 1 or 2, and $N_{\mathrm{tra,ext}} >= 3$',
                ylim=(0,1000), outstr='oneortwotraonly_', yscale='linear', outdir=outdir)

        if fourteen_thru_X:

            number_of_detections_vs_period_barchart(
                df, txt='SNR>10 and $N_{\mathrm{tra,ext}} >= 3$',
                ylim=(1,5000), outdir=outdir, yscale='log')

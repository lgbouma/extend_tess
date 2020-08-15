"""
Quick sanity checks of TOI catalog vs Barclay 2018 predictions.
"""

import numpy as np, pandas as pd, matplotlib.pyplot as plt
from astropy.io import fits
from astropy.table import Table

from cdips.utils.catalogs import get_exofop_toi_catalog
from aesthetic.plot import set_style, set_style_grid, savefig
from matplotlib.lines import Line2D

from os.path import join
from numpy import array as arr

RESULTSDIR = '/Users/luke/Dropbox/proj/extend_tess/results/grr_20200814'

def autolabel(rects, ax):
    """
    Attach a text label above each bar displaying its height
    """
    for rect in rects:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2., height+(height*0.02),
                '%d' % int(height),
                ha='center', va='bottom')



def plot_cdf(b18_key, b18_df, toi_key, toi_df, **kwargs):
    """
    e.g., "SNR", b18_df, "Planet SNR", toi_df
    """

    N_b18 = len(b18_df)
    N_toi = len(toi_df)

    #
    # make the plot!
    #
    set_style()

    f, ax = plt.subplots(figsize=(4,3))

    if b18_key == 'SNR':
        bins = np.arange(1,200,0.1)
        savekey = 'planetsnr'
    else:
        raise NotImplementedError

    cnt, bin_edges = np.histogram(arr(b18_df[b18_key]), bins=bins, normed=True)
    cdf = np.cumsum(cnt)
    ax.plot(bin_edges[:-1], cdf/cdf[-1], label=f'Barclay+18')

    cnt, bin_edges = np.histogram(arr(toi_df[toi_key]), bins=bins, normed=True)
    cdf = np.cumsum(cnt)
    ax.plot(bin_edges[:-1], cdf/cdf[-1], label=f'TOI 2020.08.14')

    if b18_key == 'SNR':
        ax.set_xlim([5,200])
        ax.set_xscale('log')

    ax.legend(fontsize='small', loc='best')

    ax.set_xlabel('Planet SNR')
    ax.set_ylabel('Fraction of total')

    for k,v in kwargs.items():
        if k == 'snrabove10' and v:
            savekey += '_snrabove10'
            titlestr = 'SNR > 10'
        if k == 'snrabove10' and not v:
            savekey += '_allsnr'
            titlestr = 'All SNR'
        if k == 'rpbelow4' and v:
            savekey += '_rpbelow4'
            titlestr += ', R$_\mathrm{p}$ < 4R$_\oplus$'
        if k == 'rpbelow4' and not v:
            savekey += '_allrp'
            titlestr += ', all R$_\mathrm{p}$'

    ax.set_title(titlestr)

    outpath = join(
        RESULTSDIR, f'cdf_{savekey}.png'
    )
    savefig(f, outpath, writepdf=0)


def plot_hist(b18_key, b18_df, toi_key, toi_df, **kwargs):
    """
    e.g., "SNR", b18_df, "Planet SNR", toi_df
    """

    N_b18 = len(b18_df)
    N_toi = len(toi_df)

    #
    # make the plot!
    #
    set_style()

    f, axs = plt.subplots(nrows=2, ncols=1, figsize=(4,6), sharex=True)

    ax = axs[0]

    ##########
    if b18_key == 'Tmag':
        diff = 1
        xmax = 16+diff
        bins = np.arange(5,xmax,diff)
        divbins = np.arange(5+diff/2,xmax-diff,diff)
        savekey = 'Tmag'
        xlabel = 'Tmag'

    else:
        raise NotImplementedError

    h0 = ax.hist(b18_df[b18_key], bins=bins, cumulative=False, fill=False,
                 density=False, histtype='step', label=f'Barclay+18')

    h1 = ax.hist(toi_df[toi_key], bins=bins, cumulative=False, fill=False,
                 density=False, histtype='step', label=f'TOI 2020.08.14')

    # Create new legend handles but use the colors from the existing ones
    handles, labels = ax.get_legend_handles_labels()
    new_handles = [Line2D([], [], c=h.get_edgecolor()) for h in handles]

    ax.legend(handles=new_handles, labels=labels, fontsize='x-small', loc='upper left')

    for k,v in kwargs.items():
        if k == 'snrabove10' and v:
            savekey += '_snrabove10'
            titlestr = 'SNR > 10'
        if k == 'snrabove10' and not v:
            savekey += '_allsnr'
            titlestr = 'All SNR'
        if k == 'rpbelow4' and v:
            savekey += '_rpbelow4'
            titlestr += ', R$_\mathrm{p}$ < 4R$_\oplus$'
        if k == 'rpbelow4' and not v:
            savekey += '_allrp'
            titlestr += ', all R$_\mathrm{p}$'


    ax.set_title(titlestr)

    ax.set_ylabel('Number per bin')
    ##########

    ax = axs[1]

    ax.plot(divbins, h1[0]/h0[0], c='k', marker='o', ls='--', lw=1)

    ax.set_ylabel('TOI / Barclay')

    ax.set_xlabel(xlabel)

    f.tight_layout()

    outpath = join(
        RESULTSDIR, f'hist_{savekey}.png'
    )
    savefig(f, outpath, writepdf=0)


def plot_twobar(b18_key, b18_df, toi_key, toi_df, **kwargs):
    """
    e.g., 'Rp' and 'Planet Radius (R_Earth)'
    """

    if b18_key in ['Rp', 'Planet Radius (R_Earth)']:
        bins = [0,1.25,2,4,25]
        tl = ['$<1.25$', '$1.25-2.0$', '$2.0-4.0$', '$>4.0$']
        xlabel = 'Planet radius ($R_\oplus$)'
        yscale = 'log'
        descriptor = 'twobar_rp'
    elif b18_key in ['Perp', 'Period (days)']:
        bins = [0,10,20,40,360]
        tl = ['$<10$', '$10-20$', '$20-40$', '$>40$']
        xlabel = 'Orbital period [days]'
        yscale = 'log'
        descriptor = 'twobar_period'
    elif b18_key in ['Teff', 'Stellar Eff Temp (K)']:
        bins = [2285,3905,5310,5980,7330.,10050]
        tl = ['M', 'K', 'G', 'F', 'A']
        xlabel = 'SpType'
        yscale = 'linear'
        descriptor = 'twobar_SpType'
    else:
        raise NotImplementedError

    set_style()

    plt.close("all")
    fig, ax = plt.subplots(1,1,figsize=[4,4])

    width = 0.35 # width of bars

    counts = np.histogram(b18_df[b18_key], bins=bins)
    yval = counts[0]
    h = ax.bar(np.arange(len(tl)), yval, width, tick_label=tl,
               label='Barclay+18', zorder=2)
    autolabel(h, ax)

    counts = np.histogram(toi_df[toi_key], bins=bins)
    yval = counts[0]
    h = ax.bar(np.arange(len(tl))+width, yval, width, tick_label=tl,
               label='TOI 2020.08.14', zorder=2)
    autolabel(h, ax)

    ax.set_xticks(np.arange(len(tl)) + width/2)

    ax.set_xlabel(xlabel)
    ax.set_ylabel('Number')

    ax.set_yscale(yscale)

    ax.set_ylim((ax.get_ylim()[0], ax.get_ylim()[1]+300))

    loc = 'upper left'
    if b18_key == 'Perp':
        loc = 'upper right'
    ax.legend(fontsize='small', loc=loc)

    savekey = ''
    for k,v in kwargs.items():
        if k == 'snrabove10' and v:
            savekey += '_snrabove10'
            titlestr = 'SNR > 10'
        if k == 'snrabove10' and not v:
            savekey += '_allsnr'
            titlestr = 'All SNR'
        if k == 'rpbelow4' and v:
            savekey += '_rpbelow4'
            titlestr += ', R$_\mathrm{p}$ < 4R$_\oplus$'
        if k == 'rpbelow4' and not v:
            savekey += '_allrp'
            titlestr += ', all R$_\mathrm{p}$'


    ax.set_title(titlestr)

    ax.minorticks_off()

    outpath = join(
        RESULTSDIR, f'{descriptor}{savekey}.png'
    )
    savefig(fig, outpath, writepdf=0)


B18KEYDICT = {
    'Rp': 'Rp',
    'Period': 'Perp',
    'Dur': 'Dur',
    'Tmag': 'Tmag'
}

TOIKEYDICT = {
    'Rp': 'Planet Radius (R_Earth)',
    'Period': 'Period (days)',
    'Dur': 'Duration (hours)',
    'Tmag': 'TESS Mag'
}

LABELDICT = {
    'Rp': 'R$_\mathrm{p}$ [R$_\oplus$]',
    'Period': 'Orbital period [days]',
    'Dur': 'Transit duration [hours]',
    'Tmag': 'T [mag]'
}


def plot_scatter(ykey, xkey, b18_df, toi_df, **kwargs):

    set_style()

    f, ax = plt.subplots(figsize=(4,3))

    ax.scatter(
        b18_df[B18KEYDICT[xkey]], b18_df[B18KEYDICT[ykey]], s=2, zorder=1,
        label='Barclay+18'
    )
    ax.scatter(
        toi_df[TOIKEYDICT[xkey]], toi_df[TOIKEYDICT[ykey]], s=2, zorder=2,
        label='TOI 2020.08.14'
    )

    ax.set_xlabel(LABELDICT[xkey])
    ax.set_ylabel(LABELDICT[ykey])

    ax.legend(fontsize='x-small', loc='best')

    ax.set_xscale('log')
    ax.set_yscale('log')

    savekey = ''
    for k,v in kwargs.items():
        if k == 'snrabove10' and v:
            savekey += '_snrabove10'
            titlestr = 'SNR > 10'
        if k == 'snrabove10' and not v:
            savekey += '_allsnr'
            titlestr = 'All SNR'
        if k == 'rpbelow4' and v:
            savekey += '_rpbelow4'
            titlestr += ', R$_\mathrm{p}$ < 4R$_\oplus$'
            ax.set_yscale('linear')
            if xkey == 'Rp':
                ax.set_xscale('linear')
        if k == 'rpbelow4' and not v:
            savekey += '_allrp'
            titlestr += ', all R$_\mathrm{p}$'
            ax.set_ylim([0.5, 102])

    if ykey == 'Tmag':
        ax.set_ylim(ax.get_ylim()[::-1])

    ax.set_title(titlestr)

    outpath = join(RESULTSDIR, f'scatter_{ykey}_vs_{xkey}{savekey}.png')
    savefig(f, outpath, writepdf=0)


def _get_toi_df(**kwargs):

    toi_df = get_exofop_toi_catalog()
    sel = (toi_df['TFOPWG Disposition'] != 'FP')

    for k,v in kwargs.items():
        if k == 'snrabove10' and bool(v):
            sel &= (toi_df['Planet SNR'] > 10)
        if k == 'rpbelow4' and bool(v):
            sel &= (toi_df['Planet Radius (R_Earth)'] < 4)

    return toi_df[sel]


def _get_barclay_df(**kwargs):

    # from vizier
    hdul = fits.open('../data/Barclay_2018_table2.fits')
    t = Table(hdul[1].data)
    b18_df = t.to_pandas()

    sel = np.ones(len(b18_df)).astype(bool)

    for k,v in kwargs.items():
        if k == 'snrabove10' and bool(v):
            sel &= ( (b18_df['SNR'] > 10) & (b18_df['Nt'] > 2) )
        if k == 'rpbelow4' and bool(v):
            sel &= (b18_df['Rp'] < 4)

    return b18_df[sel]



def main():

    for s in [0,1]:
        for r in [0,1]:

            kwargs = {
                'snrabove10': s, 'rpbelow4': r
            }

            toi_df = _get_toi_df(**kwargs)
            b18_df = _get_barclay_df(**kwargs)

            plot_scatter('Tmag', 'Rp', b18_df, toi_df, **kwargs)
            plot_scatter('Dur', 'Rp', b18_df, toi_df, **kwargs)

            plot_cdf("SNR", b18_df, "Planet SNR", toi_df, **kwargs)

            plot_hist('Tmag', b18_df, "TESS Mag", toi_df, **kwargs)

            plot_scatter('Rp', 'Period', b18_df, toi_df, **kwargs)

            plot_twobar('Teff', b18_df, 'Stellar Eff Temp (K)', toi_df, **kwargs)
            plot_twobar('Perp', b18_df, 'Period (days)', toi_df, **kwargs)

            if not kwargs['rpbelow4']:
                plot_twobar('Rp', b18_df, 'Planet Radius (R_Earth)', toi_df, **kwargs)

            del toi_df, b18_df, kwargs


if __name__ == "__main__":
    main()

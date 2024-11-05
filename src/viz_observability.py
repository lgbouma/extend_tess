import tesswcs
from tesswcs.locate import get_observability_mask
import numpy as np
import pandas as pd
from astropy.coordinates import SkyCoord
from tqdm import tqdm
import os
from os.path import join

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from matplotlib.ticker import FixedLocator
from matplotlib import rcParams
from matplotlib.projections import get_projection_class

from convert_sc_orientation import ecliptic_to_equatorial_orientation

def _shift_lon(lon, origin):
    x = np.array(np.remainder(lon+360-origin,360)) # shift lon values
    ind = x>180
    x[ind] -=360    # scale conversion to [-180, 180]
    x=-x    # reverse the scale: East to the left
    return x

def viz_obs(namestr, projection=None, showtitles=0, showticklabels=0):

    if namestr.startswith('em3'):
        sector_min, sector_max = 97, 134
    elif namestr.startswith('cumul_em3'):
        sector_min, sector_max = 1, 134
    elif namestr.startswith('cumul_pm'):
        sector_min, sector_max = 1, 26
    elif namestr.startswith('cumul_em1'):
        sector_min, sector_max = 1, 55
    elif namestr.startswith('cumul_em2'):
        sector_min, sector_max = 1, 96
    else:
        sector_min, sector_max = 1, 96

    if namestr.startswith('cumul_em3'):
        # iteration below is only over em3
        iter_sector_range = list(range(97, 134))
    else:
        iter_sector_range = list(range(sector_min, sector_max))
    sstr = f'S{str(sector_min).zfill(4)}_S{str(sector_max).zfill(4)}'

    namedict = {
        'em3_v00': 'luke_S2_S2+4_54_14_e2_G2_C11_G3_e2.csv',
        'em3_v01': 'luke_S2_S2+4_54_14_e2_G2_C12_G2_e2.csv',
        'em3_v02': 'luke_S2_S2+4_C2_G3_C9_e2_G2_54_14_e2.csv',
        'em3_v03': 'luke_S2_S2+4_C3_G2_C9_e2_G2_54_14_e2.csv',
        'em3_v04': 'luke_S2+4_70_16_e2_G2_C11_G3_e2_roll1.csv'
    }
    namekey = namestr.lstrip("cumul_")
    if namekey in namedict:
        csvname = namedict[namekey].rstrip(".csv")
    else:
        csvname = namestr

    csvpath = f'/Users/luke/Dropbox/proj/extend_tess/data/{csvname}.csv'
    if os.path.exists(csvpath):
        df = pd.read_csv(csvpath)

        ras, decs, rolls = [],[],[]
        for sectornum in iter_sector_range:
            sc_roll_key = 'sc_eroll' if 'sc_eroll' in df else 'sc_roll'
            t = tuple(df.loc[df.S==sectornum, ['sc_elon', 'sc_elat', sc_roll_key]].iloc[0])
            eq = ecliptic_to_equatorial_orientation(t[0], t[1], t[2])
            ras.append(eq[0])
            decs.append(eq[1])
            rolls.append(eq[2])

        # All the TESS pointings
        hypot_pointings = np.vstack([ras, decs, rolls]).T
    else:
        print('WRN! Did not load in any CSV file for pointings.')
        pass

    real_pointings = tesswcs.pointings[['RA', "Dec", "Roll"]].to_pandas().values

    if namestr.startswith("cumul_em3"):
        # concatenate hypothetical and real pointings for cumulative view
        hypot_pointings = np.vstack((real_pointings, hypot_pointings))
    elif namestr.startswith("cumul"):
        hypot_pointings = real_pointings[sector_min-1:sector_max,:]

    outdir = '/Users/luke/Dropbox/proj/extend_tess/results/visualize_survey_designs/EM3_BRAINSTORM'
    cachepath = join(outdir, f'cache_{namestr}_{sstr}.npz')
    if os.path.exists(cachepath):
        print(f'loading {cachepath}')
        # Load arrays from the .npz file
        with np.load(cachepath) as data:
            RA = data['RA']
            Dec = data['Dec']
            nobs = data['nobs']

    else:

        # If you increase the resolution this will take longer to calculate.
        #RA, Dec = np.mgrid[:360:2000j, -90:90:1201j] # default takes 10 minutes

        # about 1.7M points ; kind of slow
        ra = np.linspace(0, 360, 2*824)   # 824 points from 0 to 360 degrees
        dec = np.linspace(-90, 90, 2*513)  # 513 points from -90 to 90 degrees
        RA, Dec = np.meshgrid(ra, dec, indexing='ij')

        #RA, Dec = np.mgrid[:360:824j, -90:90:513j]  # takes 30 sec laptop
        #RA, Dec = np.mgrid[:360:100j, -90:90:100j]  # takes 30 sec laptop

        # Array to accumulate number of observations
        nobs = np.zeros(RA.shape, dtype=float)
        # Loop through all the ra, dec and roll of the pointings
        for ra, dec, roll in tqdm(hypot_pointings, desc='Pointing', leave=True, position=0):
            # Loop through each camera
            for camera in np.arange(1, 5):
                # Loop through each CCD
                for ccd in np.arange(1, 5):
                    wcs = tesswcs.WCS.predict(ra, dec, roll, camera, ccd)
                    mask = get_observability_mask(wcs, SkyCoord(RA, Dec, unit='deg')).astype(int)
                    nobs += mask

        np.savez_compressed(cachepath, RA=RA, Dec=Dec, nobs=nobs)
        print(f'saved {cachepath}')

    rcParams['text.usetex'] = True
    rcParams['font.family'] = 'sans-serif'
    rcParams['font.sans-serif'] = 'Arial'
    rcParams['text.latex.preamble'] = r'\usepackage{helvet}\renewcommand{\familydefault}{\sfdefault}'

    plt.close('all')

    if projection == 'rect':
        fig, ax = plt.subplots(figsize=(4,3), dpi=300)

    elif projection == 'mollweide':
        fig = plt.figure(figsize=(6, 3))
        ax = fig.add_subplot(111, projection='mollweide', facecolor='White')

    #
    # Make the colormap
    #
    if sector_min >= 97:
        colors = [#"#ffffff", # N=0 white
                  "#84ccff", # N=1 pale blue
                  "#35aaff", # N=2 a little more saturated
                  "#126199", # N=3-12 a little more saturated
                  "#0a3a5c", # little more sat
                  "#000000" # N>=13 black
                 ]
        bounds = [0.5, 1.5, 2.5, 3.5, 4.5, 39.5 ]
        ylabels = ['0', '1', '2', '3', '4', '$\geq \!$ 5']
        ticks = [0, 1, 2, 3, 4, (4.5+39.5)/2]
    else:
        colors = [#"#ffffff", # N=0 white
                  "#84ccff", # N=1 pale blue
                  "#35aaff", # N=2 a little more saturated
                  "#126199", # N=3-12 a little more saturated
                  "#0a3a5c", # little more sat
                  "#000000" # N>=13 black
                 ]
        bounds = [0.5, 1.5, 2.5, 3.5, 12.5, 39.5 ]
        ylabels = ['0', '1', '2', '3', '4-12', '$\geq \!$ 13']
        ticks = [0, 1, 2, 3, (3.5+12.5)/2, (12.5+39.5)/2]

    cmap = LinearSegmentedColormap.from_list('my_cmap', colors, N=len(colors))
    cmap.set_extremes(under='white')
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    #
    # Plot the points & instantiate colorbar
    #
    if projection == 'rect':

        im = ax.pcolormesh(RA, Dec, nobs, cmap=cmap, norm=norm, shading='nearest')

        cbar = fig.colorbar(
            im, ax=ax, cmap=cmap, norm=norm, boundaries=bounds, fraction=0.037,
            pad=0.03, ticks=ticks, orientation='vertical'
        )
        cbar.ax.set_yticklabels(ylabels)
        cbar.ax.tick_params(direction='in')
        cbar.set_label('Months observed', rotation=270, labelpad=10)

        tick_labels = np.arange(0, 361, 60)
        tick_locs = tick_labels
        ax.set_xticks(tick_locs)
        ax.set_xticklabels(tick_labels, fontsize='small')
        ax.set_yticks([-60, -30, 0, 30, 60])
        ax.set_yticklabels(['-60', '-30', '0', '30', '60'], fontsize='small')

    elif projection == 'mollweide':

        # I tried pcolormesh, contourf, etc-similar style methods, and struggled.
        # Part of the issue is the non-rectangular projection.  But even using
        # cartopy's methods that are built around this, I would get segfaults!
        # So, fall back to the old method: ax.scatter layers.

        # Flatten the arrays for plotting with scatter
        RA_flat = RA.flatten(order='F')
        Dec_flat = Dec.flatten(order='F')
        nobs_flat = nobs.flatten(order='F')

        # Convert RA and Dec to radians for plotting
        # Shift RA by 180 degrees to center the projection
        RA_rad = np.deg2rad(RA_flat - 180)
        Dec_rad = np.deg2rad(Dec_flat)

        # Plot the data using scatter
        do_naive_scatter = 0
        if do_naive_scatter:
            sc = ax.scatter(RA_rad, Dec_rad, c=nobs_flat,
                            s=0.1,  # Adjust marker size for desired effect
                            cmap=cmap, norm=norm, marker='o', rasterized=True)
        else:

            fn = lambda ra_rad: -ra_rad

            # dummy to instantiate colorbar
            sc = ax.scatter(fn(RA_rad[::100]), Dec_rad[::100], c=nobs_flat[::100],
                            s=0.001, cmap=cmap, norm=norm, marker='o', lw=0,
                            rasterized=True)

            # main plotter: successive scatter layer
            size = 0.08
            max_cv = np.max(nobs_flat)
            for ix, cv in enumerate(np.sort(np.unique(nobs_flat))):
                if cv == 0:
                    continue
                sel = nobs_flat == cv
                zorder = int(- max_cv - 1 + ix)
                _ = ax.scatter(fn(RA_rad[sel]), Dec_rad[sel],
                               c=nobs_flat[sel], s=size,
                               lw=0, zorder=zorder, cmap=cmap, norm=norm,
                               marker='o',
                               rasterized=True)

            sel = nobs_flat > 0
            _ = ax.scatter(fn(RA_rad[~sel]), Dec_rad[~sel],
                           c=nobs_flat[~sel],
                           s=size/4, lw=0, zorder=-999, cmap=cmap, norm=norm,
                           marker='s', rasterized=True)

        # Roland's convention
        #ax.set_longitude_offset(180)

        cbar = fig.colorbar(
            sc, ax=ax, cmap=cmap, norm=norm, boundaries=bounds, fraction=0.05,
            pad=0.12, ticks=ticks, orientation='horizontal'
        )
        cbar.ax.set_xticklabels(ylabels, fontsize='small')
        cbar.ax.tick_params(direction='in', length=0)
        cbar.set_label('Months observed', rotation=0, labelpad=5)

        tick_label_locs = np.arange(-150, 181, 60)
        tick_locs = np.deg2rad(tick_label_locs)
        ax.set_xticks(tick_locs)
        ax.set_yticks(np.deg2rad([-60, -30, 0, 30, 60]))
        if showticklabels:
            tick_labels = ['210', '270', '330', '30', '90', '150']
            ax.set_xticklabels(tick_labels, fontsize='small')
            ax.set_yticklabels(['-60', '-30', '0', '30', '60'], fontsize='small')
        else:
            ax.set_xticklabels([], fontsize='small')
            ax.set_yticklabels([], fontsize='small')

        ax.grid(color='lightgray', linestyle='--', linewidth=0.25, zorder=-3,
                alpha=0.3)

    if showtitles:
        ax.set(xlabel='RA [deg]', ylabel='Dec [deg]',
               title=f"TESS Sectors {sector_min}-{sector_max}")

    outpath = join(
        outdir,
        f'{projection}_{namestr}_{sstr}.png'
    )
    fig.savefig(outpath, dpi=600, bbox_inches='tight')
    print(f"Wrote {outpath}")

    # get sky coverage
    skyfrac = 100*(nobs!=0).sum()/np.prod(nobs.shape)
    print(f'{namestr}_{sstr}: {skyfrac:.3f}% of sky observed')

if __name__ == "__main__":

    #names = ['em3_v00', 'em3_v02', 'cumul_em3_v00', 'em3_v01', 'em3_v03', 'em2_v09c']
    names = ['cumul_pm', 'cumul_em1', 'cumul_em2', 'cumul_em3_v00']
    names = ['em3_v04', 'cumul_em2', 'cumul_em3_v04'] # jnw's requests

    for name in names:
        viz_obs(name, projection='mollweide')
        #viz_obs(name, projection='lambert')
        #viz_obs(name, projection='rect')

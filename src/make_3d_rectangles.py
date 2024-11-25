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
from matplotlib.colors import LightSource

from convert_sc_orientation import ecliptic_to_equatorial_orientation

def _shift_lon(lon, origin):
    x = np.array(np.remainder(lon+360-origin,360)) # shift lon values
    ind = x>180
    x[ind] -=360    # scale conversion to [-180, 180]
    x=-x    # reverse the scale: East to the left
    return x

def plot_rectangles_on_sphere(
    pointing, projection, view_kwargs={}
    ):

    assert pointing in ['standard', 'c3p0', '70-40', 'everyother']

    if pointing == 'standard':
        sector_min, sector_max, dsector = 20, 23, 1
    elif pointing == 'c3p0':
        sector_min, sector_max, dsector = 74, 77, 1
    elif pointing == '70-40':
        sector_min, sector_max, dsector = 101, 104, 1
    elif pointing == 'everyother':
        sector_min, sector_max, dsector = 1, 4, 2

    sstr = f'S{str(sector_min).zfill(4)}_S{str(sector_max).zfill(4)}'

    if pointing in ['standard', 'c3p0', 'everyother']:

        real_pointings = tesswcs.pointings[['RA', "Dec", "Roll"]].to_pandas().values

        sel_pointings = real_pointings[
            slice(sector_min-1, sector_max, dsector),:
        ]

        if pointing == 'everyother':
            sel_pointings = np.vstack([sel_pointings, sel_pointings])

    elif pointing == '70-40':

        csvname = 'luke_S2+4_70_16_e2_G2_C11_G3_e2_roll1_editdup'
        csvpath = f'/Users/luke/Dropbox/proj/extend_tess/data/{csvname}.csv'
        df = pd.read_csv(csvpath)

        ras, decs, rolls = [],[],[]
        for sectornum in range(sector_min, sector_max+1, dsector):
            sc_roll_key = 'sc_eroll' if 'sc_eroll' in df else 'sc_roll'
            t = tuple(df.loc[df.S==sectornum, ['sc_elon', 'sc_elat', sc_roll_key]].iloc[0])
            eq = ecliptic_to_equatorial_orientation(t[0], t[1], t[2])
            ras.append(eq[0])
            decs.append(eq[1])
            rolls.append(eq[2])

            # All the hypothetical TESS pointings
            sel_pointings = np.vstack([ras, decs, rolls]).T

    outdir = (
        '/Users/luke/Dropbox/proj/extend_tess/results/'+
        'EM3_3d_rectangle_pitch'
    )
    if not os.path.exists(outdir): os.mkdir(outdir)
    cachepath = join(outdir, f'cache_{pointing}_{sstr}.npz')
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

        DEBUG = False

        if not DEBUG:
            # about 1.7M points ; kind of slow
            ra = np.linspace(0, 360, 3*824)   # 824 points from 0 to 360 degrees
            dec = np.linspace(-90, 90, 3*513)  # 513 points from -90 to 90 degrees
            RA, Dec = np.meshgrid(ra, dec, indexing='ij')
        else:
            ra = np.linspace(0, 360, 100)   # 824 points from 0 to 360 degrees
            dec = np.linspace(-90, 90, 101)  # 513 points from -90 to 90 degrees
            RA, Dec = np.meshgrid(ra, dec, indexing='ij')

        # Array to accumulate number of observations
        nobs = np.zeros(RA.shape, dtype=float)
        # Loop through all the ra, dec and roll of the pointings
        for ra, dec, roll in tqdm(sel_pointings, desc='Pointing', leave=True, position=0):
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

    elif projection == '3d':
        fig = plt.figure(figsize=(4, 5), dpi=300)
        ax = fig.add_subplot(111, projection='3d')
        fig.patch.set_facecolor('white')
        ax.set_facecolor('white')

    #
    # Make the colormap
    #
    colors = [#"#ffffff", # N=0 white
              "#84ccff", # N=1 pale blue
              "#35aaff", # N=2 a little more saturated
              "#126199", # N=3-12 a little more saturated
              "#0a3a5c", # little more sat
              #"#000000" # N>=13 black
             ]
    bounds = [0.5, 1.5, 2.5, 3.5, 4.5]#, 39.5 ]
    ylabels = ['0', '1', '2', '3', '4']#-12', '$\geq \!$ 13']
    ticks = [0, 1, 2, 3, 4]#(3.5+12.5)/2, (12.5+39.5)/2]

    cmap = LinearSegmentedColormap.from_list('my_cmap', colors, N=len(colors))
    cmap.set_extremes(under='white')
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

    ###########
    # Process #
    ###########
    # Flatten the arrays for plotting with scatter
    RA_flat = RA.flatten(order='F')
    Dec_flat = Dec.flatten(order='F')
    nobs_flat = nobs.flatten(order='F')

    coords = SkyCoord(RA_flat, Dec_flat, unit='deg')
    elon_flat = coords.barycentrictrueecliptic.lon.value
    elat_flat = coords.barycentrictrueecliptic.lat.value

    ECLIPTIC_COORDS = True

    # Convert RA and Dec to radians for plotting
    # Shift RA by 180 degrees to center the projection
    if not ECLIPTIC_COORDS:
        RA_rad = np.deg2rad(RA_flat - 180)
        Dec_rad = np.deg2rad(Dec_flat)
    else:
        RA_rad = np.deg2rad(elon_flat- 180)
        Dec_rad = np.deg2rad(elat_flat)

        if pointing in ['70-40', 'everyother']:
            Dec_rad *= -1.  # sign flip based on available pointings

    #
    # Plot the points & instantiate colorbar
    #
    sstr = ''

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

        # Plot the data using layered scatter
        LAYERED = 1
        if LAYERED:
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
        showticklabels = False
        if showticklabels:
            tick_labels = ['210', '270', '330', '30', '90', '150']
            ax.set_xticklabels(tick_labels, fontsize='small')
            ax.set_yticklabels(['-60', '-30', '0', '30', '60'], fontsize='small')
        else:
            ax.set_xticklabels([], fontsize='small')
            ax.set_yticklabels([], fontsize='small')

        ax.grid(color='lightgray', linestyle='--', linewidth=0.25, zorder=-3,
                alpha=0.3)


    elif projection == '3d':

        # Plot the unit sphere
        u, v = np.mgrid[0:2*np.pi:500j, 0:np.pi:250j]
        r = 0.99
        x = r*np.cos(u)*np.sin(v)
        y = r*np.sin(u)*np.sin(v)
        z = r*np.cos(v)
        #ax.plot_surface(x, y, z, color='white', alpha=0.03, linewidth=0)

        # Get points
        x = np.cos(Dec_rad) * np.cos(RA_rad)
        y = np.cos(Dec_rad) * np.sin(RA_rad)
        z = np.sin(Dec_rad)

        # dummy for colorbar
        sc = ax.scatter(x[::100], y[::100], z[::100], c=nobs_flat[::100],
                        s=0.0, cmap=cmap, norm=norm, marker='o', lw=0,
                        rasterized=True)

        # main plotter: successive scatter layer
        size = 0.05
        sel = nobs_flat != 0
        _ = ax.scatter(x[sel], y[sel], z[sel], c=nobs_flat[sel], s=size,
                       lw=0, cmap=cmap, norm=norm,
                       marker='o', rasterized=True)

        #max_cv = np.max(nobs_flat)
        #for ix, cv in enumerate(np.sort(np.unique(nobs_flat))):
        #    if cv == 0:
        #        continue
        #    sel = nobs_flat == cv
        #    #zorder = int(- max_cv - 1 + ix)
        #    #_ = ax.scatter(x[sel], y[sel], z[sel], c=nobs_flat[sel], s=size,
        #    #               lw=0, zorder=zorder, cmap=cmap, norm=norm,
        #    #               marker='o', rasterized=True)

        #    xs, ys, zs, cs = x[sel], y[sel], z[sel], nobs_flat[sel]

        #    _ = ax.scatter(xs, ys, zs, c=cs, s=size,
        #                   lw=0, cmap=cmap, norm=norm,
        #                   marker='o', rasterized=True)


        #sel = nobs_flat > 0
        #_ = ax.scatter(fn(RA_rad[~sel]), Dec_rad[~sel],
        #               c=nobs_flat[~sel],
        #               s=size/4, lw=0, zorder=-999, cmap=cmap, norm=norm,
        #               marker='s', rasterized=True)

        # Plot equator
        azim = view_kwargs['azim']
        elev = view_kwargs['elev']
        theta = np.linspace(
            -np.pi/2, +np.pi/2, 1000
        ) + np.deg2rad(azim)
        if elev == 90:
            theta = np.linspace(0, 2*np.pi, 1000)
        r = 1.01
        x_eq = r*np.cos(theta)
        y_eq = r*np.sin(theta)
        z_eq = r*np.zeros_like(theta)
        ax.plot(x_eq, y_eq, z_eq, color='lightgray', linewidth=0.5, ls=':',
                zorder=1)

        # Plot rim
        if elev != 90:
            latitudes = np.linspace(-np.pi/2, np.pi/2, 1000)
            longitudes = [
                np.deg2rad(azim) + np.pi/2,
                np.deg2rad(azim) - np.pi/2,
            ]
            for longitude in longitudes:
                r = 1.01
                xr = r * np.cos(latitudes) * np.cos(longitude)
                yr = r * np.cos(latitudes) * np.sin(longitude)
                zr = r * np.sin(latitudes)
                ax.plot(xr, yr, zr, linewidth=0.5, ls='-', color='lightgray',
                        zorder=0)

        ax.set_xlim([-1.05, 1.05])
        ax.set_ylim([-1.05, 1.05])
        ax.set_zlim([-1.05, 1.05])

        # Adjust the view
        ax.set_box_aspect([1, 1, 544/515])  # Equal aspect ratio
        #ax.set_box_aspect([1, 1, 1])  # Equal aspect ratio
        ax.axis('off')

        # Apply view kwargs
        if 'elev' in view_kwargs and 'azim' in view_kwargs:
            ax.view_init(elev=view_kwargs['elev'], azim=view_kwargs['azim'])
            elev = view_kwargs['elev']
        else:
            ax.view_init(elev=20, azim=90)
            elev = 20


        cbar = fig.colorbar(
            sc, ax=ax, cmap=cmap, norm=norm, boundaries=bounds, fraction=0.05,
            pad=0.12, ticks=ticks, orientation='horizontal'
        )
        cbar.ax.set_xticklabels(ylabels, fontsize='small')
        cbar.ax.tick_params(direction='in', length=0)
        cbar.set_label('Months observed', rotation=0, labelpad=5)

        sstr = f'_elev{elev}'

        fig.tight_layout()

    showtitles = False
    if showtitles:
        ax.set(xlabel='RA [deg]', ylabel='Dec [deg]',
               title=f"TESS Sectors {sector_min}-{sector_max}")

    outpath = join(
        outdir,
        f'{projection}_{pointing}{sstr}.png'
    )

    fig.savefig(outpath, dpi=600, bbox_inches='tight')
    print(f"Wrote {outpath}")


if __name__ == "__main__":

    pointings = ['standard']#, 'c3p0', '70-40', 'everyother']
    pointings = ['standard', 'c3p0', '70-40', 'everyother']
    projection = '3d' # or 'mollweide'

    elevs = [20, 90]

    azim_dict = {
        'standard': 310,
        'c3p0': 290,
        '70-40': 0,
        'everyother': 160,
    }

    for pointing in pointings:
        for elev in elevs:
            view_kwargs = {'elev':elev, 'azim':azim_dict[pointing]}
            plot_rectangles_on_sphere(pointing, projection, view_kwargs=view_kwargs)

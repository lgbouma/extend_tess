# -*- coding: utf-8 -*-

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

def _shift_lon_get_x(lon, origin):
    x = np.array(np.remainder(lon+360-origin,360)) # shift lon values
    ind = x>180
    x[ind] -=360    # scale conversion to [-180, 180]
    x=-x    # reverse the scale: East to the left
    return x


def plot_mwd(lon,dec,color_val,origin=0,size=3,title='Mollweide projection',
             projection='mollweide',savdir='../results/',savname='mwd_0.pdf',
             overplot_galactic_plane=True, is_tess=False, is_radec=None,
             cbarbounds=None):

    '''
    args, kwargs:

        lon, lat are arrays of same length. they can be (RA,dec), or (ecliptic
            long, ecliptic lat). lon takes values in [0,360), lat in [-90,90],

        is_radec: mandatory. True if (RA,dec), else elong/elat.

        title is the title of the figure.

        projection is the kind of projection: 'mollweide', 'aitoff', ...

    comments: see
    http://balbuceosastropy.blogspot.com/2013/09/the-mollweide-projection.html.
    '''
    if is_radec == None:
        raise AssertionError

    # for matplotlib mollweide projection, x and y values (usually RA/dec, or
    # lat/lon) must be in -pi<x<pi, -pi/2<y<pi/2.
    # In astronomical coords, RA increases east (left on celestial charts).
    # Here, the default horizontal scale has positive to the right.

    x = _shift_lon_get_x(lon, origin)

    fig = plt.figure(figsize=(6, 4.5))
    ax = fig.add_subplot(111, projection=projection, facecolor='White')

    if is_tess:
        # set up colormap
        import seaborn as sns
        #rgbs = sns.color_palette("cubehelix_r", n_colors=17)
        #rgbs = sns.color_palette('Paired', n_colors=17, desat=0.9)
        #rgbs = sns.color_palette('Set2', n_colors=17)
        rgbs = sns.color_palette('viridis', n_colors=17)
        cmap = mpl.colors.ListedColormap(rgbs)
        if isinstance(cbarbounds,np.ndarray):
            bounds=cbarbounds
        else:
            bounds = np.arange(-27.32/2, np.max(df['obs_duration'])+1/2*27.32, 27.32)

        norm = mpl.colors.BoundaryNorm(bounds, cmap.N)

        # plot the stars
        cax = ax.scatter(np.radians(x),np.radians(dec), c=color_val, s=size,
                         lw=0, zorder=2, cmap=cmap, norm=norm, rasterized=True)

        # set up colorbar
        cbar = fig.colorbar(cax, cmap=cmap, norm=norm, boundaries=bounds,
                            fraction=0.025, pad=0.03, #ticks=np.arange(13)+1,
                            orientation='vertical')

        # ylabels = np.arange(1,18,1)*27.32

        # cbar.ax.set_yticklabels(map(str, ylabels))
        cbar.set_label('days observed', rotation=270, labelpad=10)
        cbar.ax.tick_params(direction='in')

        # label each sector. this involves computing the positions first...
        #TODO

    else:
        ax.scatter(np.radians(x),np.radians(dec), c=color_val, s=size,
                   zorder=2)


    if overplot_galactic_plane:

        ##########
        # make many points, and also label the galactic center. ideally you
        # will never need to follow these coordinate transformations.
        glons = np.arange(0,360)
        glats = np.zeros_like(glons)
        coords = SkyCoord(glons*u.degree, glats*u.degree, frame='galactic')
        gplane_ra, gplane_dec = coords.icrs.ra.value, coords.icrs.dec.value
        gplane_elon = coords.barycentrictrueecliptic.lon.value
        gplane_elat = coords.barycentrictrueecliptic.lat.value
        if is_radec:
            gplane_x = _shift_lon_get_x(gplane_ra, origin)
        else:
            gplane_x = _shift_lon_get_x(gplane_elon, origin)
            gplane_dec = gplane_elat
        ax.scatter(np.radians(gplane_x),np.radians(gplane_dec),
                   c='lightgray', s=2, zorder=3)
        gcenter = SkyCoord('17h45m40.04s', '-29d00m28.1s', frame='icrs')
        gcenter_ra, gcenter_dec = gcenter.icrs.ra.value, gcenter.icrs.dec.value
        gcenter_elon = gcenter.barycentrictrueecliptic.lon.value
        gcenter_elat = gcenter.barycentrictrueecliptic.lat.value
        if is_radec:
            gcenter_x = _shift_lon_get_x(np.array(gcenter_ra), origin)
        else:
            gcenter_x = _shift_lon_get_x(np.array(gcenter_elon), origin)
            gcenter_dec = gcenter_elat
        ax.scatter(np.radians(gcenter_x),np.radians(gcenter_dec),
                   c='black', s=2, zorder=4, marker='X')
        ax.text(np.radians(gcenter_x), np.radians(gcenter_dec), 'GC',
                fontsize='x-small', ha='left', va='top')
        ##########


    tick_labels = np.array([150, 120, 90, 60, 30, 0, 330, 300, 270, 240, 210])
    tick_labels = np.remainder(tick_labels+360+origin,360)
    ax.set_xticklabels(tick_labels, fontsize='x-small')
    ax.set_yticklabels(np.arange(-75,75+15,15), fontsize='x-small')

    ax.set_title(title, y=1.05, fontsize='small')
    if is_radec:
        ax.set_xlabel('ra', fontsize='x-small')
        ax.set_ylabel('dec', fontsize='x-small')
    else:
        ax.set_xlabel('ecl lon', fontsize='x-small')
        ax.set_ylabel('ecl lat', fontsize='x-small')
    ax.grid(color='lightgray', linestyle='--', linewidth=0.5, zorder=-1)

    ax.text(0.99,0.01,'github.com/lgbouma/extend_tess',
            fontsize='xx-small',transform=ax.transAxes,
            ha='right',va='bottom')
    fig.tight_layout()
    fig.savefig(os.path.join(savdir,savname),dpi=350, bbox_inches='tight')
    print('saved {}'.format(os.path.join(savdir,savname)))



def get_n_observations(dirnfile, outpath, n_stars):

    np.random.seed(42)

    # pick points uniformly on celestial sphere. they don't strictly need to be
    # random.  takes >~1 minute to draw random numbers after ~2*10^5. faster to
    # just do it on an appropriate grid.

    # e.g., http://mathworld.wolfram.com/SpherePointPicking.html
    uniform0 = np.linspace(0,1,n_stars)
    uniform1 = np.linspace(0,1,n_stars)
    #rand0 = np.random.uniform(low=0,high=1,size=n_stars)
    #rand1 = np.random.uniform(low=0,high=1,size=n_stars)

    theta = (2*np.pi*uniform0 * u.rad).to(u.deg).value
    phi = (np.arccos(2*uniform1 - 1) * u.rad).to(u.deg).value - 90

    ras = theta*u.deg
    decs = phi*u.deg

    coords = SkyCoord(ra=ras, dec=decs, frame='icrs')

    df = pd.read_csv(dirnfile, sep=';')

    lats = nparr([
        nparr(df['cam1_elat']),
        nparr(df['cam2_elat']),
        nparr(df['cam3_elat']),
        nparr(df['cam4_elat'])]).T
    lons = nparr([
        nparr(df['cam1_elon']),
        nparr(df['cam2_elon']),
        nparr(df['cam3_elon']),
        nparr(df['cam4_elon'])]).T

    cam_directions = []
    for lat, lon in zip(lats, lons):

        c1lat,c2lat,c3lat,c4lat = lat[0],lat[1],lat[2],lat[3]
        c1lon,c2lon,c3lon,c4lon = lon[0],lon[1],lon[2],lon[3]

        this_cam_dirn = [(c1lat, c1lon),
                         (c2lat, c2lon),
                         (c3lat, c3lon),
                         (c4lat, c4lon)]

        cam_directions.append(this_cam_dirn)

    df['camdirection'] = cam_directions

    n_observations = np.zeros_like(coords)

    for ix, row in df.iterrows():

        print(row['start'])
        cam_direction = row['camdirection']

        onchip = gcgss(coords, cam_direction, verbose=False)

        n_observations += onchip

    outdf = pd.DataFrame({'ra':coords.ra.value,
                          'dec':coords.dec.value,
                          'elat':coords.barycentrictrueecliptic.lat.value,
                          'elon':coords.barycentrictrueecliptic.lon.value,
                          'n_observations': n_observations })
    outdf[['ra','dec','elon','elat','n_observations']].to_csv(
        outpath, index=False, sep=';')
    print('saved {}'.format(outpath))


if __name__=="__main__":

    datadir = '../data/'
    savdir = '../results/visualize_survey_designs/'
    orbit_duration_days = 27.32 / 2

    # things to change
    filenames = ['idea_1_SN_ecliptic.csv',
                 'idea_2_SNSNS_hemi.csv',
                 'idea_3_SNNSN_hemi.csv']

    eclsavnames = ['idea_1_SN_ecliptic_eclmap.png',
                   'idea_2_SNSNS_hemi_eclmap.png',
                   'idea_3_SNNSN_eclmap.png']

    icrssavnames = ['idea_1_SN_ecliptic_icrsmap.png',
                    'idea_2_SNSNS_hemi_icrsmap.png',
                    'idea_3_SNNSN_icrsmap.png']

    titles = ['idea 1 SN->N(6)->ecliptic(10)->S(26)->N(remain)',
              'idea 2 SN->S(26)->N(28)->S(remain)',
              'idea 3 SN->N(26)->S(26)->N(remain)'  ]

    dirnfiles = [ os.path.join(datadir,fname) for fname in filenames]

    for dirnfile, eclsavname, icrssavname, title in zip(
        dirnfiles, eclsavnames, icrssavnames, titles):

        obsdpath = dirnfile.replace('.csv', '_coords_observed.csv')

        if not os.path.exists(obsdpath):
            get_n_observations(dirnfile, obsdpath, int(2e5))

        df = pd.read_csv(obsdpath, sep=';')
        df['obs_duration'] = orbit_duration_days*df['n_observations']

        plot_mwd(nparr(df['elon']),
                 nparr(df['elat']),
                 nparr(df['obs_duration']),
                 origin=0, size=2, title=title,
                 projection='mollweide', savdir=savdir,
                 savname=eclsavname,
                 overplot_galactic_plane=True, is_tess=True, is_radec=False,
                 cbarbounds=None)

        plot_mwd(nparr(df['ra']),
                 nparr(df['dec']),
                 nparr(df['obs_duration']),
                 origin=0, size=2, title=title,
                 projection='mollweide', savdir=savdir,
                 savname=icrssavname,
                 overplot_galactic_plane=True, is_tess=True, is_radec=True,
                 cbarbounds=None)

        df.to_csv(obsdpath, sep=';', index=False)
        print('rewrote {}'.format(obsdpath))

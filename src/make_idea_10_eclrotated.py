"""
per GRR request, by way of JNW.

"what if we did an observing strategy in which we aligned the TESS field
parallel to lines of celestial latitude?"
"""

import numpy as np, pandas as pd
from astropy.coordinates import SkyCoord

def elat_to_colatitude(elat):

    return -elat + 90

def colatitude_to_elat(colatitude):

    return 90 - colatitude

def rotate_about_x_axis(theta, phi, alpha):
    # following
    # https://math.stackexchange.com/questions/906095/rotating-a-point-in-spherical-coordinates-around-cartesian-axis
    #
    # takes: theta, phi, alpha in radians. theta and phi are the colatitude and
    # longitude. alpha is the rotation angle.
    #
    # returns: thetaprime, phiprime in radians, rotated by alpha radians about
    # the x axis.

    def Q(theta, phi):

        return np.array([ np.cos(theta/2) + 0j, np.exp(1j*phi) * np.sin(theta/2) ])

    def Qinv(z_1, z_2):

        arg_z1 = np.arctan( z_1.imag / z_1.real)
        arg_z2 = np.arctan( z_2.imag / z_2.real)

        mod_z1 = np.abs(z_1)
        mod_z2 = np.abs(z_2)

        return np.array([ 2*np.arctan(mod_z2/mod_z1), arg_z2 - arg_z1])

    R_x = np.array( [
        [np.cos(alpha/2)+0j, -1j*np.sin(alpha/2)],
        [-1j*np.sin(alpha/2), np.cos(alpha/2)+0j]
    ])

    output_1 = np.matmul( R_x, Q(theta,phi) )
    output = Qinv(output_1[0], output_1[1])

    return output


def main():

    df = pd.read_csv('../data/idea_2_SNSNS_hemi.csv', sep=';')

    # start with "idea 2 SNSNS" to get basic coordinates.
    inds = np.array((df['orbit'] >=61) & (df['orbit'] < 61+26)).astype(bool)

    edf = df.iloc[inds]
    restdf = df.iloc[~inds]

    # populate easy things. sc boresight latitudes all go to zero. sc boresight
    # longitudes stay the same. sc eroll to 23.4 degrees.
    obliq_of_ecliptic = 23.4

    edf['spacecraft_boresight_elat'] = 0
    edf['spacecraft_eroll'] = obliq_of_ecliptic

    rotelats, rotelons = [], []

    for cam in [1,2,3,4]:

        # get the first orbit's elat and elon, for each camera.
        elat = -36 + (cam-1)*24
        elon = 0

        # for math to work, convert ecliptic latitude to a colatitude.
        theta = np.deg2rad(elat_to_colatitude(elat))
        phi = np.deg2rad(elon)

        # do rotation convert between radians and degrees.
        alpha = np.deg2rad(obliq_of_ecliptic)

        thetaprime,phiprime = rotate_about_x_axis(theta, phi, alpha)

        thetaprime,phiprime = np.rad2deg(thetaprime), np.rad2deg(phiprime)

        elatprime = colatitude_to_elat(thetaprime)
        elonprime = phiprime

        rotelons.append(elonprime)
        rotelats.append(elatprime)

    # the offset in ecliptic longitude is equal to whatever the starting
    # ecliptic longitude is
    elon_offset = edf['spacecraft_boresight_elon'].iloc[0]

    initial_elon = rotelons + elon_offset

    # NOTE this line says: "advance by 360/13 ~= 27 degrees in ecliptic
    # longitude every 2 orbits".
    elonarr = (
        initial_elon[:,None] + (360/13.)*np.repeat(np.arange(0,13),2)[None,:]
    )

    # finally, write to edf
    edf['cam1_elat'] = rotelats[0]
    edf['cam2_elat'] = rotelats[1]
    edf['cam3_elat'] = rotelats[2]
    edf['cam4_elat'] = rotelats[3]

    edf['cam1_elon'] = np.mod(elonarr[0,:],360)
    edf['cam2_elon'] = np.mod(elonarr[1,:],360)
    edf['cam3_elon'] = np.mod(elonarr[2,:],360)
    edf['cam4_elon'] = np.mod(elonarr[3,:],360)

    # merge edf and restdf
    outdf = pd.concat((edf, restdf))

    outpath = '../data/idea_10_eclrotated.csv'
    outdf.to_csv(outpath, index=False, sep=';')
    print('made {}'.format(outpath))


if __name__=="__main__":
    main()

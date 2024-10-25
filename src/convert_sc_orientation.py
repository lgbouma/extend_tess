import numpy as np
from astropy.coordinates import SkyCoord, Angle, Longitude, Latitude
import astropy.units as u

def ecliptic_to_equatorial_orientation(sc_elon_deg, sc_elat_deg, sc_eroll_deg):
    # Constants
    epsilon_deg = 23.43929111  # Obliquity of the ecliptic
    epsilon = np.radians(epsilon_deg)

    # Convert input angles to radians
    sc_elon_rad = np.radians(sc_elon_deg)
    sc_elat_rad = np.radians(sc_elat_deg)
    sc_eroll_rad = np.radians(sc_eroll_deg)

    # Boresight vector in ecliptic Cartesian coordinates
    x_ecl = np.cos(sc_elon_rad) * np.cos(sc_elat_rad)
    y_ecl = np.sin(sc_elon_rad) * np.cos(sc_elat_rad)
    z_ecl = np.sin(sc_elat_rad)
    Z_ecl = np.array([x_ecl, y_ecl, z_ecl])

    # Vector towards increasing ecliptic latitude (north)
    Nx_ecl = -np.cos(sc_elon_rad) * np.sin(sc_elat_rad)
    Ny_ecl = -np.sin(sc_elon_rad) * np.sin(sc_elat_rad)
    Nz_ecl = np.cos(sc_elat_rad)
    N_ecl = np.array([Nx_ecl, Ny_ecl, Nz_ecl])
    N_ecl_perp = N_ecl - np.dot(N_ecl, Z_ecl) * Z_ecl
    N_ecl_perp /= np.linalg.norm(N_ecl_perp)

    # Vector towards increasing ecliptic longitude (east)
    Ex_ecl = -np.sin(sc_elon_rad) * np.cos(sc_elat_rad)
    Ey_ecl = np.cos(sc_elon_rad) * np.cos(sc_elat_rad)
    Ez_ecl = 0.0
    E_ecl = np.array([Ex_ecl, Ey_ecl, Ez_ecl])
    E_ecl_perp = E_ecl - np.dot(E_ecl, Z_ecl) * Z_ecl
    E_ecl_perp /= np.linalg.norm(E_ecl_perp)

    # Spacecraft's +Y axis in ecliptic coordinates
    Y_ecl = np.cos(sc_eroll_rad) * N_ecl_perp - np.sin(sc_eroll_rad) * E_ecl_perp

    # Rotation matrix from ecliptic to equatorial coordinates
    R = np.array([
        [1, 0, 0],
        [0, np.cos(epsilon), -np.sin(epsilon)],
        [0, np.sin(epsilon), np.cos(epsilon)]
    ])

    # Transform vectors to equatorial coordinates
    Z_eq = R @ Z_ecl
    Y_eq = R @ Y_ecl

    # Offset along +Y axis
    delta = 1e-6  # Small angle offset
    P_eq = Z_eq + delta * Y_eq
    P_eq /= np.linalg.norm(P_eq)

    # Convert boresight vector to RA, Dec
    sc_ra_rad = np.arctan2(Z_eq[1], Z_eq[0])
    sc_dec_rad = np.arcsin(Z_eq[2])
    sc_ra_deg = np.degrees(sc_ra_rad) % 360
    sc_dec_deg = np.degrees(sc_dec_rad)

    # Convert offset vector to RA, Dec
    offset_ra_rad = np.arctan2(P_eq[1], P_eq[0])
    offset_dec_rad = np.arcsin(P_eq[2])

    # Create SkyCoord objects
    boresight_coord = SkyCoord(ra=sc_ra_deg * u.deg,
                               dec=sc_dec_deg * u.deg,
                               frame='icrs')

    offset_coord = SkyCoord(ra=np.degrees(offset_ra_rad) * u.deg,
                            dec=np.degrees(offset_dec_rad) * u.deg,
                            frame='icrs')

    # Compute roll angle
    sc_roll_deg = 360 - boresight_coord.position_angle(offset_coord).deg

    return sc_ra_deg, sc_dec_deg, sc_roll_deg

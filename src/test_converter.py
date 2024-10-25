"""
test your ability to convert from sc boresight elon elat and eroll to
equatorial analog.

harder than it sounds.
"""
import os
from os.path import join
import numpy as np, pandas as pd
from convert_sc_orientation import ecliptic_to_equatorial_orientation
import tesswcs

def test_converted(elon, elat, eroll, sectornum):

    eq = ecliptic_to_equatorial_orientation(elon, elat, eroll)
    print(f'Passed in sc_elon_deg, sc_elat_deg, sc_eroll_deg: {elon} {elat} {eroll}')
    print(f'Got out: {eq[0]:.4f}, {eq[1]:.4f}, {eq[2]:.4f}')

    v = tuple(tesswcs.pointings.to_pandas()[['RA','Dec','Roll']].iloc[sectornum-1])
    print(f'Truth: {v[0]}, {v[1]}, {v[2]}')

    diff = v[2] - eq[2]
    print(f'diff roll: {diff:.2f}')

    print(42*'-')


csvpath = '/Users/luke/Dropbox/proj/extend_tess/data/em2_v09c.csv'
df = pd.read_csv(csvpath)

for sectornum in [1, 14, 70, 42, 52]:#, 1, 2, 3, 10, 14, 16, 17, 70]:

    print(f"S: {sectornum}")

    t = tuple(df.loc[df.S==sectornum, ['sc_elon', 'sc_elat', 'sc_eroll']].iloc[0])

    test_converted(t[0], t[1], t[2], sectornum)

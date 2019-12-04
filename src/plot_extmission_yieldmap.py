import numpy as np, pandas as pd, matplotlib.pyplot as plt
import os

from planet_yield_plots import get_data


names = ['idea_1_SNE-v5']
outdirs = [os.path.join('../results/planet_yield_plots/',dirname)
       for dirname in names]
datanames = [n.split('_')[-1] for n in names]
datapaths = [os.path.join( '../data/tommyb',
                          'detected_planet_catalog_{:s}.csv.bz2'.format(dn))
             for dn in datanames]

datapath = datapaths[0]
outdir = outdirs[0]

df = get_data(datapath)

#     ['Unnamed: 0', 'Unnamed: 0.1', 'TICID', 'RA', 'DEC', 'PLX', 'ECLONG',
#     'ECLAT', 'V', 'J', 'Ks', 'TESSMAG', 'TEFF', 'RADIUS', 'MASS', 'CONTRATIO',
#     'PRIORITY', 'isMdwarf', 'isGiant', 'isSubgiant', 'cosi', 'noise_level',
#     'Nplanets', 'planetRadius', 'planetPeriod', 'starID', 'T0', 'ars', 'ecc',
#     'omega', 'rprs', 'impact', 'duration', 'duration_correction',
#     'transit_depth', 'transit_depth_diluted', 'has_transits', '1', '2', '3',
#     '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16',
#     '17', '18', '19', '20', '21', '22', '23', '24', '25', '26', '27', '28',
#     '29', '30', '31', '32', '33', '34', '35', '36', '37', '38', '39', '40',
#     '41', '42', '43', '44', '45', '46', '47', '48', '49', '50', '51', '52',
#     '53', '54', '55', '56', '57', '58', '59', '60', '61', '62', '63', '64',
#     '65', 'isObserved', 'Ntransits', 'Ntransits_primary', 'SNR', 'SNR_primary',
#     'needed_for_detection', 'detected', 'needed_for_detection_primary',
#     'detected_primary', 'insol', 'inOptimisticHZ']

f,ax = plt.subplots()

pdf = df[df.detected_primary]
edf = df[df.detected]

ax.scatter(
    pdf.RA, pdf.DEC, color='k', label='detected in primary',
    zorder=2, s=3
)

ax.scatter(
    edf.RA, edf.DEC, color='C0', label='detected after extended',
    zorder=1, s=3
)

ax.legend(loc='best')

ax.set_xlabel('right ascension [degrees]')
ax.set_ylabel('declination [degrees]')

outpath = os.path.join(outdir, 'extmission_yieldmap.png')
f.savefig(outpath, dpi=300, bbox_inches='tight')
print('made {}'.format(outpath))

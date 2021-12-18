"""
We want to convert Vanderspek's output to my csv file format:
    "start";"end";"orbit";"sector";"spacecraft_boresight_elat";"spacecraft_boresight_elon";"spacecraft_eroll";
    "cam4_elat";"cam4_elon";"cam3_elat";"cam3_elon";"cam2_elat";"cam2_elon";"cam1_elat";"cam1_elon";"comment"
"""
import numpy as np, pandas as pd

# NOTE: only content that is pulled from this table is sector numbers, orbit
# numbers, and times; doesn't change with new pointings.
df0 = (
    pd.read_csv('/Users/luke/Dropbox/proj/extend_tess/data/20211013_vanderspek_EM2/Y1-7_ops.tbl',
                delim_whitespace=True, comment='#')
)
df0 = df0.rename(
    columns={'spacecraft_boresight_elat': 'sc_elat',
             'RA':'sc_ra',
             'Dec':'sc_dec',
             'Roll':'sc_roll'}
)
selcols = [
    'Start(UTC)', 'End(UTC)', 'S', 'O1', 'O2',
]
sdf = df0[selcols]

foo = pd.DataFrame(
    {'Start(UTC)':np.repeat('',23),'End(UTC)':np.repeat('',23),'S':np.arange(101,124,1),
     'O1':np.arange(209,209+23*2,2), 'O2':np.arange(210,210+23*2,2)}
)
sdf = pd.concat((sdf,foo)).reset_index(drop=True)

# OK. now, for each sector get the camera centers... and rolls[?]
vnumstr = '09'
vanderpath = (
    '/Users/luke/Dropbox/proj/extend_tess/data/20211013_vanderspek_EM2/'
    f'luke_Scenario_{vnumstr}.out'
)
with open(vanderpath) as f:
    lines = f.readlines()

def get_sc_cam_coords(sectornumber, lines):

    rnd = lambda x: np.round(x, 3)

    assert isinstance(sectornumber, int)

    d = {}

    for ix, l in enumerate(lines):

        head = l.split(' ')[0]

        if head == str(sectornumber):

            d['sc_elon'] = rnd(float(lines[ix+2].split(' ')[-3]))
            d['sc_elat'] = rnd(float(lines[ix+2].split(' ')[-2]))
            d['sc_eroll'] = rnd(float(lines[ix+2].split(' ')[-1]))

            for camera_index in range(1,5):
                d[f'cam{camera_index}_elon'] = (
                    rnd(float(lines[ix+2+camera_index].split(' ')[-3]))
                )
                d[f'cam{camera_index}_elat'] = (
                    rnd(float(lines[ix+2+camera_index].split(' ')[-2]))
                )

        continue

    return d

dictlist = [get_sc_cam_coords(SECTOR, lines) for SECTOR in range(1,124)]

outdf = pd.concat([sdf, pd.DataFrame(dictlist)], axis=1)

outdf['S'] = outdf['S'].astype(int)
outdf['O1'] = outdf['O1'].astype(int)
outdf['O2'] = outdf['O2'].astype(int)

outpath = f'../data/em2_v{vnumstr}.csv'
outdf.to_csv(outpath, index=False)
print(f"made {outpath}")

#
# standard imports
#
import numpy as np, pandas as pd, matplotlib.pyplot as plt
from astropy.time import Time

#
# custom imports:
# plotting imports from https://github.com/lgbouma/aesthetic
# download latest TOI catalog from TEV per https://github.com/lgbouma/cdips
#
from aesthetic.plot import savefig
from cdips.utils.catalogs import get_toi_catalog

# get the data, and clean the timestamps.
df = get_toi_catalog()

t = np.array(df.Alerted).astype(str)
t = [str(_t).replace('+0000','') for _t in t]
t = Time(t, format='iso')

# First day of Sector 27 (2020-07-04), https://tess.mit.edu/observations/.
# Allot 26 days + 1 month for acquisition and vetting.

t_EXT_START = Time('2020-07-04', format='iso')
t_TO_VET = 56 # days (1 month + 26 day sector processing)

# "EXT" targets are those from at least 56 days after EM started.
sel_EXT = ( t > (t_EXT_START + t_TO_VET) )

# The total time observed in "EXT"
t_interval = max(t) - (t_EXT_START + t_TO_VET)

# from https://tess.mit.edu/observations/, S1 started on 07/25, 2018. Add the
# same interval 
t_PRI_START = Time('2018-07-25', format='iso')

sel_PRI = (
    ( t > (t_PRI_START + t_TO_VET) )
    &
    ( t <= (t_PRI_START + t_TO_VET + t_interval) )
)

print(f'PRI: {min(t[sel_PRI])}  to  {max(t[sel_PRI])}  '
      f'({max(t[sel_PRI]) - min(t[sel_PRI])} days) ')
print(f'EXT: {min(t[sel_EXT])}  to  {max(t[sel_EXT])}  '
      f'({max(t[sel_EXT]) - min(t[sel_EXT])} days) ')

#
# make Tmag histograms
#
sdict = {
    'pri':sel_PRI,
    'ext':sel_EXT
}

fig, ax = plt.subplots(figsize=(4,3))

bins = np.arange(0,20,1)

i = 0
for k,s in sdict.items():
    ax.hist(df[s]['TMag Value'], bins=bins, cumulative=False,
            color=f'C{i}', fill=False,  histtype='step',
            linewidth=2, alpha=0.8, label=k)
    print(f'{k}: {len(df[s])} TOIs alerted, stats are:')
    print(df[s]['TMag Value'].describe())
    i += 1
ax.legend()
ax.set_xlabel('T [mag]')
ax.set_ylabel('Count')

outpath = '../results/year1_year3_comp/tmag_hist.pdf'
savefig(fig, outpath)


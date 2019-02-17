import matplotlib.pyplot as plt, numpy as np, pandas as pd

df_merged = pd.read_csv(
    '../data/idea_9_SNEshifted_coords_observed_merged_forproposal.csv', sep=';'
)
df_primary = pd.read_csv(
    '../data/primary_mission_coords_observed_forproposal.csv', sep=';'
)

# after EM1
plt.hist(df_merged['n_observations'],
         bins=np.arange(0,np.max(df_merged['n_observations'])))
plt.yscale('log')

plt.xlabel('n orbits observed')
plt.ylabel('n stars (of {})'.format(len(df_merged)))
plt.ylim((50,4e5))

plt.savefig('after_EM1.png',dpi=300)

# after primary
plt.close('all')
plt.hist(df_primary['n_observations'],
         bins=np.arange(0,np.max(df_merged['n_observations'])))
plt.yscale('log')

plt.xlabel('n orbits observed')
plt.ylabel('n stars (of {})'.format(len(df_merged)))
plt.ylim((50,4e5))

plt.savefig('after_primary.png',dpi=300)



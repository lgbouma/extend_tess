# -*- coding: utf-8 -*-

from __future__ import division, print_function

import numpy as np, pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt

def make_dataframe():
    raw_data = {
        'period_bin_days': ['<20','>50'],
        'C3PO*3': [2235, 291],
        'hemiS*3': [2573, 192],
        'allsky*3': [2842, 123]
    }

    df = pd.DataFrame(
        raw_data, columns = ['period_bin_days','C3PO*3','hemiS*3','allsky*3']
    )

    return df

def make_period_barchart(df, linscale=True, darkcolor=False):

    if darkcolor:
        plt.style.use('dark_background')
        # c3po, hemis, allsky
        colors = ['#FBE792','#91A187','#75CAF7']
    else:
        colors = ['#EE3224','#F78F1E','#FFC222']

    # Setting the positions and width for the bars
    pos = [0, 1]
    width = 0.25

    # Plotting the bars
    fig, ax = plt.subplots(figsize=(4, 3))

    ax.bar(pos, df['C3PO*3'], width, alpha=0.5, color=colors[0],
            label='C3PO*3')
    if not linscale:
        for p, ix in enumerate(pos):
            ax.text(p, 15,
                    list(map(str,np.array(df['C3PO*3'])))[ix],
                    fontsize='small',
                    ha='center', va='center')

    ax.bar([p+width for p in pos], df['hemiS*3'], width, alpha=0.5,
            color=colors[1], label='hemiS*3')
    if not linscale:
        for p, ix in enumerate(pos):
            ax.text(p+width , 15,
                    list(map(str,np.array(df['hemiS*3'])))[ix],
                    fontsize='small',
                    ha='center', va='center')

    ax.bar([p+width*2 for p in pos], df['allsky*3'], width, alpha=0.5,
            color=colors[2], label='allsky*3')
    if not linscale:
        for p, ix in enumerate(pos):
            ax.text(p+width*2, 15,
                    list(map(str,np.array(df['allsky*3'])))[ix],
                    fontsize='small',
                    ha='center', va='center')

    ax.set_xticks([-0.125+p+1.5*width for p in pos])

    ax.set_xticklabels(df['period_bin_days'])

    ax.set_ylabel('new detected planets (T<12)')
    ax.set_xlabel('period [days]')

    if linscale:
        pass
    else:
        ax.set_yscale('log')
        ax.set_ylim([1e1, 3.2e3])

    ax.legend(loc='best', fontsize='small')

    fig.tight_layout()
    cstr = 'darkcolor' if darkcolor else 'lightcolor'
    savname = (
        '../results/period_barchart_lineary_{:s}.png'.format(cstr)
        if linscale else
        '../results/period_barchart_logy_{:s}.png'.format(cstr)
    )
    fig.savefig(savname, dpi=350)
    print('saved {:s}'.format(savname))

if __name__ == '__main__':

    df = make_dataframe()
    make_period_barchart(df, linscale=True)
    make_period_barchart(df, linscale=False)
    make_period_barchart(df, linscale=True, darkcolor=True)
    make_period_barchart(df, linscale=False, darkcolor=True)

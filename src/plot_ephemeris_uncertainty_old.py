def plot_ephemeris_uncertainty_old(seed=42):

    N_tra = 3 # in primary...
    sigma_tc_single_tra = (4*u.minute).to(u.hr).value
    period_yr = 10/365.25 # days
    ext_transit_epoch = 73 # 73 transits later... roughly 2 years
    _epoch = np.arange(0, 400, 1)

    primary_epochs = np.arange(0,N_tra)
    epochs = np.concatenate((np.atleast_1d(primary_epochs),
                             np.atleast_1d(np.array(ext_transit_epoch)))
                           )

    primarystart = Time('2018-07-18')
    emstart = Time('2020-06-27')
    emend = Time('2022-10-08')
    jwstlaunch = Time('2021-03-30') # https://www.nasa.gov/press-release/nasa-completes-webb-telescope-review-commits-to-launch-in-early-2021
    gmteltstart = Time('2024-06-01') # roughly simultaneous (https://www.gmto.org/overview/, https://www.eso.org/sci/facilities/eelt/)
    tmtstart = Time('2027-07-20') # https://www.tmt.org/page/timeline?category=Observatory+Construction

    true_transit_times = (
        primarystart.byear + period_yr*epochs
    )

    np.random.seed(seed)
    err_tc_hr = np.random.normal(0., sigma_tc_single_tra, size=N_tra+1)

    # add errors to true times. get observed transit times (in years)
    observed_transit_times = true_transit_times + err_tc_hr/(24*365.25)

    # fit t_mid_observed vs epoch. tmid = t_0 + P*E. "P" is in years.
    (slope,
     intercept,
     r_value,
     p_value,
     std_err) = linregress(epochs[:-1], observed_transit_times[:-1])

    (slope_extra,
     intercept_extra,
     r_value_extra,
     p_value_extra,
     std_err_extra) = linregress(epochs, observed_transit_times)

    ##########################################

    f, ax = plt.subplots(figsize=(6,4))

    # epochs in units of time
    t_interp_byear = primarystart.byear + _epoch*period_yr
    t_epochs_byear = primarystart.byear + epochs*period_yr

    tmid_calc = intercept + slope*epochs

    # plot observed - calculated
    #FIXME
    # ax.errorbar(
    #     t_epochs_byear, (observed_transit_times - tmid_calc)*(365.25*24),
    #     yerr=sigma_tc_single_tra, fmt='.', c='k', markersize=5
    # )

    # tmid = t_0 + period*epoch
    tmid_interp = intercept + slope*_epoch
    tmid_interp_plus_1_sig = intercept + (slope+1.96*std_err)*_epoch
    tmid_interp_minus_1_sig = intercept + (slope-1.96*std_err)*_epoch
    tmid_interp_extra = intercept_extra + slope_extra*_epoch
    tmid_interp_plus_1_sig_extra = (
        intercept_extra + (slope_extra+1.96*std_err_extra)*_epoch
    )
    tmid_interp_minus_1_sig_extra = (
        intercept_extra + (slope_extra-1.96*std_err_extra)*_epoch
    )

    prime = (
        (tmid_interp_plus_1_sig-tmid_interp_extra)
        -
        (tmid_interp_minus_1_sig-tmid_interp_extra)
    )*(365.25*24)

    extend = (
        (tmid_interp_plus_1_sig_extra-tmid_interp_extra)
        -
        (tmid_interp_minus_1_sig_extra-tmid_interp_extra)
    )*(365.25*24)

    ax.plot(t_interp_byear[prime>0], prime[prime>0], color='#1f77b4',
            zorder=-3, alpha=1, lw=3)
    sel = (extend>0) & (t_interp_byear > emstart.byear)
    ax.plot(t_interp_byear[sel], extend[sel], color='#ff7f0e',
            zorder=-2, alpha=1, lw=3)

    # ax.hlines(0.5, 2018, 2029, linestyles='-', color='k', zorder=5, alpha=0.7,
    #           lw=1)
    # ax.hlines(-0.5, 2018, 2029, linestyles='-', color='k', zorder=5, alpha=0.7,
    #           lw=1)
    # arrowprops = dict(facecolor='black', edgecolor='black', arrowstyle='->',
    #                   linewidth=0.5, connectionstyle='arc3,rad=-0.05')
    # bbox = dict(facecolor='white',edgecolor='none',
    #             alpha=0.95,linewidth=0.5,pad=0.2)
    # ax.annotate('$\pm 30\,{\mathrm{minutes}}$', xy=(2019.5, -0.5),
    #             xycoords='data', xytext=(2019.8,-10), ha='center', va='center',
    #             arrowprops=arrowprops, bbox=bbox, zorder=4)

    # ax.text(0.96, 0.96-0.25, 'Prime only ($2\sigma$)', color='#1f77b4',
    #         ha='right', va='top', transform=ax.transAxes,
    #         bbox=dict(facecolor='white', lw=0))

    # ax.text(0.96, 0.89-0.25, 'Prime and Extended ($2\sigma$)', color='#ff7f0e',
    #         ha='right', va='top', transform=ax.transAxes,
    #         bbox=dict(facecolor='white', lw=0))

    ax.set_xlabel('Year', fontsize='large')
    ax.set_ylabel('Uncertainty in transit epoch [hours]', fontsize='large')
    ax.set_yscale('log')

    ax.xaxis.set_major_locator(MaxNLocator(integer=True))

    ylim = ax.get_ylim()
    # ax.vlines(primarystart.byear, min(ylim), max(ylim), color='gray',
    #           zorder=-3, alpha=0.8, linestyles='--', lw=1)
    # ax.vlines(emstart.byear, min(ylim), max(ylim), color='gray', zorder=-3,
    #           alpha=0.8, linestyles='--', lw=1)
    # ax.vlines(emend.byear, min(ylim), max(ylim), color='gray', zorder=-3,
    #           alpha=0.8, linestyles='--', lw=1)

    # ax.vlines(jwstlaunch.byear, min(ylim), max(ylim), color='goldenrod',
    #           zorder=-3, alpha=0.7, linestyles=':', lw=1)
    # ax.vlines(gmteltstart.byear, min(ylim), max(ylim), color='goldenrod',
    #           zorder=-3, alpha=0.7, linestyles=':', lw=1)
    # ax.vlines(tmtstart.byear, min(ylim), max(ylim), color='goldenrod',
    #           zorder=-3, alpha=0.7, linestyles=':', lw=1)

    # x is data units, y is axes-units
    trans = transforms.blended_transform_factory(
            ax.transData, ax.transAxes)

    # ax.text(
    #     primarystart.byear + (emstart.byear - primarystart.byear)/2,
    #     0.97,
    #     'TESS\nPrime',
    #     ha='center',va='top',
    #     transform=trans, bbox=dict(facecolor='white', lw=0)
    # )
    # ax.text(
    #     emstart.byear + (emend.byear - emstart.byear)/2,
    #     0.97,
    #     'TESS\nExtended',
    #     ha='center',va='top',
    #     transform=trans, bbox=dict(facecolor='white', lw=0)
    # )

    # ax.text(
    #     jwstlaunch.byear+0.05,
    #     0.3,
    #     'JWST launch',
    #     ha='left',va='center',
    #     rotation=90,
    #     fontsize='small',
    #     transform=trans
    # )
    # ax.text(
    #     gmteltstart.byear+0.05,
    #     0.23,
    #     'GMT/ELT first light',
    #     ha='left',va='center',
    #     rotation=90,
    #     fontsize='small',
    #     transform=trans
    # )
    # ax.text(
    #     tmtstart.byear+0.05,
    #     0.25,
    #     'TMT first light',
    #     ha='left',va='center',
    #     rotation=90,
    #     fontsize='small',
    #     transform=trans
    # )

    ax.get_yaxis().set_tick_params(which='both', direction='in')
    ax.get_xaxis().set_tick_params(which='both', direction='in')

    ax.set_ylim((0.08,108))
    ax.set_xlim([2018,2029])

    if np.array_equal(
        ax.get_yticks(),
        np.array([1.e-03, 1.e-02, 1.e-01, 1.e+00, 1.e+01, 1.e+02, 1.e+03,
                  1.e+04])
    ):
        ax.set_yticklabels('0.001,0.01,0.1,0,10,100,1000,10000'.split(','))

    f.tight_layout()
    outpath = '../results/ephemeris_uncertainty.pdf'
    f.savefig(outpath, bbox_inches='tight')
    print('made {}'.format(outpath))




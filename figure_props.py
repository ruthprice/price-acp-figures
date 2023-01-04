# figure_props.py
# Contains settings to control figures
# =====================================================================
import datetime as dt

dpi = 300      # dots per sq. inch
cm = 1/2.54    # converts inches to cm for setting figure sizes

# UM suite codes of the simulations plotted in fig.2
fig3_suites = [
    'u-ci572',
    # 'u-ch103',
    'u-ch625',
    'u-cm612',
    # 'u-cm614'
]

fig4_suites = [
    'u-ci572',
    'u-ch625',
    'u-cj340'
]

fig5_suites = [
    'u-ci572',
    'u-ch625',
    # 'u-ch627',
    'u-ci364'
]

fig7_suites = [
    'u-ci572',
    'u-ch625',
    'u-ci364',
    # 'u-ci384',
    'u-cp166'
]

# Simulation names for legends
suite_labels = {
    # CONTROL
    'u-ci515': 'UKESM_DEFAULT',
    'u-ci572': 'CONTROL',
    # K06_BL
    'u-cf509': 'SA_BL_DEFAULT',
    'u-ch103': 'SA_BL',
    # IA_BL(_SecOrg)
    'u-cm611': 'IA_BL_DEFAULT',
    'u-cm612': 'IA_BL',
    'u-cm614': 'IA_BL_OXID',
    # M10_ALL
    'u-ce705': 'SOA_DEFAULT',
    'u-ch621': 'SOA',
    'u-ch627': 'SOA_OXID',
    # M10_BL
    'u-ce703': 'SOA_BL_DEFAULT',
    'u-ch625': 'SOA_BL',
    # M10_FT
    'u-ch785': 'SOA_FT',
    # M10_Prsc(_SecOrg)
    'u-cg237': 'SOA_PRSC_DEFAULT',
    'u-ci364': 'SOA_PRSC',
    'u-ci384': 'SOA_PRSC_OXID',
    # IA_BL_M10
    'u-cp166': 'IA_BL_SOA',
    # M10_BL_70N, 80N, 85N
    'u-cj701': 'SOA_BL_75N',
    'u-cj702': 'SOA_BL_80N',
    'u-cj340': 'SOA_BL_85N',
}
# Line colours of simulations
colours = {
    # CONTROL
    'u-ci515': '#CC78BC',
    'u-ci572': '#CC78BC',
    # K06_BL
    'u-cf509': '#949494',
    'u-ch103': '#949494',
    # IA_BL(_SecOrg)
    'u-cm611': '#DE8F05',
    'u-cm612': '#DE8F05',
    'u-cm614': '#DE8F05',
    # M10_ALL
    'u-ce705': '#D65D00',
    'u-ch621': '#D65D00',
    'u-ch627': '#D65D00',
    # M10_BL
    'u-ce703': '#0173B2',
    'u-ch625': '#0173B2',
    # M10_FT
    'u-ch785': '#CA9161',
    # M10_Prsc(_SecOrg)
    'u-ci364': '#029E73',
    'u-ci384': '#029E73',
    'u-cg237': '#029E73',
    # IA_BL_M10
    'u-cp166': 'crimson',
    # M10_BL_70N, 80N, 85N
    'u-cj701': 'darkblue',
    'u-cj702': 'blue',
    'u-cj340': 'lightblue',
}
# Line thickness of simulations (thick or thin)
thick_line = 1.2
thin_line = 0.5
linewidths = {
    # CONTROL
    'u-ci515': thin_line,
    'u-ci572': thick_line,
    # K06_BL
    'u-cf509': thin_line,
    'u-ch103': thick_line,
    'u-cp135': thin_line,
    # IA_BL(_SecOrg)
    'u-cm611': thin_line,
    'u-cm612': thick_line,
    'u-cm614': thin_line,
    # M10_ALL
    'u-ce705': thin_line,
    'u-ch621': thick_line,
    'u-ch627': thin_line,
    # M10_BL
    'u-ce703': thin_line,
    'u-ch625': thick_line,
    # Metzger FT
    'u-ch785': thick_line,
    # M10_Prsc
    'u-ci364': thick_line,
    'u-ci384': thin_line,
    'u-cg237': thin_line,
    # IA_BL_M10
    'u-cp166': thick_line,
    # M10_BL_70N, 80N, 85N
    'u-cj701': thick_line,
    'u-cj702': thick_line,
    'u-cj340': thick_line,
}

# fontsizes
ax_fs = 5
ax_label_fs = 6
label_fs = 8
legend_fs = 5

# for marking AO2018 NPF events on time series
ao2018_npf_events = [
        [dt.datetime(2018,8,27,13),dt.datetime(2018,8,27,19)],
        [dt.datetime(2018,8,30,12),dt.datetime(2018,8,31,20)],
        [dt.datetime(2018,9,5,0),dt.datetime(2018,9,6,18)],
        [dt.datetime(2018,9,9,0),dt.datetime(2018,9,11,12)],
        [dt.datetime(2018,9,15,15),dt.datetime(2018,9,16,9)],
        [dt.datetime(2018,9,17,10),dt.datetime(2018,9,17,22)],
        [dt.datetime(2018,9,19,10),dt.datetime(2018,9,19,19)],
        ]
npf_event_y = [2.5e4,1.5e4,5e2]

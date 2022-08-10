# price-acp-figures.py 
# Main script
# RP Aug 2022
# =====================================================================

# main body
import datetime as dt
import aerosols as aero
from obs_map import plot_obs_map

import matplotlib.pyplot as plt
# =====================================================================

# ---------------------------------------------------------------------
start = dt.datetime.now()
verbose = True

# ---------------------------------------------------------------------
# PRE-PROCESSING
# ---------------------------------------------------------------------
model_output_path = '/gws/nopw/j04/asci/rprice/ukca_output/'
suite = 'u-ci572'

if verbose:
    print('\nGetting the aerosol concs from model output..')

N, bin_edges, times = aero.get_ao2018_aerosol_conc(model_output_path, suite)

# ---------------------------------------------------------------------
# PLOTTING
# ---------------------------------------------------------------------

# ---------------------------------------------------------------------
# FIGURE 1: map of ATom, ASCOS and AO2018 coordinates
if verbose:
    print('\nMaking figure 1..')
atom_path = '/home/users/eersp/atom_data/ATom_nav_1613/data/'
ascos_path = "/home/users/eersp/ascos_flights/data_files/"
fig1_filename = 'figures/fig01.pdf'
plot_obs_map(atom_path, ascos_path, fig1_filename)
if verbose:
    print('Done.')

# =====================================================================
end = dt.datetime.now()
print('Finished script at {} in {}'.format(end, end-start))
# =====================================================================
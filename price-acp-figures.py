# price-acp-figures.py 
# Main script
# RP Aug 2022
# =====================================================================
import datetime as dt
from obs_map import plot_obs_map
# =====================================================================
start = dt.datetime.now()
verbose = True

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
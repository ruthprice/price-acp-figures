# price-acp-figures.py 
# Main script
# RP Aug 2022
# =====================================================================

# main body
import datetime as dt
import figure_props as fprops
import aerosols as aero
import time_series_pdf as ts
from obs_map import plot_obs_map

# =====================================================================

# ---------------------------------------------------------------------
start = dt.datetime.now()
verbose = True

# ---------------------------------------------------------------------
# PRE-PROCESSING
# ---------------------------------------------------------------------
model_output_path = '/gws/nopw/j04/asci/rprice/ukca_output/'
N = {}
bin_edges = {}
times = {}
for suite in fprops.fig2_suites:
    if verbose:
        print('\nLoading output from {}..'.format(suite))
        N[suite], bin_edges[suite], times[suite] = aero.get_ao2018_aerosol_conc(model_output_path, suite)

if verbose:
    print('\nLoading observations..')
# Load CPC and DMPS observations
# get AB's ultrafine particle concentration
ao2018_data_file_dir = "/home/users/eersp/ao2018_observations/"
ufp_data_in = ao2018_data_file_dir + "ao2018-aerosol-ufp.csv"
if verbose:
    print('\nTaking running mean of UFP data..')
# take running mean of CPC dataset
ufp_running_means, ufp_std_devs, ufp_running_mean_times = aero.running_mean(*aero.load_ao2018_ufp(ufp_data_in), 3, 1)

# get DMPS dN/dlogD
# dmps_data_in = ao2018_data_file_dir + "ao2018-aerosol-dmps.csv"
dmps_data_in = ao2018_data_file_dir + 'DMPS_for_Ruth.csv' # with pollution flag
dmps_dNdlogD, dmps_diam_mid, dmps_times = aero.load_ao2018_dmps(dmps_data_in)
# integrate DMPS dN/dlogD
dmps_N = aero.total_number_from_dN(dmps_dNdlogD, dmps_diam_mid, [15,100,500])

if verbose:
    print('\nTaking running mean of DMPS data..')
# take running mean of DMPS datasets
dmps_N15_100_running_means, dmps_N15_100_std_devs, dmps_running_mean_times = aero.running_mean(dmps_N[:,0], dmps_times, 3, 1)
dmps_N100_500_running_means, dmps_N100_500_std_devs, dmps_running_mean_times = aero.running_mean(dmps_N[:,1], dmps_times, 3, 1)


if verbose:
    print('\nMaking PDFs..')
n_pdf_bins = 25
hist_obs, hist_model, pdf_bins_mid = aero.ao2018_melt_freeze_pdfs(n_pdf_bins, fprops.fig2_suites, N, bin_edges, times)

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

# ---------------------------------------------------------------------
# FIGURE 2: time series of aerosol concentrations during AO2018
if verbose:
    print('\nMaking figure 2..')

obs_N = [ufp_running_means, dmps_N15_100_running_means, dmps_N100_500_running_means]
obs_T = [ufp_running_mean_times, dmps_running_mean_times, dmps_running_mean_times]
obs_stdev = [ufp_std_devs, dmps_N15_100_std_devs, dmps_N100_500_std_devs]
ts.plot_time_series_pdfs(fprops.fig2_suites, N, times, 
                         obs_N, obs_T, obs_stdev, 
                         hist_obs, hist_model, pdf_bins_mid,
                         'figures/fig02.pdf')

# =====================================================================
end = dt.datetime.now()
print('Finished script at {} in {}'.format(end, end-start))
# =====================================================================
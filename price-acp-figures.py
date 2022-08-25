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
# figure 3
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

# load hio3 data
import files
import numpy as np
import datetime as dt
# import aerosols as aero
# load hio3 from model
from glob import glob
import iris
import campaign_routes as cr
import constants
import numpy.ma as ma


# =====================================================================

# ---------------------------------------------------------------------
start = dt.datetime.now()
verbose = True

# ---------------------------------------------------------------------
# PRE-PROCESSING
# ---------------------------------------------------------------------

# ---------------------------------------------------------------------
# LOAD MODEL DATA FOR FIGURE 2
model_output_path = '/gws/nopw/j04/asci/rprice/ukca_output/'
N = {}
bin_edges = {}
times = {}
for suite in fprops.fig2_suites:
    if verbose:
        print('\nLoading output from {}..'.format(suite))
        N[suite], bin_edges[suite], times[suite] = aero.get_ao2018_aerosol_conc(model_output_path, suite)

# ---------------------------------------------------------------------
# LOAD MEASUREMENTS FOR FIGURE 2
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
# LOAD HIO3 MEASUREMENTS FOR FIGURE 3
# Load observations and Baccarini model data
if verbose:
    print('\nLoading HIO3..')
hio3_path = "/home/users/eersp/ao2018_observations/"
hio3_file = "ao2018-aerosol-cims.csv"
hio3_file_contents, n_rows, n_cols = files.get_csv_contents(hio3_path + hio3_file)
obs_hio3_n_time = n_rows - 1    # -1 bc there is a header
fmt = '%Y-%m-%d %H:%M:%S'
obs_hio3_times = np.array([dt.datetime.strptime(str(T),fmt) for T in hio3_file_contents[1:,0]])
obs_hio3 = hio3_file_contents[1:,2].astype(float)

obs_hio3_means, obs_hio3_stdev, obs_hio3_mean_times = aero.running_mean(obs_hio3, obs_hio3_times, 3, 1)

# ---------------------------------------------------------------------
# LOAD MODEL OUTPUT FOR FIGURE 3
model_output_path = '/gws/nopw/j04/asci/rprice/ukca_output/u-cm612/All_time_steps/pk_files/'
hio3_stashcode = "m01s34i064"
hio3_file = glob(model_output_path+'*'+hio3_stashcode+'*')
print(hio3_file)
model_hio3 = iris.load(hio3_file)[0]
model_hio3.data
print(model_hio3)
air_density = 1.33 # kg m-3
model_hio3 = model_hio3 * air_density
colocated_hio3 = cr.colocate_with_ao2018_drift(model_hio3, constants.model_res)
colocated_hio3_number = iris.cube.CubeList([])
for t,cube in enumerate(colocated_hio3):
    cube.long_name = "mass_concentration_of_hio3"
    number_conc = cube * constants.avo / constants.mm_hio3
    number_conc.units = "m-3"
    number_conc.convert_units('cm-3')
    colocated_hio3_number.append(number_conc.data)
model_times = aero.get_cube_times(model_hio3, ao_drift=True)
colocated_hio3_number = np.array(colocated_hio3_number)

if verbose:
    print('\nMaking HIO3 PDFs..')
n_hio3_pdf_bins = 20

freeze_up = dt.datetime(2018,8,27)
hio3_melt_times, hio3_freeze_times = cr.ao2018_melt_freeze_times(freeze_up, obs_hio3_mean_times)
model_melt_times, model_freeze_times = cr.ao2018_melt_freeze_times(freeze_up, model_times)

# define bins
max_hio3 = ma.amax(obs_hio3_means)
hio3_pdf_bins = np.linspace(0, max_hio3, n_hio3_pdf_bins+1)
hio3_pdf_bins_mid = 0.5*(hio3_pdf_bins[1:] + hio3_pdf_bins[:-1])

hio3_hist_obs = np.zeros((2, n_hio3_pdf_bins))
hio3_hist_obs[0] = np.histogram(obs_hio3_means[hio3_melt_times], density=True, bins=hio3_pdf_bins)[0]
hio3_hist_obs[1] = np.histogram(obs_hio3_means[hio3_freeze_times], density=True, bins=hio3_pdf_bins)[0]
hio3_hist_model = np.zeros((2, n_hio3_pdf_bins))
hio3_hist_model[0] = np.histogram(colocated_hio3_number[model_melt_times], density=True, bins=hio3_pdf_bins)[0]
hio3_hist_model[1] = np.histogram(colocated_hio3_number[model_freeze_times], density=True, bins=hio3_pdf_bins)[0]

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

# ---------------------------------------------------------------------
# FIGURE 3: time series and PDF of surface HIO3 concentration during AO2018
# plot time series and PDF on same figure with subplots
if verbose:
    print('\nMaking figure 3..')
npf_event_marker_y = 8e6

fig = plt.figure(figsize=(16*fprops.cm,3*fprops.cm), dpi=300)
widths = [4,1]
spec = fig.add_gridspec(ncols=2, nrows=1,
                        width_ratios=widths,
                        hspace=0.4, wspace=0.1)
ax1 = fig.add_subplot(spec[0])
# obs
ax1.plot(obs_hio3_mean_times, obs_hio3_means, color='black', label="Observations (3hr mean)", linewidth=fprops.thick_line)
# model
ax1.plot(model_times, colocated_hio3_number, label=fprops.suite_labels['u-cm612'],
         color=fprops.colours['u-cm612'], linewidth=fprops.linewidths['u-cm612'])
# mark NPF events
for j,event in enumerate(fprops.ao2018_npf_events):
            ax1.hlines(npf_event_marker_y,event[0],event[1],
                       linewidth=0.8, color='red', linestyle=(0,(1.8,0.9)))
            ax1.plot(event[0], npf_event_marker_y, marker='|', ms=2.5, color='red')
            ax1.plot(event[1], npf_event_marker_y, marker='|', ms=2.5,color='red')
ax1.set_xlabel('DoY', fontsize=fprops.ax_fs)
ax1.set_xlim(right=dt.datetime(2018,9,20))
ax1.set_ylabel('Surface IA conc [cm$^{-3}$]', fontsize=fprops.ax_fs)
ax1.tick_params(axis='y', which='major', labelsize=fprops.ax_fs)
ax1.tick_params(axis='y', which='minor', labelsize=fprops.ax_fs)
ax1.tick_params(axis='y', which='major', pad=1)

ax1.legend(fontsize=fprops.legend_fs)
ax1.ticklabel_format(style='sci', axis='y', scilimits=(0,0))

# use DOY on x axis
left_lim_date = dt.datetime(2018, 8, 1)
right_lim_date = dt.datetime(2018, 9, 19)
ax1.set_xlim(left=left_lim_date, right=right_lim_date)
left,right = plt.xlim()
doy_left = left_lim_date.timetuple().tm_yday
doy_right = right_lim_date.timetuple().tm_yday
x_axis_array = np.arange(left,right+1)
doy_array = np.arange(doy_left, doy_right+1)
ax1.set_xticks(x_axis_array[::3])
ax1.set_xticklabels(doy_array[::3])
ax1.tick_params(axis='x', which='major', labelsize=fprops.ax_fs)
ax1.tick_params(axis='x', which='minor', labelsize=fprops.ax_fs)

ax1.set_title('(a)', loc='left', fontsize=fprops.label_fs, y=1.1)

ax2 = fig.add_subplot(spec[1])
ax2.plot(hio3_pdf_bins_mid, hio3_hist_obs[1], color='k',
         label='Observations', linewidth=fprops.thick_line)
ax2.plot(hio3_pdf_bins_mid, hio3_hist_model[1], label=fprops.suite_labels['u-cm612'],
         color=fprops.colours['u-cm612'], linewidth=fprops.linewidths['u-cm612'])
ax2.set_xlabel("Surface IA conc [cm$^{-3}$]", labelpad=12, fontsize=fprops.ax_fs)
ax2.tick_params(axis='both', which='major', labelsize=fprops.ax_fs)
ax2.tick_params(axis='both', which='minor', labelsize=fprops.ax_fs)
ax2.tick_params(axis='y', which='major', pad=1)
ax2.set_title('Freeze season', fontsize=fprops.label_fs, y=1.1)
ax2.set_title('(b)', loc='left', fontsize=fprops.label_fs, y=1.1)

ax1.grid()
ax2.grid()
filename = "figures/fig03"
plt.savefig(filename+".pdf", bbox_inches='tight', facecolor='white', format='pdf')
# =====================================================================
end = dt.datetime.now()
print('Finished script at {} in {}'.format(end, end-start))
# =====================================================================
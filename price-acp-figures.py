# price-acp-figures.py 
# Main script
# RP Aug 2022
# =====================================================================

# ascos loading
# import aerosols as aero

# main body
import datetime as dt
import figure_props as fprops
import aerosols as aero
import time_series_pdf as ts
from obs_map import plot_obs_map
# figure 4
import matplotlib.pyplot as plt

# load hio3 data
import files
import numpy as np
import datetime as dt
# import aerosols as aero
# load hio3 from model
from glob import glob
import iris
import iris.coord_categorisation
import campaign_routes as cr
import constants
import numpy.ma as ma


# =====================================================================
def return_gridpoint_lon(X, model_resolution_x):
    # for a given longitude X, return the longitude bounds
    # of the gridbox the point is in
    # given a resolution of model_resolution_x
    if X < 0:
        X += 360
    grid_lon_bounds = np.linspace(0,360,(360/model_resolution_x)+1)
    
    for i,line_of_lon in enumerate(grid_lon_bounds[:-1]):
        if X >= line_of_lon and X < grid_lon_bounds[i+1]:
            return line_of_lon, grid_lon_bounds[i+1]

def return_gridpoint_lat(Y, model_resolution_y):
    # as above for lats
    if Y < 0:
        Y += 360
    grid_lat_bounds = np.linspace(-90,90,(180/model_resolution_y)+1)
    
    for i,line_of_lat in enumerate(grid_lat_bounds[:-1]):
        if Y >= line_of_lat and Y < grid_lat_bounds[i+1]:
            return line_of_lat, grid_lat_bounds[i+1]
# ---------------------------------------------------------------------
start = dt.datetime.now()
verbose = True

# ---------------------------------------------------------------------
# PRE-PROCESSING
# ---------------------------------------------------------------------

# # ---------------------------------------------------------------------
# # LOAD MODEL DATA FOR FIGURE 2
# model_output_path = '/gws/nopw/j04/asci/rprice/ukca_output/'
# N = {}
# bin_edges = {}
# times = {}
# for suite in fprops.fig2_suites:
#     if verbose:
#         print('\nLoading output from {}..'.format(suite))
#     N[suite], bin_edges[suite], times[suite] = aero.get_ao2018_aerosol_conc(model_output_path, suite)

# # ---------------------------------------------------------------------
# # LOAD MEASUREMENTS FOR FIGURE 2
# if verbose:
#     print('\nLoading observations..')
# # Load CPC and DMPS observations
# # get AB's ultrafine particle concentration
# ao2018_data_file_dir = "/home/users/eersp/ao2018_observations/"
# ufp_data_in = ao2018_data_file_dir + "ao2018-aerosol-ufp.csv"
# if verbose:
#     print('\nTaking running mean of UFP data..')
# # take running mean of CPC dataset
# ufp_running_means, ufp_std_devs, ufp_running_mean_times = aero.running_mean(*aero.load_ao2018_ufp(ufp_data_in), 3, 1)

# # get DMPS dN/dlogD
# # dmps_data_in = ao2018_data_file_dir + "ao2018-aerosol-dmps.csv"
# dmps_data_in = ao2018_data_file_dir + 'DMPS_for_Ruth.csv' # with pollution flag
# dmps_dNdlogD, dmps_diam_mid, dmps_times = aero.load_ao2018_dmps(dmps_data_in)
# # integrate DMPS dN/dlogD
# dmps_N = aero.total_number_from_dN(dmps_dNdlogD, dmps_diam_mid, [15,100,500])

# if verbose:
#     print('\nTaking running mean of DMPS data..')
# # take running mean of DMPS datasets
# dmps_N15_100_running_means, dmps_N15_100_std_devs, dmps_running_mean_times = aero.running_mean(dmps_N[:,0], dmps_times, 3, 1)
# dmps_N100_500_running_means, dmps_N100_500_std_devs, dmps_running_mean_times = aero.running_mean(dmps_N[:,1], dmps_times, 3, 1)


# if verbose:
#     print('\nMaking PDFs..')
# n_pdf_bins = 25
# hist_obs, hist_model, pdf_bins_mid = aero.ao2018_melt_freeze_pdfs(n_pdf_bins, fprops.fig2_suites, N, bin_edges, times)

# # ---------------------------------------------------------------------
# # LOAD HIO3 MEASUREMENTS FOR FIGURE 3
# # Load observations and Baccarini model data
# if verbose:
#     print('\nLoading HIO3..')
# hio3_path = "/home/users/eersp/ao2018_observations/"
# hio3_file = "ao2018-aerosol-cims.csv"
# hio3_file_contents, n_rows, n_cols = files.get_csv_contents(hio3_path + hio3_file)
# obs_hio3_n_time = n_rows - 1    # -1 bc there is a header
# fmt = '%Y-%m-%d %H:%M:%S'
# obs_hio3_times = np.array([dt.datetime.strptime(str(T),fmt) for T in hio3_file_contents[1:,0]])
# obs_hio3 = hio3_file_contents[1:,2].astype(float)

# obs_hio3_means, obs_hio3_stdev, obs_hio3_mean_times = aero.running_mean(obs_hio3, obs_hio3_times, 3, 1)

# # ---------------------------------------------------------------------
# # LOAD MODEL OUTPUT FOR FIGURE 3
# model_output_path = '/gws/nopw/j04/asci/rprice/ukca_output/u-cm612/All_time_steps/pk_files/'
# hio3_stashcode = "m01s34i064"
# hio3_file = glob(model_output_path+'*'+hio3_stashcode+'*')
# if verbose:
#     print(hio3_file)
# model_hio3 = iris.load(hio3_file)[0]
# model_hio3.data
# if verbose:
#     print(model_hio3)
# air_density = 1.33 # kg m-3
# model_hio3 = model_hio3 * air_density
# colocated_hio3 = cr.colocate_with_ao2018_drift(model_hio3, constants.model_res)
# colocated_hio3_number = iris.cube.CubeList([])
# for t,cube in enumerate(colocated_hio3):
#     cube.long_name = "mass_concentration_of_hio3"
#     number_conc = cube * constants.avo / constants.mm_hio3
#     number_conc.units = "m-3"
#     number_conc.convert_units('cm-3')
#     colocated_hio3_number.append(number_conc.data)
# model_times = aero.get_cube_times(model_hio3, ao_drift=True)
# colocated_hio3_number = np.array(colocated_hio3_number)

# if verbose:
#     print('\nMaking HIO3 PDFs..')
# n_hio3_pdf_bins = 20

# freeze_up = dt.datetime(2018,8,27)
# hio3_melt_times, hio3_freeze_times = cr.ao2018_melt_freeze_times(freeze_up, obs_hio3_mean_times)
# model_melt_times, model_freeze_times = cr.ao2018_melt_freeze_times(freeze_up, model_times)

# # define bins
# max_hio3 = ma.amax(obs_hio3_means)
# hio3_pdf_bins = np.linspace(0, max_hio3, n_hio3_pdf_bins+1)
# hio3_pdf_bins_mid = 0.5*(hio3_pdf_bins[1:] + hio3_pdf_bins[:-1])

# hio3_hist_obs = np.zeros((2, n_hio3_pdf_bins))
# hio3_hist_obs[0] = np.histogram(obs_hio3_means[hio3_melt_times], density=True, bins=hio3_pdf_bins)[0]
# hio3_hist_obs[1] = np.histogram(obs_hio3_means[hio3_freeze_times], density=True, bins=hio3_pdf_bins)[0]
# hio3_hist_model = np.zeros((2, n_hio3_pdf_bins))
# hio3_hist_model[0] = np.histogram(colocated_hio3_number[model_melt_times], density=True, bins=hio3_pdf_bins)[0]
# hio3_hist_model[1] = np.histogram(colocated_hio3_number[model_freeze_times], density=True, bins=hio3_pdf_bins)[0]

# # ---------------------------------------------------------------------
# # LOAD MODEL OUTPUT FOR FIG 4
# if verbose:
#     print('\nLoading output for fig 4..')
# model_output_path = '/gws/nopw/j04/asci/rprice/ukca_output/'
# nuc_stash = 'm01s34i101'
# sol_nuc_N = {}
# zbl = {}
# for s,suite in enumerate(fprops.fig4_suites):
#     density_file = "{}{}/L1/daily_3d/L1_air_density_Density_of_air.nc".format(model_output_path,suite)
#     air_density = iris.load(density_file)[0]
#     particle_density_of_air = (air_density/constants.mm_air)*constants.avo
#     try:
#         iris.coord_categorisation.add_day_of_year(particle_density_of_air, 'time', name='day_of_year')
#         iris.coord_categorisation.add_year(particle_density_of_air, 'time', name='year')
#     except ValueError:
#         # ValueError raised if coord already exists
#         # don't know if anything else would trigger it so be careful
#         pass
#     nuc_file = glob('{}{}/All_time_steps/pl_files/*{}*'.format(model_output_path,suite,nuc_stash))
#     nuc_per_air_mol = iris.load(nuc_file)[0]
#     nuc_per_air_mol.data
#     if verbose:
#         print(nuc_per_air_mol)
#     n_conc = nuc_per_air_mol*particle_density_of_air
#     n_conc.long_name = "number_concentration_of_soluble_nucleation_mode_aerosol"
#     n_conc.units = "m-3"
#     n_conc.convert_units('cm-3')
#     n_conc_coloc = cr.colocate_with_ao2018_drift(n_conc, constants.model_res)
#     sol_nuc_N[suite] = np.array([cube.data for cube in n_conc_coloc])
#     fig4_times = aero.get_cube_times(n_conc, ao_drift=True)
#     fig4_heights = n_conc.coord('level_height').points/1000

#     zbl_stash = 'm01s00i025'
#     zbl_file = glob('{}{}/All_time_steps/pl_files/*{}*'.format(model_output_path,suite,zbl_stash))
#     z = iris.load(zbl_file)[0]
#     z.data
#     if verbose:
#         print(z)
#     zbl_coloc = cr.colocate_with_ao2018_drift(z, constants.model_res)
#     zbl[suite] = np.array([cube.data for cube in zbl_coloc])

# ---------------------------------------------------------------------
# LOAD DATA FOR FIGURE 5
if verbose:
    print('\nLoading data for figure 5..')
# load observations
ascos_obs_dir = "/home/users/eersp/ascos_flights/"
ascos_flight_data = aero.load_ascos_data(ascos_obs_dir)
zmax = np.amax([ascos_flight_data[flight]['z_bins'][-1] for flight in ascos_flight_data])
bin_width = 200
top_bin = bin_width * int(np.ceil(zmax/bin_width))
z_bins = np.arange(0,top_bin+bin_width,bin_width)
z_bin_centres = z_bins[1:] - (bin_width/2)
n_z_bins = len(z_bins)
ascos_N3_14_mean = ma.zeros((n_z_bins))
ascos_N14_300_mean = ma.zeros((n_z_bins))
ascos_N3_14_median = ma.zeros((n_z_bins))
ascos_N14_300_median = ma.zeros((n_z_bins))

for z, bin_lower_bound in enumerate(z_bins[:-1]):
    ascos_N3_14_binned = []
    ascos_N14_300_binned = []
    for flight in ascos_flight_data:
        flight_z = ma.masked_invalid(ascos_flight_data[flight]['altitude'])
        N3_14 = ma.masked_invalid(ascos_flight_data[flight]['N3_14'])
        N14_300 = ma.masked_invalid(ascos_flight_data[flight]['N14_300'])
        ascos_N3_14_binned.extend([x for i,x in enumerate(N3_14) if flight_z[i] >= bin_lower_bound and flight_z[i] < z_bins[z+1]])
        ascos_N14_300_binned.extend([x for i,x in enumerate(N14_300) if flight_z[i] >= bin_lower_bound and flight_z[i] < z_bins[z+1]])
    ascos_N3_14_mean[z] = ma.mean(ma.array(ascos_N3_14_binned))
    ascos_N14_300_mean[z] = ma.mean(ma.array(ascos_N14_300_binned))
    ascos_N3_14_median[z] = ma.median(ma.array(ascos_N3_14_binned))
    ascos_N14_300_median[z] = ma.median(ma.array(ascos_N14_300_binned))

# ---------------------------------------------------------------------
# LOAD MODEL OUTPUT FOR FIG 5
if verbose:
    print('\nLoading output for fig 5..')
model_output_path = '/gws/nopw/j04/asci/rprice/ukca_output/'

model_N = {}
for s,suite in enumerate(fprops.fig5_suites):
    if verbose:
        print(suite)
    model_N[suite], heights = aero.calculate_ascos_aerosol_conc(model_output_path, suite, ascos_flight_data)
# ---------------------------------------------------------------------
# PLOTTING
# ---------------------------------------------------------------------

# # ---------------------------------------------------------------------
# # FIGURE 1: map of ATom, ASCOS and AO2018 coordinates
# if verbose:
#     print('\nMaking figure 1..')
# atom_path = '/home/users/eersp/atom_data/ATom_nav_1613/data/'
# ascos_path = "/home/users/eersp/ascos_flights/data_files/"
# fig1_filename = 'figures/fig01.pdf'
# plot_obs_map(atom_path, ascos_path, fig1_filename)
# if verbose:
#     print('Done.')

# # ---------------------------------------------------------------------
# # FIGURE 2: time series of aerosol concentrations during AO2018
# if verbose:
#     print('\nMaking figure 2..')

# obs_N = [ufp_running_means, dmps_N15_100_running_means, dmps_N100_500_running_means]
# obs_T = [ufp_running_mean_times, dmps_running_mean_times, dmps_running_mean_times]
# obs_stdev = [ufp_std_devs, dmps_N15_100_std_devs, dmps_N100_500_std_devs]
# ts.plot_time_series_pdfs(fprops.fig2_suites, N, times, 
#                          obs_N, obs_T, obs_stdev, 
#                          hist_obs, hist_model, pdf_bins_mid,
#                          'figures/fig02.pdf')
# if verbose:
#     print('Done.')
# # ---------------------------------------------------------------------
# # FIGURE 3: time series and PDF of surface HIO3 concentration during AO2018
# # plot time series and PDF on same figure with subplots
# if verbose:
#     print('\nMaking figure 3..')

# ts.plot_IA_time_series_pdf(obs_hio3_means, obs_hio3_mean_times,
#                            colocated_hio3_number, model_times,
#                            hio3_hist_obs[1], hio3_pdf_bins_mid,
#                            hio3_hist_model[1], hio3_pdf_bins_mid,
#                            "figures/fig03")
# if verbose:
#     print('Done.')

# # ---------------------------------------------------------------------
# # FIGURE 4: colocated aerosol height profiles for nucleation mode
# if verbose:
#     print('\nMaking figure 4..')

# ts.time_series_3d(sol_nuc_N, zbl, fig4_times, fig4_heights, 'figures/fig04')

# if verbose:
#     print('Done.')

# ---------------------------------------------------------------------
# FIGURE 5: ASCOS profiles
fig = plt.figure(figsize=(16*fprops.cm,6*fprops.cm), dpi=300)
gs = fig.add_gridspec(ncols=2, nrows=1, wspace=0.1)

# 3-14 nm
ax1 = fig.add_subplot(gs[0])
for f,flight in enumerate(ascos_flight_data):
    ax1.plot(ascos_flight_data[flight]['N3_14_mean'],
             ascos_flight_data[flight]['z_bins'][:-1]/1000,
             zorder=1, alpha=0.4, color='grey', linewidth=fprops.thick_line)
ax1.plot(ascos_N3_14_mean[:-1], z_bins[:-1]/1000, color='k',
         label='ASCOS mean', zorder=2, linewidth=fprops.thick_line)
ax1.plot(ascos_N3_14_median[:-1], z_bins[:-1]/1000, color='k',
         linestyle='dashed', label='ASCOS median', zorder=2, linewidth=fprops.thick_line)
for s,suite in enumerate(fprops.fig5_suites):
    N = model_N[suite][0]
    colour = fprops.colours[suite]
    linewidth = fprops.linewidths[suite]
    ax1.plot(N, heights, label=fprops.suite_labels[suite],
             zorder=4, color=colour, linewidth=linewidth)

ax1.grid()
ax1.set_xscale('log')
ax1.set_ylim(0,5)
ax1.set_xlim(1e-2,1e4)
ax1.set_title('3\u201314 nm', fontsize=fprops.label_fs)
ax1.set_title('(a)', loc='left', fontsize=fprops.label_fs)
ax1.set_ylabel('Altitude [km]', fontsize=fprops.label_fs)
ax1.set_xlabel('Particle concentration [cm$^{-3}$]', fontsize=fprops.ax_label_fs)
ax1.tick_params(axis='x', labelsize=fprops.ax_fs)
ax1.tick_params(axis='y', labelsize=fprops.ax_fs)

# 14-300 nm
ax2 = fig.add_subplot(gs[1])
for f,flight in enumerate(ascos_flight_data):
    if f == 0 :
        ax2.plot(ascos_flight_data[flight]['N14_300_mean'],ascos_flight_data[flight]['z_bins'][:-1]/1000,
                 zorder=1, alpha=0.4, color='grey', linewidth=fprops.thick_line,
                 label='Flight means')
    else:
        ax2.plot(ascos_flight_data[flight]['N14_300_mean'],ascos_flight_data[flight]['z_bins'][:-1]/1000,
                 zorder=1, alpha=0.4, color='grey', linewidth=fprops.thick_line)
ax2.plot(ascos_N14_300_mean[:-1], z_bins[:-1]/1000, color='k',
         label='ASCOS mean', zorder=2, linewidth=fprops.thick_line)
ax2.plot(ascos_N14_300_median[:-1], z_bins[:-1]/1000, color='k',
         linestyle='dashed', label='ASCOS median', zorder=2, linewidth=fprops.thick_line)
for s,suite in enumerate(fprops.fig5_suites):
    N = model_N[suite][1]
    colour = fprops.colours[suite]
    linewidth = fprops.linewidths[suite]
    ax2.plot(N, heights, label=fprops.suite_labels[suite],
             zorder=4, color=colour,
             linewidth=linewidth)
ax2.grid()
ax2.set_xscale('log')
ax2.set_ylim(0,5)
ax2.set_xlim(1e-2,1e4)
ax2.set_title('14\u2013300 nm', fontsize=fprops.label_fs)
ax2.set_title('(b)', loc='left', fontsize=fprops.label_fs)
ax2.set_xlabel('Particle concentration [cm$^{-3}$]', fontsize=fprops.ax_label_fs)
ax2.tick_params(axis='x', labelsize=fprops.ax_fs)
ax2.tick_params(axis='y', labelsize=fprops.ax_fs)

plt.legend(fontsize=fprops.legend_fs, ncol=2,
           columnspacing=1, borderpad=0.2,
           handletextpad=0.3, labelspacing=0.5,
           handlelength=1.5, loc='upper center')
filename = 'figures/fig05.pdf'
plt.savefig(filename, bbox_inches="tight", facecolor='white', format='pdf')
if verbose:
    print('Done')

# =====================================================================
end = dt.datetime.now()
print('Finished script at {} in {}'.format(end, end-start))
# =====================================================================
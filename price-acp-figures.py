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
# figure 4
import matplotlib.cm as cm
from matplotlib.colors import BoundaryNorm
import string
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

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

---------------------------------------------------------------------
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
# LOAD MODEL OUTPUT FOR FIG 4
if verbose:
    print('\nLoading output for fig 4..')
model_output_path = '/gws/nopw/j04/asci/rprice/ukca_output/'
nuc_stash = 'm01s34i101'
sol_nuc_N = {}
zbl = {}
for s,suite in enumerate(fprops.fig4_suites):
    density_file = "{}{}/L1/daily_3d/L1_air_density_Density_of_air.nc".format(model_output_path,suite)
    air_density = iris.load(density_file)[0]
    particle_density_of_air = (air_density/constants.mm_air)*constants.avo
    try:
        iris.coord_categorisation.add_day_of_year(particle_density_of_air, 'time', name='day_of_year')
        iris.coord_categorisation.add_year(particle_density_of_air, 'time', name='year')
    except ValueError:
        # ValueError raised if coord already exists
        # don't know if anything else would trigger it so be careful
        pass
    nuc_file = glob('{}{}/All_time_steps/pl_files/*{}*'.format(model_output_path,suite,nuc_stash))
    nuc_per_air_mol = iris.load(nuc_file)[0]
    nuc_per_air_mol.data
    print(nuc_per_air_mol)
    n_conc = nuc_per_air_mol*particle_density_of_air
    n_conc.long_name = "number_concentration_of_soluble_nucleation_mode_aerosol"
    n_conc.units = "m-3"
    n_conc.convert_units('cm-3')
    n_conc_coloc = cr.colocate_with_ao2018_drift(n_conc, constants.model_res)
    sol_nuc_N[suite] = np.array([cube.data for cube in n_conc_coloc])
    fig4_times = aero.get_cube_times(n_conc, ao_drift=True)
    fig4_heights = n_conc.coord('level_height').points/1000

    zbl_stash = 'm01s00i025'
    zbl_file = glob('{}{}/All_time_steps/pl_files/*{}*'.format(model_output_path,suite,zbl_stash))
    z = iris.load(zbl_file)[0]
    z.data
    print(z)
    zbl_coloc = cr.colocate_with_ao2018_drift(z, constants.model_res)
    zbl[suite] = np.array([cube.data for cube in zbl_coloc])


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
if verbose:
    print('Done.')
# ---------------------------------------------------------------------
# FIGURE 3: time series and PDF of surface HIO3 concentration during AO2018
# plot time series and PDF on same figure with subplots
if verbose:
    print('\nMaking figure 3..')

ts.plot_IA_time_series_pdf(obs_hio3_means, obs_hio3_mean_times,
                           colocated_hio3_number, model_times,
                           hio3_hist_obs[1], hio3_pdf_bins_mid,
                           hio3_hist_model[1], hio3_pdf_bins_mid,
                           "figures/fig03")
if verbose:
    print('Done.')

# ---------------------------------------------------------------------
# FIGURE 4: colocated aerosol height profiles for nucleation mode
if verbose:
    print('\nMaking figure 4..')

cmap = cm.get_cmap('magma')
levels = [1e-1, 1e0, 1e1, 1e2, 1e3, 1e4]
norm = BoundaryNorm(levels,ncolors=cmap.N,extend='both')
# for m,mode in enumerate(modes):
#     print(mode)
fig = plt.figure(figsize=(16*fprops.cm,12*fprops.cm),dpi=300)
spec = fig.add_gridspec(ncols=1,nrows=3,hspace=0.3)
ax = []
subfig_labels = string.ascii_lowercase
for s,suite in enumerate(fprops.fig4_suites):
    ax1 = fig.add_subplot(spec[s])
    N = sol_nuc_N[suite]
    C = ax1.pcolormesh(fig4_times, fig4_heights, N.transpose(),
                       norm=norm, cmap=cmap, shading='nearest',
                       rasterized=True)
    plt.plot(fig4_times, zbl[suite]/1000,
             color='white', label='BL height', linewidth=0.75)
    if s == 0:
        ax1.legend(fontsize=fprops.legend_fs)
    ax1.set_title(fprops.suite_labels[suite], fontsize=fprops.label_fs)
    ax1.set_title('({})'.format(subfig_labels[s]), loc='left', fontsize=fprops.label_fs)
    # make pretty axes
    left_lim_date = dt.datetime(2018, 8, 2)
    right_lim_date = dt.datetime(2018, 9, 19)
    ax1.set_xlim(left=left_lim_date, right=right_lim_date)
    left,right = ax1.get_xlim()
    doy_left = left_lim_date.timetuple().tm_yday
    doy_right = right_lim_date.timetuple().tm_yday
    x_axis_array = np.arange(left,right+1)
    doy_array = np.arange(doy_left, doy_right+1)
    ax1.set_xticks(x_axis_array[::2])
    if s == len(fprops.fig4_suites)-1:
        ax1.set_xticklabels(doy_array[::2])
    else:
        ax1.set_xticklabels([])
    ax1.tick_params(axis='x', labelsize=fprops.ax_fs)
    ax1.tick_params(axis='y', labelsize=fprops.ax_fs)
    ax.append(ax1)
fig.supxlabel('Day of year', fontsize=fprops.ax_label_fs,
          y=-0.25, x=0.5, ha='center', transform=ax[-1].transAxes)
fig.supylabel('Altitude [km]', fontsize=fprops.ax_label_fs, x=0.075)
# colorbar
fmt = lambda x, pos: '10$^{{{:.0f}}}$'.format(np.log10(x))
cbar = plt.colorbar(C, extend='both', label="N [cm$^{-3}$]",
                    ax=ax, format=FuncFormatter(fmt), shrink=0.6)
cbar.ax.tick_params(labelsize=fprops.ax_fs)
cbar.ax.set_ylabel(ylabel="Nucleation particle concentration [cm$^{-3}$]", fontsize=fprops.ax_label_fs)
filename = 'figures/fig04'
plt.savefig(filename+".pdf", bbox_inches='tight', facecolor='white', format='pdf')


if verbose:
    print('Done.')
# =====================================================================
end = dt.datetime.now()
print('Finished script at {} in {}'.format(end, end-start))
# =====================================================================
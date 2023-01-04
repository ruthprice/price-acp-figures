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

# ---------------------------------------------------------------------
# LOAD MODEL DATA FOR FIGURE 3
model_output_path = '/gws/nopw/j04/asci/rprice/ukca_output/'
N = {}
bin_edges = {}
times = {}
for suite in fprops.fig3_suites:
    if verbose:
        print('\nLoading output from {}..'.format(suite))
    N[suite], bin_edges[suite], times[suite] = aero.get_ao2018_aerosol_conc(model_output_path, suite)

# ---------------------------------------------------------------------
# LOAD MEASUREMENTS FOR FIGURE 3 and 7
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
if verbose:
    print('\nloading DMPS for time series..')
dmps_N = aero.total_number_from_dN(dmps_dNdlogD, dmps_diam_mid, [15,100,500], base='10')

if verbose:
    print('\nTaking running mean of DMPS data..')
# take running mean of DMPS datasets
dmps_N15_100_running_means, dmps_N15_100_std_devs, dmps_running_mean_times = aero.running_mean(dmps_N[:,0], dmps_times, 3, 1)
dmps_N100_500_running_means, dmps_N100_500_std_devs, dmps_running_mean_times = aero.running_mean(dmps_N[:,1], dmps_times, 3, 1)


if verbose:
    print('\nMaking PDFs..')
n_pdf_bins = 25
hist_obs, hist_model, pdf_bins = aero.ao2018_melt_freeze_pdfs(n_pdf_bins, fprops.fig2_suites, N, bin_edges, times)
pdf_bins_mid = [0.5*(X[1:] + X[:-1]) for X in pdf_bins]

# ---------------------------------------------------------------------
# LOAD HIO3 MEASUREMENTS FOR FIGURE 2
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
# LOAD MODEL OUTPUT FOR FIGURE 2
model_output_path = '/gws/nopw/j04/asci/rprice/ukca_output/u-cm612/All_time_steps/pk_files/'
hio3_stashcode = "m01s34i064"
hio3_file = glob(model_output_path+'*'+hio3_stashcode+'*')
if verbose:
    print(hio3_file)
model_hio3 = iris.load(hio3_file)[0]
model_hio3.data
if verbose:
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
    if verbose:
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
    if verbose:
        print(z)
    zbl_coloc = cr.colocate_with_ao2018_drift(z, constants.model_res)
    zbl[suite] = np.array([cube.data for cube in zbl_coloc])

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
# LOAD DATA FOR FIGURE 6
if verbose:
    print('\nLoading data for fig6..')
atom_data_path = '/home/users/eersp/atom_data/ATom_AMP_Instrument_Data_1671/data/'
atom_nav_path = '/home/users/eersp/atom_data/ATom_nav_1613/data/'
# load aerosol data
atom1_sd_files = 'SDAerosol_DC8_2016*.ict'
atom_aero_T, atom_dNdlogD, atom_Dp = {}, {}, {}
n = 1
files = glob(atom_data_path+atom1_sd_files)
atom_aero_T, atom_dNdlogD, atom_Dp = aero.read_atom_sd_aerosol_files(files)
# load navigation data
atom_nav_t1, atom_nav_t2, atom_lats, atom_lons, atom_alts = aero.read_atom_nav_file(atom_nav_path)[:5]
n_nav_time = len(atom_nav_t1)
# get common timesteps for nav and measurement data
# nav data has 10s interval, measurements 1s
aero_start_times = np.array([atom_aero_T[f][0] for f in atom_aero_T])
aero_end_times = np.array([atom_aero_T[f][-1] for f in atom_aero_T])
# get lists of indices of the nav_t1/nav_t2 arrays
# that correspond to beginning and end of aerosol files
# eg the latitude at the start of the first file would be
# atom_lats[i_nav_aero_start[<first flight code>]
i_nav_aero_start = {}
i_nav_aero_end = {}
for i_nav in np.arange(n_nav_time):
    for i_aero_file,f in enumerate(atom_aero_T):
        if aero_start_times[i_aero_file] > atom_nav_t1[i_nav]:
            if aero_start_times[i_aero_file] <= atom_nav_t2[i_nav]:
                # this file starts within this nav time interval
                i_nav_aero_start[f] = i_nav
        if aero_end_times[i_aero_file] > atom_nav_t1[i_nav]:
            if aero_end_times[i_aero_file] <= atom_nav_t2[i_nav]:
                # this file ends within this nav time interval
                i_nav_aero_end[f] = i_nav
# get lists of indices of navigation arrays
# corresponding to each aerosol time array element
# eg latitude at timestep atom_aero_T[f][x] would be
# atom_lats[i_nav_for_aero[f][x]]
i_nav_for_aero = {}
for i_aero_file,f in enumerate(atom_aero_T):
    i_nav_for_aero[f] = []
    for i_aero, aero_T in enumerate(atom_aero_T[f]):
        found_timestep_flag = 0
        for i_nav in np.arange(i_nav_aero_start[f],i_nav_aero_end[f]+1):
            if aero_T > atom_nav_t1[i_nav] and aero_T <= atom_nav_t2[i_nav]:
                found_timestep_flag = 1
                i_nav_for_aero[f].append(i_nav)
        if found_timestep_flag == 0:
            print("ERROR: couldn't find aero timestep {} from file {} in navigation file".format(i_aero,f))
    i_nav_for_aero[f] = np.array(i_nav_for_aero[f])
# get indices of timesteps that are high latitude
# eg for a file f, timesteps where latitude is > min_lat
# are given by obs_aero_T[f][i_aero_high_lats[f]]
min_lat = 60
i_aero_high_lats = {}
for i_aero_file, f in enumerate(atom_aero_T):
    i_aero_high_lats[f] = []
    for i_aero, aero_T in enumerate(atom_aero_T[f]):
        if atom_lats[int(i_nav_for_aero[f][i_aero])] >= min_lat:
            i_aero_high_lats[f].append(i_aero)
    i_aero_high_lats[f] = np.array(i_aero_high_lats[f]).astype(int)
# Sum over size distribution to match UKCA sizes
size_bins = [2.5, 15, 100, 500] # edges of bins
n_size_bins = len(size_bins) - 1
size_limit_inds = {}
for f in atom_Dp:
    size_limit_inds[f] = []
    for n in size_bins:
        if n < atom_Dp[f][0]:
            # if n is smaller than first bin
            loc = 0
        else:
            loc = np.nonzero(atom_Dp[f] < n)[0][-1]
        # loc is index of largest diameter bin that's still less than n
        size_limit_inds[f].append(loc + 1)
atom_N = {}
atom_N_mean = {}
atom_N_median = {}
atom_N_stddev = {}
atom_z_bins = {}
for i_aero_file, f in enumerate(atom_aero_T):
    if len(i_aero_high_lats[f]) > 0:
        atom_N[f] = ma.zeros((len(i_aero_high_lats[f]),n_size_bins))
        bin_bounds = aero.calc_bin_boundaries(atom_Dp[f])
        log_bin_bounds = np.log10(bin_bounds)
        dlogD = log_bin_bounds[1:] - log_bin_bounds[:-1]
        for i_high_lat, i_aero in enumerate(i_aero_high_lats[f]):
            for D in np.arange(n_size_bins):
                i1 = size_limit_inds[f][D]
                i2 = size_limit_inds[f][D+1]
                atom_N[f][i_high_lat,D] = ma.sum(dlogD[i1:i2] * atom_dNdlogD[f][i_aero][i1:i2])
            atom_N[f][i_high_lat] = ma.masked_invalid(atom_N[f][i_high_lat])

        i_high_lat_nav = i_nav_for_aero[f][i_aero_high_lats[f]]
        z_bin_size = 250 # [m]
        z_max = np.amax(atom_alts[i_high_lat_nav])
        n_z_bins = int(np.ceil(z_max/z_bin_size))
        atom_z_bins[f] = np.linspace(0,n_z_bins*z_bin_size, num=n_z_bins)

        # for each element in atom_N or atom_alts,
        # value of i_z_bins gives index of z_bins
        # for that element
        i_z_bins = np.digitize(atom_alts[i_high_lat_nav], atom_z_bins[f])

        atom_N_mean[f] = ma.zeros((n_z_bins,n_size_bins))
        atom_N_median[f] = ma.zeros((n_z_bins,n_size_bins))
        atom_N_stddev[f] = ma.zeros((n_z_bins,n_size_bins))
        for d in np.arange(n_size_bins):
            N_binned = []
            for z in np.arange(n_z_bins):
                N_binned.append([])

            for i in np.arange(len(i_aero_high_lats[f])):
                N_binned[i_z_bins[i]].append(atom_N[f][i,d])

            z_bin_sample_size = np.array([len(x) for x in N_binned])

            for z in np.arange(n_z_bins):
                N_z = ma.masked_where(np.array(N_binned[z]) < 0, N_binned[z])
                N_z = ma.masked_invalid(N_z)
                atom_N_mean[f][z,d] = ma.mean(N_z)
                atom_N_median[f][z,d] = ma.median(N_z, axis=0)
                atom_N_stddev[f][z,d] = ma.std(N_z, axis=0)

# now bin merged data
atom_alts_all = np.array([atom_alts[i] for f in i_nav_for_aero for i in i_nav_for_aero[f][i_aero_high_lats[f]]])
atom_N_all = []
for f in atom_N:
    atom_N_all.extend(atom_N[f])
atom_N_all = ma.array(atom_N_all)
z_bin_size = 250 # [m]
z_max_merged = np.amax(atom_alts_all)
n_z_bins_merged = int(np.ceil(z_max_merged/z_bin_size))
z_bins_merged = np.linspace(0,n_z_bins_merged*z_bin_size, num=n_z_bins_merged)

# for each element in obs_N[:,d] or alts,
# value of i_z_bins gives index of z_bins
# for that element
i_z_bins_merged = np.digitize(atom_alts_all, z_bins_merged)

atom_N_all_mean = ma.zeros((n_z_bins_merged,n_size_bins))
atom_N_all_median = ma.zeros((n_z_bins_merged,n_size_bins))
atom_N_all_stddev = ma.zeros((n_z_bins_merged,n_size_bins))
for d in np.arange(n_size_bins):
    N_binned_merged = []
    for z in np.arange(n_z_bins_merged):
        N_binned_merged.append([])

    for i in np.arange(len(atom_N_all)):
        N_binned_merged[i_z_bins_merged[i]].append(atom_N_all[i,d])

    for z in np.arange(n_z_bins_merged):
        N_binned_merged[z] = ma.masked_where(np.array(N_binned_merged[z]) < 0, N_binned_merged[z])
        N_binned_merged[z] = ma.masked_invalid(N_binned_merged[z])
        atom_N_all_mean[z,d] = ma.mean(N_binned_merged[z], axis=0)
        atom_N_all_median[z,d] = ma.median(N_binned_merged[z], axis=0)
        atom_N_all_stddev[z,d] = ma.std(N_binned_merged[z], axis=0)

# ---------------------------------------------------------------------
# LOAD MODEL OUTPUT FOR FIGURE 6
if verbose:
    print('\nLoading output for fig6..')
model_output_path = '/gws/nopw/j04/asci/rprice/ukca_output/'
number_conc = {}
number_conc_STP = {}
radius = {}
STP_factor = {}
number_conc_aug = {}
radius_aug = {}
n_conc_extracted_dict = {}
diams_extracted_dict = {}
atom_model_N = {}
atom_model_z = {}
for s,suite in enumerate(fprops.fig5_suites):
    if verbose:
        print(suite)
    number_conc[suite] = aero.load_monthly_mean_aerosol_number_conc(model_output_path, suite)
    STP_factor[suite] = aero.calculate_STP_factor(model_output_path, suite)
    number_conc_STP[suite] = {}
    for mode in number_conc[suite]:
        number_conc_STP[suite][mode] = number_conc[suite][mode] * STP_factor[suite]
    radius[suite] = aero.get_mode_radius(model_output_path, suite, number_conc=number_conc[suite], time_res='monthly')
    # Extract August 2018
    number_conc_aug[suite] = {}
    radius_aug[suite] = {}
    for mode in number_conc[suite]:
        Naero = number_conc_STP[suite][mode]
        R = radius[suite][mode]
        try:
            iris.coord_categorisation.add_month(Naero, 'time', name='month')
            iris.coord_categorisation.add_year(Naero, 'time', name='year')
        except ValueError:
            pass
        Naero = Naero.extract(iris.Constraint(month='Aug'))
        Naero = Naero.extract(iris.Constraint(year=2018))
        number_conc_aug[suite][mode] = Naero
        try:
            iris.coord_categorisation.add_month(R, 'time', name='month')
            iris.coord_categorisation.add_year(R, 'time', name='year')
        except ValueError:
            pass
        R = R.extract(iris.Constraint(month='Aug'))
        R = R.extract(iris.Constraint(year=2018))
        radius_aug[suite][mode] = R

    lat_coord = number_conc[suite][constants.modes[0]].coord('latitude')
    lon_coord = number_conc[suite][constants.modes[0]].coord('longitude')
    if not lat_coord.has_bounds():
        lat_coord.guess_bounds()
    if not lon_coord.has_bounds():
        lon_coord.guess_bounds()
    model_lats = np.unique(lat_coord.bounds.flatten())
    model_lons = np.unique(lon_coord.bounds.flatten())
    model_flight_coords = []
    n = 0
    for i_aero_file, f in enumerate(atom_aero_T):
        for i_atom_coord in i_nav_for_aero[f][i_aero_high_lats[f]]:
            lat = atom_lats[i_atom_coord]
            lon = atom_lons[i_atom_coord]
            if lon < 0:
                lon += 360
            i_model_lon = [x for x,model_lon in enumerate(model_lons[:-1]) if lon >= model_lon and lon < model_lons[x+1]]
            i_model_lat = [y for y,model_lat in enumerate(model_lats[:-1]) if lat >= model_lat and lat < model_lats[y+1]]
            if len(i_model_lon) > 1 or len(i_model_lat) > 1:
                print("ERROR: location of model gridbox didn't work")
                print(lon, lat, model_lon, model_lat)
            if len(i_model_lon) == 1 and len(i_model_lat) == 1:
                model_flight_coords.append([i_model_lon[0],i_model_lat[0]])
    model_flight_coords = np.unique(model_flight_coords, axis=0)

#     # plot map of selected gridboxes with flight track
#     fig = plt.figure()#figsize=(12,6),dpi=300)
#     ax = plt.axes(projection=ccrs.Mercator())
#     ax.set_extent([-180,-30,55,90], crs=ccrs.PlateCarree())
#     ax.coastlines(linewidth=0.5)

#     for i_aero_file, f in enumerate(atom_aero_T):
#         if len(i_aero_high_lats[f]) > 0:
#             i_coords = i_nav_for_aero[f][i_aero_high_lats[f]]
#             plt.plot(atom_lons[i_coords], atom_lats[i_coords], transform=ccrs.PlateCarree(),
#                      linestyle='None', marker='.', ms=1, color='tab:orange')
#     for i,coord in enumerate(model_flight_coords):
#         E_coord = model_lons[coord[0]]
#         W_coord = model_lons[coord[0]+1]
#         S_coord = model_lats[coord[1]]
#         N_coord = model_lats[coord[1]+1]
#         plt.vlines(E_coord, S_coord, N_coord, transform=ccrs.PlateCarree(), color='grey')
#         plt.vlines(W_coord, S_coord, N_coord, transform=ccrs.PlateCarree(), color='grey')
#         plt.hlines(S_coord, E_coord, W_coord, transform=ccrs.PlateCarree(), color='grey')
#         plt.hlines(N_coord, E_coord, W_coord, transform=ccrs.PlateCarree(), color='grey')
#     plt.show()

#     fig = plt.figure()
#     ax = plt.axes(projection=ccrs.PlateCarree())
#     ax.set_extent([-180,-105,55,90], crs=ccrs.PlateCarree())
#     ax.coastlines(linewidth=0.5)
#     for i_aero_file, f in enumerate(atom_aero_T):
#         if len(i_aero_high_lats[f]) > 0:
#             i_coords = i_nav_for_aero[f][i_aero_high_lats[f]]
#             plt.plot(atom_lons[i_coords], atom_lats[i_coords], transform=ccrs.PlateCarree(),
#                      linestyle='None', marker='.', ms=1, color='tab:orange')
#     for i,coord in enumerate(model_flight_coords):
#         E_coord = model_lons[coord[0]]
#         W_coord = model_lons[coord[0]+1]
#         S_coord = model_lats[coord[1]]
#         N_coord = model_lats[coord[1]+1]
#         plt.vlines(E_coord, S_coord, N_coord, transform=ccrs.PlateCarree(), color='grey')
#         plt.vlines(W_coord, S_coord, N_coord, transform=ccrs.PlateCarree(), color='grey')
#         plt.hlines(S_coord, E_coord, W_coord, transform=ccrs.PlateCarree(), color='grey')
#         plt.hlines(N_coord, E_coord, W_coord, transform=ccrs.PlateCarree(), color='grey')
#     plt.show()
#     fig = plt.figure()#figsize=(12,6),dpi=300)
#     ax = plt.axes(projection=ccrs.PlateCarree())
#     ax.set_extent([-105,-30,55,90], crs=ccrs.PlateCarree())
#     ax.coastlines(linewidth=0.5)
#     for i_aero_file, f in enumerate(atom_aero_T):
#         if len(i_aero_high_lats[f]) > 0:
#             i_coords = i_nav_for_aero[f][i_aero_high_lats[f]]
#             plt.plot(atom_lons[i_coords], atom_lats[i_coords], transform=ccrs.PlateCarree(),
#                      linestyle='None', marker='.', ms=1, color='tab:orange')
#     for i,coord in enumerate(model_flight_coords):
#         E_coord = model_lons[coord[0]]
#         W_coord = model_lons[coord[0]+1]
#         S_coord = model_lats[coord[1]]
#         N_coord = model_lats[coord[1]+1]
#         plt.vlines(E_coord, S_coord, N_coord, transform=ccrs.PlateCarree(), color='grey')
#         plt.vlines(W_coord, S_coord, N_coord, transform=ccrs.PlateCarree(), color='grey')
#         plt.hlines(S_coord, E_coord, W_coord, transform=ccrs.PlateCarree(), color='grey')
#         plt.hlines(N_coord, E_coord, W_coord, transform=ccrs.PlateCarree(), color='grey')
#     plt.show()

    n_conc_extracted_dict[suite] = {}
    diams_extracted_dict[suite] = {}
    for m,mode in enumerate(number_conc_aug[suite]):
        n_debug = 0
        n_conc_extracted_dict[suite][mode] = iris.cube.CubeList([])
        diams_extracted_dict[suite][mode] = iris.cube.CubeList([])
        model_lon_points = lon_coord.points
        model_lat_points = lat_coord.points
        for coord in model_flight_coords:
            lon = model_lon_points[coord[0]]
            if lon < 0:
                lon += 360
            lat = model_lat_points[coord[1]]
            lon_constraint = iris.Constraint(longitude = lon)
            lat_constraint = iris.Constraint(latitude = lat)
            N_cube = number_conc_aug[suite][mode].extract(lon_constraint)
            N_cube = N_cube.extract(lat_constraint)
            N_cube.data
            D_cube = radius_aug[suite][mode].extract(lon_constraint) * 2
            D_cube = D_cube.extract(lat_constraint)
            D_cube.data
            n_conc_extracted_dict[suite][mode].append(N_cube)
            diams_extracted_dict[suite][mode].append(D_cube)
            n_debug += 1
    # calculate total number in size ranges from model data
    first_bound = 1    #[nm]
    last_bound = 1e3
    model_n_bins = 500
    model_bins = np.logspace(np.log(first_bound),np.log(last_bound),num=model_n_bins,base=np.exp(1))
    model_bin_bounds = aero.calc_bin_boundaries(model_bins)
    model_dlogD = np.log10(model_bin_bounds[1:]) - np.log10(model_bin_bounds[:-1])
    model_size_limit_inds = []
    for n in size_bins:
        loc = np.nonzero(model_bin_bounds < n)[0][-1]
        model_size_limit_inds.append(loc + 1)
    n_levels = 41
    n_coords = len(model_flight_coords)
    z_coord = number_conc_aug[suite]['soluble_nucleation'].coord('atmosphere_hybrid_height_coordinate')
    atom_model_z[suite] = z_coord.points[:n_levels]/1000
    atom_model_N[suite] = np.zeros((n_size_bins, n_levels))
    dist = np.zeros((model_n_bins, n_levels, n_coords))
    dist_area_mean = np.zeros((model_n_bins, n_levels))
    for m,mode in enumerate(number_conc[suite]):
        Naero = n_conc_extracted_dict[suite][mode]
        D = diams_extracted_dict[suite][mode]
        for z in np.arange(n_levels):
            for i in np.arange(n_coords):
                pdf = aero.lognormal_pdf(np.log(D[i][z].data*1e09), np.log(constants.mode_sig[mode]), model_bins)
                dist[:,z,i] += pdf * Naero[i][z].data
    for z in np.arange(n_levels):
        dist_area_mean[:,z] = np.mean(dist[:,z],axis=1)

    for d in np.arange(n_size_bins):
        for z in np.arange(n_levels):
            i1 = model_size_limit_inds[d]
            i2 = model_size_limit_inds[d+1]
            atom_model_N[suite][d,z] = ma.sum(model_dlogD[i1:i2] * dist_area_mean[i1:i2,z])

# ---------------------------------------------------------------------
# LOAD MODEL DATA FOR FIGURE 7
model_output_path = '/gws/nopw/j04/asci/rprice/ukca_output/'
for suite in fprops.fig7_suites:
    try:
        a = np.any(N[suite])
        if a and verbose:
            print('\nAlready loaded {}'.format(suite))
    except KeyError:
        if verbose:
            print('\nLoading output from {}..'.format(suite))
        N[suite], bin_edges[suite], times[suite] = aero.get_ao2018_aerosol_conc(model_output_path, suite)

if verbose:
    print('\nMaking PDFs..')
n_pdf_bins = 25
hist_obs, hist_model, pdf_bins = aero.ao2018_melt_freeze_pdfs(n_pdf_bins, fprops.fig7_suites, N, bin_edges, times, hist_model=hist_model, pdf_bins=pdf_bins)
pdf_bins_mid = [0.5*(X[1:] + X[:-1]) for X in pdf_bins]

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
# FIGURE 2: time series and PDF of surface HIO3 concentration during AO2018
# plot time series and PDF on same figure with subplots
if verbose:
    print('\nMaking figure 3..')

ts.plot_IA_time_series_pdf(obs_hio3_means, obs_hio3_mean_times,
                           colocated_hio3_number, model_times,
                           hio3_hist_obs[1], hio3_pdf_bins_mid,
                           hio3_hist_model[1], hio3_pdf_bins_mid,
                           "figures/fig02")
if verbose:
    print('Done.')

# ---------------------------------------------------------------------
# FIGURE 3: time series of aerosol concentrations during AO2018
if verbose:
    print('\nMaking figure 2..')

obs_N = [ufp_running_means, dmps_N15_100_running_means, dmps_N100_500_running_means]
obs_T = [ufp_running_mean_times, dmps_running_mean_times, dmps_running_mean_times]
obs_stdev = [ufp_std_devs, dmps_N15_100_std_devs, dmps_N100_500_std_devs]
ts.plot_time_series_pdfs(fprops.fig3_suites, N, times, 
                         obs_N, obs_T, obs_stdev, 
                         hist_obs, hist_model, pdf_bins_mid,
                         'figures/fig03.pdf')
if verbose:
    print('Done.')
# ---------------------------------------------------------------------
# FIGURE 4: colocated aerosol height profiles for nucleation mode
if verbose:
    print('\nMaking figure 4..')

ts.time_series_3d(sol_nuc_N, zbl, fig4_times, fig4_heights, 'figures/fig04')

if verbose:
    print('Done.')

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
    Naero = model_N[suite][0]
    colour = fprops.colours[suite]
    linewidth = fprops.linewidths[suite]
    ax1.plot(Naero, heights, label=fprops.suite_labels[suite],
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
    Naero = model_N[suite][1]
    colour = fprops.colours[suite]
    linewidth = fprops.linewidths[suite]
    ax2.plot(Naero, heights, label=fprops.suite_labels[suite],
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

# ---------------------------------------------------------------------
# FIGURE 6: ATom profiles
# Make plot
fig = plt.figure(figsize=(16*fprops.cm,6*fprops.cm), dpi=300)
gs = fig.add_gridspec(ncols=3, nrows=1, wspace=0.15)

subfiglabels = ['(a)','(b)','(c)']

for d in np.arange(n_size_bins):
    ax = fig.add_subplot(gs[d])
    for i_aero_file, f in enumerate(atom_aero_T):
        if len(i_aero_high_lats[f]) > 0:
            ax.plot(atom_N_mean[f][:,d], atom_z_bins[f]/1000,
                     alpha=0.4, color='grey', zorder=1, linewidth=fprops.thick_line)
    ax.plot(atom_N_all_mean[:,d], z_bins_merged/1000, color='k', label='ATom mean', linewidth=fprops.thick_line)
    ax.plot(atom_N_all_median[:,d], z_bins_merged/1000, color='k', linestyle='dashed', label='ATom median', linewidth=fprops.thick_line)
    for s,suite in enumerate(fprops.fig5_suites):
        Naero = atom_model_N[suite][d]
        z = atom_model_z[suite]
        ax.plot(Naero, z, label=fprops.suite_labels[suite],
                zorder=4, color=fprops.colours[suite],
                linewidth=fprops.linewidths[suite])
    ax.grid()
    ax.set_xscale('log')
    ax.set_ylim(0,16)
    ax.set_xlim(1e-1,1e4)
    ax.set_title('{}\u2013{} nm'.format(size_bins[d],size_bins[d+1]), fontsize=fprops.label_fs)
    ax.set_title(subfiglabels[d], loc='left', fontsize=fprops.label_fs)
    ax.tick_params(axis='x', labelsize=fprops.ax_fs)
    ax.tick_params(axis='y', labelsize=fprops.ax_fs)
    if d == 0:
        ax.set_ylabel('Altitude [km]', fontsize=fprops.ax_label_fs)
    if d == 1:
        ax.set_xlabel('Particle concentration [cm$^{-3}$]', fontsize=fprops.ax_label_fs)
plt.legend(fontsize=fprops.legend_fs*0.8, ncol=2,
           columnspacing=1, borderpad=0.2,
           handletextpad=0.3, labelspacing=0.5,
           handlelength=1.5)
filename = 'figures/fig06.pdf'
plt.savefig(filename, bbox_inches="tight", facecolor='white', format='pdf')
plt.close()
if verbose:
    print('Done')
# ---------------------------------------------------------------------
# FIGURE 7: time series of aerosol concentrations during AO2018
if verbose:
    print('\nMaking figure 7..')

obs_N = [ufp_running_means, dmps_N15_100_running_means, dmps_N100_500_running_means]
obs_T = [ufp_running_mean_times, dmps_running_mean_times, dmps_running_mean_times]
obs_stdev = [ufp_std_devs, dmps_N15_100_std_devs, dmps_N100_500_std_devs]
ts.plot_time_series_pdfs(fprops.fig7_suites, N, times, 
                         obs_N, obs_T, obs_stdev, 
                         hist_obs, hist_model, pdf_bins_mid,
                         'figures/fig07.pdf', alpha_pale=0.6)
if verbose:
    print('Done.')

# =====================================================================
end = dt.datetime.now()
print('Finished script at {} in {}'.format(end, end-start))
# =====================================================================
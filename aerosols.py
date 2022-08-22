# aerosols.py
# Functions for loading, manipulating aerosol data
# =====================================================================
import numpy as np
import iris
import campaign_routes as cr
import constants
from glob import glob
import re
from files import save_cube
import datetime as dt
from files import get_csv_contents
import numpy.ma as ma
# =====================================================================
def lognormal_pdf(mean, sigma, arr):
    # returns pdf of specified lognormal dist at points given by arr
    # formula adapted from exmaple on numpy.random.lognormal webpage
    # first line (probably commented out) returns pdf as function of arr
    # second line returns pdf as a function of log(arr)
    #pdf = (np.exp(-(np.log(arr) - mean)**2 / (2 * sigma**2)) / (arr * sigma * np.sqrt(2 * np.pi)))
    pdf = (np.exp(-(np.log(arr) - mean)**2 / (2 * sigma**2)) / (sigma * np.sqrt(2 * np.pi)))

    return pdf 

def calc_bin_boundaries(bin_midpoints):
    # Function that calculates bin boundaries
    # from array of bin midpoints
    n_bins = len(bin_midpoints)
    bin_bounds = np.zeros((n_bins + 1))

    for j in np.arange(1, n_bins):
        # bin_bound[j] will be lower limit of bin with midpoint bin_midpoints[j]
        bin_bounds[j] = (bin_midpoints[j] + bin_midpoints[j-1]) / 2

    bin_bounds[0] = bin_midpoints[0] - (bin_bounds[1] - bin_midpoints[0])
    bin_bounds[-1] = bin_midpoints[-1] + (bin_midpoints[-1] - bin_bounds[-2])

    return bin_bounds

def get_cube_times(cube, ao_drift=False):
    '''
    Return cube timesteps as dt objects
    '''
    model_times = cube.coord('time').units.num2date(cube.coord('time').points)
    model_times = np.array([dt.datetime(T.year,T.month,T.day,T.hour,T.minute,T.second) for T in model_times])
    if ao_drift:
        time_inds = cr.take_drift_timesteps_from_cube(cube)[1]
        model_times = model_times[time_inds]
    return model_times

def integrate_modes(number_conc, radius, bins):
    '''
    Integrate over the UKCA modes to get a size distribution
    on specified diameter bins (dN/dlogD)
    '''
    n_time = len(number_conc['soluble_nucleation'])
    n_bins = len(bins)
    dNdlogD = np.zeros((n_time, n_bins))
    for t in np.arange(n_time):
        for i,mode in enumerate(radius):
            R = radius[mode][t]
            N = number_conc[mode][t]
            pdf = lognormal_pdf(np.log(R.data*1e09*2), np.log(constants.mode_sig[mode]),bins*1e09)
            dist = pdf*N.data
            dNdlogD[t] += dist
    return dNdlogD

def total_number_from_dN(dNdlogD, bins, target_bins, base='e'):
    '''
    From an array of dN/dlogD in bins, sum to get N
    in specified target_bins
    '''
    if base=='e':
        dlogD = np.log(bins[1:]) - np.log(bins[:-1])
        logD = np.log(bins)
    elif base=='10':
        dlogD = np.log10(bins[1:]) - np.log10(bins[:-1])
        logD = np.log10(bins)
    else:
        print("[N_from_dN] ERROR: base should be 'e' or '10'")
    n_time = len(dNdlogD)
    n_target_bins = len(target_bins) - 1
    # find indices of diam list corresponding to the target bins
    # start with 0 - first index of smallest size range
    target_inds = []
    for n in target_bins:
        if n == 0:
            target_inds.append(0)
        else:
            loc = np.nonzero(bins < n)[0][-1]
            # loc is index of largest diameter bin that's still less than n
            target_inds.append(loc + 1)
#     N = np.zeros((n_time, n_target_bins))
#     for t in np.arange(n_time):
#         for d,D in enumerate(target_inds[:-1]):
#             dN = ma.sum(dlogD[D:target_inds[d+1]] * dNdlogD[t,D:target_inds[d+1]])
#             N[t,d] += dN
#     return N

    N = ma.zeros((n_time, n_target_bins))
    for d,D in enumerate(target_inds[:-1]):
        N[:,d] = np.trapz(dNdlogD[:,D:target_inds[d+1]], x=logD[D:target_inds[d+1]], axis=1)
    return N

def load_3hr_surface_aerosol_number_conc(model_output_path, suite):
    '''From model output files of number mixing ratio,
       calculate number concentration (ie go from number/air molecule
       -> number/cm3). This is for 3hrly means at surface only.
    '''
    # Calculate number of air molecules per unit volume:
    # using an estimate for air density because I don't have it at 3 hrly res
    # estimate is from taking an average of daily air density values at ship
    # about 2.5% variation from this average
    air_density = 1.33 # kg m-3
    particle_density_of_air = (air_density/constants.mm_air)*constants.avo
    # Load aerosol data as [# / molecules of air] (number mixing ratio)
    nmr_stashcodes = [ # model output filenames are based on UM stashcodes
        'm01s34i101',    # sol nucleation
        'm01s34i103',    # sol aitken
        'm01s34i107',    # sol accumulation
        'm01s34i113',    # sol coarse
        'm01s34i119'     # insol aitken
        ]
    nmr_path = "{}/{}/All_time_steps/pk_files/".format(model_output_path, suite)
    nmr_cubes_list = [iris.load(glob(nmr_path+'*{}*'.format(stash)))[0] for stash in nmr_stashcodes]
    ncon_cubes = {}
    for cube in nmr_cubes_list:
        label = [mode for mode in constants.modes if "_"+mode in cube.long_name.lower()][0]
        if cube.shape[1] == 2:
            cube = cube[:,0]
        cube.data
        n_conc = cube*particle_density_of_air
        n_conc.long_name = "number_concentration_of_{}_mode_aerosol".format(label)
        n_conc.units = "m-3"
        n_conc.convert_units('cm-3')
        ncon_cubes[label] = n_conc
    return ncon_cubes

def load_3hr_surface_aerosol_mass_conc(model_output_path, suite):
    '''From model output files of mass mixing ratio,
       calculate mass concentration (ie go from mass/mass of air
       -> kg/cm3). This is for 3hrly means at surface only.
    '''
    # using an estimate for air density because I don't have it at 3 hrly res
    # estimate is from taking an average of daily air density values at ship
    # about 2.5% variation from this average
    air_density = 1.33 # kg m-3
    # Load aerosol data as [# / molecules of air] (number mixing ratio)
    mmr_stashcodes = {
        'sulfuric_acid_in_soluble_nucleation': 'm01s34i102',
        'sulfuric_acid_in_soluble_aitken': 'm01s34i104',
        'black_carbon_in_soluble_aitken': 'm01s34i105',
        'particulate_organic_matter_in_soluble_aitken': 'm01s34i106',
        'sulfuric_acid_in_soluble_accumulation': 'm01s34i108',
        'black_carbon_in_soluble_accumulation': 'm01s34i109',
        'particulate_organic_matter_in_soluble_accumulation': 'm01s34i110',
        'seasalt_in_soluble_accumulation': 'm01s34i111',
        'sulfuric_acid_in_soluble_coarse': 'm01s34i114',
        'black_carbon_in_soluble_coarse': 'm01s34i115',
        'particulate_organic_matter_in_soluble_coarse': 'm01s34i116',
        'seasalt_in_soluble_coarse': 'm01s34i117',
        'black_carbon_in_insoluble_aitken': 'm01s34i120',
        'particulate_organic_matter_in_insoluble_aitken': 'm01s34i121',
        'particulate_organic_matter_in_soluble_nucleation': 'm01s34i126'
        }
    mmr_path = "{}/{}/All_time_steps/pk_files/".format(model_output_path, suite)
    mmr_cubes_list = [iris.load(glob(mmr_path+'*{}*'.format(stash)))[0] for stash in mmr_stashcodes]
    mcon_cubes = {}
    for cube in mmr_cubes_list:
        label = [cpt for cpt in mmr_stashcodes.keys() if "_"+cpt in cube.long_name.lower()][0]
        if cube.shape[1] == 2:
            cube = cube[:,0]
        cube.data
        m_conc = cube*air_density
        m_conc.units = "kg m-3"
        m_conc.convert_units('kg cm-3')
        mcon_cubes[label] = m_conc
    return mcon_cubes

def radius_from_number_and_mass_conc(number_conc, mass_conc):
    '''Use aerosol number in [cm-3] and mass in [kg cm-3]
       to calculate the mean radius of each mode.
       Inputs should be dictionaries of cubes where
       number_conc keys are the modes and mass_conc keys are
       sulfuric_acid_in_soluble_nucleation, etc etc
    '''
    # Add up mass per mode
    mass_conc_per_mode = {}
    for mode in constants.modes:
        cubes_to_sum = [mass_conc[label] for label in mass_conc if "_"+mode in label]
        mass_conc_per_mode[mode] = np.sum(cubes_to_sum)
    # Calculate mean volume of each mode
    volume_cubes_list = iris.cube.CubeList([])
    for label in mass_conc:
        cube = mass_conc[label]
        mode = [mode for mode in constants.modes if "_"+mode in label][0]
        species = [species for species in constants.species_rho.keys() if species in label][0]
        volume = cube / constants.species_rho[species]
        volume.long_name = "volume_fraction_of_{}_in_{}_mode".format(species,mode)
        volume_cubes_list.append(volume)
    volume_per_mode = {}
    for mode in constants.modes:
        cubes_to_sum = [cube for cube in volume_cubes_list if "_"+mode in cube.long_name]
        volume_per_mode[mode] = np.sum(cubes_to_sum)
    # calculate radius
    mode_radius_cubes = {}
    for i,mode in enumerate(constants.modes):
        volume = volume_per_mode[mode]
        number = number_conc[mode]
        sigma = constants.mode_sig[mode]
        radius = 0.5*iris.analysis.maths.exponentiate((6*volume/number)/(
                                     np.pi*np.exp(4.5*(np.log(sigma))**2)),1./3.)
        radius.long_name = "radius_of_{}_mode".format(mode)
        radius.convert_units('m')
        mode_radius_cubes[mode] = radius
    return mode_radius_cubes

def get_mode_radius(model_output_path, suite, number_conc=None, mass_conc=None):
    '''
    Function loads saved radius files if they exist
    calculates radius from number and mass if not
    '''
    calculate_modes = []
    radius_folder = 'data/processed/'
    radius_files = []
    # check which modes have been saved
    for mode in constants.modes:
        file = glob('{}*{}*radius*{}*'.format(radius_folder, suite, '_'+mode))
        if len(file) == 1:
            radius_files.append(file[0])
        else:
            calculate_modes.append(mode)

    if calculate_modes:
        # if some modes need calculating
        # calculate..
        if number_conc==None:
            number_conc = load_3hr_surface_aerosol_number_conc(model_output_path, suite)
        if mass_conc==None:
            mass_conc = load_3hr_surface_aerosol_mass_conc(model_output_path, suite)
        radius = radius_from_number_and_mass_conc(number_conc, mass_conc)
        # ..then save
        saving_folder = 'data/processed/'
        for mode in radius:
            cube = radius[mode]
            save_name = 'L1_{}_{}.nc'.format(suite, cube.long_name)
            save_name = re.sub("\s+", "_", save_name)    
            save_cube(cube, saving_folder, save_name)
    else:
        # can all be loaded
        radius = {}
        for m,mode in enumerate(constants.modes):
            radius[mode] = iris.load(radius_files[m])[0]
            if radius[mode].shape[1] == 2:
                radius[mode] = radius[mode][:,0]
            radius[mode].data
    return radius

def save_ao2018_aerosol_conc(suite, N, bin_edges, times):
    '''
    Save aerosol concentrations from model output
    colocated with AO2018 as a .csv file.
    File will look like this:

    Timestep, D[0], D[1],..
    YYYY-MM-DD HH:MM:SS, N[0,0], N[0,1]..0,
    YYYY-MM-DD HH:MM:SS, N[1,0], N[1,1]..0,

    where D is the lower limit of that bin
    and N is the data.
    '''
    output = 'Timestep,'
    for Dp in bin_edges:
        output += "{},".format(str(Dp))
    output = output[:-1]
    output += "\n"
    fmt = '%Y-%m-%d %H:%M:%S'
    for t in np.arange(len(times)):
        output += "{},".format(dt.datetime.strftime(times[t],fmt))
        for d in np.arange(len(bin_edges)-1):
            output += "{},".format(N[t,d])
        # add a 0 at end of column because len(bin_edges) = len(N[0])+1
        output += "0\n"

    save_path = 'data/processed/'
    filename = "{}_ao2018_colocated_psd.txt".format(suite)

    file_writer = open(save_path+filename, "w")
    file_writer.write(output)
    file_writer.close()
    return

def calculate_ao2018_aerosol_conc(model_output_path, suite):
    '''
    Get aerosol concentration at surface for 2.5-15, 15-100,
    100-500 nm size ranges from model output,
    colocated with Oden. Returns np arrays with concentration
    and corresponding timesteps (as dt objects)
    '''
    # Get the output
    number_conc = load_3hr_surface_aerosol_number_conc(model_output_path, suite)
    radius = get_mode_radius(model_output_path, suite, number_conc=number_conc)
    times = get_cube_times(radius['soluble_nucleation'], ao_drift=True)
    # Do colocation
    colocated_number_conc = {}
    colocated_radius = {}
    for mode in constants.modes:
        cube = number_conc[mode]
        colocated_number_conc[mode] = cr.colocate_with_ao2018_drift(cube, constants.model_res)
        cube = radius[mode]
        colocated_radius[mode] = cr.colocate_with_ao2018_drift(cube, constants.model_res)
    # Integrate modal output for particle concentrations
    first_bound = 1e-9    #[m]
    last_bound = 10e-6
    n_bins = 500
    bins = np.logspace(np.log(first_bound),np.log(last_bound),num=n_bins,base=np.exp(1))
    dNdlogD = integrate_modes(colocated_number_conc, colocated_radius, bins)
    N_limits = [2.5e-9, 15e-9, 100e-9, 500e-9]  # nm
    N = total_number_from_dN(dNdlogD, bins, N_limits)
    save_ao2018_aerosol_conc(suite, N, N_limits, times)
    return N, N_limits, times

def load_ao2018_aerosol_conc(suite):
    '''
    Load time series of aerosol concentrations during AO2018
    as written to file by save_ao2018_aerosol_conc
    '''
    filename = 'data/processed/{}_ao2018_colocated_psd.txt'.format(suite)
    file_contents, n_rows, n_cols = get_csv_contents(filename)
    fmt = '%Y-%m-%d %H:%M:%S'
    times = np.array([dt.datetime.strptime(T, fmt) for T in file_contents[1:,0]])
    bin_edges = file_contents[0,1:].astype(float)
    N = file_contents[1:,1:-1].astype(float)
    return N, bin_edges, times

def get_ao2018_aerosol_conc(model_output_path, suite):
    '''
    Load the time series if it's been saved,
    calculate if not
    '''
    filename = glob('data/processed/{}_ao2018_colocated_psd.txt'.format(suite))
    if filename:
        print('Loading.')
        N, bin_edges, times = load_ao2018_aerosol_conc(suite)
    else:
        print('Calculating.')
        N, bin_edges, times = calculate_ao2018_aerosol_conc(model_output_path, suite)
    return N, bin_edges, times

def load_ao2018_ufp(data_path):
    '''
    Load the ultrafine particle (N_2.5) dataset from the AO2018 campaign
    '''
    file_contents, n_rows, n_cols = get_csv_contents(data_path)
    n_time_ufp = n_rows - 1
    timestamps_str = file_contents[1:,0]
    fmt = '%Y-%m-%d %H:%M:%S'
    ufp_timestamps = np.array([dt.datetime.strptime(T, fmt) for T in timestamps_str])
    ufp_conc = file_contents[1:,1].astype(float)
    # set minimum threshold
    ufp_conc[(ufp_conc<0.5)] = 0.
    ufp_conc = ma.masked_where(ufp_conc==0, ufp_conc)
    return ufp_conc, ufp_timestamps

def load_ao2018_dmps(data_path):
    '''
    Load the DMPS dataset from the AO2018 campaign
    '''
    file_contents, n_rows, n_cols = get_csv_contents(data_path)
    n_time_dmps = n_rows - 1
    n_bins_dmps = n_cols - 2
    timestamps_str = file_contents[1:,0]
    fmt = '%Y-%m-%d %H:%M:%S'
    dmps_timestamps = np.array([dt.datetime.strptime(T, fmt) for T in timestamps_str])
    dmps_diam_mid = file_contents[0,1:-1].astype(float)
    dmps_dNdlogD = file_contents[1:,1:-1].astype(float)
    qc_flag = file_contents[1:,-1].astype(int)
    mask = np.zeros(dmps_dNdlogD.shape)
    mask[(qc_flag==2)] = 1
    dmps_dNdlogD = ma.masked_array(dmps_dNdlogD, mask=mask)
    return dmps_dNdlogD, dmps_diam_mid, dmps_timestamps

def running_mean(X, times, n_hours, spacing):
    '''
    Calculates running means of array X for n_hours window
    e.g. if n_hours=3, 3 hourly means.
    spacing is number of timesteps to move each mean
    e.g. if spacing=1, give running mean for each timestep
    '''
    running_means = []
    std_devs = []
    running_mean_times = []

    lower = 0
    upper = 0
    while upper < len(times) - 1:
        # get index of T+n hrs
        upper = np.nonzero(times < times[lower] + dt.timedelta(hours=n_hours))[0][-1]
        subset = np.arange(lower, upper+1) # indices of 3hr window
        subset_times = times[subset]
        subset_range = subset_times[-1] - subset_times[0]
        subset_median_time = subset_times[0] + subset_range/2
        running_mean_times.append(subset_median_time)
        running_means.append(ma.mean(X[subset], axis=0))
        std_devs.append(ma.std(X[subset], axis=0))
        lower += spacing
    running_means = ma.array(running_means)
    running_mean_times = np.array(running_mean_times)
    std_devs = ma.array(std_devs)
    return running_means, std_devs, running_mean_times

def ao2018_melt_freeze_pdfs(n_pdf_bins, suites, N, Dp_bin_edges, times,
                            obs=None, freeze_up=dt.datetime(2018,8,27)):
    '''
    Creates PDFs of observations and model output (for suites)
    of aerosol concentrations in 3 sizes during AO2018
    split into the melt and freeze periods
    returns PDFS for observations, model and the bins
    '''
    if obs == None:
        # Load data if need be
        ao2018_data_file_dir = "/home/users/eersp/ao2018_observations/"
        ufp_data_in = ao2018_data_file_dir + "ao2018-aerosol-ufp.csv"
        dmps_data_in = ao2018_data_file_dir + 'DMPS_for_Ruth.csv'
        ufp_conc, ufp_times = load_ao2018_ufp(ufp_data_in)
        dmps_dNdlogD, dmps_diam_mid, dmps_times = load_ao2018_dmps(dmps_data_in)
        # integrate dmps data
        i_15_dmps = np.nonzero(dmps_diam_mid < 15.0)[0][-1]
        i_100_dmps = np.nonzero(dmps_diam_mid < 100.0)[0][-1]
        i_500_dmps = np.nonzero(dmps_diam_mid < 500.0)[0][-1]
        dmps_N_15_100 = np.trapz(dmps_dNdlogD[:,i_15_dmps:i_100_dmps], x=np.log10(dmps_diam_mid[i_15_dmps:i_100_dmps]), axis=1)
        dmps_N_100_500 = np.trapz(dmps_dNdlogD[:,i_100_dmps:i_500_dmps], x=np.log10(dmps_diam_mid[i_100_dmps:i_500_dmps]), axis=1)
        # take running means
        ufp_running_means, stdev, ufp_running_mean_times = running_mean(ufp_conc, ufp_times, 3, 1)
        dmps_N15_100_running_means, stdev, dmps_N100_500_running_means = running_mean(dmps_N_15_100, dmps_times, 3, 1)
        dmps_N100_500_running_means, stdev, dmps_running_mean_times = running_mean(dmps_N_100_500, dmps_times, 3, 1)

    # separate obs into melt and freeze
    ufp_melt_times, ufp_freeze_times = cr.ao2018_melt_freeze_times(freeze_up, ufp_running_mean_times)
    ufp_melt = ufp_running_means[ufp_melt_times]
    ufp_freeze = ufp_running_means[ufp_freeze_times]
    dmps_melt_times, dmps_freeze_times = cr.ao2018_melt_freeze_times(freeze_up, dmps_running_mean_times)
    dmps_N15_100_melt = dmps_N15_100_running_means[dmps_melt_times]
    dmps_N15_100_freeze = dmps_N15_100_running_means[dmps_freeze_times]
    dmps_N100_500_melt = dmps_N100_500_running_means[dmps_melt_times]
    dmps_N100_500_freeze = dmps_N100_500_running_means[dmps_freeze_times]
    # now the model output
    model_melt_times = {}
    model_freeze_times = {}
    model_melt = {}
    model_freeze = {}
    for i,suite in enumerate(suites):
        model_melt_times[suite], model_freeze_times[suite] = cr.ao2018_melt_freeze_times(freeze_up, times[suite])
        model_melt[suite] = N[suite][model_melt_times[suite]]
        model_freeze[suite] = N[suite][model_freeze_times[suite]]
    n_Dp_bins = len(Dp_bin_edges[suites[0]]) - 1
    # get max values
    obs_data_all = [ufp_running_means,dmps_N15_100_running_means,dmps_N100_500_running_means]
    obs_data_melt = [ufp_melt,dmps_N15_100_melt,dmps_N100_500_melt]
    obs_data_freeze = [ufp_freeze,dmps_N15_100_freeze,dmps_N100_500_freeze]
    max_N = [np.amax(X) for X in obs_data_all]
    max_N_log = np.log10(max_N) 
    # define bins
    pdf_bins = [np.logspace(-4, N, n_pdf_bins+1) for N in max_N_log]
    pdf_bins_mid = [0.5*(X[1:] + X[:-1]) for X in pdf_bins]

    # Observations
    hist_obs = np.zeros((n_Dp_bins, 2, n_pdf_bins))
    for n in np.arange(n_Dp_bins):
        hist_obs[n,0] = np.histogram(obs_data_melt[n], density=True, bins=pdf_bins[n])[0]
        hist_obs[n,1] = np.histogram(obs_data_freeze[n], density=True, bins=pdf_bins[n])[0]
    # Model
    hist_model = {}
    for s,suite in enumerate(suites):
        hist_model[suite] = np.zeros((n_Dp_bins, 2, n_pdf_bins))
        for n in np.arange(n_Dp_bins):
            hist_model[suite][n,0] = np.histogram(model_melt[suite][:,n], density=True, bins=pdf_bins[n])[0]
            hist_model[suite][n,1] = np.histogram(model_freeze[suite][:,n], density=True, bins=pdf_bins[n])[0]
    return hist_obs, hist_model, pdf_bins_mid

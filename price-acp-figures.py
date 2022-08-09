# price-acp-figures.py 
# Main script
# RP Aug 2022
# =====================================================================
# calc_bin_boundaries
# import numpy as np

# lognormal_pdf
# import numpy as np

# load_3hr_surface_aerosol_number_conc
import constants
from glob import glob
import iris

# radius_from_numer_and_mass
import numpy as np
# import constants

# get_radius
# import constants
# from glob import glob
# import iris
import re
# from files import save_cube

# main body
import datetime as dt
import campaign_routes as cr
from files import save_cube
from obs_map import plot_obs_map

import matplotlib.pyplot as plt
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
            radius[mode].data
    return radius

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

def total_number_from_dN(dNdlogD, bins, target_bins):
    '''
    From an array of dN/dlogD in bins, sum to get N
    in specified target_bins
    '''
    dlogD = np.log(bins[1:]) - np.log(bins[:-1])
    n_time = len(dNdlogD)
    n_target_bins = len(target_bins) - 1
    # find indices of diam list corresponding to the target bins
    # start with 0 - first index of smallest size range
    target_inds = []
    for n in target_bins:
        loc = np.nonzero(bins < n)[0][-1]
        # loc is index of largest diameter bin that's still less than n
        target_inds.append(loc + 1)
    N = np.zeros((n_time, n_target_bins))
    for t in np.arange(n_time):
        for d,D in enumerate(target_inds[:-1]):
            dN = np.sum(dlogD[D:target_inds[d+1]] * dNdlogD[t,D:target_inds[d+1]])
            N[t,d] += dN
    return N

def get_ao2018_aerosol_conc(model_output_path, suite):
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
        if verbose:
            print(mode)

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
    return N, times

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

N, times = get_ao2018_aerosol_conc(model_output_path, suite)

if verbose:
    print('\nPlotting a time series to check functions')

fig = plt.figure(figsize=(25,5), dpi=300)
plt.plot(times, N[:,0])
plt.yscale('log')
plt.grid()
plt.show()
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
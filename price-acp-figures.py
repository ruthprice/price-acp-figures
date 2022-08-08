# price-acp-figures.py 
# Main script
# RP Aug 2022
# =====================================================================
# load_3hr_surface_aerosol_number_conc
import constants
from glob import glob
import iris
# radius_from_numer_and_mass
# import constants
# colocate_with_ao2018_drift
import campaign_routes as cr
import numpy as np
# main body
import datetime as dt
from obs_map import plot_obs_map
# =====================================================================
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
        number = ncon_cubes[mode]
        sigma = constants.mode_sig[mode]
        radius = 0.5*iris.analysis.maths.exponentiate((6*volume/number)/(
                                     np.pi*np.exp(4.5*(np.log(sigma))**2)),1./3.)
        radius.long_name = "radius_of_{}_mode".format(mode)
        radius.convert_units('m')
        mode_radius_cubes[mode] = radius
    return mode_radius_cubes

def take_drift_timesteps_from_cube(cube):
    # function takes an iris cube and returns only the cube timesteps
    # within AO2018 drift period. Also returns the indices
    # of the cube corresponding to these timesteps
    drift_begin = dt.datetime(year=2018,month=8,day=2)
    drift_end = dt.datetime(year=2018,month=9,day=19)

    drift_time_inds = []
    cube_datestamps = cube.coord('time').units.num2date(cube.coord('time').points)
    for t,datestamp in enumerate(cube_datestamps):
        if datestamp >= drift_begin and datestamp <= drift_end:
            drift_time_inds.append(t)

    if drift_time_inds == []:
        # cube timesteps don't cover drift period
        print("[take_drift_timesteps_from_cube] no drift times found, returning original cube")
        return [cube, drift_time_inds]
    else:
        return [cube[drift_time_inds], drift_time_inds]

def colocate_with_ao2018_drift(cube, model_resolution):
    # function takes iris cube of whatever time and space
    # and restricts to location of AO2018 ship for each timestep during drift
    ship_coords = cr.get_ao2018_track()

    if take_drift_timesteps_from_cube(cube)[1] == []:
        # cube doesn't contain data from drift time
        print("[colocate_with_ao2018_drift] can't colocate, no data from drift period")
        return []

    cube = take_drift_timesteps_from_cube(cube)[0]
    timestamps = cube.coord('time').units.num2date(cube.coord('time').points)
    timestamps = [dt.datetime(T.year,T.month,T.day,T.hour,T.minute,T.second) for T in timestamps]
    n_drift_time = len(timestamps)

    print("[colocate_with_ao2018] constraining data to ship location")
    colocated_cubes = iris.cube.CubeList([])
    for t in np.arange(n_drift_time):
        T = timestamps[t]
        ship_coord_ind = cr.find_ship_at_time(T)
        if ship_coord_ind == None:
            print("ERROR: can't find ship at {}".format(timestamps[t]))
            colocated_cubes.append(None)
        else:
            ship_lon = ship_coords['lon'][ship_coord_ind]
            ship_lat = ship_coords['lat'][ship_coord_ind]
            lon_constraint1 = iris.Constraint(longitude=lambda cell: cell > ship_lon-(0.5*model_resolution[0]))
            lon_constraint2 = iris.Constraint(longitude=lambda cell: cell < ship_lon+(0.5*model_resolution[0]))
            lat_constraint1 = iris.Constraint(latitude =lambda cell: cell > ship_lat-(0.5*model_resolution[1]))
            lat_constraint2 = iris.Constraint(latitude =lambda cell: cell < ship_lat+(0.5*model_resolution[1]))
            colocated_cube = cube[t].extract(lon_constraint1 & lon_constraint2 & lat_constraint1 & lat_constraint2)
            colocated_cubes.append(colocated_cube)

    return colocated_cubes


start = dt.datetime.now()
verbose = True

# ---------------------------------------------------------------------
# PRE-PROCESSING
# ---------------------------------------------------------------------
model_output_path = '/gws/nopw/j04/asci/rprice/ukca_output/'
suite = 'u-ci572'

if verbose:
    print('Loading aerosol number concentrations..')

ncon_cubes = load_3hr_surface_aerosol_number_conc(model_output_path, suite)

if verbose:
    print(ncon_cubes)
    print('\nLoading aerosol mass concentrations..')

mcon_cubes = load_3hr_surface_aerosol_mass_conc(model_output_path, suite)

if verbose:
    print(mcon_cubes)
    print('\nCalculating radius..')

radius = radius_from_number_and_mass_conc(ncon_cubes, mcon_cubes)

if verbose:
    print(radius)

# if verbose:
#     print('\nColocating cubes with AO2018..')

# colocated_ncon_cubes = {}
# for mode in ncon_cubes:

#     if verbose:
#         print(mode)

#     cube = ncon_cubes[mode]
#     colocated_ncon_cubes[mode] = colocate_with_ao2018_drift(cube, constants.model_res)

# if verbose:
#     print(colocated_ncon_cubes)
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
# constants.py
# Thermodynamic constants etc
# =====================================================================
import iris.coords
# Thermodynamcs
mm_air = 28.991e-3  # molar mass of air, J kg-1 K-1
mm_hio3 = 175.91e-3 # kg mole-1
avo = 6.022e23      # Avogadro's number, J kg-1 K-1
p0 = iris.coords.AuxCoord(1000.0,
                          long_name='reference_pressure',
                          units='hPa')
p0.convert_units('Pa')
Rd=287.05 # J/kg/K
cp=1005.46 # J/kg/K
Rd_cp=Rd/cp
R_specific=iris.coords.AuxCoord(287.058,
            long_name='R_specific',
            units='J-kilogram^-1-kelvin^-1')#J/(kg K)
mm_air_coord=iris.coords.AuxCoord(28.991e-3,
            long_name='Molar mass of air',
            units='kilogram-mole^-1')#J/(kg K)
avo_coord=iris.coords.AuxCoord(6.022e23,
            long_name='Avogadros number - particles per mol',
            units='mole^-1')#J/(kg K)

# UKCA
mode_sig = {
    'soluble_nucleation': 1.59,
    'soluble_aitken':1.59,
    'soluble_accumulation': 1.4,
    'soluble_coarse': 2.0,
    'insoluble_aitken': 1.59,
    }
modes = [mode for mode in mode_sig]
species_rho = {
    # these are given as iris coords because they're used with cubes
	'black_carbon': iris.coords.AuxCoord(1500,units='kg m-3'),
	'particulate_organic_matter': iris.coords.AuxCoord(1500,units='kg m-3'),
	'sulfuric_acid': iris.coords.AuxCoord(1769,units='kg m-3'),
	'seasalt': iris.coords.AuxCoord(2650,units='kg m-3')
	}
model_res = [1.875, 1.25]     # [lon,lat] in degrees

# constants.py
# Thermodynamic constants etc
# =====================================================================
import iris.coords
# Thermodynamcs
mm_air = 28.991e-3 # molar mass of air, J kg-1 K-1
avo = 6.022e23     # Avogadro's number, J kg-1 K-1

# UKCA
mode_sig = {
    'soluble_nucleation': 1.59,
    'soluble_aitken':1.59,
    'soluble_accumulation': 1.4,
    'soluble_coarse': 2.0,
    'insoluble_aitken': 1.59,
    }
modes = mode_sig.keys()
species_rho = {
    # these are given as iris coords because they're used with cubes
	'black_carbon': iris.coords.AuxCoord(1500,units='kg m-3'),
	'particulate_organic_matter': iris.coords.AuxCoord(1500,units='kg m-3'),
	'sulfuric_acid': iris.coords.AuxCoord(1769,units='kg m-3'),
	'seasalt': iris.coords.AuxCoord(2650,units='kg m-3')
	}
model_res = [1.875, 1.25]     # [lon,lat] in degrees

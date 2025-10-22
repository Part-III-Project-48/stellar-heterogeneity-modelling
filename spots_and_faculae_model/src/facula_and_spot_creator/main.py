# read in data

# for a given FeH; log_g; sum all the spectra according to w_i F_i where i \member of T_eff values
# and we want 1 of the w_i to be >> than the rest (as that will represent the star)
# the rest will be small and I guess randomised
# store FeH log_g {w_i's for each T_eff in some structured way}

from matplotlib import pyplot as plt
import numpy as np

from phoenix_grid_creator.fits_to_hdf5 import TEFF_COLUMN, WAVELENGTH_COLUMN, FLUX_COLUMN
from phoenix_grid_creator.basic_plotter import get_hdf5_data
 
table = get_hdf5_data()

# debug variables
FeH : float = 0.0
log_g : float = 4.0

star_T_eff_Kelvin : float = 3000

# max dT of facula, spots : inclusive of endpoints
delta_T_max_Kelvin : float = 500

spectra_grid_temperature_resolution_Kelvin : float = 50

min_spot_temperature : float = star_T_eff_Kelvin - delta_T_max_Kelvin
max_spot_temperature : float = star_T_eff_Kelvin + delta_T_max_Kelvin
valid_spot_temperatures : np.array = np.arange(min_spot_temperature,  max_spot_temperature + spectra_grid_temperature_resolution_Kelvin, spectra_grid_temperature_resolution_Kelvin)

# remove the star temperature itself; as that is already accounted for
valid_spot_temperatures = valid_spot_temperatures[valid_spot_temperatures != star_T_eff_Kelvin]

# get spectra for each temperature; sum them according to some weights

def generate_weights(star_temperature : float, spot_temperatures : np.array) -> dict:
	"""
	creates a dictionary with kvp of the form temperature : weight such that the sum of all weights in the dictionary = 1
	"""
	spot_weights = np.random.random((len(spot_temperatures)))
	# normalise
	spot_weights = spot_weights / np.sum(spot_weights)
	
	return dict(zip(spot_temperatures, spot_weights))

t_effs_and_weights = generate_weights(star_T_eff_Kelvin, valid_spot_temperatures)

# if you want to double check its 1
# print(sum(t_effs_and_weights.values()))

# print(t_effs_and_weights)
 
# this will be our total spectrum; combined from the star + all spots / faculae sources
# keep columns and metadata; but no rows // data
spectrum = table[:0]

from astropy import units as u
from astropy.table import vstack, QTable


for (temperature, weight) in t_effs_and_weights.items():
	# read in wavelength, flux value for temperature
	# if this doesnt exist then probably just print a msg and skip
	# or could interpolate on the fly; might be a nicer method
	
	# astropy masks require units (if the relevant column is unit-ed)
	sub_spectrum = table[table[TEFF_COLUMN] == temperature * u.K]
	sub_spectrum[FLUX_COLUMN] *= weight
	
	# same logic as in fits_to_hdf5.py, but just for an astropy qtable instead of for a pandas dataframe
	if len(spectrum) != 0:
		spectrum = vstack([spectrum, sub_spectrum])
	else:
		spectrum = sub_spectrum
		
plt.plot(spectrum[WAVELENGTH_COLUMN], spectrum[FLUX_COLUMN], linewidth=1)
plt.show()
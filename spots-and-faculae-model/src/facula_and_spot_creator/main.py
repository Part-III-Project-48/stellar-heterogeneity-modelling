# read in data

# for a given FeH; log_g; sum all the spectra according to w_i F_i where i \member of T_eff values
# and we want 1 of the w_i to be >> than the rest (as that will represent the star)
# the rest will be small and I guess randomised
# store FeH log_g {w_i's for each T_eff in some structured way}

import numpy as np
import pandas as pd

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

test = np.random.random((5))
print(test)

def generate_weights(star_temperature : float, spot_temperatures : np.array) -> dict:
	"""
	creates a dictionary with kvp of the form temperature : weight such that the sum of all weights in the dictionary = 1
	"""
	spot_weights = np.random.random((len(spot_temperatures)))
	# normalise
	spot_weights = spot_weights / np.sum(spot_weights)
	
	return dict(zip(spot_temperatures, spot_weights))

t_effs_and_weights = generate_weights(star_T_eff_Kelvin, valid_spot_temperatures)

# check its 1
print(sum(t_effs_and_weights.values()))

WAVELENGTH_COLUMN : str = "Wavelength / Angstroms"
FLUX_COLUMN : str = "Total Flux / Counts"

# this will be our total spectrum; combined from the star + all spots / faculae sources
spectrum = pd.DataFrame(columns=[WAVELENGTH_COLUMN, FLUX_COLUMN])

print(spectrum)
# for a given FeH; log_g; sum all the spectra according to w_i F_i where i \member of T_eff values
# and we want 1 of the w_i to be >> than the rest (as that will represent the star)
# the rest will be small and I guess randomised
# store FeH log_g {w_i's for each T_eff in some structured way}

# external
from matplotlib import pyplot as plt
import numpy as np
import scipy as sp
from astropy import units as u

# internal
from phoenix_grid_creator.fits_to_hdf5 import TEFF_COLUMN, WAVELENGTH_COLUMN, FLUX_COLUMN, LOGG_COLUMN, FEH_COLUMN
from phoenix_grid_creator.basic_plotter import get_hdf5_data

# # # # # helper functions # # # # #

def generate_weights(star_temperature : float, spot_temperatures : np.array) -> dict:
	"""
	creates a dictionary with kvp of the form temperature : weight such that the sum of all weights in the dictionary = 1
	
	eventually this would maybe not be 100% random and have some physical basis
	"""
	spot_weights = np.random.random((len(spot_temperatures)))
	
	STAR_TEMPERATURE_WEIGHT : float = 5
	
	# add these in seperately - we might want to form the weight for the star differently
	spot_temperatures = np.append(spot_temperatures, star_temperature)
	spot_weights = np.append(spot_weights, STAR_TEMPERATURE_WEIGHT)
	
	# normalise
	
	spot_weights = spot_weights / np.sum(spot_weights)
	
	# bodge to make the debug graph for the input and found weights to have a line that doesn't go back on itself due to dictionaries being unsorted
	result = dict((zip(spot_temperatures, spot_weights)))
	result = {i[0]:i[1] for i in sorted(result.items(), key=lambda x: x[0])}
	return result

from astropy.constants import h, c

def normalise_flux(wavelength : np.array, counts : np.array) -> np.array:
	energy = h * c / wavelength
	
	total_energy = energy * counts
	
	# these units have no physical meaning, but a consistent between graphs, so I guess its okay?
	normalising_constant = np.sum(total_energy)
	
	total_energy /= normalising_constant
	
	# return counts again
	return (total_energy / energy).to_value()
	
	
	
# # # # # # # # # #

def get_composite_spectrum(star_T_eff_Kelvin : float,
						   delta_T_max_Kelvin : float,
						   spectra_grid_temperature_resolution_Kelvin : float,
						   FeH : float,
						   log_g : float):
	
	table = get_hdf5_data()
	
	table = table[(table[FEH_COLUMN] == FeH) & (table[LOGG_COLUMN] == log_g)]

	min_spot_temperature : float = star_T_eff_Kelvin - delta_T_max_Kelvin
	max_spot_temperature : float = star_T_eff_Kelvin + delta_T_max_Kelvin
	valid_spot_temperatures : np.array = np.arange(min_spot_temperature,  max_spot_temperature + spectra_grid_temperature_resolution_Kelvin, spectra_grid_temperature_resolution_Kelvin)
	# remove the star temperature itself; as that is already accounted for
	valid_spot_temperatures = valid_spot_temperatures[valid_spot_temperatures != star_T_eff_Kelvin]

	t_effs_and_weights = generate_weights(star_T_eff_Kelvin, valid_spot_temperatures)
	
	for temperature, weight in sorted(t_effs_and_weights.items()):
		print(f"temperature {temperature} : {weight}")
	
	plt.plot(t_effs_and_weights.keys(), t_effs_and_weights.values(), color="red")
	
	# this will be our total spectrum; combined from the star + all spots / faculae sources. keep columns and metadata from table; but no rows // data
	spectrum = table[:0]

	# get spectra for each temperature; sum them according to some weights
	for (temperature, weight) in t_effs_and_weights.items():
		sub_spectrum = table[table[TEFF_COLUMN] == temperature * u.K]

		if (len(sub_spectrum) == 0):
			# if the subset is empty / doesnt exist, could instead interpolate on the fly; might be a nicer method
			print(f"[FACULA AND SPOT CREATOR] : no spectrum data found for temperature {temperature}, FeH {FeH} and log_g {log_g}. Maybe check if the hdf5 contains the data you are trying to access. Skipping to next temperature value.")
			continue
		
		# normalise the sub spectrum somehow? rn this isn't physical; I assume we need to convert to energy per metre^2 or some analogy to that. this is just a placeholder. ig we want the property that the final energy contained in the spectrum corresponds to some fittable energy / the actual energy of the star ?
		normalising_constant = sp.integrate.simpson(sub_spectrum[FLUX_COLUMN], x=sub_spectrum[WAVELENGTH_COLUMN])
		# sub_spectrum[FLUX_COLUMN] /= normalising_constant ## including this breaks b in component_analyzer: aka we need to normalise everything or nothing to maintain consistency
		sub_spectrum[FLUX_COLUMN] = normalise_flux(sub_spectrum[WAVELENGTH_COLUMN], sub_spectrum[FLUX_COLUMN])
		sub_spectrum[FLUX_COLUMN] *= weight
		
		if len(spectrum) != 0:
			spectrum[FLUX_COLUMN] += sub_spectrum[FLUX_COLUMN]
		else:
			spectrum = sub_spectrum

	return spectrum

# just want these exposed to other python files for debugging, but there's no intrinsic need for them to be global / out of a helper function scope
FeH : float = 0.0
log_g : float = 4.0

def get_example_spectrum():
	# debug variables

	star_T_eff_Kelvin : float = 2650

	# max dT of facula, spots : inclusive of endpoints
	delta_T_max_Kelvin : float = 350

	spectra_grid_temperature_resolution_Kelvin : float = 50
	
	spectrum = get_composite_spectrum(star_T_eff_Kelvin, delta_T_max_Kelvin, spectra_grid_temperature_resolution_Kelvin, FeH, log_g)
	
	return spectrum
	
if __name__ == "__main__":
	spectrum = get_example_spectrum()
	
	plt.plot(spectrum[WAVELENGTH_COLUMN], spectrum[FLUX_COLUMN], linewidth=1)
	plt.show()

# for a given FeH; log_g; sum all the spectra according to w_i F_i where i \member of T_eff values
# and we want 1 of the w_i to be >> than the rest (as that will represent the star)
# the rest will be small and I guess randomised
# store FeH log_g {w_i's for each T_eff in some structured way}

# external
import astropy
from matplotlib import pyplot as plt
import numpy as np
import scipy as sp
from astropy import units as u
import specutils
from specutils.spectra import Spectrum

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
	
	STAR_TEMPERATURE_WEIGHT : float = 2
	
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

def old_normalise_flux(wavelength : np.array, counts : np.array) -> np.array:
	energy = h * c / wavelength
	
	total_energy = energy * counts
	
	# these units have no physical meaning, but a consistent between graphs, so I guess its okay?
	normalising_constant = np.sum(total_energy)
	
	total_energy /= normalising_constant
	
	# return counts again
	return (total_energy / energy).to_value()
 
from specutils.fitting import fit_generic_continuum, fit_continuum
from specutils import SpectralRegion
from astropy.visualization import quantity_support

def normalise_Janskys(wavelengths : np.array, counts : np.array, normalised_point = 2.2 * u.um, smoothing_range = 0.5 * u.um) -> np.array:
	"""
	this will fail if wavelengths does not span at least 0.5um
	
	the canonical way to normalise spectra is to choose a portion that's continuum and make that be at a consistent scale
	
	inputs should both be lists of astropy quantities (aka values with astropy units)
	
	I tried doing this but the continuum fit was terrible: https://specutils.readthedocs.io/en/stable/fitting.html#continuum-fitting
	"""
	
	# kernel size of about 501 with 9999 points between 5 and 15 um seemed good - this range corresponds (roughly) to that 
	wavelengths_in_range = wavelengths[(wavelengths[0] <= wavelengths) & (wavelengths <= wavelengths[0] + smoothing_range)]
	kernel_size = len(wavelengths_in_range)
	if kernel_size % 2 == 0:
		kernel_size +=1
		
	# smooth to make sure there's no spikes
	counts = np.array(counts, dtype=np.float64)
	smoothed_counts = sp.signal.medfilt(counts, kernel_size=[kernel_size])
	# normalise the counts at normalised_point (or next nearest value) to be 1
	counts /= smoothed_counts[(normalised_point <= wavelengths)][0]
	
	return counts
	
	
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
	
	# for temperature, weight in sorted(t_effs_and_weights.items()):
		# print(f"temperature {temperature} : {weight}")
	
	plt.plot(t_effs_and_weights.keys(), t_effs_and_weights.values(), color="red", linestyle="-", alpha=0.5, label="input weights", linewidth=4, marker="x", markersize=12)
	
	# this will be our total spectrum; combined from the star + all spots / faculae sources. keep columns and metadata from table; but no rows // data
	spectrum = table[:0]

	# get spectra for each temperature; sum them according to some weights
	for (temperature, weight) in t_effs_and_weights.items():
		sub_spectrum = table[table[TEFF_COLUMN] == temperature * u.K]

		if (len(sub_spectrum) == 0):
			# if the subset is empty / doesnt exist, could instead interpolate on the fly; might be a nicer method
			print(f"[FACULA AND SPOT CREATOR] : no spectrum data found for temperature {temperature}, FeH {FeH} and log_g {log_g}. Maybe check if the hdf5 contains the data you are trying to access. Skipping to next temperature value.")
			continue
		
		print(f"[FACULA AND SPOT CREATOR] : spectrum data found for temperature {temperature}, FeH {FeH} and log_g {log_g}.")
		
		# normalise the sub spectrum somehow? rn this isn't physical; I assume we need to convert to energy per metre^2 or some analogy to that. this is just a placeholder. ig we want the property that the final energy contained in the spectrum corresponds to some fittable energy / the actual energy of the star ?
		normalising_constant = sp.integrate.simpson(sub_spectrum[FLUX_COLUMN], x=sub_spectrum[WAVELENGTH_COLUMN])
		# sub_spectrum[FLUX_COLUMN] /= normalising_constant ## including this breaks b in component_analyzer: aka we need to normalise everything or nothing to maintain consistency
		sub_spectrum[FLUX_COLUMN] = normalise_Janskys(sub_spectrum[WAVELENGTH_COLUMN], sub_spectrum[FLUX_COLUMN])
		sub_spectrum[FLUX_COLUMN] *= weight
		
		if len(spectrum) != 0:
			spectrum[FLUX_COLUMN] += sub_spectrum[FLUX_COLUMN]
		else:
			spectrum = sub_spectrum
	return spectrum

# just want these exposed to other python files for debugging, but there's no intrinsic need for them to be global / out of a helper function scope
FeH : float = 0.0
log_g : float = 4.5

def get_example_spectrum():
	# debug variables

	star_T_eff_Kelvin : float = 3200

	# max dT of facula, spots : inclusive of endpoints
	delta_T_max_Kelvin : float = 500

	spectra_grid_temperature_resolution_Kelvin : float = 50
	
	spectrum = get_composite_spectrum(star_T_eff_Kelvin, delta_T_max_Kelvin, spectra_grid_temperature_resolution_Kelvin, FeH, log_g)
	return spectrum
	
if __name__ == "__main__":
	spectrum = get_example_spectrum()
	
	plt.plot(spectrum[WAVELENGTH_COLUMN], spectrum[FLUX_COLUMN], linewidth=1)
	plt.show()

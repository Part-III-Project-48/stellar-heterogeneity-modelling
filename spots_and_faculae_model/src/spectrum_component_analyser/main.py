# for now, just request a random composite spectrum from facula_and_spot_creator
# and try to decompose it - aka can we regenerate the w's

# eventually can read in external data or some training data from a large hdf5 file etc

from pathlib import Path
import astropy
from astropy.table import QTable
import numpy as np
from matplotlib import pyplot as plt
import scipy as sp
from astropy.visualization import quantity_support
from tqdm import tqdm
quantity_support()
import astropy.units as u
from scipy.interpolate import interp1d

from facula_and_spot_creator.main import normalise_counts, get_example_spectrum, FeH, log_g, old_normalise_flux
from phoenix_grid_creator.fits_to_hdf5 import FLUX_COLUMN, TEFF_COLUMN, FEH_COLUMN, LOGG_COLUMN, WAVELENGTH_COLUMN
from phoenix_grid_creator.basic_plotter import get_hdf5_data
from spectrum_component_analyser.external_spectrum_reader import get_external_spectra, read_JWST_fits

# # # # # helper functions # # # # #

from astropy.constants import h, c

def convert_counts_to_Janskys(wavelength : np.array, counts : np.array):
	"""
	wavelength should be an array of astropy quantities (with units of length)
	"""

	return (counts * [a.value**2 for a in wavelength]) * u.Jy
	energy = h * c / wavelength
	
	frequency = c / wavelength

	total_energy = energy * counts
	janskys = 10e-26 * (total_energy / u.second) / (u.m**2 * frequency)
	return janskys

# spectrum_to_decompose = get_example_spectrum()
# # add some random noise
# units = spectrum_to_decompose[FLUX_COLUMN][0].unit
# spectrum_to_decompose[FLUX_COLUMN] += (np.random.rand(len(spectrum_to_decompose[FLUX_COLUMN])) - 0.5) * units * (1/10) * np.average(spectrum_to_decompose[FLUX_COLUMN])
if __name__ == "__main__":
	external_spectrum_path = Path("spots_and_faculae_model/assets/MAST_2025-10-26T08_10_09.071Z/MAST_2025-10-26T08_10_09.071Z/JWST/jw02722003001_04101_00001-seg001_nis_x1dints.fits")

	absolute_external_spectrum = external_spectrum_path.resolve()

	spectrum_to_decompose = read_JWST_fits(absolute_external_spectrum)
	spectrum_to_decompose = spectrum_to_decompose[np.isfinite(spectrum_to_decompose[FLUX_COLUMN])]
	spectrum_to_decompose = spectrum_to_decompose[spectrum_to_decompose[FLUX_COLUMN] != np.nan]
	spectrum_to_decompose = spectrum_to_decompose[(.5 * u.um <= spectrum_to_decompose[WAVELENGTH_COLUMN])]
	# plt.plot(spectrum_to_decompose[WAVELENGTH_COLUMN], spectrum_to_decompose[FLUX_COLUMN])
	# plt.show()
	# read in all the available spectra we have (we are assuming FeH and log_g, so this is an easier problem for now)
	print("reading in hdf5")
	all_data : QTable = get_hdf5_data()
	print("finished reading in hdf5")
	# print(all_data)

	# we are effectively now carrying out the minimisation
	# ((\sum w_i f_i) - f_known) to find w_i for a ton of f_i's (our component spectra) and f_known being the total spectrum
	# this might be doable just with classical fitting but I guess there might be too much data ? we'll find out ...
	print("normalising counts")
	# spectrum_to_decompose[FLUX_COLUMN] *= 10e10

	plt.plot(spectrum_to_decompose[WAVELENGTH_COLUMN], spectrum_to_decompose[FLUX_COLUMN], label="unnormalised real spectrum to componentise")
	spectrum_to_decompose[FLUX_COLUMN] = normalise_counts(spectrum_to_decompose[WAVELENGTH_COLUMN], spectrum_to_decompose[FLUX_COLUMN])
	# spectrum_to_decompose[FLUX_COLUMN] = sp.signal.medfilt(spectrum_to_decompose[FLUX_COLUMN], kernel_size=[5])
	plt.plot(spectrum_to_decompose[WAVELENGTH_COLUMN], spectrum_to_decompose[FLUX_COLUMN], label="normalised real spectrum to componentise")
	plt.legend()
	plt.show()

	print("finished normalising counts")
	# print(spectrum_to_decompose)
	# plt.show()

	A = np.empty((0, 0))
	print("start filter")
	# print(all_data)
	all_data_subset = all_data[(all_data[FEH_COLUMN] == FeH) & (all_data[LOGG_COLUMN] == log_g)]
	print("finished filter")

	# this is slow but gives us all the T_effs saved to the spectral grid
	available_T_effs = astropy.table.unique(all_data_subset, keys=[TEFF_COLUMN])[TEFF_COLUMN] 

	# columns to be f_column number's and rows to increment the x value
	for T_eff in tqdm(available_T_effs, desc="appending fluxes to A matrix"):
		subset = all_data_subset[all_data_subset[TEFF_COLUMN] == T_eff]
		index_per_wavelength = len(subset[WAVELENGTH_COLUMN]) / (np.max(subset[WAVELENGTH_COLUMN]) - np.min(subset[WAVELENGTH_COLUMN]))
		convolution_range = int(index_per_wavelength * .001 * u.um) # should be something like the resolution of the jwst spectrum
		tqdm.write(f"convolution_range={convolution_range}")
		# plt.plot(subset[WAVELENGTH_COLUMN], subset[FLUX_COLUMN], label="unnormalised synthetic data")
		# plt.plot(spectrum_to_decompose[WAVELENGTH_COLUMN], spectrum_to_decompose[FLUX_COLUMN], label="normalised real spectrum to componentise")
		flux = sp.ndimage.gaussian_filter(subset[FLUX_COLUMN], convolution_range)
		interp = interp1d(subset[WAVELENGTH_COLUMN], flux, kind='linear')
		flux = interp(spectrum_to_decompose[WAVELENGTH_COLUMN].to(u.Angstrom))
		flux = convert_counts_to_Janskys(spectrum_to_decompose[WAVELENGTH_COLUMN], flux)
		# flux = sp.signal.medfilt(flux, kernel_size=[5])
		# plt.show()
		flux = normalise_counts(spectrum_to_decompose[WAVELENGTH_COLUMN], flux)
		plt.plot(spectrum_to_decompose[WAVELENGTH_COLUMN], flux)
		# plt.plot(spectrum_to_decompose[WAVELENGTH_COLUMN], flux, label="interpolated synthetic data onto the real spectral wavelengths")
		# plt.legend()
		# plt.show()
		if A.size == 0:
			A = np.empty((len(flux), 0))
		A = np.column_stack((A, flux))
	plt.show()
	print("minimising")
	# assume that w \in [0,1] : but I think this will only be true for real data if normalisation has been done correctly (???)
	result = sp.optimize.lsq_linear(A, spectrum_to_decompose[FLUX_COLUMN], bounds = (0, 100))#, tol=1e-15, lsmr_tol=1e-10, max_iter=2000)
	# alternative fortran method
	# result = sp.optimize.nnls(A, spectrum_to_decompose[FLUX_COLUMN])
	print(result)
	print(f"sum of weights={np.sum(result.x)}")
	# plt.plot(result[0])
	# plt.show()

	# print("found weights")
	# for i, T_eff in enumerate(available_T_effs):
	# 	weight = result.x[i]
	# 	print(f"temperature {T_eff} : {weight}")

	# # # plot some data # # #
	plt.plot(available_T_effs, result.x, color="blue", linestyle="dashed", marker="o", alpha=0.7, label="found weights")

	plt.xlabel("Temperature / K")
	plt.ylabel("weight / unitless")

	plt.legend()
	plt.show()

	plt.plot(spectrum_to_decompose[WAVELENGTH_COLUMN], spectrum_to_decompose[FLUX_COLUMN], label="experimental spectrum")
	plt.plot(spectrum_to_decompose[WAVELENGTH_COLUMN], A @ result.x, label="numerically found solution (sum of spectral components)")

	plt.legend()
	plt.show()
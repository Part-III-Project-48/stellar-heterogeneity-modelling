"""
this is going to download our spectrum data that meets some specified data ranges
this file should just be run once when you want to create the data grid, which is stored in the hdf5 format. Then you can load this hdf5 in quickly into a jupyter notebook etc to do the science. See basic_plotter.py for a way to read in the data
"""

# external imports
import numpy as np
from astropy import units as u
from pathlib import Path

# internal imports
from spectrum_component_analyser.internals.spectral_grid import spectral_grid
from spectrum_component_analyser.internals.spectrum import spectrum
from spectrum_component_analyser.internals.readers import read_JWST_fits, read_HARPS_fits, JWST_RESOLUTION, HARPS_resolution, JWST_NORMALISING_POINT, HARPS_normalising_point

SPECTRAL_GRID_FILENAME : Path = Path("test_JWST_not_oversmoothed.hdf5")

# data to request (these numbers have to be included in the PHOENIX dataset; view PHOENIX_filename_conventions.py for which are allowed)
T_effs = np.arange(2300, 4001, 100) * u.K
FeHs = np.array([-4, -3, -2, -1.5, -1, -0.5, 0, 0.5, 1])
log_gs = np.arange(0, 6.1, 0.5)

# test values
# T_effs = [2300] * u.K
# FeHs = np.array([0])

# # # flags # # #

REGULARISE_WAVELENGTH_GRID : bool = True
# the wavelength in the df starts out in angstroms (we add units to an astropy QTable later)
MIN_WAVELENGTH_ANGSTROMS : float = 0.5 * 10**(-6) * 10**(10)
MAX_WAVELENGTH_ANGSTROMS : float = 5.5 * 10**(-6) * 10**(10) # phoenix only goes up to 5.5?
WAVELENGTH_NUMBER_OF_POINTS : int = 5
regularised_wavelengths : np.array = np.linspace(MIN_WAVELENGTH_ANGSTROMS, MAX_WAVELENGTH_ANGSTROMS, WAVELENGTH_NUMBER_OF_POINTS) * u.Angstrom

# temperature interpolation
REGULARISE_TEMPERATURE_GRID : bool = False
MIN_TEMPERATURE_KELVIN = 2300
MAX_TEMPERATURE_KELVIN = 4500
TEMPERATURE_RESOLUTION_KELVIN = 50
regularised_temperatures = np.arange(MIN_TEMPERATURE_KELVIN, MAX_TEMPERATURE_KELVIN + TEMPERATURE_RESOLUTION_KELVIN, TEMPERATURE_RESOLUTION_KELVIN)


# unsupported atm
# set to np.inf to ignore
# DEBUG_MAX_NUMBER_OF_SPECTRA_TO_DOWNLOAD : int = np.inf
# if REGULARISE_TEMPERATURE_GRID:
# 	grid.regularise_temperatures(regularised_temperatures)

if __name__ == "__main__":

	# load in a spectrum and use that as the regularised wavelengths
	external_spectrum_path = Path("../../observed_spectra/MAST_2025-10-26T08_10_09.071Z - K2-18/MAST_2025-10-26T08_10_09.071Z/JWST/jw02722003001_04101_00001-seg001_nis_x1dints.fits")
	script_dir = Path(__file__).resolve().parent
	wavelength_grid_absolute_path = (script_dir / external_spectrum_path).resolve()

	spectrum_to_decompose : spectrum = read_JWST_fits(wavelength_grid_absolute_path)

	spec_grid : spectral_grid = spectral_grid.from_internet(T_effs=T_effs,
														 FeHs=FeHs,
														 log_gs=log_gs,
														 normalising_point=JWST_NORMALISING_POINT,
														 observational_resolution=JWST_RESOLUTION,
														 observational_wavelengths=spectrum_to_decompose.Wavelengths,
														 name="phoenix_data")

	spec_grid.save(absolute_path=SPECTRAL_GRID_FILENAME, overwrite=True)

	test_read : spectral_grid = spectral_grid.from_hdf5(SPECTRAL_GRID_FILENAME)

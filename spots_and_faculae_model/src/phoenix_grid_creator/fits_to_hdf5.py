"""
this is going to download our spectrum data that meets some specified data ranges
this file should just be run once when you want to create the data grid, which is stored in the hdf5 format. Then you can load this hdf5 in quickly into a jupyter notebook etc to do the science. See basic_plotter.py for a way to read in the data
"""

# external imports
from io import BytesIO
from joblib import Parallel, delayed
import pandas as pd
import requests
from tqdm import tqdm
import numpy as np
import astropy
from astropy.io import fits
from astropy import units as u
from astropy.table import Table
from itertools import product
from pathlib import Path
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import datetime
from astropy.visualization import quantity_support
from astropy.table import QTable
import scipy as sp
from astropy.table import Table, vstack
import h5py
from typing import Sequence

# internal imports
from phoenix_grid_creator.PHOENIX_filename_conventions import *
from spots_and_faculae_model.spectrum_grid import spectrum_grid
from spots_and_faculae_model.phoenix_spectrum import phoenix_spectrum
from spots_and_faculae_model.simpler_spectral_grid import simpler_spectral_grid

SAVE_TO_HDF : bool = False
SPECTRAL_GRID_FILENAME : str = 'spectral_grid.hdf5'

# data to request (these numbers have to be included in the PHOENIX dataset; view PHOENIX_filename_conventions.py for which are allowed)
T_effs = np.arange(2300, 4001, 100) * u.K
FeHs = np.array([-4, -3, -2, -1.5, -1, -0.5, 0, 0.5, 1])
log_gs = np.arange(0, 6.1, 0.5)
alphaM = 0

T_effs = np.array([2500, 2700]) * u.K
FeHs = [-4, -3, -2]
log_gs = [0]

PHOENIX_FLUX_UNITS = u.erg / (u.s * u.cm**2 * u.cm)

lte : bool = True
# # # flags # # #

REGULARISE_WAVELENGTH_GRID : bool = True
# the wavelength in the df starts out in angstroms (we add units to an astropy QTable later)
MIN_WAVELENGTH_ANGSTROMS : float = 0.5 * 10**(-6) * 10**(10)
MAX_WAVELENGTH_ANGSTROMS : float = 5.5 * 10**(-6) * 10**(10) # phoenix only goes up to 5.5?
WAVELENGTH_NUMBER_OF_POINTS : int = 5
regularised_wavelengths : np.array = np.linspace(MIN_WAVELENGTH_ANGSTROMS, MAX_WAVELENGTH_ANGSTROMS, WAVELENGTH_NUMBER_OF_POINTS) * u.Angstrom

# # ipynb files complain about this otherwise
# if __name__ == "__main__":
# 	# get the wavelengths from a jwst file and use that to convolve against
# 	path = Path("spots_and_faculae_model/assets/MAST_2025-10-26T08_10_09.071Z/MAST_2025-10-26T08_10_09.071Z/JWST/jw02722003001_04101_00001-seg001_nis_x1dints.fits")
# 	# regularised_wavelengths = read_JWST_fits(path).Wavelengths

# temperature interpolation
REGULARISE_TEMPERATURE_GRID : bool = False
MIN_TEMPERATURE_KELVIN = 2300
MAX_TEMPERATURE_KELVIN = 4500
TEMPERATURE_RESOLUTION_KELVIN = 50
regularised_temperatures = np.arange(MIN_TEMPERATURE_KELVIN, MAX_TEMPERATURE_KELVIN + TEMPERATURE_RESOLUTION_KELVIN, TEMPERATURE_RESOLUTION_KELVIN)

CONVERT_WAVELENGTHS_TO_AIR : bool = False

# set to np.inf to ignore
DEBUG_MAX_NUMBER_OF_SPECTRA_TO_DOWNLOAD : int = np.inf

# # # # # #

# # # # # helper functions # # # # # 

if __name__ == "__main__":

	# spec_grid : simpler_spectral_grid = simpler_spectral_grid.from_internet(T_effs, FeHs, log_gs, regularised_wavelengths=regularised_wavelengths)

	path : Path = "test2.hdf5"
	# spec_grid.save(absolute_path=Path(path), overwrite=True)
	
	# if REGULARISE_TEMPERATURE_GRID:
	# 	grid.regularise_temperatures(regularised_temperatures)
	test_read : simpler_spectral_grid = simpler_spectral_grid.from_hdf5(path)

	print(test_read)




# if __name__ == "__main__":
#     import cProfile, pstats
#     with cProfile.Profile() as pr:
#         main()   # or whatever entry point you want to run
#     stats = pstats.Stats(pr)
#     # stats.sort_stats("tottime").print_stats(20)

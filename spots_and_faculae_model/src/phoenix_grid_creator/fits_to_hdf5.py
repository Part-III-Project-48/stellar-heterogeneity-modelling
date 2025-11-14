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
import specutils
from scipy.interpolate import interp1d
import datetime
from astropy.visualization import quantity_support
from astropy.table import QTable
import scipy as sp
from astropy.table import Table, vstack

# internal imports
from phoenix_grid_creator.PHOENIX_filename_conventions import *
from spots_and_faculae_model.spectrum_grid import spectrum_grid
from spots_and_faculae_model.external_spectrum_reader import read_JWST_fits

# column names


SAVE_TO_HDF : bool = False
SPECTRAL_GRID_FILENAME : str = 'spectral_grid.hdf5'

# data to request (these numbers have to be included in the PHOENIX dataset; view PHOENIX_filename_conventions.py for which are allowed)
T_effs = np.arange(2300, 4001, 100) * u.K
FeHs = np.array([-4, -3, -2, -1.5, -1, -0.5, 0, 0.5, 1])
log_gs = np.arange(0, 6.1, 0.5)
alphaM = 0

T_effs = [2500 * u.K]

# seems like the only data I can find is LTE data (?)
lte : bool = True
# # # flags # # #

REGULARISE_WAVELENGTH_GRID : bool = True
# the wavelength in the df starts out in angstroms (we add units to an astropy QTable later)
MIN_WAVELENGTH_ANGSTROMS : float = 0.5 * 10**(-6) * 10**(10)
MAX_WAVELENGTH_ANGSTROMS : float = 5.5 * 10**(-6) * 10**(10) # phoenix only goes up to 5.5?
WAVELENGTH_NUMBER_OF_POINTS : int = 5
regularised_wavelengths : np.array = np.linspace(MIN_WAVELENGTH_ANGSTROMS, MAX_WAVELENGTH_ANGSTROMS, WAVELENGTH_NUMBER_OF_POINTS) * u.Angstrom

# ipynb files complain about this otherwise
if __name__ == "__main__":
	# get the wavelengths from a jwst file and use that to convolve against
	path = Path("spots_and_faculae_model/assets/MAST_2025-10-26T08_10_09.071Z/MAST_2025-10-26T08_10_09.071Z/JWST/jw02722003001_04101_00001-seg001_nis_x1dints.fits")
	# regularised_wavelengths = read_JWST_fits(path).Wavelengths

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
	
def get_wavelength_grid() -> np.ndarray:
	"""
	returns 1D array consisting of astropy Quantities of dimension length
	"""
	# read in the wavelength (1D) grid so we can save this into our mega-grid correctly

	script_dir = Path(__file__).resolve().parent
	WAVELENGTH_GRID_RELATIVE_PATH = Path("../../assets/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits")
	wavelength_grid_absolute_path = (script_dir / WAVELENGTH_GRID_RELATIVE_PATH).resolve()

	if not wavelength_grid_absolute_path.exists():
		raise FileNotFoundError(f"wavelength grid file not found at : {wavelength_grid_absolute_path}.")

	with fits.open(wavelength_grid_absolute_path) as hdul:
		# there is only 1 HDU in this wavelength grid file
		WAVELENGTH_GRID_HDU_INDEX = 0
		wavelengths = hdul[WAVELENGTH_GRID_HDU_INDEX].data
		# the fits file is big endian; pandas requires little endian. this swaps between them
		wavelengths = wavelengths.byteswap().view(wavelengths.dtype.newbyteorder())
		wavelengths *= u.Angstrom
		print("[PHOENIX GRID CREATOR] : wavelength grid found & loaded in")

	return wavelengths

def download_spectrum(T_eff, FeH, log_g, wavelengths : np.array) -> spectrum_grid:
	"""
	well use this function to parallelise getting the spectra
	"""
	try:
		file = get_file_name(lte, T_eff, log_g, FeH, alphaM)
		url = get_url(file)
	except ValueError as e:
		tqdm.write(f"[PHOENIX GRID CREATOR] : filename or urlname error: {e}. continuing onto next requested spectrum")
		return
	try:
		response = requests.get(url)
		response.raise_for_status()
	except requests.exceptions.HTTPError as e:
		tqdm.write(f"[PHOENIX GRID CREATOR] : HTTPError raised with the following parameters.\nlte: {lte}\nT_eff={T_eff}\nlog_g={log_g}\nFeH={FeH}\nalphaM={alphaM}")
		tqdm.write(f"url = {url}")
		tqdm.write("\n continuing with the next file...")
		return
	
	# the index of the header data unit the data we want is in (looks to be 0 being the spectra, and 1 being the abundances, and those are the only 2 HDUs in the .fits files)
	SPECTRA_HDU_INDEX = 0

	with fits.open(BytesIO(response.content)) as hdul:
		# hdul.info()
		
		fluxes = hdul[SPECTRA_HDU_INDEX].data

		# for some reason, the fits file is big-endian; pandas required little-endian
		fluxes = fluxes.byteswap().view(fluxes.dtype.newbyteorder())
		fluxes *= u.erg / (u.s * u.cm**2 * u.cm)
		JWST_resolution = .001 * u.um

		if REGULARISE_WAVELENGTH_GRID:
			index_per_wavelength = len(wavelengths) / (np.max(wavelengths) - np.min(wavelengths))
			index_per_wavelength = index_per_wavelength.to(u.um**-1)
			convolution_range = int((index_per_wavelength * JWST_resolution)) # in number of adjacent points to consider
			fluxes = sp.ndimage.gaussian_filter(fluxes, convolution_range) * fluxes.unit # gaussian_filter seems to remove units
			t = spectrum_grid.from_arrays(
				T_effs=[T_eff]*len(regularised_wavelengths),
				FeHs=[FeH]*len(regularised_wavelengths),
				log_gs=[log_g]*len(regularised_wavelengths),
				wavelengths=regularised_wavelengths,
				fluxes=np.interp(regularised_wavelengths, wavelengths, fluxes)
				)
		else:
			t = spectrum_grid.from_arrays(
					T_effs=[T_eff]*len(wavelengths),
					FeHs=[FeH]*len(wavelengths),
					log_gs=[log_g]*len(wavelengths),
					wavelengths=wavelengths,
					fluxes=fluxes)
		# this might be quicker to stream data to disc rather than creating a massive df
		# temp_df.write(HDF5_FILENAME_TO_SAVE, path = "data", serialize_meta=True, overwrite=True, append=True)
		# continue

		return t

def main():

	wavelengths = get_wavelength_grid()

	# now use defined ranges for the data we want, process it and save this to a hdf5 fil
	
	grids = Parallel(n_jobs=-1, prefer="threads")(
		delayed(download_spectrum)(T_eff, FeH, log_g, wavelengths) for T_eff, FeH, log_g in tqdm(product(T_effs, FeHs, log_gs), total=len(T_effs) * len(FeHs) * len(log_gs), desc="Downloading .fits spectra files")
		)
	
	grid = spectrum_grid(vstack([grid.Table for grid in grids]))
	
	if REGULARISE_TEMPERATURE_GRID:
		grid.regularise_temperatures(regularised_temperatures)

	if CONVERT_WAVELENGTHS_TO_AIR:
		grid.convert_vacuum_to_air()

	# add some metadata to the QTable e.g. (wavelength medium = air, source, date)
	grid.Table.meta = {"wavelength medium" : "air" if CONVERT_WAVELENGTHS_TO_AIR else "vacuum",
				"source" : "https://phoenix.astro.physik.uni-goettingen.de/data/",
				"date this hdf5 file was created" : datetime.datetime.now(),
				"description" : "if it includes interpolated values, then a specified list of wavelengths and/or temperatures were given, and the simulated data was (linearly) interpolated onto those values (aka; information was removed and accuracy is not guaranteed)",
				"includes interpolated wavelengths?" : REGULARISE_WAVELENGTH_GRID,
				"regularised wavelengths (Angstroms):" : f"np.linspace({MIN_WAVELENGTH_ANGSTROMS}, {MAX_WAVELENGTH_ANGSTROMS}, {WAVELENGTH_NUMBER_OF_POINTS})" if REGULARISE_WAVELENGTH_GRID else "not applicable",
				"includes interpolated temperatures?" : REGULARISE_TEMPERATURE_GRID,
				"regularised temperatures (Kelvin)" : f"np.arange({MIN_TEMPERATURE_KELVIN}, {MAX_TEMPERATURE_KELVIN + TEMPERATURE_RESOLUTION_KELVIN}, {TEMPERATURE_RESOLUTION_KELVIN})" if REGULARISE_TEMPERATURE_GRID else "not applicable",
				"Teff (original data)" : T_effs,
				"FeH (original data)" : FeHs,
				"log_gs (original data)" : log_gs}
	
	if SAVE_TO_HDF:
		grid.save(absolute_path="data", name=SPECTRAL_GRID_FILENAME)



if __name__ == "__main__":
    import cProfile, pstats
    with cProfile.Profile() as pr:
        main()   # or whatever entry point you want to run
    stats = pstats.Stats(pr)
    # stats.sort_stats("tottime").print_stats(20)

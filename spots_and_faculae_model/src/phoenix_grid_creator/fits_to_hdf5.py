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
# temp
def read_JWST_fits(fits_absolute_path : Path) -> QTable:
	if not fits_absolute_path.exists():
		raise FileNotFoundError(f"JWST spectrum file not found at : {fits_absolute_path}.")

	with fits.open(fits_absolute_path) as hdul:
		hdul.info()
		HDU_INDEX = 3 # aka EXTRACT1D
		data = hdul[HDU_INDEX].data
		hdr = hdul[HDU_INDEX].header
		print(repr(hdr))
		INTEGRATION_INDEX : int = 100
		# print("[SPECTRUM COMPONENT ANALYSER] : external spectrum found & loaded in")
		spectrum = QTable()
		spectrum[WAVELENGTH_COLUMN] = data["WAVELENGTH"][INTEGRATION_INDEX] * u.um
		spectrum[FLUX_COLUMN] = data["FLUX"][INTEGRATION_INDEX] * u.Jy
	
	return spectrum

# column names
TEFF_COLUMN = "T_eff / K"
FEH_COLUMN = "Fe/H / relative to solar"
LOGG_COLUMN = "log_g / log(cm s^(-2))"
WAVELENGTH_COLUMN = "wavelength / angstroms"
FLUX_COLUMN = "flux / erg / (s * cm**2 * cm)"

SAVE_TO_HDF : bool = True
HDF5_FILENAME_TO_SAVE : str = 'spectral_grid.hdf5'

# data to request (these numbers have to be included in the PHOENIX dataset; view PHOENIX_filename_conventions.py for which are allowed)
T_effs = np.arange(2300, 4001, 100) * u.K
FeHs = np.array([-4, -3, -2, -1.5, -1, -0.5, 0, 0.5, 1])
log_gs = np.arange(0, 6.1, 0.5)
alphaM = 0

# seems like the only data I can find is LTE data (?)
lte : bool = True

#debug override for testing - this is the data we collect; the data we save is specified below (if interpolating // regularising)
# T_effs = np.array([3500, 3600]) * u.K
# FeHs = np.array([0, 1])
# log_gs = np.array([4.5, 5])

# # # flags # # #

REGULARISE_WAVELENGTH_GRID : bool = True
# the wavelength in the df starts out in angstroms (we add units to an astropy QTable later)
MIN_WAVELENGTH_ANGSTROMS : float = 0.5 * 10**(-6) * 10**(10)
MAX_WAVELENGTH_ANGSTROMS : float = 5.5 * 10**(-6) * 10**(10) # phoenix only goes up to 5.5?
WAVELENGTH_NUMBER_OF_POINTS : int = 10_000
regularised_wavelengths = np.linspace(MIN_WAVELENGTH_ANGSTROMS, MAX_WAVELENGTH_ANGSTROMS, WAVELENGTH_NUMBER_OF_POINTS)

# ipynb files complain about this otherwise
if __name__ == "__main__":
	path = Path("spots_and_faculae_model/assets/MAST_2025-10-26T08_10_09.071Z/MAST_2025-10-26T08_10_09.071Z/JWST/jw02722003001_04101_00001-seg001_nis_x1dints.fits")
	regularised_wavelengths = read_JWST_fits(path)[WAVELENGTH_COLUMN]

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

def debug_plot_interpolated_temperatures():
	DEBUG_MIN_GRAPH_TEMPERATURE = 2300
	DEBUG_MAX_GRAPH_TEMPERATURE = 3000
	
	DEBUG_MIN_GRAPH_WAVELENGTH = 5000
	DEBUG_MAX_GRAPH_WAVELENGTH = 6000
	
	plt.figure(figsize=(10, 5))
	
	example_FeH = FeHs[0]
	example_log_g = log_gs[0]
	
	# new (interpolated) T_effs
	for T_eff in regularised_temperatures[(DEBUG_MIN_GRAPH_TEMPERATURE <= regularised_temperatures) & (regularised_temperatures <= DEBUG_MAX_GRAPH_TEMPERATURE)]:
		subset_df = regularised_wavelength_df[(regularised_wavelength_df[TEFF_COLUMN] == T_eff) & 
												(regularised_wavelength_df[FEH_COLUMN] == example_FeH) &
												(regularised_wavelength_df[LOGG_COLUMN] == example_log_g)]
		
		subset_df = subset_df[(DEBUG_MIN_GRAPH_WAVELENGTH <= subset_df[WAVELENGTH_COLUMN]) & (subset_df[WAVELENGTH_COLUMN] <= DEBUG_MAX_GRAPH_WAVELENGTH)]
		
		plt.plot(subset_df[WAVELENGTH_COLUMN], subset_df[FLUX_COLUMN], linestyle="dashed", label=f"(interpolated) {T_eff}K")

	# old T_effs
	for T_eff in T_effs[(DEBUG_MIN_GRAPH_TEMPERATURE <= T_effs) & (T_effs <= DEBUG_MAX_GRAPH_TEMPERATURE)]:
		subset_df = df[(df[TEFF_COLUMN] == T_eff) &
						(df[FEH_COLUMN] == example_FeH) &
						(df[LOGG_COLUMN] == example_log_g)]
		
		subset_df = subset_df[(DEBUG_MIN_GRAPH_WAVELENGTH <= subset_df[WAVELENGTH_COLUMN]) & (subset_df[WAVELENGTH_COLUMN] <= DEBUG_MAX_GRAPH_WAVELENGTH)]
		
		plt.plot(subset_df[WAVELENGTH_COLUMN], subset_df[FLUX_COLUMN], label=f"(actual) {T_eff}K")
	
	
	plt.title(f"subset of interpolated data\nat [Fe/H] = {example_FeH} and log_g = {example_log_g}")
	plt.xlabel("Wavelength / Angstroms")
	plt.ylabel("Flux / counts")
	plt.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
	# plt.tight_layout()
	plt.subplots_adjust(right=0.75)
	
	plt.show()
	
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
	

# # # # # # # # # #

# just want this exposed to other python files for debugging (temporarily)
wavelengths = get_wavelength_grid()

main_table = QTable()

if __name__ == "__main__":

	# now use defined ranges for the data we want, process it and save this to a hdf5 file

	def download_spectrum(T_eff, FeH, log_g):		
		
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
				t = QTable([
					[T_eff]*len(regularised_wavelengths),
					[FeH]*len(regularised_wavelengths),
					[log_g]*len(regularised_wavelengths),
					regularised_wavelengths,
					np.interp(regularised_wavelengths, wavelengths, fluxes)
					],
				names=(TEFF_COLUMN, FEH_COLUMN, LOGG_COLUMN, WAVELENGTH_COLUMN, FLUX_COLUMN))
			else:
				t = QTable([
						[T_eff]*len(wavelengths),
						[FeH]*len(wavelengths),
						[log_g]*len(wavelengths),
						wavelengths,
						fluxes],
					names=(TEFF_COLUMN, FEH_COLUMN, LOGG_COLUMN, WAVELENGTH_COLUMN, FLUX_COLUMN))
			# this might be quicker to stream data to disc rather than creating a massive df
			# temp_df.write(HDF5_FILENAME_TO_SAVE, path = "data", serialize_meta=True, overwrite=True, append=True)
			# continue

			return t

	tables = Parallel(n_jobs=-1, prefer="threads")(
		delayed(download_spectrum)(T_eff, FeH, log_g) for T_eff, FeH, log_g in tqdm(product(T_effs, FeHs, log_gs), total=len(T_effs) * len(FeHs) * len(log_gs), desc="Downloading .fits spectra files")
		)
	
	main_table = vstack(tables)
	
	if REGULARISE_TEMPERATURE_GRID:
		
		regularised_wavelength_df = pd.DataFrame(columns=[TEFF_COLUMN, FEH_COLUMN, LOGG_COLUMN, WAVELENGTH_COLUMN, FLUX_COLUMN])
		
		for FeH, log_g, new_T_eff in tqdm(product(FeHs, log_gs, regularised_temperatures), total= len(FeHs) * len(log_gs) * len(regularised_temperatures), desc="Regularising temperature points"):
			
			# this df is at all temperatures
			subset_df = df[(df[FEH_COLUMN] == FeH) & (df[LOGG_COLUMN] == log_g)]

			# aka subset df represents 3D data which maps (wavelength to flux) over a set of temperatures
			# we want the wavelength to flux map at a different temperature; namely at new_T_eff
			
			# so we want to linearly interpolate every flux between T_1 and T_2 at all wavelengths
			# aka. np.interp(new_T_eff, subset_df[TEMPERATURE_COLUMN], subset_df[FLUX_COLUMN]) # assuming this vectorises and returns a map from all wavelengths to fluxes
			
			pivoted = subset_df.pivot(index=TEFF_COLUMN, columns=WAVELENGTH_COLUMN, values=FLUX_COLUMN)
			x = pivoted.index.to_numpy()               # shape (n_temperatures,)
			wavelengths = pivoted.columns.to_numpy()   # shape (n_wavelengths,)
			y = pivoted.values
			
			f = interp1d(x, y, axis=0, kind="cubic")
			
			wavelength_to_flux_map_at_new_T_eff = f(new_T_eff)
			
			temp_df = pd.DataFrame({
				TEFF_COLUMN : new_T_eff,
				FEH_COLUMN : FeH,
				LOGG_COLUMN : log_g,
				WAVELENGTH_COLUMN : wavelengths,
				FLUX_COLUMN : wavelength_to_flux_map_at_new_T_eff # interpolated flux function between previous and next temperatures 
			})
			
			# avoid warning about concat-ing an empty df
			if not regularised_wavelength_df.empty:
				# our df index has no meaningful meaning, and sort I think just ensures the columns are in the correct order or something?
				regularised_wavelength_df = pd.concat([regularised_wavelength_df, temp_df], ignore_index=True)#, sort=True)
			else:
				regularised_wavelength_df = temp_df
		
		### debug plotting to double check the interpolation worked ###

		# just for plotting
		debug_plot_interpolated_temperatures()
		
		### end of debugging plotting ###
		
		# this must be set after the debug plotting, as the debug plotting function depends on df and regularised_wavelength_df being distinct
		df = regularised_wavelength_df
	
	# pandas tables can't save their metadata into a HDF5 directly (can use HDFStore or smthn) - but astropy tables can have metadata, units etc. so lets convert to an astropy table

	# add astropy units to columns (this will be stored in metadata and can be read back out into an astropy QTable)
	# main_table[TEFF_COLUMN].unit = u.Kelvin
	main_table[TEFF_COLUMN].desc = "effective surface temperature"

	# main_table[FEH_COLUMN].unit = u.dimensionless_unscaled
	main_table[FEH_COLUMN].desc = "relative to solar metallacity"

	# astropy seems to have a hard time reading in log quantities from hdf5 files. so lets just save this as unitless
	# main_table[LOGG_COLUMN].unit = u.dimensionless_unscaled
	main_table[LOGG_COLUMN].desc = "log_10(u.cm * u.second**(-2)) of the surface gravity"

	# main_table[WAVELENGTH_COLUMN].unit = u.Angstrom

	# main_table[FLUX_COLUMN].unit = u.dimensionless_unscaled
	main_table[FLUX_COLUMN].desc = "in counts"

	if CONVERT_WAVELENGTHS_TO_AIR:
		t[WAVELENGTH_COLUMN] = specutils.utils.wcs_utils.vac_to_air(t[WAVELENGTH_COLUMN])

	# add some metadata to the QTable e.g. (wavelength medium = air, source, date)
	main_table.meta = {"wavelength medium" : "air" if CONVERT_WAVELENGTHS_TO_AIR else "vacuum",
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
		print("[PHOENIX GRID CREATOR] : writing dataframe to hdf5...")
		main_table.write(HDF5_FILENAME_TO_SAVE, path = "data", serialize_meta=True, overwrite=True)
		print("[PHOENIX GRID CREATOR] : hdf5 saving complete")





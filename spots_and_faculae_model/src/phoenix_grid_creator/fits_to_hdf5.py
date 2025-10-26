"""
this is going to download our spectrum data that meets some specified data ranges
this file should just be run once when you want to create the data grid, which is stored in the hdf5 format. Then you can load this hdf5 in quickly into a jupyter notebook etc to do the science. See basic_plotter.py for a way to read in the data
"""

# external imports
from io import BytesIO
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
import scipy
import matplotlib.pyplot as plt
import specutils
from scipy.interpolate import interp1d
import datetime
from astropy.visualization import quantity_support
from astropy.table import QTable

# internal imports
from phoenix_grid_creator.PHOENIX_filename_conventions import *

# column names
TEFF_COLUMN = "T_eff / K"
FEH_COLUMN = "Fe/H / relative to solar"
LOGG_COLUMN = "log_g / log(cm s^(-2))"
WAVELENGTH_COLUMN = "wavelength / angstroms"
FLUX_COLUMN = "flux / counts"

SAVE_TO_HDF : bool = True
HDF5_FILENAME_TO_SAVE : str = 'spectral_grid.hdf5'

# data to request (these numbers have to be included in the PHOENIX dataset; view PHOENIX_filename_conventions.py for which are allowed)
T_effs = np.arange(2300, 7001, 100)
FeHs = np.arange(-0.5, 0.6, 0.5)
log_gs = np.arange(3.5, 5.6, 0.5)
alphaM = 0

# seems like the only data I can find is LTE data (?)
lte : bool = True

#debug override for testing - this is the data we collect; the data we save is specified below (if interpolating // regularising)
# T_effs = np.array([3500])
FeHs = np.array([0])
log_gs = np.array([4.5])

# # # flags # # #

REGULARISE_WAVELENGTH_GRID : bool = False
# the wavelength in the df starts out in angstroms (we add units to an astropy QTable later)
MIN_WAVELENGTH_ANGSTROMS : float = 0.5 * 10**(-6) * 10**(10)
MAX_WAVELENGTH_ANGSTROMS : float = 5.5 * 10**(-6) * 10**(10) # phoenix only goes up to 5.5?
WAVELENGTH_NUMBER_OF_POINTS : int = 10_000
regularised_wavelengths = np.linspace(MIN_WAVELENGTH_ANGSTROMS, MAX_WAVELENGTH_ANGSTROMS, WAVELENGTH_NUMBER_OF_POINTS)

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
print(1 % 0.2)
plt.show()



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
	
def get_wavelength_grid():
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
		
		print("[PHOENIX GRID CREATOR] : wavelength grid found & loaded in")
	
	return wavelengths
	

# # # # # # # # # #

# just want this exposed to other python files for debugging (temporarily)
wavelengths = get_wavelength_grid()

if __name__ == "__main__":

	# now use defined ranges for the data we want, process it and save this to a hdf5 file

	total_number_of_files : int = len(T_effs) * len(FeHs) * len(log_gs)

	# we will save our grid to this df
	df = pd.DataFrame(columns=[TEFF_COLUMN, FEH_COLUMN, LOGG_COLUMN, WAVELENGTH_COLUMN, FLUX_COLUMN])

	file_number = 0

	for T_eff, FeH, log_g in tqdm(product(T_effs, FeHs, log_gs), total=total_number_of_files, desc="Downloading .fits spectra files"):
		
		if file_number >= DEBUG_MAX_NUMBER_OF_SPECTRA_TO_DOWNLOAD:
			break
		
		file_number += 1
		
		file = get_file_name(lte, T_eff, log_g, FeH, alphaM)
		url = get_url(file)
		
		try:
			response = requests.get(url)
			response.raise_for_status()
		except requests.exceptions.HTTPError as e:
			tqdm.write(f"[PHOENIX GRID CREATOR] : HTTPError raised with the following parameters.\nlte: {lte}\nT_eff={T_eff}\nlog_g={log_g}\nFeH={FeH}\nalphaM={alphaM}")
			tqdm.write(f"url = {url}")
			tqdm.write("\n continuing with the next file...")
			continue
		
		
		# if you want to write the .fits file somewhere then use this
		# temp_file_name : str = "example.fits"
		# with open(temp_file_name, "wb") as f:
		# 	f.write(response.content)
		
		# the index of the header data unit the data we want is in (looks to be 0 being the spectra, and 1 being the abundances, and those are the only 2 HDUs in the .fits files)
		SPECTRA_HDU_INDEX = 0

		with fits.open(BytesIO(response.content)) as hdul:
			# hdul.info()
			
			fluxes = hdul[SPECTRA_HDU_INDEX].data

			# min = -np.inf * u.um
			# max = np.inf * u.um
			# subset = wavelengths[(min / u.Angstrom <= wavelengths) & (wavelengths <= max / u.Angstrom)]
			# print(len(subset))
			# print(len(regularised_wavelengths))

			# for some reason, the fits file is big-endian; pandas required little-endian
			fluxes = fluxes.byteswap().view(fluxes.dtype.newbyteorder())
			# pandas will repeat the constant values len(fluxes) times for us
			temp_df = pd.DataFrame({
				TEFF_COLUMN : T_eff,
				FEH_COLUMN : FeH,
				LOGG_COLUMN : log_g,
				# need to use the row index <-> wavelength map provided to us by PHOENIX
				WAVELENGTH_COLUMN : wavelengths,
				FLUX_COLUMN : fluxes
			})

			# quantity_support()
			# fig, ax = plt.subplots()
			# t = QTable.from_pandas(temp_df)
			# t[FLUX_COLUMN] *= t[WAVELENGTH_COLUMN]**2
			# t[WAVELENGTH_COLUMN] *= u.Angstrom
			# ax.plot(t[WAVELENGTH_COLUMN], t[FLUX_COLUMN])
			# ax.xaxis.set_units(u.um)
			# ax.set_xlim(0.6 * u.um, 2.8 * u.um)
			# plt.show()
			
			if REGULARISE_WAVELENGTH_GRID:
				# linear interpolate the fluxes onto a specified grid (this can easily be changed for cubic splines etc if needed)
				temp_df = temp_df[(MIN_WAVELENGTH_ANGSTROMS <= temp_df[WAVELENGTH_COLUMN]) & (temp_df[WAVELENGTH_COLUMN] <= MAX_WAVELENGTH_ANGSTROMS)]
				
				temp_df = pd.DataFrame({
					TEFF_COLUMN : T_eff,
					FEH_COLUMN : FeH,
					LOGG_COLUMN : log_g,
					# need to use the row index <-> wavelength map provided to us by PHOENIX
					WAVELENGTH_COLUMN : regularised_wavelengths,
					FLUX_COLUMN : np.interp(regularised_wavelengths, temp_df[WAVELENGTH_COLUMN], temp_df[FLUX_COLUMN])
				})

				plt.plot(regularised_wavelengths, temp_df[FLUX_COLUMN])
			
			# this might be quicker to stream data to disc rather than creating a massive df
			# temp_df.write(HDF5_FILENAME_TO_SAVE, path = "data", serialize_meta=True, overwrite=True, append=True)
			# continue
			
			# avoid warning about concat-ing an empty df
			if not df.empty:
				# our df index has no meaningful meaning, and sort I think just ensures the columns are in the correct order or something?
				df = pd.concat([df, temp_df], ignore_index=True)#, sort=True)
			else:
				df = temp_df
	
	plt.show()
	
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
	table = QTable.from_pandas(df)

	# add astropy units to columns (this will be stored in metadata and can be read back out into an astropy QTable)
	from astropy.units import imperial
	imperial.enable()
	table[TEFF_COLUMN].unit = u.Kelvin
	table[TEFF_COLUMN].desc = "effective surface temperature"

	table[FEH_COLUMN].unit = u.dimensionless_unscaled
	table[FEH_COLUMN].desc = "relative to solar metallacity"

	# astropy seems to have a hard time reading in log quantities from hdf5 files. so lets just save this as unitless
	table[LOGG_COLUMN].unit = u.dimensionless_unscaled
	table[LOGG_COLUMN].desc = "log_10(u.cm * u.second**(-2)) of the surface gravity"

	table[WAVELENGTH_COLUMN].unit = u.Angstrom

	table[FLUX_COLUMN].unit = u.dimensionless_unscaled
	table[FLUX_COLUMN].desc = "in counts"

	# remove the wavelength ranges we don't want

	MIN_WAVELENGTH = 0 * u.micron
	MAX_WAVELENGTH = 15 * u.micron
	table = table[(MIN_WAVELENGTH <= table[WAVELENGTH_COLUMN]) & (table[WAVELENGTH_COLUMN] <= MAX_WAVELENGTH)]

	if CONVERT_WAVELENGTHS_TO_AIR:
		table[WAVELENGTH_COLUMN] = specutils.utils.wcs_utils.vac_to_air(table[WAVELENGTH_COLUMN])

	# add some metadata to the QTable e.g. (wavelength medium = air, source, date)
	table.meta = {"wavelength medium" : "air" if CONVERT_WAVELENGTHS_TO_AIR else "vacuum",
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
		table.write(HDF5_FILENAME_TO_SAVE, path = "data", serialize_meta=True, overwrite=True)
		print("[PHOENIX GRID CREATOR] : hdf5 saving complete")





"""
this is going to download our spectrum data that meets some specified data ranges
this file should just be run once when you want to create the data grid, which is stored in the hdf5 format. Then you can load this hdf5 in quickly into a jupyter notebook etc to do the science 
"""

# external imports
from io import BytesIO
import pandas as pd
import requests
from tqdm import tqdm
import numpy as np
from astropy.io import fits
# from astropy import units as u
from itertools import product
from pathlib import Path

# internal imports
from PHOENIX_filename_conventions import *

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
	
	print("wavelength grid found & loaded in")


# now define the ranges for the data we want (and save this to a hdf5 file)

T_effs = np.arange(2300, 4000, 100)
FeHs = np.arange(-0.5, 0.5, 0.5)
log_gs = np.arange(3.5, 5.5, 0.5)
alphaM = 0
lte : bool = True

total_number_of_files : int = len(T_effs) * len(FeHs) * len(log_gs)

TEFF_COLUMN = "T_eff / K"
FEH_COLUMN = "Fe/H / relative to solar"
LOGG_COLUMN = "log_g / log(cm s^(-2))"
WAVELENGTH_COLUMN = "wavelength / angstroms"
FLUX_COLUMN = "flux / counts"


# we will save our grid to this df
df = pd.DataFrame(columns=[TEFF_COLUMN, FEH_COLUMN, LOGG_COLUMN, WAVELENGTH_COLUMN, FLUX_COLUMN])

i = 0

for T_eff, FeH, log_g in tqdm(product(T_effs, FeHs, log_gs), total=total_number_of_files, desc="Downloading spectra"):
	i += 1
	if i >= 5:
		break
	file = get_file_name(lte, T_eff, log_g, FeH, alphaM)
	url = get_url(file)
	
	try:
		response = requests.get(url)
		response.raise_for_status()
	except requests.exceptions.HTTPError as e:
		print("---")
		print(f"HTTPError raised with the following parameters.\nlte: {lte}\nT_eff={T_eff}\nlog_g={log_g}\nFeH={FeH}\nalphaM={alphaM}\n continuing with the next file")
		print(f"url = {url}")
		print("---")
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
		
		# pandas will repeat the constant values len(fluxes) times for us
		temp_df = pd.DataFrame({
			TEFF_COLUMN : T_eff,
			FEH_COLUMN : FeH,
			LOGG_COLUMN : log_g,
			# need to use the row index <-> wavelength map provided to us by PHOENIX
			WAVELENGTH_COLUMN : wavelengths,
			FLUX_COLUMN : fluxes
		})
		
		# avoid warning about concat-ing an empty df
		if not df.empty:
			# our df index has no meaningful meaning, and sort I think just ensures the columns are in the correct order or something?
			df = pd.concat([df, temp_df], ignore_index=True)#, sort=True)
		else:
			df = temp_df
		
		# tqdm.write("df length = " + str(df.shape[0]))
		# print(df.tail())

import specutils

df[WAVELENGTH_COLUMN] = specutils.utils.wcs_utils.vac_to_air(df[WAVELENGTH_COLUMN])

# interpolation

## [SOME INTERPOLATION HERE] ##

df.to_hdf('data.h5', key='df', mode='w')
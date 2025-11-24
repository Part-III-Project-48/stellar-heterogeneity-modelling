"""file for reading in data from different sources e.g. JWST or HST"""

from pathlib import Path
from astropy.table import QTable
import astropy.units as u
from astropy.io import fits

import numpy as np
from spots_and_faculae_model.spectrum import spectrum

def read_JWST_fits(fits_absolute_path : Path, verbose : bool = False, name : str = None, INTEGRATION_INDEX : int = 0) -> QTable:
	"""
	Attributes
	----------
	verbose : bool (default False)
		prints a summary of all the headers found in the fits file, as well as the string representation of the header with header index HDU_INDEX
	"""
	if not fits_absolute_path.exists():
		raise FileNotFoundError(f"JWST spectrum file not found at : {fits_absolute_path}.")

	with fits.open(fits_absolute_path) as hdul:
		HDU_INDEX = 3 # aka EXTRACT1D
		data = hdul[HDU_INDEX].data
		hdr = hdul[HDU_INDEX].header
		
		# these column name strings are unique to JWST 1D 
		spec : spectrum = spectrum(wavelengths = data["WAVELENGTH"][INTEGRATION_INDEX] * u.um,
				  fluxes = data["FLUX"][INTEGRATION_INDEX] * u.MJy, name=name)

		if verbose:
			hdul.info()
			print(repr(hdr))
	
	return spec

def read_JWST_fits_all_spectra(fits_absolute_path : Path, verbose : bool = False, name : str = None) -> list[spectrum]:
	"""
	Attributes
	----------
	verbose : bool (default False)
		prints a summary of all the headers found in the fits file, as well as the string representation of the header with header index HDU_INDEX
	"""
	if not fits_absolute_path.exists():
		raise FileNotFoundError(f"JWST spectrum file not found at : {fits_absolute_path}.")

	spectra : list[spectrum] = []

	with fits.open(fits_absolute_path) as hdul:
		HDU_INDEX = 3 # aka EXTRACT1D
		data = hdul[HDU_INDEX].data
		hdr = hdul[HDU_INDEX].header
		
		for integration_index in range(len(data["WAVELENGTH"])):
			if (np.all(data["FLUX"][integration_index] == np.nan)):
				print(f"[JWST READER] : integration index {integration_index} contains all nan fluxes")
			# print("[SPECTRUM COMPONENT ANALYSER] : external spectrum found & loaded in")
			# these column name strings are unique to JWST 1D 
			spec : spectrum = spectrum(wavelengths = data["WAVELENGTH"][integration_index] * u.um,
					fluxes = data["FLUX"][integration_index] * u.Jy, name=name)

			if verbose:
				hdul.info()
				print(repr(hdr))
			
			spectra.append(spec)
			
	return spectra
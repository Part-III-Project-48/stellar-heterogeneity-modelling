"""file for reading in data from different sources e.g. JWST or HST"""

from pathlib import Path
import astropy.units as u
from astropy.io import fits
from astropy.units import Unit

import numpy as np
from spots_and_faculae_model.spectrum import spectrum

JWST_WAVELENGTH_UNITS : Unit = u.um
JWST_FLUX_UNITS : Unit = u.MJy

def read_JWST_fits(fits_absolute_path : Path, verbose : bool = False, name : str = None, INTEGRATION_INDEX : int = 0) -> spectrum:
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
		spec : spectrum = spectrum(wavelengths = data["WAVELENGTH"][INTEGRATION_INDEX] * JWST_WAVELENGTH_UNITS,
				  fluxes = data["FLUX"][INTEGRATION_INDEX] * JWST_FLUX_UNITS, name=name)

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


HARPS_WAVELENGTH_UNITS : Unit = u.Angstrom
HARPS_FLUX_UNITS : Unit = u.MJy

def read_HARPS_fits(fits_absolute_path : Path, verbose : bool = False, name : str = None, INTEGRATION_INDEX : int = 0) -> spectrum:
	"""
	Attributes
	----------
	verbose : bool (default False)
		prints a summary of all the headers found in the fits file, as well as the string representation of the header with header index HDU_INDEX
	"""
	if not fits_absolute_path.exists():
		raise FileNotFoundError(f"JWST spectrum file not found at : {fits_absolute_path}.")

	with fits.open(fits_absolute_path) as hdul:
		HDU_INDEX = 1 # aka EXTRACT1D
		

		data = hdul[HDU_INDEX].data
		hdr = hdul[HDU_INDEX].header
		
		if verbose:
			hdul.info()
			print(repr(hdr))

		# these column name strings are unique to JWST 1D 
		spec : spectrum = spectrum(wavelengths = data["WAVE"][INTEGRATION_INDEX] * HARPS_WAVELENGTH_UNITS,
				  fluxes = data["FLUX"][INTEGRATION_INDEX] * HARPS_FLUX_UNITS, name=name, normalised_point=5000*u.Angstrom, normalise_flux=True)
		
	return spec

			

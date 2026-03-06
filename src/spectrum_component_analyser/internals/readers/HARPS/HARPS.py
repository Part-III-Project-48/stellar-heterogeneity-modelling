from pathlib import Path
import astropy.units as u
from astropy.io import fits
from astropy.units import Unit
from astropy.units import Quantity

import numpy as np
from spectrum_component_analyser.internals.spectrum import spectrum

HARPS_WAVELENGTH_UNITS : Unit = u.Angstrom
HARPS_FLUX_UNITS : Unit = u.MJy

HARPS_resolution = 50 * u.Angstrom

HARPS_normalising_point = 5000 * u.Angstrom
# HARPS_smoothing_range = 50 * u.um

from astropy.units import Quantity

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
				  fluxes = data["FLUX"][INTEGRATION_INDEX] * HARPS_FLUX_UNITS, name=name, normalised_point=HARPS_normalising_point, observational_resolution=HARPS_resolution, normalise_flux=True)
		
	return spec

			

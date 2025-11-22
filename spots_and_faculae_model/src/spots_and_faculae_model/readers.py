"""file for reading in data from different sources e.g. JWST or HST"""

from pathlib import Path
from astropy.table import QTable
import astropy.units as u
from astropy.io import fits

from spots_and_faculae_model.spectrum import spectrum

def read_JWST_fits(fits_absolute_path : Path, verbose : bool = False, name : str = None) -> QTable:
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
		if verbose:
			hdul.info()
			print(repr(hdr))
		INTEGRATION_INDEX : int = 0
		# print("[SPECTRUM COMPONENT ANALYSER] : external spectrum found & loaded in")
		# these column name strings are unique to JWST 1D 
		spec : spectrum = spectrum(wavelengths = data["WAVELENGTH"][INTEGRATION_INDEX] * u.um,
				  fluxes = data["FLUX"][INTEGRATION_INDEX] * u.Jy, name=name)
	
	return spec
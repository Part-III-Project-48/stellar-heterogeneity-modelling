"""file for reading in data from different sources e.g. JWST or HST"""

from pathlib import Path
from astropy.table import QTable
import astropy.units as u
from astropy.io import fits

from spots_and_faculae_model.spectrum import spectrum

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
		spec = spectrum(wavelengths = data["WAVELENGTH"][INTEGRATION_INDEX] * u.um,
				  fluxes = data["FLUX"][INTEGRATION_INDEX] * u.Jy)
	
	return spec
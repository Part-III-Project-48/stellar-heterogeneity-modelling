from pathlib import Path
from astropy.io import fits
import astropy.units as u
from matplotlib import pyplot as plt
import numpy as np
from astropy.table import QTable
# from astropy.visualization import quantity_support

from phoenix_grid_creator.fits_to_hdf5 import WAVELENGTH_COLUMN, FLUX_COLUMN

def read_JWST_fits(fits_absolute_path : Path) -> QTable:
	if not fits_absolute_path.exists():
		raise FileNotFoundError(f"JWST spectrum file not found at : {fits_absolute_path}.")

	with fits.open(fits_absolute_path) as hdul:
		# hdul.info()
		HDU_INDEX = 3 # aka EXTRACT1D
		data = hdul[HDU_INDEX].data
		hdr = hdul[HDU_INDEX].header
		# print(repr(hdr))
		INTEGRATION_INDEX : int = 0
		# print("[SPECTRUM COMPONENT ANALYSER] : external spectrum found & loaded in")
		spectrum = QTable()
		spectrum[WAVELENGTH_COLUMN] = data["WAVELENGTH"][INTEGRATION_INDEX] * u.um
		spectrum[FLUX_COLUMN] = data["FLUX"][INTEGRATION_INDEX] * u.Jy
	
	return spectrum














# literally no clue what this was for

# def get_external_spectra(external_spectrum_absolute_path : Path) -> QTable:
# 	if not external_spectrum_absolute_path.exists():
# 		raise FileNotFoundError(f"wavelength grid file not found at : {external_spectrum_absolute_path}.")

# 	with fits.open(external_spectrum_absolute_path) as hdul:
# 		# there is only 1 HDU in this wavelength grid file
# 		WAVELENGTH_GRID_HDU_INDEX = 0
# 		data = hdul[WAVELENGTH_GRID_HDU_INDEX].data

# 		# the fits file is big endian; pandas requires little endian. this swaps between them
# 		# wavelengths = wavelengths.byteswap().view(wavelengths.dtype.newbyteorder())
		
# 		print("[SPECTRUM COMPONENT ANALYSER] : external spectrum found & loaded in")

# 	wavelengths = np.arange(0, len(data), 1)
# 	wavelengths *= u.Angstrom

# 	fluxes = data

# 	t = QTable()
# 	t[WAVELENGTH_COLUMN] = wavelengths
# 	t[FLUX_COLUMN] = fluxes

# 	t = t[t[FLUX_COLUMN] != 0]
# 	print(t)
# 	return t
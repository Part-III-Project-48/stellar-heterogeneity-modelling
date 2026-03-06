"""
reads in spectra (one or many) from NIRISS-formatted *.fits files 
"""

from enum import Enum
from pathlib import Path
import astropy.units as u
from astropy.io import fits
from astropy.units import Unit, Quantity
import numpy as np

from spectrum_component_analyser.internals.spectrum import spectrum

JWST_NORMALISING_POINT = 1.1 * u.um # believe this is unused

# ef
class NIRISSReader():
	WAVELENGTH_UNITS : Unit = u.um
	FLUX_UNITS : Unit = u.MJy
	RESOLUTION : Quantity[u.um] = .001 * u.um

	verbose = False

	def get_spectrum(self, file_path : Path, INTEGRATION_INDEX : int = 0, name="Observational JWST NIRISS spectrum", verbose : bool = False) -> spectrum:
		with fits.open(file_path) as hdul:
			HDU_INDEX = 3 	# aka EXTRACT1D
			hdr = hdul[HDU_INDEX].header
			data = hdul[HDU_INDEX].data
			
			# these column name strings are unique to JWST 1D 
			spec : spectrum = spectrum(wavelengths = data["WAVELENGTH"][INTEGRATION_INDEX] * self.WAVELENGTH_UNITS,
					fluxes = data["FLUX"][INTEGRATION_INDEX] * self.FLUX_UNITS,
					name=name,
					normalised_point = JWST_NORMALISING_POINT,
					temperature=None,
					observational_resolution=None, # this is an observational spectrum (as we are reading in a JWST fits file) - so no convolution or resampling is necessary
					observational_wavelengths=None)

			if verbose:
				hdul.info()
				print(repr(hdr))
		
		return spec

	def get_all_spectra(self, file_path : Path, name : str = "Observational JWST NIRISS spectrum", verbose : bool = False) -> list[spectrum]:
		"""
		Attributes
		----------
		verbose : bool (default False)
			prints a summary of all the headers found in the fits file, as well as the string representation of the header with header index HDU_INDEX
		"""

		spectra : list[spectrum] = []

		with fits.open(file_path) as hdul:
			HDU_INDEX = 3 # aka EXTRACT1D
			data = hdul[HDU_INDEX].data
			hdr = hdul[HDU_INDEX].header
			
			for integration_index in range(len(data["WAVELENGTH"])):
				if (np.all(data["FLUX"][integration_index] == np.nan)):
					print(f"[JWST READER] : integration index {integration_index} contains all nan fluxes")
				
				# these column name strings are unique to JWST 1D 
				spec : spectrum = spectrum(wavelengths = data["WAVELENGTH"][integration_index] * self.WAVELENGTH_UNITS,
						fluxes = data["FLUX"][integration_index] * self.FLUX_UNITS,
						normalised_point=JWST_NORMALISING_POINT, # this is an observational spectrum: no normalising or interpolation should be done on it
						temperature=None,
						observational_resolution=None,
						observational_wavelengths=None,
						name=name)

				if verbose:
					hdul.info()
					print(repr(hdr))
				
				spectra.append(spec)
		
		return spectra

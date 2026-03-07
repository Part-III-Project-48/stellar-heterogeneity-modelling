"""
reads in spectra (one or many) from NIRISS-formatted *.fits files 
"""

from enum import Enum
from pathlib import Path
import astropy.units as u
from astropy.io import fits
from astropy.units import Unit, Quantity
import numpy as np

from spectrum_component_analyser.internals.readers.JWST.instruments import Instrument
from spectrum_component_analyser.internals.spectrum import spectrum

JWST_NORMALISING_POINT = 1.1 * u.um # believe this is unused

class JWSTFileReader():
	@staticmethod
	def get_spectrum(file_path : Path, instrument : Instrument, INTEGRATION_INDEX : int = 0, name="Observational JWST NIRISS spectrum", verbose : bool = False) -> spectrum:
		"""
		returns the spectrum at the given INTEGRATION_INDEX from the specified .fits file
		"""
		with fits.open(file_path) as hdul:
			hdr = hdul[instrument.HDUIndex].header
			data = hdul[instrument.HDUIndex].data
			
			if verbose:
				hdul.info()
				print(f"Printing header for hdu index {instrument.HDUIndex}")
				print(repr(hdr))
			
			# these column name strings are unique to JWST 1D NIRISS data 
			spec : spectrum = spectrum(wavelengths = data["WAVELENGTH"][INTEGRATION_INDEX] * instrument.WavelengthUnits,
					fluxes = data["FLUX"][INTEGRATION_INDEX] * instrument.FluxUnits,
					name=name,
					normalised_point = JWST_NORMALISING_POINT,
					temperature=None,
					observational_resolution=None, # this is an observational spectrum (as we are reading in a JWST fits file) - so no convolution or resampling is necessary
					observational_wavelengths=None)

		
		return spec

	@staticmethod
	def get_all_spectra(file_path : Path, instrument : Instrument, name : str = "Observational JWST NIRISS spectrum", verbose : bool = False) -> list[spectrum]:
		"""
		returns all spectra from the specified .fits file
		
		Attributes
		----------
		verbose : bool (default False)
			prints a summary of all the headers found in the fits file, as well as the string representation of the header with header index HDU_INDEX
		"""

		spectra : list[spectrum] = []

		with fits.open(file_path) as hdul:
			data = hdul[instrument.HDUIndex].data
			hdr = hdul[instrument.HDUIndex].header
			
			if verbose:
				print(f"/// .fits file summary")
				hdul.info()
				print(f"/// header for hdu index {instrument.HDUIndex}")
				print(repr(hdr))
			
			for integration_index in range(len(data["WAVELENGTH"])):
				if (np.all(data["FLUX"][integration_index] == np.nan)):
					print(f"[JWST READER] : integration index {integration_index} contains all nan fluxes")
				
				# these column name strings are unique to JWST 1D 
				spec : spectrum = spectrum(wavelengths = data["WAVELENGTH"][integration_index] * instrument.WavelengthUnits,
						fluxes = data["FLUX"][integration_index] * instrument.FluxUnits,
						normalised_point=JWST_NORMALISING_POINT, # this is an observational spectrum: no normalising or interpolation should be done on it
						temperature=None,
						observational_resolution=None,
						observational_wavelengths=None,
						name=name)
				
				spectra.append(spec)
		
		return spectra

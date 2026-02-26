from typing import Sequence
from astropy.units import Quantity
import numpy as np

from spectrum_component_analyser.internals.phoenix_spectrum import phoenix_spectrum
from spectrum_component_analyser.internals.readers import JWST_NORMALISING_POINT, JWST_RESOLUTION
from spectrum_component_analyser.internals.spectral_component import spectral_component
from spectrum_component_analyser.internals.spectral_grid import download_spectrum, get_wavelength_grid

class spectral_list():
	"""
	This class is designed for storing e.g. 5 or so spectral components so I'm not going to implement save / loading to hdf5 for it. If you need more spectra, use a spectral grid instead.
	"""
	def __init__(
			self,
			phoenix_spectra : list[phoenix_spectrum],
			_internal=False
			):
		"""
		Don't use this init: use the wrappers instead :)
		"""
		if (not _internal):
			raise RuntimeError("Spectral List's __init__ should only be used by internal methods: use factory methods instead")
		
		self.PhoenixSpectra : list[phoenix_spectrum] = phoenix_spectra
	
	@classmethod
	def from_internet(cls,
				   spectral_components : list[spectral_component],
				   normalising_point : Quantity,
				   observational_resolution : Quantity,
				   observational_wavelengths : np.ndarray,
				   name : str = "spectral_list",
				   alphaM = 0,
				   lte = True):
		"""
		Analagous to spectral_grid.from_internet, but an explicit list of spectral_components are the only ones downloaded (instead of a cartesian product over 3 distinct lists, which is what spectral_grid does).

		Use this if you want a small list of phoenix spectra whose parameters which you know beforehands.
		"""
		phoenix_spectra : list[phoenix_spectrum] = []

		phoenix_wavelengths = get_wavelength_grid()
		
		component : spectral_component # yay for python
		for component in spectral_components:
			phoenix_spec : phoenix_spectrum = download_spectrum(
				component.T_eff,
				component.FeH,
				component.Log_g,
				lte=lte,
				alphaM=alphaM,
				phoenix_wavelengths=phoenix_wavelengths,
				normalising_point=normalising_point,
				observational_resolution=observational_resolution,
				observational_wavelengths=observational_wavelengths,
				name=name
			)

			phoenix_spectra.append(phoenix_spec)

		return cls(phoenix_spectra=phoenix_spectra, _internal=True)
		




"""
super simple class that just wraps 2 arrays for arbitrary x, y axes for a spectrum

also has some util functions for working with spectra in general
"""

import numpy as np
import astropy.units as u
import scipy as sp
from astropy.visualization import quantity_support
quantity_support()
from matplotlib import pyplot as plt
import specutils
from astropy.units import Quantity

## can do : update normalise jansksys to act on the spectrum class self and then update main.ipynb to use that

DEFAULT_FLUX_UNIT = u.Jy

class spectrum:
	def __init__(self, wavelengths : np.array, fluxes : np.array, name : str = None, normalise_flux : bool = True, normalised_point : Quantity = 2.2 * u.um):
		"""
		Flux is going to be stored in Janskys from now on

		This initialiser will normalise janskys for us; any reference to normalise_janskys outside of this class is redundant (if only python had private functions :/)

		Attributes
		----------

		wavelengths : np.array
			an array of astropy quantities

		phoenix_fluxes : np.array
			an array of astropy quantities (with some units that are convertible to Janskys by u.spectral_density equivalencies; e.g. Janskys themself or [erg / (s * cm**2 * cm)])
		"""

		wavelengths = np.atleast_1d(wavelengths)
		fluxes = np.atleast_1d(fluxes)
		
		if len(wavelengths) != len(fluxes):
			raise ValueError("wavelengths and fluxes must have the same length")
		
		fluxes_janskys = fluxes.to(DEFAULT_FLUX_UNIT, equivalencies=u.spectral_density(wavelengths))

		# make sure the wavelengths are in ascending order, so that normalising_janskys doesn't break
		indices = np.argsort(wavelengths) # get the indices that would sort the wavelengths np.array

		self.Wavelengths : np.array = wavelengths[indices]
		self.Fluxes : np.array = fluxes_janskys[indices]
		self.Name : str = name

		if normalise_flux:
			self.normalise_flux(normalised_point=normalised_point)
		
		self.Normalised_Point = normalised_point

	def __getitem__(self, idx):
		"""
		Allow slicing, indexing, and boolean masks.
		Returns a new spectrum with sliced wavelength and flux arrays.
		"""
		return spectrum(self.Wavelengths[idx], self.Fluxes[idx], name=self.Name, normalised_point=self.Normalised_Point)

	def __len__(self):
		if len(self.Wavelengths) != len(self.Fluxes):
			raise ValueError("wavelengths and fluxes must have the same length. (Length is not well defined)")
		
		return len(self.Fluxes)
	
	def __iter__(self):
		"""
		pandas complains very hard when trying to print dataframes containing spectrum objects if we dont define our own iterator
		"""
		return iter(self.Fluxes)
	
	def __repr__(self):
		unit = getattr(self.Fluxes, "unit", None)
		return f"<spectrum name={self.Name} len={len(self)} unit={unit}>"

	def __str__(self):
		return self.__repr__()
	
	def plot(self):
		plt.clf()
		plt.title(f"Observational Spectrum for {self.Name}")
		plt.plot(self.Wavelengths, self.Fluxes)
		plt.show()

	def normalise_flux(self, normalised_point, smoothing_range = 0.5 * u.um) -> np.array:
		"""
		this will fail if wavelengths does not span at least smoothing_range
		
		the canonical way to normalise spectra is to choose a portion that's continuum and make that be at a consistent scale
		
		inputs should both be lists of astropy quantities (aka values with astropy units)
		
		I tried doing this but the continuum fit was terrible: https://specutils.readthedocs.io/en/stable/fitting.html#continuum-fitting
		"""

		if (u.get_physical_type(self.Fluxes[0].unit) != u.get_physical_type(u.Jy)):
			raise ValueError(f"fluxes are in units of {self.Fluxes.unit}. this is not in a unit convertible to janskys. no normalisation will be carried out.")

		# kernel size of about 501 with 9999 points between 5 and 15 um seemed good - this range corresponds (roughly) to that 
		wavelengths_in_range = self.Wavelengths[(self.Wavelengths[0] <= self.Wavelengths) & (self.Wavelengths <= self.Wavelengths[0] + smoothing_range)]
		kernel_size = len(wavelengths_in_range)
		if kernel_size % 2 == 0:
			kernel_size +=1
		
		# smooth to make sure there's no spikes
		unit = self.Fluxes.unit # doesn't matter what this is: it just has to be the same for the dividing out and timesing (medfilt seems to silently remove units)
		counts = [(i / unit).value for i in self.Fluxes]
		counts = np.array(counts, dtype=np.float64)
		
		smoothed_counts = sp.signal.medfilt(counts, kernel_size=[kernel_size])
		# normalise the counts at normalised_point (or next nearest value) to be 1
		counts /= smoothed_counts[(normalised_point <= self.Wavelengths)][0]

		self.Fluxes = counts * unit
	
	@property
	def air_wavelengths(self):
		"""
		this assumes that the hdf5 is in vacuum units; can easily check metadata of hdf5 file
		"""
		return specutils.utils.wcs_utils.vac_to_air(self.Wavelengths)

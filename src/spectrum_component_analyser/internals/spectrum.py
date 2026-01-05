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
import warnings
from scipy.ndimage import gaussian_filter1d

## can do : update normalise jansksys to act on the spectrum class self and then update main.ipynb to use that

DEFAULT_FLUX_UNIT = u.Jy

class spectrum:
	def __init__(self, wavelengths : np.array, fluxes : np.array, normalised_point : Quantity, desired_resolution : Quantity, output_wavelengths : np.ndarray, name : str = None):
		"""
		Flux is going to be stored in Janskys from now on

		This initialiser will normalise janskys for us; any reference to normalise_janskys outside of this class is redundant (if only python had private functions :/)

		Leave desired_resolution or normalised_point both to ignore regridding to a given resolution and/or normalising respectively.

		output_wavelengths must have a resolution of (at least approximately) the input desired resolution.

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

		if desired_resolution != None:
			self.regrid_flux(desired_resolution=desired_resolution)

		if normalised_point != None:
			self.normalise_flux(normalised_point)

		if output_wavelengths != None:
			self.Wavelengths = output_wavelengths
			self.Fluxes = np.interp(output_wavelengths, self.Wavelengths, self.Fluxes) # new = np.interp(new | old | old)
		
		self.Normalised_Point = normalised_point
		self.Desired_Resolution = desired_resolution

	# def __getitem__(self, idx):
	# 	"""
	# 	Allow slicing, indexing, and boolean masks.
	# 	Returns a new spectrum with sliced wavelength and flux arrays.
	# 	"""
	# 	return spectrum(self.Wavelengths[idx], self.Fluxes[idx], name=self.Name, normalised_point=self.Normalised_Point, desired_resolution=self.Desired_Resolution)

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
	
	def normalise_flux(self, normalised_point : Quantity):
		"""
		normalise the counts at normalised_point (or next nearest value) to be 1

		normalised_point : an astropy quantity with dimension of length
		"""

		self.Fluxes /= self.Fluxes[(normalised_point <= self.Wavelengths)][0].value
	
	def regrid_flux(self, desired_resolution : Quantity) -> np.array:
		"""
		Regrid the spectrum onto a uniform wavelength array of the input resolution. Uses a gaussian to simulate how real data would be recorded.

		This method assumes that the desired resolution >> the current resolution of the spectrum when this function is called.
		"""

		if (u.get_physical_type(self.Fluxes[0].unit) != u.get_physical_type(u.Jy)):
			raise ValueError(f"fluxes are in units of {self.Fluxes.unit}. this is not in a unit convertible to janskys. no normalisation will be carried out.")
		
		# interpolate ontoa uniform wavelength grid - self.Wavelengths & self.Fluxes are assumed to currently be very high res
		wave_uniform = np.linspace(self.Wavelengths.min(), self.Wavelengths.max(), len(self.Wavelengths))
		flux_uniform = np.interp(wave_uniform, self.Wavelengths, self.Fluxes)

		# convert resolution to sigma (in array indices)
		delta_lambda = wave_uniform[1] - wave_uniform[0]
		sigma_pix = desired_resolution / (2 * np.sqrt(2 * np.log(2)) * delta_lambda)

		convolved_flux = gaussian_filter1d(flux_uniform, sigma_pix, mode="nearest")

		# resample onto desired wavelengths
		desired_number_of_wavelength_points = (self.Wavelengths.max() - self.Wavelengths.min()) / desired_resolution
		wave_desired_resolution = np.linspace(self.Wavelengths.min(), self.Wavelengths.max(), desired_number_of_wavelength_points)

		new_flux = np.interp(wave_desired_resolution, wave_uniform, convolved_flux)

		self.Wavelengths = wave_desired_resolution
		self.Fluxes = new_flux

	@property
	def air_wavelengths(self):
		"""
		this assumes that the hdf5 is in vacuum units; can easily check metadata of hdf5 file
		"""
		return specutils.utils.wcs_utils.vac_to_air(self.Wavelengths)

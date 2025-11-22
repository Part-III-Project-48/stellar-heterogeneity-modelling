"""
super simple class that just wraps 2 arrays for arbitrary x, y axes for a spectrum

also has some util functions for working with spectra in general
"""

from astropy.table import QTable
import numpy as np
import astropy.units as u
import scipy as sp
import astropy
from astropy.visualization import quantity_support
quantity_support()
from matplotlib import pyplot as plt

## can do : update normalise jansksys to act on the spectrum class self and then update main.ipynb to use that

class spectrum:
    def __init__(self, wavelengths : np.array, fluxes : np.array, name : str = None):
        """
        both arrays should be astropy quantities with units
        """
        self.Wavelengths = wavelengths
        self.Fluxes = fluxes
        self.Name : str = name

        if len(wavelengths) != len(fluxes):
            raise ValueError("wavelengths and fluxes must have the same length")

    def __getitem__(self, idx):
        """
        Allow slicing, indexing, and boolean masks.
        Returns a new spectrum with sliced wavelength and flux arrays.
        """
        return spectrum(self.Wavelengths[idx], self.Fluxes[idx], name=self.Name)

    def __len__(self):
        if len(self.Wavelengths) != len(self.Fluxes):
            raise ValueError("wavelengths and fluxes must have the same length. (Length is not well defined)")
        
        return len(self.Fluxes)
    
    def plot(self):
        plt.title(f"Observational Spectrum for {self.Name}")
        plt.plot(self.Wavelengths, self.Fluxes)
        plt.show()

    def normalise_Janskys(self, normalised_point = 2.2 * u.um, smoothing_range = 0.5 * u.um) -> np.array:
        """
        this will fail if wavelengths does not span at least smoothing_range
        
        the canonical way to normalise spectra is to choose a portion that's continuum and make that be at a consistent scale
        
        inputs should both be lists of astropy quantities (aka values with astropy units)
        
        I tried doing this but the continuum fit was terrible: https://specutils.readthedocs.io/en/stable/fitting.html#continuum-fitting
        """

        # check for correct units
        # print(f"self.Fluxes[0].unit = {self.Fluxes[0].unit}")
        # self.Fluxes.to(u.Jy)
        # return

        if (u.get_physical_type(self.Fluxes[0].unit) != u.get_physical_type(u.Jy)):
            raise ValueError(f"fluxes are in units of {self.Fluxes.unit}. this is not in a unit convertible to janskys. no normalisation will be carried out.")
            return

        # kernel size of about 501 with 9999 points between 5 and 15 um seemed good - this range corresponds (roughly) to that 
        wavelengths_in_range = self.Wavelengths[(self.Wavelengths[0] <= self.Wavelengths) & (self.Wavelengths <= self.Wavelengths[0] + smoothing_range)]
        kernel_size = len(wavelengths_in_range)
        if kernel_size % 2 == 0:
            kernel_size +=1
        
        # smooth to make sure there's no spikes
        counts = [i.value for i in self.Fluxes]
        counts = np.array(counts, dtype=np.float64)
        
        smoothed_counts = sp.signal.medfilt(counts, kernel_size=[kernel_size])
        # normalise the counts at normalised_point (or next nearest value) to be 1
        counts /= smoothed_counts[(normalised_point <= self.Wavelengths)][0]

        self.Fluxes = counts * u.Jy

    @classmethod
    def from_phoenix_units(cls, wavelengths : np.array, phoenix_fluxes : np.array):
        """
        construct a spectrum with janskys as the flux units, from phoenix units

        wavelengths should be an np.array of astropy quantities with units of length

        PHOENIX flux is in power density // spectral flux density (per unit wavelength)
        aka erg / (s * cm**2 * cm)

        Jy = erg / (s cm² Hz)

        so 1 Jy = PHOENIX units * wavelength / (frequency * c)

        e.g. see here https://physics.stackexchange.com/questions/725928/converting-between-f-nu-and-f-lambda-spectral-density
        """

        converted_fluxes = phoenix_fluxes * wavelengths**2 / astropy.constants.c

        return spectrum(wavelengths=wavelengths, fluxes=converted_fluxes)

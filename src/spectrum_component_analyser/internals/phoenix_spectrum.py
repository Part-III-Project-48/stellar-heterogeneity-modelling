import numpy as np
from astropy.units import Quantity
import astropy.units as u

from spectrum_component_analyser.internals.spectrum import spectrum

class phoenix_spectrum(spectrum):
   def __init__(self, wavelengths : np.array, fluxes : np.array, t_eff : Quantity[u.K], feh, log_g, normalising_point : Quantity, desired_resolution : Quantity, name : str = None, normalise_flux = True):
        super().__init__(wavelengths, fluxes, normalised_point=normalising_point, desired_resolution=desired_resolution, normalise_flux=normalise_flux, name=name)
        self.T_eff = t_eff
        self.FeH = feh
        self.Log_g = log_g
import numpy as np
from astropy.units import Quantity
import astropy.units as u

from spectrum_component_analyser.internals.spectrum import spectrum

class phoenix_spectrum(spectrum):
   def __init__(self, wavelengths : np.array, fluxes : np.array, t_eff : Quantity[u.K], feh, log_g, normalising_point : Quantity, smoothing_range : Quantity, name : str = None, normalise_flux = True):
        super().__init__(wavelengths, fluxes, normalised_point=normalising_point, smoothing_range=smoothing_range, normalise_flux=normalise_flux, name=name)
        self.T_eff = t_eff
        self.FeH = feh
        self.Log_g = log_g
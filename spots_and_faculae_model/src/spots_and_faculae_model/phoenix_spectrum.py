import numpy as np
from astropy.units import Quantity
import astropy.units as u

from spots_and_faculae_model.spectrum import spectrum

class phoenix_spectrum(spectrum):
   def __init__(self, wavelengths : np.array, fluxes : np.array, t_eff : Quantity[u.K], feh, log_g, name : str = None):
        super().__init__(wavelengths, fluxes, name)
        self.T_eff = t_eff
        self.FeH = feh
        self.Log_g = log_g
import numpy as np
from astropy.units import Quantity
import astropy.units as u

from spectrum_component_analyser.internals.spectrum import spectrum

class phoenix_spectrum(spectrum):
   def __init__(
           self,
           wavelengths : np.array,
           fluxes : np.array,
           t_eff : Quantity[u.K],
           feh,
           log_g,
           normalising_point : Quantity,
           observational_resolution : Quantity,
           observational_wavelengths : np.ndarray,
           name : str
         ):
        super().__init__(
            wavelengths,
            fluxes,
            normalised_point=normalising_point,
            observational_resolution=observational_resolution,
            observational_wavelengths=observational_wavelengths,
            temperature=t_eff,
            name=name
         )
        self.T_eff = t_eff
        self.FeH = feh
        self.Log_g = log_g
"""
used for comparing our approximation to the currently stated values from the exoplanet archive
"""

from pathlib import Path
from typing import Sequence, Tuple
from astropy.units import Quantity

import astropy.units as u
from matplotlib import pyplot as plt
import numpy as np
from scipy.interpolate import interpn

from spectrum_component_analyser.internals.phoenix_spectrum import phoenix_spectrum
from spectrum_component_analyser.internals.readers.JWST import JWST_NORMALISING_POINT
from spectrum_component_analyser.internals.spectral_grid import spectral_grid

def get_interpolated_phoenix_spectrum(
        T_eff : Quantity[u.K],
        FeH : Quantity[u.dimensionless_unscaled],
        Log_g : Quantity[u.dimensionless_unscaled],
        star_name : str,
        spec_grid : spectral_grid
    ) -> phoenix_spectrum:

    parameter_space : Tuple[
        Sequence[Quantity[u.K]], Sequence[Quantity[u.K]], Sequence[Quantity[u.K]]
        ] = (spec_grid.T_effs, spec_grid.FeHs, spec_grid.Log_gs, spec_grid.Wavelengths)

    interpolated_fluxes = []

    # kinda a bodge but its okay for now
    for w in spec_grid.Wavelengths:
        desired_spectrum_parameters : Tuple[float, float, float, float] = [T_eff.value, FeH.value, Log_g.value, w.value] # if only numpy supported units :(

        f = spec_grid.Fluxes

        v : np.ndarray = interpn(parameter_space, f, desired_spectrum_parameters)

        interpolated_fluxes.append(*v) # unpack v as it is returned as a 1-element np.ndarray

    # reintroduce unitsW
    interpolated_fluxes *= spec_grid.Fluxes.unit
    
    # plt.plot(spec_grid.Wavelengths, interpolated_fluxes)
    # plt.show()

    interpolated_spectrum : phoenix_spectrum = phoenix_spectrum(
        wavelengths=spec_grid.Wavelengths,
        fluxes=interpolated_fluxes,
        t_eff=T_eff,
        feh=FeH,
        log_g=Log_g,
        normalising_point=JWST_NORMALISING_POINT, 
        observational_resolution=None, # no convolution / regridding - as the spectrum should already be convoluted to JWST
        observational_wavelengths=None,
        name=star_name
    )

    return interpolated_spectrum

if __name__ == "__main__":
    # example values
    T_eff = 3358 * u.K
    FeH = 0.06 * u.dimensionless_unscaled
    Log_g = 4.85 * u.dimensionless_unscaled
    star_name = "LTT-3780-from-Literature"

    get_interpolated_phoenix_spectrum(
        T_eff=T_eff,
        FeH=FeH,
        Log_g=Log_g,
        star_name=star_name
    ).plot()
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
from spectrum_component_analyser.internals.spectral_grid import spectral_grid

# example values
# T_eff = 3358 * u.K
# FeH = 0.06 * u.dimensionless_unscaled
# Log_g = 4.85 * u.dimensionless_unscaled
# star_name = "LTT-3780-from-Literature"

def get_interpolated_phoenix_spectrum(
        T_eff : Quantity[u.K],
        FeH : Quantity[u.dimensionless_unscaled],
        Log_g : Quantity[u.dimensionless_unscaled],
        star_name : str
    ) -> phoenix_spectrum:

    spectral_grid_relative_path = Path("../../../spectral_grids/JWST_convolved_not_oversmoothed.hdf5")
    spectral_grid_absolute_path = (__file__ / spectral_grid_relative_path).resolve()
    spec_grid : spectral_grid = spectral_grid.from_hdf5(absolute_path=spectral_grid_absolute_path)

    parameter_space : Tuple[
        Sequence[Quantity[u.K]], Sequence[Quantity[u.K]], Sequence[Quantity[u.K]]
        ] = (spec_grid.T_effs, spec_grid.FeHs, spec_grid.Log_gs, spec_grid.Wavelengths)

    interpolated_fluxes = []

    # kinda a bodge but its okay for now
    for w in spec_grid.Wavelengths:
        desired_spectrum_parameters : Tuple[float, float, float, float] = [T_eff.value, FeH.value, Log_g.value, w.value] # if only numpy supported units :(

        v = interpn(parameter_space, spec_grid.Fluxes, desired_spectrum_parameters)

        interpolated_fluxes.append(v)
    
    # plt.plot(spec_grid.Wavelengths, interpolated_fluxes)
    # plt.show()

    interpolated_spectrum : phoenix_spectrum = phoenix_spectrum(
        wavelengths=spec_grid.Wavelengths,
        fluxes=interpolated_fluxes,
        t_eff=T_eff,
        feh=FeH,
        log_g=Log_g,
        normalising_point=None, # no normalisation / convolution
        observational_resolution=None,
        observational_wavelengths=None,
        name=star_name
    )

    return interpolated_spectrum
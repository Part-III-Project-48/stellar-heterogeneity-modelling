"""
used for comparing our approximation to the currently stated values from the exoplanet archive
"""

from pathlib import Path
from typing import Sequence, Tuple
from astropy.units import Quantity

import astropy.units as u
from matplotlib import pyplot as plt
import numpy as np
from scipy.interpolate import RegularGridInterpolator, interpn

from spectrum_component_analyser.internals.phoenix_spectrum import phoenix_spectrum
from spectrum_component_analyser.internals.readers.JWST.file_reader import JWST_NORMALISING_POINT
from spectrum_component_analyser.internals.spectral_grid import spectral_grid

def get_interpolated_phoenix_spectrum(
        T_eff : Quantity[u.K],
        FeH : Quantity[u.dimensionless_unscaled],
        Log_g : Quantity[u.dimensionless_unscaled],
        star_name : str,
        spec_grid : spectral_grid
    ) -> phoenix_spectrum:

    interp = RegularGridInterpolator(
        (
            spec_grid.T_effs.value,
            spec_grid.FeHs.value,
            spec_grid.Log_gs.value
        ),
        spec_grid.Fluxes,
        bounds_error=False,
        fill_value=True
    )

    params = np.array(
        [
            [T_eff.value, FeH.value, Log_g.value]
        ]
    )

    interpolated_fluxes = interp(params)[0] * spec_grid.Fluxes.unit

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

def old_get_interpolated_phoenix_spectrum(
        T_eff : Quantity[u.K],
        FeH : Quantity[u.dimensionless_unscaled],
        Log_g : Quantity[u.dimensionless_unscaled],
        star_name : str,
        spec_grid : spectral_grid
    ) -> phoenix_spectrum:
    """
    old method
    """
    parameter_space : Tuple[
        Sequence[Quantity[u.K]], Sequence[Quantity[u.K]], Sequence[Quantity[u.K]]
        ] = (spec_grid.T_effs, spec_grid.FeHs, spec_grid.Log_gs, spec_grid.Wavelengths)

    w = spec_grid.Wavelengths.value

    interpolation_points = np.column_stack([
        np.full_like(w, T_eff.value),
        np.full_like(w, FeH.value),
        np.full_like(w, Log_g.value),
        w
    ])

    # remember to reintroduce units
    interpolated_fluxes = interpn(parameter_space, spec_grid.Fluxes, interpolation_points) * spec_grid.Fluxes.unit
    
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
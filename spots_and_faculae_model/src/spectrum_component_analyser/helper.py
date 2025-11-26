
from itertools import product
from pathlib import Path
import astropy
from astropy.table import QTable
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import scipy as sp
from astropy.visualization import quantity_support
quantity_support()
from tqdm import tqdm
import astropy.units as u
from scipy.interpolate import interp1d
from astropy.units import Quantity
from joblib import Parallel, delayed
import os
from typing import Sequence, Tuple
from scipy.optimize._optimize import OptimizeResult

from spots_and_faculae_model.spectrum import spectrum
from spots_and_faculae_model.readers import read_JWST_fits
from spots_and_faculae_model.simpler_spectral_grid import simpler_spectral_grid

def calc_result(parameter_space, lookup_table, spec_grid : simpler_spectral_grid, mask, spectrum_to_decompose, total_number_of_components : int = None, verbose : bool = True) -> Tuple[np.ndarray, OptimizeResult]:
    A = np.empty((0, 0))

    def force_to_janskys(T_eff : Quantity, FeH : Quantity, log_g : Quantity, wavelengths : Sequence[Quantity], mask):
        fluxes = lookup_table[T_eff, FeH, log_g]
        return fluxes.to(u.Jy, equivalencies=u.spectral_density(wavelengths))[mask]

    normalised_and_converted_spectral_components : list[list[Quantity]] = Parallel(n_jobs=-1, prefer="threads")(
        delayed(force_to_janskys)(T_eff, FeH, log_g, spec_grid.Wavelengths, mask) for T_eff, FeH, log_g in tqdm(parameter_space, total=total_number_of_components, desc="Appending values to A matrix...", disable=not verbose)
    )

    A = np.column_stack(normalised_and_converted_spectral_components)

    if verbose:
        print("minimising")
    
    # assume that w \in [0,1] : but I think this will only be true for real data if normalisation has been done correctly (???)
    result : OptimizeResult = sp.optimize.lsq_linear(A, [i.value for i in spectrum_to_decompose.Fluxes], bounds = (0, 1), verbose = 2 if verbose else 0, max_iter=30)#, tol=1e-10, lsmr_tol=1e-5)
    
    if verbose:
        print(result)
        print(f"sum of weights={np.sum(result.x)}")

    return A, result

def get_optimality(A, result, spectrum_to_decompose : spectrum):
    determined_spectrum = spectrum(spectrum_to_decompose.Wavelengths, A @ result.x, normalised_point=None, smoothing_range=None, normalise_flux=False)
    residual = (determined_spectrum.Fluxes - spectrum_to_decompose.Fluxes) / spectrum_to_decompose.Fluxes
    rmse = np.sqrt(np.mean(residual**2))
    rss  = np.sum(residual**2)

    return rmse, rss
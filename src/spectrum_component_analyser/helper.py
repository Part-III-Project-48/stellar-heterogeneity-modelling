
import numpy as np
from matplotlib import pyplot as plt
import pandas as pd
import scipy as sp
from astropy.visualization import quantity_support
quantity_support()
from tqdm import tqdm
import astropy.units as u
from astropy.units import Quantity
from joblib import Parallel, delayed
from typing import Sequence, Tuple
from scipy.optimize._optimize import OptimizeResult

from spectrum_component_analyser.internals.spectrum import spectrum
from spectrum_component_analyser.internals.spectral_grid import spectral_grid

# units should be stored in the astropy quantity anyway
# changing these is fine, as long as a new spectral grid is created which uses these column names
TEFF_COLUMN = "T_eff / K"
FEH_COLUMN = "Fe/H / relative to solar"
LOGG_COLUMN = "log_g / log(cm s^(-2))"
WAVELENGTH_COLUMN = "wavelength / angstroms"
FLUX_COLUMN = "flux / erg / (s * cm**2 * cm)"
SPECTRUM_COLUMN : str = "spectrum object"

def calc_result(parameter_space, lookup_table, spec_grid : spectral_grid, mask, spectrum_to_decompose, total_number_of_components : int = None, verbose : bool = True) -> Tuple[np.ndarray, OptimizeResult]:
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
    from scipy import sparse

    # this seems to break the fitting
    # A = sparse.csr_matrix(A)

    # this change seems to remove units from result.x?
    result : OptimizeResult = sp.optimize.lsq_linear(A, spectrum_to_decompose.Fluxes.value, bounds = (0, 1), verbose = 2 if verbose else 0, max_iter=100)#, tol=1e-10, lsmr_tol=1e-5)

    if verbose:
        print(result)
        print(f"sum of weights={np.sum(result.x)}")

    return A, result

def get_optimality(A, result, spectrum_to_decompose : spectrum):
    determined_spectrum = spectrum(spectrum_to_decompose.Wavelengths, A @ result.x, normalised_point=None, observational_resolution=None, normalise_flux=False)
    residual = (determined_spectrum.Fluxes - spectrum_to_decompose.Fluxes) / spectrum_to_decompose.Fluxes
    residual_mean_squared_error = np.sqrt(np.mean(residual**2))
    residual_sum_of_squares  = np.sum(residual**2)

    return residual_mean_squared_error, residual_sum_of_squares


# # # plot some data # # #
# dependent on the old spectrum_grid class, but its fine for now (and its just dependent on some arbitrary strings anyway)
from spectrum_component_analyser.helper import TEFF_COLUMN, FEH_COLUMN, LOGG_COLUMN
WEIGHT_COLUMN : str = "weight"

def plot_nicely(A, result, parameter_space, spec_grid : spectral_grid, spectrum_to_decompose : spectrum):
    result_map = {}
    i = 0
    for (T_eff, FeH, log_g) in parameter_space:
        key = (T_eff, FeH, log_g)
        result_map[key] = i
        i += 1

    hash_map = pd.DataFrame(columns=[TEFF_COLUMN, FEH_COLUMN, LOGG_COLUMN, WEIGHT_COLUMN])
    
    for (T_eff, FeH, log_g) in parameter_space:
        new_row = {TEFF_COLUMN: T_eff, FEH_COLUMN: FeH, LOGG_COLUMN: log_g, WEIGHT_COLUMN: result.x[result_map[(T_eff, FeH, log_g)]]}
        hash_map = pd.concat([hash_map, pd.DataFrame([new_row])], ignore_index=True)

    print(hash_map.sort_values(WEIGHT_COLUMN, ascending=False).head(10).round(3))

    fig, axes = plt.subplots(4, 4, figsize=(15, 15), sharex=True, sharey=True)
    axes = axes.ravel()
    for i, log_g in enumerate(spec_grid.Log_gs):
        subset = hash_map[hash_map[LOGG_COLUMN] == log_g]
        x_vals = [a.value for a in subset[TEFF_COLUMN]]
        y_vals = subset[FEH_COLUMN]
        z_vals = subset[WEIGHT_COLUMN]

        sc = axes[i].scatter(x_vals, y_vals, c=z_vals**.2, cmap='plasma', vmin=0, vmax=1)

        axes[i].set_title(f"log_g={log_g}")
        axes[i].set_xlabel("Temperature / K")
        axes[i].set_ylabel("FeHs / relative to solar")
        # axes[i].set_xticks(np.arange(np.min(T_effs) / u.K, np.max(T_effs) / u.K + 1, 50) * u.K)
        # axes[i].grid()

    STAR_NAME : str = "LTT 3780"
    cbar = fig.colorbar(sc, ax=axes, orientation='vertical', fraction=0.05, pad=0.04)
    cbar.set_label("Weights")
    fig.suptitle(STAR_NAME)
    plt.show()

    plt.clf()

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 8), sharex=True)

    # --- First subplot: spectrum comparison ---
    ax1.set_title(STAR_NAME)

    determined_spectrum = spectrum(
        spectrum_to_decompose.Wavelengths,
        A @ result.x,
        normalised_point=None,
        smoothing_range=None,
        normalise_flux=False
    )

    ax1.plot(
        spectrum_to_decompose.Wavelengths,
        spectrum_to_decompose.Fluxes,
        label="Observational JWST spectrum"
    )
    ax1.plot(
        determined_spectrum.Wavelengths,
        determined_spectrum.Fluxes,
        label="Fitted Spectrum"
    )

    ax1.legend()

    # --- Second subplot: residuals ---
    residual = (determined_spectrum.Fluxes - spectrum_to_decompose.Fluxes) / spectrum_to_decompose.Fluxes

    ax2.plot(spectrum_to_decompose.Wavelengths, residual)
    ax2.set_ylabel(r"Residual = $\frac{\mathrm{Fitted\ Flux}-\mathrm{Observed\ Flux}}{\mathrm{Observed\ Flux}}$")
    ax2.set_xlabel("Wavelength / $\mu$m")

    plt.tight_layout()
    plt.show()

    return hash_map

def get_main_components(hash_map, number_of_components_to_keep : int) -> list[Tuple[Quantity, Quantity, Quantity]]:
    main_components : list[(Quantity, Quantity, Quantity)] = []

    for _, row in hash_map.sort_values(WEIGHT_COLUMN, ascending=False)[0:number_of_components_to_keep].iterrows():
        main_components.append((row[TEFF_COLUMN], row[FEH_COLUMN], row[LOGG_COLUMN]))
    
    return main_components
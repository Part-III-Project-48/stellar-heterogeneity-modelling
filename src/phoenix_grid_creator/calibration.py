"""
there looks to be some offset between the phoenix wavelengths and JWST wavelengths; lets investigate

we are not minimising against the grid: we just pick 1 roughly correct PHOENIX spectra, and np.roll the axis left/right and see where it aligns the best. this offset can then be used when creating the data grid.
"""

from pathlib import Path
from matplotlib import pyplot as plt
import numpy as np
from astropy import units as u
from astropy.visualization import quantity_support
quantity_support()

from spectrum_component_analyser.internals.phoenix_spectrum import phoenix_spectrum
from spectrum_component_analyser.internals.readers import JWST_NORMALISING_POINT, JWST_RESOLUTION, read_JWST_fits
from spectrum_component_analyser.internals.spectral_grid import download_spectrum, get_wavelength_grid
from spectrum_component_analyser.internals.spectrum import spectrum

jwst_spectrum_path : Path = (
    Path(__file__).resolve().parent / Path("../../observed_spectra/MAST_2025-10-26T11_57_04.058Z - LTT-3780/MAST_2025-10-26T11_57_04.058Z/JWST/jw03557004001_04101_00001-seg001_nis_x1dints.fits")
    ).resolve()

jwst_spectrum : spectrum = read_JWST_fits(jwst_spectrum_path, INTEGRATION_INDEX=10, name="LTT-3780")

phoenix_wavelengths = get_wavelength_grid()
# nearby as in: nearby in parameter space to the observed spectrum
nearby_phoenix_spectrum : phoenix_spectrum = download_spectrum(
    3800 * u.K,
    0.0,
    5.0,
    lte=True,
    alphaM=0,
    phoenix_wavelengths=phoenix_wavelengths,
    normalising_point= JWST_NORMALISING_POINT,
    observational_resolution= JWST_RESOLUTION,
    # observational_resolution= JWST_RESOLUTION / 10,
    observational_wavelengths = None,
    name="nearby phoenix spectrum"
)

# maybe its a vacuum vs air wavelength issue? seems reasonable
# nearby_phoenix_spectrum.plot(clear=True, show=False)
# nearby_phoenix_spectrum.Wavelengths = nearby_phoenix_spectrum.air_wavelengths(conversion_method="Griesen2006")
# nearby_phoenix_spectrum.plot(clear=False, show=True)

# plot if needed
# nearby_phoenix_spectrum.plot(clear = True, show = False)
# plt.plot([JWST_NORMALISING_POINT, JWST_NORMALISING_POINT]  * u.um, [np.min(nearby_phoenix_spectrum.Fluxes), np.max(nearby_phoenix_spectrum.Fluxes)] * nearby_phoenix_spectrum.Fluxes.unit)
# jwst_spectrum.plot(clear = False, show = True)

mask = np.isfinite(jwst_spectrum.Fluxes)

roll_resolution = 0.00001 * u.um

def get_MSE(roll : int = 0) -> float:
    # rolled_phoenix_fluxes = np.roll(nearby_phoenix_spectrum.Fluxes, roll)
    # rolled_phoenix_fluxes.Wavelengths = rolled_phoenix_fluxes.Wavelengths.value + 0.01

    # roll to find minimum
    placed_onto_fluxes = np.interp(jwst_spectrum.Wavelengths, nearby_phoenix_spectrum.Wavelengths + roll * roll_resolution, nearby_phoenix_spectrum.Fluxes) # new y = np.interp(new x | old x | old y)

    # mask _after_ rolling & interpolating
    # the jwst data always seems to have some NaN; let's remove them
    placed_onto_fluxes = placed_onto_fluxes[mask]

    residual = (placed_onto_fluxes - jwst_spectrum[mask].Fluxes) / jwst_spectrum[mask].Fluxes
    residual_mean_squared_error = np.sqrt(np.mean(residual**2))

    return residual_mean_squared_error

# how far to move the jwst spectrum to either side to check
max_roll_delta : int = 1000

MSEs = []
rolls = [i for i in reversed(range(-max_roll_delta, max_roll_delta))]

for i in rolls:
    MSEs.append(get_MSE(i))

print(f"10 rolls corresponds to: {roll_resolution * 10}")

plt.clf()
plt.cla()
plt.plot(rolls, MSEs)
plt.show()

# looks to be 0.0042 um off - but fits from main.ipynb suggest moore like 0.0010 (just from eyeing it)
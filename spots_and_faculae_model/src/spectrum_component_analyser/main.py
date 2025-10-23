# for now, just request a random composite spectrum from facula_and_spot_creator
# and try to decompose it - aka can we regenerate the w's

# eventually can read in external data or some training data from a large hdf5 file etc

#external
from astropy.table import QTable
import numpy as np
from matplotlib.lines import Line2D
from matplotlib import pyplot as plt
import scipy as sp

#internal
from facula_and_spot_creator.main import get_example_spectrum, FeH, log_g, normalise_flux
from phoenix_grid_creator.fits_to_hdf5 import FLUX_COLUMN, TEFF_COLUMN, FEH_COLUMN, LOGG_COLUMN, WAVELENGTH_COLUMN, wavelengths
from phoenix_grid_creator.basic_plotter import get_hdf5_data

example_spectrum = get_example_spectrum()

# read in all the available spectra we have (we are assuming FeH and log_g, so this is an easier problem for now)


all_data : QTable = get_hdf5_data()

# print(all_data)

# we are effectively now carrying out the minimisation
# ((\sum w_i f_i) - f_known) to find w_i for a ton of f_i's (our component spectra) and f_known being the total spectrum
# this might be doable just with classical fitting but I guess there might be too much data ? we'll find out ...


b : np.array = np.array(example_spectrum[FLUX_COLUMN])
normalising_constant = sp.integrate.simpson(example_spectrum[FLUX_COLUMN], x=example_spectrum[WAVELENGTH_COLUMN])
# b /= normalising_constant
b = normalise_flux(example_spectrum[WAVELENGTH_COLUMN], b)


# print(b)

A = np.empty((0, 0))

all_data_subset = all_data[(all_data[FEH_COLUMN] == FeH) & (all_data[LOGG_COLUMN] == log_g)]

available_T_effs = sorted(set(all_data_subset[TEFF_COLUMN]))
# columns to be f_column number's and rows to increment the x value
for T_eff in available_T_effs:
	subset = all_data_subset[all_data_subset[TEFF_COLUMN] == T_eff]
	
	flux = subset[FLUX_COLUMN]
	flux = normalise_flux(subset[WAVELENGTH_COLUMN], flux)
	# flux /= sp.integrate.simpson(subset[FLUX_COLUMN], x=subset[WAVELENGTH_COLUMN])
	
	if A.size == 0:
		A = np.empty((len(flux), 0))
	
	A = np.column_stack((A, flux))
 
# print(A)
# assume that w \in [0,1] : but I think this will only be true for real data if normalisation has been done correctly (???)

result = sp.optimize.lsq_linear(A, b, bounds = (0, 1), tol =1e-10, lsmr_tol=1e-15)


print("found weights")
for i, T_eff in enumerate(available_T_effs):
	weight = result.x[i]
	print(f"temperature {T_eff} : {weight}")
	
plt.plot([a.to_value() for a in available_T_effs] , result.x, color="blue", linestyle="-", marker="o")

# alternative fortran method	
# result = sp.optimize.nnls(A, b)


plt.xlabel("Temperature / K")
plt.ylabel("weight / unitless")

# add a fake legend
blue = Line2D([0], [0], label='found weights', marker='s', color="blue", linestyle='')
red = Line2D([0], [0], label='found weights', marker='s',  color="red", linestyle='')
points = [blue, red]
labels = ['found weights', 'input weights']
plt.legend(points, labels)

plt.show()
# this reads in a hdf5 spectral grid of an assumed format and then plots it

from pathlib import Path
import matplotlib.pyplot as plt

from astropy.table import QTable
from astropy import units as u
import numpy as np



# import column names
from fits_to_hdf5 import TEFF_COLUMN, FEH_COLUMN, LOGG_COLUMN, WAVELENGTH_COLUMN, FLUX_COLUMN

script_dir = Path(__file__).resolve().parent
SPECTRAL_GRID_RELATIVE_PATH = Path("../..")
spectral_grid_absolute_path = (script_dir / SPECTRAL_GRID_RELATIVE_PATH).resolve()
print(spectral_grid_absolute_path)
t = QTable.read('data.hdf5', path=f"./")
# print(t.tail(5))

# check for negative wavelength
# we want this to be empty
print(t[t[WAVELENGTH_COLUMN] <= 0])

# # check for nans
bad_nans = np.logical_or.reduce([np.isnan(col) for col in t.itercols()])
# we want this to be 0
print(bad_nans[bad_nans == True].shape[0])

# check some things

print(t.meta)

plt.scatter(t[LOGG_COLUMN], t[FLUX_COLUMN])
plt.show()

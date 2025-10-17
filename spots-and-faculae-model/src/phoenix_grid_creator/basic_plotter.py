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
filename = Path("data.hdf5")
spectral_grid_absolute_path = (script_dir / SPECTRAL_GRID_RELATIVE_PATH / filename).resolve()
# QTable.read.help()
t = QTable.read(spectral_grid_absolute_path, format="hdf5")
# print(t.tail(5))

# some sanity checks - we expect these to return 0#
negative_wavelength_rows = t[t[WAVELENGTH_COLUMN] <= 0]
print(f"number of rows with a negative wavelength: {len(negative_wavelength_rows)}")

bad_nans = np.logical_or.reduce([np.isnan(col) for col in t.itercols()])
print(f"number of rows with a nan: {bad_nans[bad_nans == True].shape[0]}")

print(t.meta)
print(t[WAVELENGTH_COLUMN].unit)

plt.scatter(t[LOGG_COLUMN], t[FLUX_COLUMN])
plt.show()

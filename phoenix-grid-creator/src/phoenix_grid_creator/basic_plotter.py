# this reads in a hdf5 spectral grid of an assumed format and then plots it

from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

TEFF_COLUMN = "T_eff / K"
FEH_COLUMN = "Fe/H / relative to solar"
LOGG_COLUMN = "log_g / log(cm s^(-2))"
WAVELENGTH_COLUMN = "wavelength / angstroms"
FLUX_COLUMN = "flux / counts"

script_dir = Path(__file__).resolve().parent
SPECTRAL_GRID_RELATIVE_PATH = Path("../../data.h5")
spectral_grid_absolute_path = (script_dir / SPECTRAL_GRID_RELATIVE_PATH).resolve()
df = pd.read_hdf(spectral_grid_absolute_path)
print(df.tail(5))

plt.scatter(df[LOGG_COLUMN], df[FLUX_COLUMN])
plt.show()

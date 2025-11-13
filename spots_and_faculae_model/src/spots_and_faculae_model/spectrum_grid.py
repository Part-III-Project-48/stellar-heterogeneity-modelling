from pathlib import Path
from astropy.table import QTable

from phoenix_grid_creator.fits_to_hdf5 import TEFF_COLUMN, FEH_COLUMN, LOGG_COLUMN, WAVELENGTH_COLUMN, FLUX_COLUMN
    
class spectrum_grid:
    """
    a wrapper for an astropy qtable which stores PHOENIX spectrum data.
    
    typically will be convolved so that its resolution matches some known telescope
    """

    def __init__(self, table : QTable = None):
        """default initialiser from a qtable object"""
        self.Table = table

    @classmethod
    def from_hdf5_file(cls, absolute_hdf5_path : Path) -> None:
        """Alternative constructor from a hdf5 path"""
        table = QTable.read(absolute_hdf5_path, format="hdf5")
        return cls(table)

    def save(self, absolute_path : Path = "default_path", name : str = "spectrum_grid"):
        self.Table.write(name, path = absolute_path, serialize_meta=True, overwrite=True)
    
    # all columns should have units
    @property
    def T_effs(self):
        return self.Table[TEFF_COLUMN]

    @property
    def FeHs(self):
        return self.Table[FEH_COLUMN]

    @property
    def log_gs(self):
        return self.Table[LOGG_COLUMN]

    @property
    def wavelengths(self):
        return self.Table[WAVELENGTH_COLUMN]

    @property
    def fluxes(self):
        return self.Table[FLUX_COLUMN]
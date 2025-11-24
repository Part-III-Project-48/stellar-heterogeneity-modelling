from itertools import product
from pathlib import Path
from astropy.table import QTable
import numpy as np
import pandas as pd
from tqdm import tqdm
from scipy.interpolate import interp1d
import astropy
from astropy.units import Quantity
import astropy.units as u

from spots_and_faculae_model.spectrum import spectrum

# units should be stored in the astropy quantity anyway
# changing these is fine, as long as a new spectral grid is created which uses these column names
TEFF_COLUMN = "T_eff / K"
FEH_COLUMN = "Fe/H / relative to solar"
LOGG_COLUMN = "log_g / log(cm s^(-2))"
WAVELENGTH_COLUMN = "wavelength / angstroms"
FLUX_COLUMN = "flux / erg / (s * cm**2 * cm)"
SPECTRUM_COLUMN : str = "spectrum object"
	
class new_spectrum_grid:
	"""
	a wrapper for a pandas dataframe which stores PHOENIX spectrum data in the format
	T_eff | FeH | Log_g | spectrum object containing wavelength & fluxes
	
	currently, the PHOENIX HDF5 data is convolved so that its resolution matches some known telescope
	"""

	def __init__(self, table : QTable = None, name : str = None):
		"""default initialiser from a qtable object"""
		self.Name = name
		self.Table = table
		self.FancyTable = pd.DataFrame(columns=[TEFF_COLUMN, FEH_COLUMN, LOGG_COLUMN, SPECTRUM_COLUMN])
		self.update_fancy_table()
			
	def update_fancy_table(self):
		"""
		updates the fancy table to align to the information in self.Table
		WARNING: this deletes any alterations made to the spectra (although changing spectra is not preferred anyway)
		"""

		self.FancyTable = pd.DataFrame(columns=[TEFF_COLUMN, FEH_COLUMN, LOGG_COLUMN, SPECTRUM_COLUMN])
		
		for T_eff, FeH, log_g in tqdm(product(self.T_effs, self.FeHs, self.log_gs), total= len(self.T_effs) * len(self.FeHs) * len(self.log_gs), desc="Creating fancier spectral grid..."):
			subset = self.Table[(self.Table[TEFF_COLUMN] == T_eff) &
								(self.Table[FEH_COLUMN] == FeH) &
								(self.Table[LOGG_COLUMN] == log_g)]
			
			spectrum_name : str = f"simulated phoenix spectrum for {self.Name}" if self.Name != None else f"simulated phoenix spectrum"
			spec : spectrum = spectrum(wavelengths = subset[WAVELENGTH_COLUMN], fluxes = subset[FLUX_COLUMN], name=spectrum_name)
			
			new_row = pd.DataFrame({
				TEFF_COLUMN : [T_eff],
				FEH_COLUMN : [FeH],
				LOGG_COLUMN : [log_g],
				SPECTRUM_COLUMN : [spec]
			})

			self.FancyTable = pd.concat([self.FancyTable, new_row])
			
		# pandas tables can't save their metadata into a HDF5 directly (can use HDFStore or smthn) - but astropy tables can have metadata, units etc. so lets convert to an astropy table

		# add astropy units to columns (this will be stored in metadata and can be read back out into an astropy QTable)
		# self.Table[TEFF_COLUMN].unit = u.Kelvin
		self.Table[TEFF_COLUMN].desc = "effective surface temperature"

		# self.Table[FEH_COLUMN].unit = u.dimensionless_unscaled
		self.Table[FEH_COLUMN].desc = "relative to solar metallacity"

		# astropy seems to have a hard time reading in log quantities from hdf5 files. so lets just save this as unitless
		# self.Table[LOGG_COLUMN].unit = u.dimensionless_unscaled
		self.Table[LOGG_COLUMN].desc = "log_10(u.cm * u.second**(-2)) of the surface gravity"

		# self.Table[WAVELENGTH_COLUMN].unit = u.Angstrom

		# self.Table[FLUX_COLUMN].unit = u.dimensionless_unscaled
		self.Table[FLUX_COLUMN].desc = "in erg s^-1 cm^-2 cm^-1"

	@classmethod
	def from_hdf5_file(cls, absolute_hdf5_path : Path):
		"""Alternative constructor from a hdf5 path"""
		table = QTable.read(absolute_hdf5_path, format="hdf5")
		return cls(table)

	@classmethod
	def from_arrays(cls, T_effs : np.array, FeHs : np.array, log_gs : np.array, wavelengths : np.array, fluxes : np.array):
		table = QTable([
			T_effs,
			FeHs,
			log_gs,
			wavelengths,
			fluxes
			],
		names=(TEFF_COLUMN, FEH_COLUMN, LOGG_COLUMN, WAVELENGTH_COLUMN, FLUX_COLUMN))
		return cls(table)

	# def save(self, absolute_path : Path = "default_path", name : str = "spectrum_grid", Overwrite : bool = False):
	# 	print("[PHOENIX GRID CREATOR] : writing dataframe to hdf5...")
	# 	self.Table.write(name, path = absolute_path, serialize_meta=True, overwrite=Overwrite)
	# 	print("[PHOENIX GRID CREATOR] : hdf5 saving complete")

	def regularise_temperatures(self, new_T_effs : np.array) -> None:
		"""
		interpolates the data onto a given T_eff array.

		Attributes:
		new_T_effs (np.array): temperatures to interpolate the data onto
		"""
		regularised_wavelength_df = pd.DataFrame(columns=[TEFF_COLUMN, FEH_COLUMN, LOGG_COLUMN, WAVELENGTH_COLUMN, FLUX_COLUMN])
		
		for FeH, log_g, new_T_eff in tqdm(product(self.FeHs, self.log_gs, new_T_effs), total= len(self.FeHs) * len(self.log_gs) * len(new_T_effs), desc="Regularising temperature points"):
			
			# this df is at all temperatures
			subset = self.Table[(self.Table[FEH_COLUMN] == FeH) & (self.Table[LOGG_COLUMN] == log_g)]

			# aka subset df represents 3D data which maps (wavelength to flux) over a set of temperatures
			# we want the wavelength to flux map at a different temperature; namely at new_T_eff
			
			# so we want to linearly interpolate every flux between T_1 and T_2 at all wavelengths
			# aka. np.interp(new_T_eff, subset_df[TEMPERATURE_COLUMN], subset_df[FLUX_COLUMN]) # assuming this vectorises and returns a map from all wavelengths to fluxes
			
			pivoted = subset.pivot(index=TEFF_COLUMN, columns=WAVELENGTH_COLUMN, values=FLUX_COLUMN)
			x = pivoted.index.to_numpy()               # shape (n_temperatures,)
			wavelengths = pivoted.columns.to_numpy()   # shape (n_wavelengths,)
			y = pivoted.values
			
			f = interp1d(x, y, axis=0, kind="cubic")
			
			wavelength_to_flux_map_at_new_T_eff = f(new_T_eff)
			
			temp_df = pd.DataFrame({
				TEFF_COLUMN : new_T_eff,
				FEH_COLUMN : FeH,
				LOGG_COLUMN : log_g,
				WAVELENGTH_COLUMN : wavelengths,
				FLUX_COLUMN : wavelength_to_flux_map_at_new_T_eff # interpolated flux function between previous and next temperatures 
			})
			
			# avoid warning about concat-ing an empty df
			if not regularised_wavelength_df.empty:
				# our df index has no meaningful meaning, and sort I think just ensures the columns are in the correct order or something?
				regularised_wavelength_df = pd.concat([regularised_wavelength_df, temp_df], ignore_index=True)#, sort=True)
			else:
				regularised_wavelength_df = temp_df
		
		self.Table = regularised_wavelength_df
		self.update_fancy_table()
	


	# spectrum_to_decompose is temporary here; it can probably be changed to be self.wavelengths but that isn't tested
	def process_single_spectral_component(self, T_eff : Quantity[u.K], FeH : float, log_g : float, mask : np.array, spectrum_to_decompose : spectrum) -> np.array:
		"""
		returns a np array of astropy quantities with units of Janskys

		Attributes
		----------
		mask : np.array
			something of length of spectra stored in the spectral grid. contains False where infinities were found in the observed spectra to fit


		"""
		subset_table = self.Table[(self.Table[TEFF_COLUMN] == T_eff) &
						  (self.Table[FEH_COLUMN] == FeH) &
						  (self.Table[LOGG_COLUMN] == log_g)]
		
		subset : new_spectrum_grid = new_spectrum_grid(subset_table)
		
		# remove the indices that were nan in the spectrum
		# must be in the same order as we did for the spectrum_to_decompose
		subset.Table = subset.Table[mask] # should make the table's spectra have the same x axis (wavelengths) as the spectrum_to_decompose 
		subset_spectrum = spectrum.from_phoenix_units(wavelengths=spectrum_to_decompose.Wavelengths, phoenix_fluxes=subset.Table[FLUX_COLUMN])
		subset_spectrum.normalise_Janskys()
		# subset_spectrum = subset_spectrum[(1.25 * u.um <= subset_spectrum.Wavelengths) & (subset_spectrum.Wavelengths <= 2 * u.um)]

		return subset_spectrum.Fluxes
	
	def new_process_single_spectral_component(self, T_eff : Quantity[u.K], FeH : float, log_g : float, mask : np.array) -> np.array:
		subset_table = self.FancyTable[(self.FancyTable[TEFF_COLUMN] == T_eff) &
					(self.FancyTable[FEH_COLUMN] == FeH) &
					(self.FancyTable[LOGG_COLUMN] == log_g)]
		
		if len(subset_table[SPECTRUM_COLUMN] != 0):
			raise LookupError(f"new_spectrum_grid is degenerate; it contains multiple spectra for the parameters T_eff = {T_eff}, FeH = {FeH}, log_g = {log_g}")
		
		spec : spectrum = subset_table[SPECTRUM_COLUMN][0]
		return spec.Fluxes[mask]
	
	# just to expose stuff
	# these give the list of unique T_effs : if you wanna do slicing etc you'll need the whole self.Table object
	# this might be slow / a bottleneck
	@property
	def T_effs(self):
		return astropy.table.unique(self.Table, keys=[TEFF_COLUMN])[TEFF_COLUMN]

	@property
	def FeHs(self):
		return astropy.table.unique(self.Table, keys=[FEH_COLUMN])[FEH_COLUMN]

	@property
	def log_gs(self):
		return astropy.table.unique(self.Table, keys=[LOGG_COLUMN])[LOGG_COLUMN]
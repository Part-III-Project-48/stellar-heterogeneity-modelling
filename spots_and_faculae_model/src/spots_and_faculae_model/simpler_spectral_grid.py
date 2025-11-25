import datetime
from pathlib import Path
from typing import Sequence
from astropy.units import Quantity
import numpy as np
import astropy.units as u
from io import BytesIO
from joblib import Parallel, delayed
import requests
from tqdm import tqdm
import numpy as np
from astropy.io import fits
from astropy import units as u
from pathlib import Path
import datetime
import scipy as sp
import h5py
from typing import Sequence
import specutils

# internal imports
from phoenix_grid_creator.PHOENIX_filename_conventions import *
from spots_and_faculae_model.phoenix_spectrum import phoenix_spectrum
import h5py

PHOENIX_FLUX_UNITS = u.erg / (u.s * u.cm**2 * u.cm)

UNIT_METADATA_NAME : str = "units"
WAVELENGTH_DATASET_NAME : str = "wavelengths"
TEFF_DATASET_NAME : str = "Teff"
FEH_DATASET_NAME : str = "FeH"
LOGG_DATASET_NAME : str = "log_g"
FLUX_DATASET_NAME : str = "fluxes"

MAIN_GRID_NAME : str = "main_grid"

USES_REGULARISED_WAVELENGTHS_METADATA_NAME : str = "includes interpolated wavelengths?"
USES_REGULARISED_TEMPERATURES_METADATA_NAME : str = "includes interpolated temperatures?"

def get_wavelength_grid() -> Sequence[Quantity]:
	"""
	returns 1D array consisting of astropy Quantities of dimension length
	"""
	# read in the wavelength (1D) grid so we can save this into our mega-grid correctly

	script_dir = Path(__file__).resolve().parent
	WAVELENGTH_GRID_RELATIVE_PATH = Path("../../assets/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits")
	wavelength_grid_absolute_path = (script_dir / WAVELENGTH_GRID_RELATIVE_PATH).resolve()

	if not wavelength_grid_absolute_path.exists():
		raise FileNotFoundError(f"wavelength grid file not found at : {wavelength_grid_absolute_path}.")

	with fits.open(wavelength_grid_absolute_path) as hdul:
		# there is only 1 HDU in this wavelength grid file
		WAVELENGTH_GRID_HDU_INDEX = 0
		wavelengths = hdul[WAVELENGTH_GRID_HDU_INDEX].data
		# the fits file is big endian; pandas requires little endian. this swaps between them
		wavelengths = wavelengths.byteswap().view(wavelengths.dtype.newbyteorder())
		wavelengths *= u.Angstrom
		print("[PHOENIX GRID CREATOR] : wavelength grid found & loaded in")

	return wavelengths
		
def download_spectrum(T_eff, FeH, log_g, lte : bool, alphaM : float, phoenix_wavelengths : np.array, regularised_wavelengths : np.array = None) -> phoenix_spectrum:
	"""
	well use this function to parallelise getting the spectra
	"""
	try:
		file = get_file_name(lte, T_eff, log_g, FeH, alphaM)
		url = get_url(file)
	except ValueError as e:
		tqdm.write(f"[PHOENIX GRID CREATOR] : filename or urlname error: {e}. continuing onto next requested spectrum")
		return
	try:
		response = requests.get(url)
		response.raise_for_status()
	except requests.exceptions.HTTPError as e:
		tqdm.write(f"[PHOENIX GRID CREATOR] : HTTPError raised with the following parameters.\nlte: {lte}\nT_eff={T_eff}\nlog_g={log_g}\nFeH={FeH}\nalphaM={alphaM}")
		tqdm.write(f"url = {url}")
		tqdm.write("\n continuing with the next file...")
		return
	
	# the index of the header data unit the data we want is in (looks to be 0 being the spectra, and 1 being the abundances, and those are the only 2 HDUs in the .fits files)
	SPECTRA_HDU_INDEX = 0

	with fits.open(BytesIO(response.content)) as hdul:
		# hdul.info()
		
		fluxes = hdul[SPECTRA_HDU_INDEX].data

		# for some reason, the fits file is big-endian; pandas required little-endian
		fluxes = fluxes.byteswap().view(fluxes.dtype.newbyteorder())
		fluxes *= PHOENIX_FLUX_UNITS
		JWST_resolution = .001 * u.um

		if regularised_wavelengths != None:
			index_per_wavelength = len(phoenix_wavelengths) / (np.max(phoenix_wavelengths) - np.min(phoenix_wavelengths))
			index_per_wavelength = index_per_wavelength.to(u.um**-1)
			convolution_range = int((index_per_wavelength * JWST_resolution)) # in number of adjacent points to consider
			fluxes = sp.ndimage.gaussian_filter(fluxes, convolution_range) * fluxes.unit # gaussian_filter seems to remove units
			spec = phoenix_spectrum(
				wavelengths=regularised_wavelengths,
				fluxes = np.interp(regularised_wavelengths, phoenix_wavelengths, fluxes),
				t_eff=T_eff,
				feh=FeH,
				log_g=log_g)
		else:
			spec = phoenix_spectrum(
					wavelengths=phoenix_wavelengths,
					fluxes=fluxes,
					t_eff=T_eff,
					feh=FeH,
					log_g=log_g)

		return spec

class simpler_spectral_grid():
	"""
	this is an internal class really: its much more human readable to just use list[spectrum]; this is just for saving to / loading from a hdf5 file
	"""
	def __init__(self, wavelengths : np.array, t_effs : np.array, fehs : np.array, log_gs : np.array, fluxes : np.array, uses_regularised_wavelengths : bool, uses_regularised_temperatures : bool):
		"""
		don't use this init: use the other wrappers that download things or load in from a hdf5 file
		(the structure of fluxes is non trivial)
		"""

		# 1D arrays
		self.Wavelengths = wavelengths
		self.T_effs = t_effs
		self.FeHs = fehs
		self.Log_gs = log_gs

		# this is a 4D array
		self.Fluxes = fluxes

		self.Uses_Regularised_Wavelengths = uses_regularised_wavelengths
		self.Uses_Regularised_Temperatures = uses_regularised_temperatures

	@classmethod
	def from_internet(cls, T_effs, FeHs, log_gs, regularised_wavelengths, alphaM = 0, lte = True, regularised_temperatures = None):
		"""
		seems like the only data I can find is LTE data (?)
		"""

		if regularised_temperatures != None:
			raise NotImplementedError("regularising temperatures is not added into simpler spectral grid atm")
		
		phoenix_wavelengths = get_wavelength_grid()
		
		# effectively redundant, as spectrum class has a built in converter now anyway. maybe delete this logic
		# if CONVERT_WAVELENGTHS_TO_AIR:
		# 	phoenix_wavelengths = specutils.utils.wcs_utils.vac_to_air(phoenix_wavelengths)
		
		def fetch_spectra_and_indices(i, j, k, T_eff, FeH, log_g):
			spec : phoenix_spectrum = download_spectrum(T_eff, FeH, log_g, lte, alphaM, phoenix_wavelengths, regularised_wavelengths)
			return i, j, k, spec
		
		# pre - allocate 4d flux array
		fluxes = np.zeros((len(T_effs), len(FeHs), len(log_gs), len(regularised_wavelengths)))
			
		tasks = [
			(i, j, k, t, f, g)
			for i, t in enumerate(T_effs)
			for j, f in enumerate(FeHs)
			for k, g in enumerate(log_gs)
		]

		results = Parallel(n_jobs=-1, prefer="threads")(
			delayed(fetch_spectra_and_indices)(*task) for task in tqdm(tasks)
		)
		
		spec : phoenix_spectrum
		for i, j, k, spec in results:
			# i think this removes units from fluxes silently - as this is some 4D array. maybe we can readd them somehow; doesnt rly matter for now though
			fluxes[i, j, k, :] = spec.Fluxes

		# assume all spectra have the same wavelengths
		return cls(spec.Wavelengths, T_effs, FeHs, log_gs, fluxes, regularised_wavelengths != None, regularised_temperatures != None)

	def save(self, absolute_path : Path, overwrite : bool):
		if absolute_path.exists() and not overwrite:
			raise FileExistsError(f"specified path already exists; and overwrite is set to false. Change the file name or turn on overwrite. (File path: {absolute_path})")
		
		with h5py.File(absolute_path, "w") as f:
			f.attrs["creator"] = "Ben Green"
			f.attrs["description"] = "Collection of synthetic spectra from PHOENIX dataset"
			f.attrs["version"] = "0.1"
			f.attrs["date"] = str(datetime.datetime.now())
			f.attrs["notes"] = "All spectra share the same wavelength grid."

			# add some metadata to the QTable e.g. (wavelength medium = vacuum, source, date)
			f.attrs["wavelength medium"] = "vacuum"
			f.attrs["source"] = "https://phoenix.astro.physik.uni-goettingen.de/data/"

			f.attrs[USES_REGULARISED_WAVELENGTHS_METADATA_NAME] = self.Uses_Regularised_Wavelengths
			f.attrs[USES_REGULARISED_TEMPERATURES_METADATA_NAME] = self.Uses_Regularised_Temperatures

			g = f.create_group(MAIN_GRID_NAME)

			# we have to remove units from our np arrays and write them to metadata to be retrieved later

			wavelength_unit = self.Wavelengths.unit
			T_eff_unit = self.T_effs.unit
			flux_unit = PHOENIX_FLUX_UNITS

			wavelength_dataset = g.create_dataset(WAVELENGTH_DATASET_NAME, data=np.array(self.Wavelengths))
			wavelength_dataset.attrs[UNIT_METADATA_NAME] = str(wavelength_unit)

			T_eff_dataset = g.create_dataset(TEFF_DATASET_NAME, data=np.array([i.value for i in self.T_effs]))
			T_eff_dataset.attrs[UNIT_METADATA_NAME] = str(T_eff_unit)

			FeH_dataset = g.create_dataset(FEH_DATASET_NAME, data=self.FeHs)
			FeH_dataset.attrs[UNIT_METADATA_NAME] = str(u.dimensionless_unscaled)

			log_g_dataset = g.create_dataset(LOGG_DATASET_NAME, data=self.Log_gs)
			log_g_dataset.attrs[UNIT_METADATA_NAME] = str(u.dimensionless_unscaled)

			flux_dataset = g.create_dataset(FLUX_DATASET_NAME, data=np.array(self.Fluxes))
			flux_dataset.attrs[UNIT_METADATA_NAME] = str(flux_unit)

		print("[PHOENIX GRID CREATOR] : hdf5 saving complete")

	

	@classmethod
	def from_hdf5(cls, absolute_path : Path):
		from astropy.units import Unit
		
		with h5py.File(absolute_path, "r") as f:
			main_grid = f[MAIN_GRID_NAME]

			wavelengths = main_grid[WAVELENGTH_DATASET_NAME]
			wavelengths = np.array(wavelengths) * Unit(wavelengths.attrs[UNIT_METADATA_NAME], parse_strict='raise')

			T_effs = main_grid[TEFF_DATASET_NAME]
			T_effs = np.array(T_effs) * Unit(T_effs.attrs[UNIT_METADATA_NAME], parse_strict='raise')

			FeHs = main_grid[FEH_DATASET_NAME]
			FeHs = np.array(FeHs) * Unit(FeHs.attrs[UNIT_METADATA_NAME], parse_strict='raise')

			log_gs = main_grid[LOGG_DATASET_NAME]
			log_gs = np.array(log_gs) * Unit(log_gs.attrs[UNIT_METADATA_NAME], parse_strict='raise')

			fluxes = main_grid[FLUX_DATASET_NAME]
			fluxes = np.array(fluxes) * Unit(fluxes.attrs[UNIT_METADATA_NAME], parse_strict='raise')

			uses_regularised_wavelengths = f.attrs[USES_REGULARISED_WAVELENGTHS_METADATA_NAME]

			uses_regularised_temperatures = f.attrs[USES_REGULARISED_TEMPERATURES_METADATA_NAME]
			
		return cls(wavelengths, T_effs, FeHs, log_gs, fluxes, uses_regularised_wavelengths, uses_regularised_temperatures)
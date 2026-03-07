from enum import Enum
from pathlib import Path

import matplotlib.pyplot as plt
from tqdm import tqdm
import typer

from spectrum_component_analyser.internals.readers.JWST.instruments import Instrument, NIRISS, NIRSPEC 
from spectrum_component_analyser.internals.readers.JWST.file_reader import JWSTFileReader
from spectrum_component_analyser.internals.spectrum import spectrum
from internal_constants import package_path

from vanity.printing_colors import *

class JWSTTargets(Enum):
	LTT3780 = "LTT 3780"
	LTT144A = "LTT 1445A"
	K218 = "K2-18"
	TRAPPIST1 = "TRAPPIST-1"

class JWSTFolderReader():
	@staticmethod
	def get_all_spectra(target : JWSTTargets, instrument : Instrument) -> list[spectrum]:
		"""
		returns all spectra from all .fits files contained in fits_directory (specification is in get_file_path() and also a markdown file)
		"""

		def get_file_path(target : JWSTTargets, instrument : Instrument) -> Path:
			return Path(f"{package_path}/neater_observed_spectra_folder/{target.value}/{instrument.FolderPath}/")

		fits_directory : Path = get_file_path(target, instrument)

		print(fits_directory)

		all_spectra : list[spectrum] = []

		files = list(fits_directory.glob("*.fits"))

		for path in tqdm(files, desc=f"Loading in JWST *.fits from {LIGHT_BLUE}{fits_directory}{RESET}"):
			all_spectra.extend(
				JWSTFileReader.get_all_spectra(
					file_path=path,
					instrument=instrument,
					verbose=False
				)
			)

		return all_spectra

def main() -> None:
	niriss_spectra = []
	nirspec_spectra = []
	niriss_spectra = JWSTFolderReader.get_all_spectra(JWSTTargets.TRAPPIST1, NIRISS)
	# nirspec_spectra = read_all_JWST_fits(JWSTTargets.LTT3780, NIRSPEC)

	all_spectra = [*nirspec_spectra, *niriss_spectra]

	s : spectrum
	for s in all_spectra:
		s.plot(clear=False)
	
	plt.show()

if __name__ == "__main__":
	typer.run(main)
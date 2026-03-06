from enum import Enum
from pathlib import Path

import matplotlib.pyplot as plt
from tqdm import tqdm
import typer

from spectrum_component_analyser.internals.readers.JWST.instruments import Instrument, NIRISS, NIRSPEC 
from spectrum_component_analyser.internals.readers.JWST.reader import JWSTReader
from spectrum_component_analyser.internals.spectrum import spectrum

from vanity.printing_colors import *

class JWSTTargets(Enum):
	LTT3780 = "LTT 3780"

def read_all_JWST_fits(target : JWSTTargets, instrument : Instrument) -> list[spectrum]:
	"""
	returns all spectra from all .fits files contained in fits_directory (specification is in a markdown file)
	"""
	fits_directory : Path = Path(f"neater_observed_spectra_folder/{target.value}/{instrument.FolderPath}/")

	all_spectra : list[spectrum] = []

	files = list(fits_directory.glob("*.fits"))

	for path in tqdm(files, desc=f"Loading in JWST *.fits from {LIGHT_BLUE}{fits_directory}{RESET}"):
		all_spectra.extend(
			JWSTReader().get_all_spectra(
				file_path=path,
				instrument=instrument,
				verbose=False
			)
		)

	return all_spectra

def main() -> None:
	niriss_spectra = read_all_JWST_fits(JWSTTargets.LTT3780, NIRISS)
	nirspec_spectra = read_all_JWST_fits(JWSTTargets.LTT3780, NIRSPEC)
	niriss_spectra = []

	all_spectra = [*nirspec_spectra, *niriss_spectra]

	s : spectrum
	for s in all_spectra:
		s.plot(clear=False)
	
	plt.show()

if __name__ == "__main__":
	typer.run(main)
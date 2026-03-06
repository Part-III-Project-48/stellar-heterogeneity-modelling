from enum import Enum
from pathlib import Path

import matplotlib.pyplot as plt
from tqdm import tqdm
import typer

from spectrum_component_analyser.internals.readers.JWST.NIRISS import NIRISSReader
from spectrum_component_analyser.internals.spectrum import spectrum

from vanity.printing_colors import *

class JWSTInstrument(Enum):
	NIRISS = "NIRISS"
	NIRSPEC = "NIRSPEC"

class JWSTTargets(Enum):
	LTT3780 = "LTT 3780"

def read_JWST_fits(target : JWSTTargets, instrument : JWSTInstrument) -> spectrum:
	fits_directory : Path = Path(f"neater_observed_spectra_folder/{target.value}/{instrument.value}/")

	all_spectra : list[spectrum] = []

	files = list(fits_directory.glob("*.fits"))

	for path in tqdm(files, desc=f"Loading in JWST *.fits from {LIGHT_BLUE}{fits_directory}{RESET}"):
		all_spectra.extend(NIRISSReader().get_all_spectra(path))
	
	s : spectrum
	for s in all_spectra:
		s.plot(clear=False)
	
	plt.show()
	

def main() -> None:
	read_JWST_fits(JWSTTargets.LTT3780, JWSTInstrument.NIRISS)

if __name__ == "__main__":
	typer.run(main)
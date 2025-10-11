# this is going to download our data that meets some desired data ranges
# this file should just be run once when you want to create the data grid
# this file defines the filename formats used by PHOENIX in a readable way

import requests
from tqdm import tqdm

import numpy as np

import astropy as ap
from astropy import units as u
from itertools import product

# working example
# url = "https://phoenix.astro.physik.uni-goettingen.de/data/HiResFITS/PHOENIX-ACES-AGSS-COND-2011/Z-0.0/lte06000-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits"
# filename = "lte06000-4.50-0.0.PHOENIX.fits"

GRID : str = "PHOENIX-ACES-AGSS-COND-2011"

def get_file_name(lte : bool,
				  T_eff : int,
				  log_g : float,
				  FeH : float,
				  alphaM : float) -> str:
	
	### --- check for step and range validity --- ###
	
	if not ((2_300 <= T_eff <= 7_000 and T_eff % 100 == 0) or (7_000 <= T_eff <= 12000 and T_eff % 200 == 0)):
		raise ValueError("T_eff value not in range or of incorrect step")
	
	if not (0.0 <= log_g <= 6.0 and log_g % 0.5 == 0):
		raise ValueError("log_g value not in range or of incorrect step")
	
	if not ((-4.0 <= FeH <= -2.0 and FeH % 1.0 == 0) or (-2.0 <= FeH <= 1.0 and FeH % 0.5 == 0)):
		raise ValueError("FeH value not in range or of incorrect step")
	
	if not (-0.2 <= alphaM <= 1.2 and alphaM % 0.2 == 0):
		raise ValueError("alphaM value not in range or of incorrect step")
	
	### --- check for data availability here --- ###
	
	if (not -3.0 <= FeH <= 0.0) and alphaM != 0:
		raise ValueError("Alpha element abundances [α/M]≠0 are available for -3.0≤[Fe/H]≤0.0. only.")
	
	### --- force data into strings of the correct format --- ###
	
	NLTE_OR_LTE : str = "lte" if lte else "nlte"

	# always 5 digits with leading 0s
	T_eff = f"{T_eff:05d}"

	# 2 d.p, unsigned
	log_g = f"{log_g:.2f}"

	# 1 d.p., signed (and 0.0 is rendered as -0.0, not +0.0)
	FeH = -0.0 if FeH == 0 else FeH
	FeH = f"{FeH:+.1f}"

	# 2d.p.; signed
	alphaM = f"{alphaM:+.2f}"
	
	### --- define the subgrid name --- ###
	
	# describes FeH and alphaM
	SUBGRID : str
	
	# spectra with alphaM = 0 have the alphaM omitted from their filename
	if (float(alphaM) != 0):
		SUBGRID = f"{FeH}.Alpha={alphaM}"
	else:
		SUBGRID = f"{FeH}"
	
	### --- define the filename --- ###
	
	file_name : str = f"{GRID}/Z{SUBGRID}/{NLTE_OR_LTE}{T_eff}-{log_g}{SUBGRID}.{GRID}-HiRes.fits"
	return file_name

def get_url(file_name : str) -> str:
	url = f"https://phoenix.astro.physik.uni-goettingen.de/data/HiResFITS/{file_name}"
	return url

T_effs = np.arange(2300, 4000, 100)
FeHs = np.arange(-0.5, 0.5, 0.5)
log_gs = np.arange(3.5, 5.5, 0.5)
alphaM = 0
lte : bool = True

file_name_to_save : str = "example.fits"

total_number_of_files : int = len(T_effs) * len(FeHs) * len(log_gs)

for T_eff, FeH, log_g in tqdm(product(T_effs, FeHs, log_gs), total=total_number_of_files, desc="Downloading spectra"):
	file = get_file_name(lte, T_eff, log_g, FeH, alphaM)
	url = get_url(file)
	
	try:
		response = requests.get(url)
		response.raise_for_status()
	except requests.exceptions.HTTPError as e:
		print("---")
		print(f"HTTPError raised with the following parameters.\nlte: {lte}\nT_eff={T_eff}\nlog_g={log_g}\nFeH={FeH}\nalphaM={alphaM}\n continuing with the next file")
		print(f"url = {url}")
		print("---")
		continue

	with open("example.fits", "wb") as f:
		f.write(response.content)

# with ap.io.fits.open(fits_image_filename) as hdul:

#     hdul.info()

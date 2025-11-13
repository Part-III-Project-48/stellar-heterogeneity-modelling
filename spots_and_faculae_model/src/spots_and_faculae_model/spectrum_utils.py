import numpy as np
import astropy.units as u
import scipy as sp

def normalise_Janskys(wavelengths : np.array, counts : np.array, normalised_point = 2.2 * u.um, smoothing_range = 0.5 * u.um) -> np.array:
	"""
	this will fail if wavelengths does not span at least smoothing_range
	
	the canonical way to normalise spectra is to choose a portion that's continuum and make that be at a consistent scale
	
	inputs should both be lists of astropy quantities (aka values with astropy units)
	
	I tried doing this but the continuum fit was terrible: https://specutils.readthedocs.io/en/stable/fitting.html#continuum-fitting
	"""
	
	# kernel size of about 501 with 9999 points between 5 and 15 um seemed good - this range corresponds (roughly) to that 
	wavelengths_in_range = wavelengths[(wavelengths[0] <= wavelengths) & (wavelengths <= wavelengths[0] + smoothing_range)]
	kernel_size = len(wavelengths_in_range)
	if kernel_size % 2 == 0:
		kernel_size +=1
		
	# smooth to make sure there's no spikes
	counts = np.array(counts, dtype=np.float64)
	smoothed_counts = sp.signal.medfilt(counts, kernel_size=[kernel_size])
	# normalise the counts at normalised_point (or next nearest value) to be 1
	counts /= smoothed_counts[(normalised_point <= wavelengths)][0]
	
	return counts
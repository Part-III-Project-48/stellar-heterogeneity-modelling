# this is going to download our data that meets some desired data ranges
# this file should just be run once when you want to create the data grid

import requests

import astropy as ap
from astropy import units as u

with ap.io.fits.open(fits_image_filename) as hdul:

    hdul.info()
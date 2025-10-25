# Setting up the dependencies

To use this folder you'll need to download poetry: https://python-poetry.org/docs/

To install python dependencies use:
```bash
poetry install
```

This creates a virtual env in some OS-dependent config file (see the poetry docs above) which you can find using `poetry config --list` and select in vscode etc. Or you can change the venv to be local to this folder by doing

```bash
poetry config virtualenvs.in-project true

# then make the venv appear:
poetry env list  # shows the name of the current environment
poetry env remove <current environment>
poetry install  # will create a new environment using your updated configuration
```

Then you should just be able to run the python file in whichever way you like.

Using poetry means that the python dependencies are fully specified and don't break in the future.

# Downloading and using the wavelength grid

To create the HDF5 file, you need to download the wavelength grid into the assets folder. The wavelength grid can be downloaded using (assuming you start in the top-level of the repo):

```bash
cd spots_and_faculae_model/assets
wget https://phoenix.astro.physik.uni-goettingen.de/data/HiResFITS/WAVE_PHOENIX-ACES-AGSS-COND-2011.fits
```

This .fits contains 1x HDU which contains 1D data. The length of this 1D array is the same as the length for all the spectral data; and the values in this array is the corresponding wavelength for the index. i.e. this file maps index -> wavelength (angstroms)

Then to create the HDF5 file, just run `spots_and_faculae_model/src/phoenix_grid_creator/fits_to_hdf5.py`.

# Downloading other .fits files

All .fits files are gitignored currently. This means that the wavelength grid is NOT stored in this repo and must be downloaded locally (e.g. into assets/) before recreating the HDF5 spectra grid.

All other spectra downloads are performed automatically by fits_to_hdf5.py (currently), and the file format conventions for the spectra can just be read directly from within PHOENIX_filename_conventions.py.

working example url format:
https://phoenix.astro.physik.uni-goettingen.de/data/HiResFITS/PHOENIX-ACES-AGSS-COND-2011/Z-0.0/lte06000-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits

# Notes

- the only wavelength grid I can find is WAVE (not the AWAV) which is specified as VACUUM wavelengths. astropy specutils is used to convert from vacuum to air wavelength (this apparently has a low error)
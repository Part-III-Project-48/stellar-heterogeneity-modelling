# Stellar Spectrum Decomposition using PHOENIX

Cambridge Part III Project 48

- PHOENIX homepage : [link](https://phoenix.astro.physik.uni-goettingen.de/)
- PHOENIX reference paper : [link](https://arxiv.org/abs/1303.5632v2)
- JWST MAST search : [link](https://mast.stsci.edu/search/ui/#/jwst)

NOTE: the ftp links on the PHOENIX website are broken. Converting the format to start with https:// works; see spots_and_faculae_model/README.md for an example

# Setting up the dependencies

To use this repo you'll need to download poetry: https://python-poetry.org/docs/

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

All .fits files are gitignored currently. This means that the wavelength grid is NOT stored in this repo and must be downloaded locally (e.g. into spectral_grids/) before recreating the HDF5 spectra grid.

All other spectra downloads are performed automatically by fits_to_hdf5.py (currently), and the file format conventions for the spectra can just be read directly from within PHOENIX_filename_conventions.py.

working example url format:
https://phoenix.astro.physik.uni-goettingen.de/data/HiResFITS/PHOENIX-ACES-AGSS-COND-2011/Z-0.0/lte06000-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits

# Notes

- the only wavelength grid I can find is WAVE (not the AWAV) which is specified as VACUUM wavelengths. astropy specutils can be used to convert from vacuum to air wavelength if needed (this apparently has a low error)
- there appears to be 1,569,128 wavelength points in total in the PHOENIX fits files I've checked
- JWST Stage 2 products are specified in vacuum wavelength [source](https://jwst-docs.stsci.edu/accessing-jwst-data/jwst-science-data-overview#gsc.tab=0:~:text=Note%20that%20spectroscopic%20data%20products%20have%20wavelengths%20given%20in%20the%20barycentric%20vacuum%20rest%20frame)

# Conventions

I've tried to use .hdf5 and not use .h5, just to keep things simple. (although before commit 6cc4ecd I used .h5, and it doesn't really matter, the files themself are all the same anyway)

## Commit Types

- **`feat`**: A new feature  
- **`fix`**: A bug fix  
- **`docs`**: Documentation-only changes
- **`notes`**: Adding dev-notes / diary entries
- **`style`**: Code style changes (whitespace, formatting, etc. — no code behavior change)  
- **`refactor`**: A code change that neither fixes a bug nor adds a feature  
- **`perf`**: A code change that improves performance  
- **`test`**: Adding or correcting tests  
- **`build`**: Changes to the build system or dependencies (e.g. npm, Makefile)  
- **`ci`**: Changes to CI configuration or scripts (e.g. GitHub Actions, Travis)
- **`chore`**: General maintenance tasks not related to features, fixes, docs, or build

## Branch Names

Branch names can follow the same names, but are formatted like
- **`feat/adding-x-from-y`**
- **`fix-enemymeshes/removing-extraneous-faces`**

(as **`:`**, **`()`** etc would have to be escaped and branch names cannot have whitespace.)

## Example Messages

```bash
feat(bvh): implement initial bounding volume hierarchy generation
fix(plot): correct axis scaling in star visualisation
docs: update README with usage instructions
build: update python_requirements.txt
chore: remove .vscode folder from git tracking (cache)
```

## Sources

This project follows the [Conventional Commits specification v1.0.0](https://www.conventionalcommits.org/en/v1.0.0/#summary).

Some of the commit types listed above are from the [Angular Commit Message Guidelines](https://github.com/angular/angular/blob/22b96b9/CONTRIBUTING.md#-commit-message-guidelines). That link also gives some general guidelines for message formatting that I think are nice.
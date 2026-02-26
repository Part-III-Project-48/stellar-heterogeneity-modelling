<p align="center">
    <a href="/LICENSE" alt="License">
        <img src="https://img.shields.io/badge/License-MIT-blue.svg" /></a>
    <picture><img src="https://img.shields.io/github/commit-activity/m/Part-III-Project-48/stellar-heterogeneity-modelling"/></picture>
    <a href="https://www.python.org/downloads/" alt="Download Python">
        <img src="https://img.shields.io/python/required-version-toml?tomlFilePath=https%3A%2F%2Fraw.githubusercontent.com%2FPart-III-Project-48%2Fstellar-heterogeneity-modelling%2Frefs%2Fheads%2Fmain%2Fpyproject.toml&logo=python&logoColor=white"/></a>
    <a href="/pyproject.toml" alt="Poetry Package Version">
        <img src="https://img.shields.io/badge/dynamic/toml?url=https%3A%2F%2Fraw.githubusercontent.com%2FPart-III-Project-48%2Fstellar-heterogeneity-modelling%2Fmain%2Fpyproject.toml&query=%24.project.version&label=version"/></a>
</p>
<p align="center">
    <picture><img src="https://img.shields.io/github/languages/code-size/Part-III-Project-48/stellar-heterogeneity-modelling"/></picture>
    <picture><img src="https://img.shields.io/github/repo-size/Part-III-Project-48/stellar-heterogeneity-modelling"/></picture>
    <a href="https://www.conventionalcommits.org/en/v1.0.0/" alt="Conventional Commits">
        <img src="https://img.shields.io/badge/Conventional%20Commits-1.0.0-%23FE5196?style=flat-square&logo=conventionalcommits"/></a>
</p>

<!-- See here for prevent images showing as links to opening the image in a new tab on Github : https://stackoverflow.com/questions/40625614/is-it-possible-to-disable-the-automatic-linking-of-images-in-github-markdown-ren -->

# Stellar Spectrum Decomposition using PHOENIX

- PHOENIX homepage : https://phoenix.astro.physik.uni-goettingen.de/
- PHOENIX reference paper : https://arxiv.org/abs/1303.5632v2
- JWST MAST search : https://mast.stsci.edu/search/ui/#/jwst

NOTE: the FTP links on the PHOENIX website are broken. Converting the format to start with https:// works; see below for an example URL.

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

The PHOENIX paper and page talk about ftp (which is now deprecated). The files are actually found (and can be browsed) here: https://phoenix.astro.physik.uni-goettingen.de/data/. (even the http server has a welcome message which calls itself an FTP server!).

I believe (rather confusingly) that the "2011" in many of the file names actually corresponds to the year of the spectral model, rather than the year of upload (which is instead just marked in the server). But the paper doesn't actually state the reasons behind their filenaming so I don't know for sure.

We are using the 'HiResFITS' files currently.

All .fits files are gitignored currently. This means that the wavelength grid is NOT stored in this repo and must be downloaded locally (e.g. into spectral_grids/) before recreating the HDF5 spectra grid.

All other spectra downloads are performed automatically & internally by fits_to_hdf5.py. If you need them, the file format conventions for the spectra can be read directly from within PHOENIX_filename_conventions.py. (although for this code you don't need to understand the PHOENIX filenames to use this code; it deals with it all for you.)

example working url format:
https://phoenix.astro.physik.uni-goettingen.de/data/HiResFITS/PHOENIX-ACES-AGSS-COND-2011/Z-0.0/lte06000-4.50-0.0.PHOENIX-ACES-AGSS-COND-2011-HiRes.fits

# Notes

- the only wavelength grid I can find is WAVE (not the AWAV) which is specified as VACUUM wavelengths. astropy specutils can be used to convert from vacuum to air wavelength if needed (this apparently has a low error)
- there appears to be 1,569,128 wavelength points in total in the PHOENIX fits files I've checked
- JWST Stage 2 products are specified in vacuum wavelength [source](https://jwst-docs.stsci.edu/accessing-jwst-data/jwst-science-data-overview#gsc.tab=0:~:text=Note%20that%20spectroscopic%20data%20products%20have%20wavelengths%20given%20in%20the%20barycentric%20vacuum%20rest%20frame)

<!-- # Conventions

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

Some of the commit types listed above are from the [Angular Commit Message Guidelines](https://github.com/angular/angular/blob/22b96b9/CONTRIBUTING.md#-commit-message-guidelines). That link also gives some general guidelines for message formatting that I think are nice. -->
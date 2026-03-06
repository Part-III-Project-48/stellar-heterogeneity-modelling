# Folder Structure Specification

The structure assumed by `readers.py`:

`neater_observed_spectra_folder/target_name/instrument_type/<any_file_name>.fits`

- `target_name`: name of stellar target
- `instrument_type`: NIRISS (lower wavelengths) or NIRSPEC (higher wavelengths)
- `<any_file_name>` is ignored; just read in all *.fits in the specified folder and stitch them together to form the transmission curve (I guess we assume they are ordered)

from astropy.units import Unit, Quantity
import astropy.units as u

class Instrument():
    def __init__(
            self,
            wavelength_units : Unit,
            flux_units : Unit,
            resolution : Quantity,
            extract_1d_hdu_index : int, # the .fits hdu index which contains the Extract1D spectra
            folder_path : str
    ):
        self.WavelengthUnits = wavelength_units
        self.FluxUnits = flux_units
        self.Resolution = resolution
        self.HDUIndex = extract_1d_hdu_index

        # the folder path convention for this instrument i.e. this instrument's data is assumed to be placed in <target_name>/<Instrument.FolderPath>/
        self.FolderPath = folder_path
        

# treated like static classes

# lower wavelengths
NIRISS : Instrument = Instrument(
	wavelength_units=u.um,
	flux_units=u.MJy,
	resolution=.001 * u.um,
	extract_1d_hdu_index=3,
    folder_path="NIRISS",
)

# higher wavelengths
NIRSPEC : Instrument = Instrument(
	wavelength_units=u.um,
	flux_units=u.MJy,
	resolution=.001 * u.um, # not sure on this
	extract_1d_hdu_index=2,
    folder_path="NIRSPEC",
)
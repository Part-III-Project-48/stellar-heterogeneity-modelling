from astropy import units as u
from astropy.units import Quantity

class spectral_component():
    """
    Represents a tuple of T_eff (in Kelvin), FeH (in dex) and Log_g (in dex)

    Provides an iterable which is unpacked in the order: T_eff, FeH, Log_g
    """
    def __init__(self,
			  t_eff : Quantity[u.K],
              feh : Quantity[u.dimensionless_unscaled],
              log_g : Quantity[u.dimensionless_unscaled]):
        self.T_eff = t_eff
        self.FeH = feh
        self.Log_g = log_g

    def __iter__(self):
        yield self.T_eff
        yield self.FeH
        yield self.Log_g

from .hcts import generate_lightcurve as lc

__all__ = ['generate_lightcurve']


def generate_lightcurve(image, R_planet_pixels=17.25, background=269):
    """
    Generate a light curve
    """
    return lc.generate_lightcurve(image, R_planet_pixels=R_planet_pixels,
                                  background=background)
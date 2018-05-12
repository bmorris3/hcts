import numpy as np
import astropy.units as u
from astropy.constants import R_sun, R_earth
from . import generate_lightcurve as lc

__all__ = ['generate_lightcurve']

from .lightcurve import LightCurve


def generate_lightcurve(image, period, a, R_planet_physical=R_earth,
                        background=269):
    """
    Generate a light curve
    """

    diff = np.diff(np.sum(image, axis=1))
    left, right = np.argmax(diff), np.argmin(diff)

    star_width_pixels = right - left
    star_width = 2 * R_sun
    distance_per_pixel = star_width / star_width_pixels

    orbital_velocity_earth = (2*np.pi*a)/period
    time_per_step = (distance_per_pixel/orbital_velocity_earth).to(u.day).value # days

    radius_planet_pixels = (R_earth / distance_per_pixel).value

    fluxes = lc.generate_lightcurve(image, R_planet_pixels=radius_planet_pixels,
                                  background=background)

    times = np.arange(0, len(fluxes)) * time_per_step

    return LightCurve(times, fluxes)

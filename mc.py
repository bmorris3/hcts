import batman
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.time import Time

from exoplanets import extract_transiting_database, parse_raw_database


exo = extract_transiting_database(parse_raw_database('exoplanets.csv'))
kepler_times = np.load('times_long_cadence.npy')


def draw_prominence_latitude():
    """
    Randomly draw a prominence from the latitude distribution
    """
    max_heliographic_latitude = 55
    latitude = max_heliographic_latitude * (2 * np.random.rand() - 1)
    return latitude


def daily_prom_no_dist(numbers):
    mean_daily_no_of_proms = 14.5
    fwhm_daily_no_of_proms = 10
    sigma = fwhm_daily_no_of_proms/2.35
    return np.exp(-0.5 * (mean_daily_no_of_proms - numbers)**2 / sigma**2)

numbers = np.linspace(0, 32, 200)
prom_no_dist = daily_prom_no_dist(numbers)
cdf = np.cumsum(prom_no_dist)
cdf /= np.max(cdf)
interp_cdf = lambda x: np.interp(x, cdf, numbers)

def draw_daily_prominence_number():
    """
    Randomly draw a number of spots to put on the sun that day
    """
    return int(interp_cdf(np.random.rand()))


def angle_subtended_by_planet(Rp_Rs):
    """
    angle subtended by the planet on the stellar surface in radians
    """
    return np.arctan(2*Rp_Rs / 1.0)


def duration_ecc_term(e, omega):
    return np.sqrt(1 - e**2) / (1 + e * np.sin(omega))


def duration_14(P, Rp_Rs, a_Rs, b, inc, e, omega):
    """
    Winn 2011 Eqn 14
    """
    d =  P / np.pi * np.arcsin(np.sqrt((1 + Rp_Rs**2)**2 - b**2) / a_Rs /
                               np.sin(inc)) * duration_ecc_term(e, omega)
    return d.value

def duration_23(P, Rp_Rs, a_Rs, b, inc, e, omega):
    """
    Winn 2011 Eqn 15
    """
    d =  P / np.pi * np.arcsin(np.sqrt((1 - Rp_Rs**2)**2 - b**2) / a_Rs /
                               np.sin(inc)) * duration_ecc_term(e, omega)
    return d.value


def duration_ingress(P, Rp_Rs, a_Rs, b, inc, e, omega):
    return (duration_14(P, Rp_Rs, a_Rs, b, inc, e, omega) -
            duration_23(P, Rp_Rs, a_Rs, b, inc, e, omega)) / 2.0


class Planet(object):
    def __init__(self, planet_row):

        self.planet_row = planet_row

        # Initialize transit model
        self.params = batman.TransitParams()
        self.params.t0 = planet_row['TT'][0]
        self.params.inc = planet_row['I'][0].value
        self.params.a = planet_row['AR'][0].value
        self.params.rp = planet_row['DEPTH'][0].value**0.5 #planet_row['RR'][0].value
        self.params.ecc = planet_row['ECC'][0].value
        self.params.w = planet_row['OM'][0].value
        self.params.u = [0.0]
        self.params.limb_dark = 'linear'
        self.params.per = planet_row['PER'][0].value

        self.b = planet_row['B'][0]

        self.transit_duration = duration_14(self.params.per, self.params.rp,
                                            self.params.a, self.b,
                                            np.radians(self.params.inc),
                                            self.params.ecc,
                                            np.radians(self.params.w))

        # Initialize some attributes to fill later
        self.transit_model = None
        self.phase = None

    def transit_photosphere(self):

        phase_units_days = (kepler_times - self.params.t0) % self.params.per
        phase_units_days[phase_units_days > 0.5] -= self.params.per
        near_transit = np.abs(phase_units_days) < self.transit_duration

        times = kepler_times[near_transit]
        photosphere_transit_model = batman.TransitModel(self.params, times)
        phased_times = (((times - self.params.t0) % self.params.per) /
                        self.params.per)
        phased_times[phased_times > 0.5] -= 1.0
        self.transit_model = photosphere_transit_model.light_curve(self.params)
        self.phase = phased_times



""" In most cases (particularly for quiet region structures) the initial rise
phase of the filament is quite slow (1–15 km s–1) and can last a few hours,
while this phase can be much shorter (10 min) for active region structures
(Sterling and Moore, 2004; Williams et al., 2005). Parenti 2014"""

filament_duration = 3 * u.hour

"""
If there are ~14 prominences/24 hours, and they last ~3 hours each, then you'd
expect ~2 prom's to be visible at any instant

instantaneous number ~ prom generation rate * prom duration
"""

planet_row = exo[exo['NAME'] == 'HAT-P-11 b']

hat11 = Planet(planet_row)
hat11.transit_photosphere()
p, f = hat11.phase, hat11.transit_model

plt.plot(p, f, '.')

plt.show()
# plt.hist([draw_prominence_latitude() for i in range(1000)])
# #plt.show()
#
# plt.figure()
# plt.hist([draw_daily_prominence_number() for i in range(1000)])
# plt.show()
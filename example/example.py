import numpy as np
from hcts import simulate_lightcurve
import astropy.units as u
from sunpy.net import Fido, attrs as a
from sunpy.map import Map
import matplotlib.pyplot as plt

result = Fido.search(a.Time('2013/5/13 15:55', '2013/5/13 15:56'),
                     a.Instrument('HMI'), a.vso.Physobs('intensity'))

downloaded_files = Fido.fetch(result, path='data/.')

m = Map(downloaded_files[0])

image = m.data

background = np.percentile(image[(image < 500) & (image > 100)], 10)
image[np.isnan(image)] = background

period = 365.25 * u.day
a = 1 * u.AU

control_lc = simulate_lightcurve(image/image.sum(), period, a)

control_lc.plot()
plt.show()
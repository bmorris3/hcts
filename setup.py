from distutils.core import setup, Extension
import numpy as np

module1 = Extension('hcts.generate_lightcurve',
                    sources = ['hcts/generate_lightcurve.c'],
                    include_dirs=[np.get_include()])

setup(name = 'hcts',
      version = '1.0',
      description = 'package',
      packages = ['hcts'],
      include_dirs = [np.get_include()],
      ext_modules = [module1])
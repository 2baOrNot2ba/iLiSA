#!/usr/bin/env python

from setuptools import setup, find_packages

setup(name = 'LiSA',
      version = '2.0',
      description = 'Package for international LOFAR station in stand-alone mode',
      author = 'Tobia D. Carozzi',
      author_email = 'tobia.carozzi@chalmers.se',
      packages = find_packages(),
      package_data = {
            'ilisa.observations': ['access_config.yml'],
            'ilisa.antennameta': ['share/CalTables/*.dat',
                                  'share/StaticMetaData/*.conf']},
      license = 'ISC',
      classifiers = [
          'Development Status :: 1 - Planning',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: ISC License',
          'Programming Language :: Python :: 2.7',
          'Topic :: Scientific/Engineering :: Astronomy',
          'Topic :: Scientific/Engineering :: Mathematics',
          'Topic :: Scientific/Engineering :: Visualization'
      ],
      install_requires=[
          'numpy>=1.10',
          'python-casacore',
          'matplotlib>=1.5'
      ],
      scripts = ['scripts/']
     )

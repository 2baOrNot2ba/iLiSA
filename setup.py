#!/usr/bin/env python

from setuptools import setup, find_packages
import ilisa

setup(name='iLiSA',
      version=ilisa.__version__,
      description="""Package for handling international LOFAR station in
                  stand-alone mode""",
      author='Tobia D. Carozzi',
      author_email='tobia.carozzi@chalmers.se',
      packages=find_packages(),
      package_data={
            'ilisa.observations': ['config/access_lclstn_DIST.conf',
                                   'config/project_0_DIST.conf',
                                   'scansess_examples/*.yaml'],
            'ilisa.antennameta': ['share/CalTables/*/data/*.dat',
                                  'share/StaticMetaData/*.conf']},
      license='ISC',
      classifiers=[
          'Development Status :: 3 - Alpha',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: ISC License',
          'Programming Language :: Python :: 2.7',
          'Topic :: Scientific/Engineering :: Astronomy',
          'Topic :: Scientific/Engineering :: Mathematics',
          'Topic :: Scientific/Engineering :: Visualization'
      ],
      install_requires=[
          'numpy',
          'python-casacore',
          'matplotlib',
          'PyYAML',
          'h5py'
      ],
      entry_points={
          'console_scripts': [
              'ilisa_rec = scripts.ilisa_rec:main',
              'ilisa_cmd = scripts.ilisa_cmd:main',
              'ilisa_view = scripts.ilisa_view:main',
          ]
      }
      #,scripts=['scripts/ilisa_cmd.py', 'scripts/ilisa_rec.py', 'scripts/ilisa_view.py']
      )

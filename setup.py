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
            'ilisa.operations': ['config/access_lclstn_DIST.conf',
                                      'config/project_0_DIST.conf',
                                      'scansess_examples/*.yaml'],
            'ilisa.antennameta': ['share/CalTables/*/data/*.dat',
                                  'share/StaticMetaData/*.conf']},
      license='ISC',
      classifiers=[
          'Development Status :: 3 - Alpha',
          'Intended Audience :: Science/Research',
          'License :: OSI Approved :: ISC License',
          'Programming Language :: Python :: 3.6',
          'Topic :: Scientific/Engineering :: Astronomy',
          'Topic :: Scientific/Engineering :: Mathematics',
          'Topic :: Scientific/Engineering :: Visualization'
      ],
      install_requires=[
          'numpy',
          'scipy',
          'python-casacore',
          'matplotlib',
          'PyYAML',
          'h5py'
      ],
      entry_points={
          'console_scripts': [
              'ilisa_adm = ilisa.operations.stationdriver:main_cli',
              'ilisa_rec = ilisa.operations.scan:main_cli',
              'ilisa_obs = ilisa.operations.scansession:main_cli',
              'ilisa_view = ilisa.operations.data_io:main',
              'ilisa_applycal = ilisa.calim.calibration:main_cli',
              'ilisa_image = ilisa.calim.imaging:main_cli',
              'ilisa_sched = ilisa.scripts.ilisa_sched:main',
              'pl_bfs = ilisa.pipelines.bfs_data:main_cli',
              'pl_rec = ilisa.pipelines.bfbackend:bfsrec_main_cli'
          ]
      }
      )

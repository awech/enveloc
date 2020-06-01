from setuptools import setup, find_packages

long_description = '''
enveloc: open-source python package built on obspy to locate seismic signals
using envelope cross correlation and a grid search
'''

setup(
  name = 'enveloc',
  packages = find_packages(),
  include_package_data=True,
  version = '1.1.2',
  license='LGPL',
  description = 'enveloc - envelope cross-correlation location method',
  long_description = long_description,
  author = 'Aaron Wech',
  author_email = 'awech@usgs.gov',
  url = 'https://github.com/awech/enveloc',
  download_url = 'https://github.com/awech/enveloc/archive/v1.1.2.tar.gz',
  keywords = ['envelope cross correlation', 'seimsology', 'tremor'],
  install_requires=[
          'numpy',
          'scipy',
          'matplotlib',
          'proj',
          'obspy',
          'cartopy',
          'scikit-learn',
          'utm',
      ],
  setup_requires=[
          'numpy',
          'proj',
      ],
  classifiers=[
    'Development Status :: 4 - Beta',
    'Intended Audience :: Science/Research',
    'Topic :: Scientific/Engineering',
    'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
    'Programming Language :: Python :: 3.7',
    'Programming Language :: Python :: 3.8',
  ],
)
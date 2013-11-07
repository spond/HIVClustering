#!/usr/bin/env python3.2

import sys

from os.path import abspath, join, split
from setuptools import setup

sys.path.insert(0, join(split(abspath(__file__))[0], 'lib'))
from hivclustering import __version__ as _hivclustering_version

setup(name='hivclustering',
      version=_hivclustering_version,
      description='HIV molecular clustering tools',
      author='Sergei Kosakovsky Pond',
      author_email='spond@ucsd.edu',
      url='http://github.com/spond/hivclustering',
      license='MIT License',
      packages=['hivclustering'],
      package_dir={'': 'lib'},
      package_data={'hivclustering': [
            'data/HBL/*.bf',
    ]},
      data_files=[('/usr/local/bin', [
            'bin/hivnetworkcsv', 'bin/networkbuild.py', 'bin/TNS'
      ])],
      requires=['hppy']
     )

#!/usr/bin/env python

from distutils.core import setup

setup(name='utils',
      version='1.0',
      description='utilities for Offshore Geotechnics',
      author='Kevin Duffy',
      author_email='k.duffy@tudelft.nl',
      url='https://github.com/triaduct/offshoregeotechnics',
      packages=['utils'],
      install_requires=['matplotlib', 'numpy']
      )
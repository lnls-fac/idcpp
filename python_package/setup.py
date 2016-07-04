#!/usr/bin/env python3

from setuptools import setup

with open('VERSION','r') as _f:
    __version__ = _f.read().strip()

setup(
    name='insertion_devices',
    version=__version__,
    author='lnls-fac',
    description='insertion_devices python package',
    url='https://github.com/lnls-fac/insertion_devices',
    download_url='https://github.com/lnls-fac/insertion_devices',
    classifiers=[
        'Intended Audience :: Science/Research',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering'
    ],
    packages=['insertion_devices'],
    package_data={'insertion_devices': ['_insertion_devices.so', 'VERSION']},
    zip_safe=False
)

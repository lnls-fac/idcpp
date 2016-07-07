#!/usr/bin/env python3

from setuptools import setup

with open('VERSION','r') as _f:
    __version__ = _f.read().strip()

setup(
    name='idcpp',
    version=__version__,
    author='lnls-fac',
    description='idcpp python package',
    url='https://github.com/lnls-fac/idcpp',
    download_url='https://github.com/lnls-fac/idcpp',
    classifiers=[
        'Intended Audience :: Science/Research',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering'
    ],
    packages=['idcpp'],
    package_data={'idcpp': ['_idcpp.so', 'VERSION']},
    zip_safe=False
)

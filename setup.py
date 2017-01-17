#!/usr/bin/env python

try:
    from setuptools import setup
    setup
except ImportError:
    from distutils.core import setup
    setup

setup(
    name="dynamics",
    version='0.0.1',
    author="Asher Wasserman",
    author_email="adwasser@ucsc.edu",
    description="Astrophysical Dynamics",
    packages=["dynamics"],
    package_data={
        "": ["LICENSE"],
    },
    include_package_data=True,
)

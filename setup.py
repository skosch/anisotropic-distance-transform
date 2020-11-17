#!/usr/bin/env python

from distutils.core import setup, Extension
import numpy as np
import pathlib

HERE = pathlib.Path(__file__).parent
README = (HERE / "README.md").read_text()

c_adt = Extension('adt', sources = ['adt.c'], include_dirs = [np.get_include()])

setup(name = 'anisotropic-distance-transform',
      version = '0.1',
      description = 'Anisotropic Euclidean Distance Transform (2D)',
      long_description=README,
      ext_modules = [c_adt],
      long_description_content_type="text/markdown",
      url="https://github.com/skosch/anisotropic-distance-transform",
      author="Jan Hosang, Sebastian Kosch",
      license="MIT",
      classifiers=[
          "License :: OSI Approved :: MIT License",
          "Programming Language :: Python :: 3",
          "Programming Language :: Python :: 3.7",
      ],
      packages=["anisotropic-distance-transform"],
      include_package_data=True,
      install_requires=["numpy"],
      extras_require={"dev": ["matplotlib"]})

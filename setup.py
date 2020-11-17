#!/usr/bin/env python

from distutils.core import setup, Extension
import numpy as np
from os import path

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()

def parse_requirements(filename):
        """Load requirements from a pip requirements file."""
        lineiter = (line.strip() for line in open(filename))
        return [line for line in lineiter if line and not line.startswith("#")]

c_adt = Extension('adt', sources=['anisotropic-distance-transform/ext/adt.c'], include_dirs=[np.get_include()])

setup(name='anisotropic-distance-transform',
      version='0.1.7',
      description='Anisotropic Euclidean Distance Transform (2D)',
      long_description=long_description,
      long_description_content_type="text/markdown",
      ext_modules=[c_adt],
      author="Sebastian Kosch",
      author_email="skosch@users.noreply.github.com",
      url="https://github.com/skosch/anisotropic-distance-transform",
      license="MIT",
      classifiers=[
          "License :: OSI Approved :: MIT License",
          "Programming Language :: Python :: 3",
          "Programming Language :: Python :: 3.7",
      ],
      packages=["anisotropic-distance-transform"],
)

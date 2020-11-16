#!/usr/bin/env python

from distutils.core import setup, Extension
import numpy as np

c_gdt = Extension('C_gdt', sources = ['C_gdt.c'],
       include_dirs = [np.get_include()],
)

setup (name = 'GDT',
       version = '1.0',
       description = 'Generalized Distance Transform',
       ext_modules = [c_gdt])

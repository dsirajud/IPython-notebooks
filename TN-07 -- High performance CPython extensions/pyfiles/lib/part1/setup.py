#!/usr/bin/env python

from distutils.core import setup, Extension

setup(name             = "numpy_c_ext_example",
      version          = "1.0",
      description      = "Example code for blog post.",
      author           = "J. David Lee",
      author_email     = "contact@crumpington.com",
      maintainer       = "contact@crumpington.com",
      url              = "https://www.crumpington.com",
      ext_modules      = [
          Extension(
              'evolve1', ['evolve1.c'],
              extra_compile_args=["-Ofast", "-march=native"]),
      ], 
      
)

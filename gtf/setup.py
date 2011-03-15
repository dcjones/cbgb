#!/usr/bin/env python

from Cython.Distutils import build_ext
from setuptools import setup, Extension

setup( name         = 'gtf',
       version      = '0.1',
       author       = 'Daniel Jones',
       author_email = 'dcjones@cs.washington.edu',
       description  = 'A library for parsing the Gene Transfer Format (GTF)',
       url          = '',
       requires     = [ 'python (>=2.6, <3.0)' ],
       ext_modules  = [ Extension( 'gtf', ['gtf.pyx'] ) ],
       cmdclass     = { 'build_ext' : build_ext }
       )


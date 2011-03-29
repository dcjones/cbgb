#!/usr/bin/env python


from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext


ext_modules = [Extension(name         = 'ushuffle',
                         sources      = ['ushuffle.pyx', 'ushuffle/ushuffle.c'],
                         include_dirs = ['ushuffle'])]


setup(
        name = 'ushuffle',
        cmdclass = {'build_ext': build_ext},
        ext_modules = ext_modules
)



#!/usr/bin/env python

from distutils.core      import setup
from distutils.extension import Extension
from Cython.Distutils    import build_ext


ext_models = [Extension(name    = 'gtf',
                        sources = ['gtf.pyx', 'common.c',
                                   'gtf_parse.c', 'str_map.c'])]


setup( name         = 'gtf',
       version      = '0.2',
       author       = 'Daniel Jones',
       author_email = 'dcjones@cs.washington.edu',
       description  = 'A library for parsing the Gene Transfer Format (GTF)',
       url          = '',
       requires     = [ 'python (>=2.6, <3.0)' ],
       ext_modules  = ext_models,
       cmdclass     = {'build_ext' : build_ext}
       )


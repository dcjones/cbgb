#!/usr/bin/env python

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

samtools_src = \
    '''
    samtools/bam_aux.c
    samtools/bam.c
    samtools/bam_import.c
    samtools/bam_index.c
    samtools/bam_pileup.c
    samtools/bgzf.c
    samtools/faidx.c
    samtools/knetfile.c
    samtools/kstring.c
    samtools/razf.c
    samtools/sam.c
    samtools/sam_header.c
    '''.split()

ext_modules = [Extension(name          = 'bam',
                         sources       = ['bam.pyx'] + samtools_src,
                         include_dirs  = ['samtoools'],
                         libraries     = ['z'],
                         define_macros = [('_USE_KNETFILE', None),
                                          ('_LARGEFILE64_SOURCE', None),
                                          ('_FILE_OFFSET_BITS', '64')])]

setup(
  name = 'bam',
  cmdclass = {'build_ext': build_ext},
  ext_modules = ext_modules
)


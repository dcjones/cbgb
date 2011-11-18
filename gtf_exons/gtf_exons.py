#!/usr/bin/env python

from gtf import gtf_file
from sys import argv, stdout, stderr, stdin

if len(argv) < 2:
    stderr.write('useage: gtf_exons.py genes.gtf')

for row in gtf_file(argv[1]):
    if row.feature != 'exon': continue

    stdout.write('{seqname}\t{start}\t{end}\t{name}\t{score}\t{strand}\n'.format(
        seqname = row.seqname,
        start   = int(row.start) - 1,
        end     = row.end,
        name    = row.attributes['gene_id'],
        score   = 0,
        strand  = row.strand))



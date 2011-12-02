#!/usr/bin/env python


'''
Given a GTF file, make a new one with the unspliced form of each gene included
as well.
'''

from gtf         import gtf_file
from sys         import argv, stdout, stderr, stdin
from collections import defaultdict

if len(argv) < 2:
    stderr.write('usage: gtf_pre_mrna.py genes.gtf\n')
    exit(1)

fn = argv[1]


class transcript:
    def __init__(self):
        self.strand  = None
        self.seqname = None
        self.start   = None
        self.end     = None
        self.gene_id = None
        self.n       = 0

    def add_row(self, row):
        if self.strand  is None: self.strand  = row.strand
        if self.seqname is None: self.seqname = row.seqname
        if self.gene_id is None: self.gene_id = row.attributes['gene_id']

        if self.start is None or row.start < self.start: self.start = row.start
        if self.end   is None or row.end   > self.end:   self.end = row.end

        self.n += 1


T = defaultdict(transcript)

for row in gtf_file(fn):
    if row.feature != 'exon': continue
    tid = row.attributes['transcript_id']
    T[tid].add_row(row)



# use only unique pre-mrnas (that are spliced)
S = defaultdict(set)
for (tid, t) in T.iteritems():
    if t.n <= 1: continue
    key = (t.gene_id, t.seqname, t.strand, t.start, t.end)
    S[key].add(tid)



gtf = '{seqname}\t{source}\t{feature}\t{start}\t{end}\t{score}\t' \
      '{strand}\t{frame}\tgene_id "{gene_id}"; transcript_id "{transcript_id}";' \
      ' spliced_transcript_ids "{spliced_transcript_ids}";\n'
gene_pre_count = defaultdict(lambda: 1)

for ((gid, seqname, strand, start, end), tids) in S.iteritems():

    transcript_id = '{gene_id}.pre-mrna.{k:03d}'.format(
                        gene_id = gid,
                        k       = gene_pre_count[gid])
    gene_pre_count[gid] += 1

    stdout.write(gtf.format(
        seqname = seqname,
        source  = 'gtf_pre_mrna',
        feature = 'exon',
        start   = start,
        end     = end,
        score   = '.',
        strand  = strand,
        frame   = '.',
        gene_id = gid,
        transcript_id = transcript_id,
        spliced_transcript_ids = ','.join(tids)))




#!/usr/bin/env python

'''
bamcov
------
Measure the aggregate read coverage across transcripts given a BAM file and a GTF file.

'''


import bam
import gtf
import argparse
import numpy as np
from sys         import stdout, stderr, stdin
from collections import defaultdict
from math        import ceil


class gene:
    '''
    An exteremely simple gene representation.
    '''

    def __init__(self):
        self.gene_id = None
        self.seqname = None
        self.strand  = None
        self.exons   = []


    def add_row(self, row):
        '''
        Add an exon to the gene.
        '''

        if self.gene_id is None: self.gene_id = row.attributes['gene_id']
        if self.seqname is None: self.seqname = row.seqname
        if self.strand  is None: self.strand  = row.strand
        self.exons.append([row.start - 1, row.end - 1])

    def flatten(self):
        '''
        Take the union of all exons, resulting in a non-overlapping set.
        '''

        xs = []
        self.exons.sort()
        i = 0
        j = 1
        while j < len(self.exons):
            if self.exons[i][1] >= self.exons[j][0]:
                self.exons[i][1] = max(self.exons[i][1], self.exons[j][1])
                self.exons[j] = None
            else:
                xs.append(self.exons[i])
                i = j
            j += 1

        xs.append(self.exons[i])

        self.exons = xs



def read_genes(genes_fn):
    '''
    For each gene in the GTF file, compute the union of all exons, in sorted
    order.
    '''

    stderr.write('parsing GTF file ... ');

    genes = defaultdict(gene)

    for row in gtf.gtf_file(genes_fn):
        if row.feature != 'exon': continue
        if 'gene_id' not in row.attributes: continue

        genes[row.attributes['gene_id']].add_row(row)

    stderr.write('done. ({0} genes)\n'.format(len(genes)))

    for g in genes.itervalues():
        g.flatten()

    return genes



if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('-n', metavar = 'n', type = int,
                    help = 'number of bins', default = 100)
    ap.add_argument('-s', '--stranded', action = 'store_true')
    ap.add_argument('--min-length', default = 1000, type = int)
    ap.add_argument('--min-count', default = 25, type = int)
    ap.add_argument('--gene-list', default = None, type = str, dest = 'gene_list_fn')
    ap.add_argument('genes_fn', metavar = 'genes.gtf')
    ap.add_argument('reads_fn', metavar = 'reads.bam')
    args = ap.parse_args()

    reads = bam.Bam(args.reads_fn)
    genes = read_genes(args.genes_fn)


    # averaged histogram
    xs = np.zeros(args.n)

    # histogram for the given gene
    ys = np.zeros(args.n)

    # record used genes
    if args.gene_list_fn is not None:
        genes_f = open(args.gene_list_fn, 'w')
    else:
        genes_f = None


    for (i, g) in enumerate(genes.itervalues()):
        stderr.write('({i}/{n}) {gene_id}\n'.format(
            i = i, n = len(genes), gene_id = g.gene_id))

        if args.stranded:
            strand = 0 if g.strand == '+' else 1
        else:
            strand = None

        # position along the transcript
        u = 0

        # gene length
        glen = float(sum(end - start + 1 for (start, end) in g.exons) + 1)
        if glen < args.min_length: continue

        ys.fill(0)

        for (start, end) in g.exons:
            cs = reads.counts(g.seqname, start, end, strand)
            if cs is None: continue

            if g.strand == '+':
                for (v, c) in enumerate(cs):
                    idx = np.floor(args.n * ((u + v) / glen))
                    ys[idx] += c
            else:
                for (v, c) in enumerate(cs):
                    idx = np.floor(args.n * (((glen - 1) - (u + v)) / glen))
                    ys[idx] += c

            u += end - start + 1

        ys_sum = sum(ys)
        if ys_sum < args.min_count: continue
        if genes_f is not None:
            genes_f.write('{0}\n'.format(g.gene_id))
        xs += ys / ys_sum


    xs = xs / sum(xs)

    for (i, x) in enumerate(xs):
        stdout.write('{0:0.2f}\t{1:e}\n'.format(
            100.0 * float(i) / args.n, x))



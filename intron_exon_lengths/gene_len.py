#!/usr/bin/python

#!/usr/bin/python

import re
from collections import namedtuple, defaultdict
from sys         import argv, stdout, stdin, stderr
from itertools   import izip


gtf_row = namedtuple( 'gtf_row', 'seqname source feature start end ' \
                                 'score strand frame attributes' )

bed_row = namedtuple( 'bed_row', 'seqname start end name score strand' )

def gtf_row_from_gtf_line( line ):
    line = line.split('\t')
    if len(line) < 9:
        raise TypeError

    seqname = line[0]
    source  = line[1]
    feature = line[2]

    # gff/gtf is 1-based, end-inclusive
    start   = int(line[3])
    end     = int(line[4])

    score   = line[5]
    strand  = line[6]
    frame      = line[7]
    attributes = {}

    for mat in re.finditer( r'\s*(\w+)\s+(([\.\w]+)|"([\.\w]+)")\s*(;|\s*$)', line[8] ):
        attributes[mat.group(1)] = mat.group(2) # .strip('"')

    return gtf_row( seqname, source, feature, \
                    start, end, score, strand, frame, \
                    attributes )




def parse_gtf( gtf_f, feature_filter = 'exon' ):
    ''' extract annotations from a gtf file '''

    stderr.write( 'parsing ... ' )

    rows = []
    i = 0
    for line in gtf_f:
        i+=1

        try:
            row = gtf_row_from_gtf_line( line )
        except TypeError:
            stderr.write( 'Only %d fields found on line %d of %s. Skipping.\n' % \
                          (len(line.split('\t')), i, gtf_fn) )
            continue

        if row.feature == feature_filter:
            rows.append( row )

    stderr.write( 'done (%d features).\n' % len(rows) )

    return rows



def get_gene_lengths( rows ):

    # organize by gene

    stderr.write( 'organizing genes ... ' )
    genes = defaultdict(list)
    for row in rows:
        if 'gene_id' in row.attributes:
            genes[row.attributes['gene_id']].append(row)

    stderr.write( 'done (%d genes)\n' % len(genes) )


    stderr.write( 'counting lengths ... ' )

    gene_lens = dict()

    for (gene_id,rows) in genes.iteritems():
        rows.sort( key = lambda row: row.start )
        N = 0
        prev_end = -1
        for row in rows:
            if prev_end == -1 or prev_end < row.start:
                N += row.end - row.start + 1
                prev_end = row.end
            elif prev_end < row.end:
                N += row.end - prev_end
                prev_end = row.end

        gene_lens[gene_id] = N


    stderr.write( 'done (%d genes)\n' % len(gene_lens) )

    return gene_lens



def print_gene_lens( gene_lens ):
    for (gene_id,N) in gene_lens.iteritems():
        stdout.write( "%s\t%d\n" % (gene_id, N) )


def main():
    if len(argv) > 1:
        rows = parse_gtf(stdin,argv[1])
    else:
        rows = parse_gtf(stdin)

    gene_lens = get_gene_lengths(rows)
    print_gene_lens(gene_lens)


if __name__ == '__main__':
    main()

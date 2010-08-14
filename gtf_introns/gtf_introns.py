#!/usr/bin/python

import re
from collections import namedtuple, defaultdict
from copy        import copy
from sys         import stdout, stdin, stderr
from itertools   import izip


gtf_row = namedtuple( 'gtf_row', 'seqname source feature start end ' \
                                 'score strand frame attributes' )

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

    for mat in re.finditer( r'\s*(\w+)\s+(([\.\w]+)|"([\.\w]+)")\s*;', line[8] ):
        attributes[mat.group(1)] = mat.group(2) # .strip('"')

    return gtf_row( seqname, source, feature, \
                    start, end, score, strand, frame, \
                    attributes )




def parse_gtf( gtf_f ):
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

        if row.feature == 'exon':
            rows.append( row )

    stderr.write( 'done (%d exons).\n' % len(rows) )

    return rows


def get_introns( rows ):

    # organize by transcript

    stderr.write( 'getting transcripts ... ' )

    transcripts = defaultdict(list)
    for row in rows:
        if 'transcript_id' in row.attributes:
            transcripts[row.attributes['transcript_id']].append( row )

    stderr.write( 'done (%d transcripts).\n' % len(transcripts) )

    # get introns

    stderr.write( 'getting introns ... ' )

    introns = []
    for (transcript_id,rows) in transcripts.iteritems():
        if len(set( (row.seqname,row.strand) for row in rows )) != 1:
            stderr.write( 'Malformed transcript "%s". Skipping.' % transcript_id)

        rows.sort( key = lambda row: row.start )
        starts = [ row.start for row in rows ]
        ends   = [ row.end   for row in rows ]

        for (i,j) in izip(ends[:-1],starts[1:]):
            assert i < j
            if i+1 <= j-1:
                intron_row = gtf_row(
                        seqname = rows[0].seqname,
                        source  = rows[0].source,
                        feature = 'intron',
                        start   = i+1,
                        end     = j-1,
                        score   = '.',
                        strand  = rows[0].strand,
                        frame   = '.',
                        attributes = rows[0].attributes )

                introns.append( intron_row )


    stderr.write( 'done (%d introns).\n' % len(introns) )

    return introns


def print_gtf( fout, rows ):

    gtf_str = '{seqname}\t{source}\t{feature}\t{start}\t{end}\t{score}\t{strand}\t{frame}\t{attr_str}\n'

    for row in rows:
        attr_str = '; '.join( [ '%s %s' % attr for attr in row.attributes.iteritems() ] )

        fout.write( gtf_str.format(
            seqname = row.seqname,
            source  = row.source,
            feature = row.feature,
            start   = row.start,
            end     = row.end,
            score   = row.score,
            strand  = row.strand,
            frame   = row.frame,
            attr_str = attr_str ) )





def main():
    rows = parse_gtf(stdin)
    rows = rows + get_introns(rows)

    stderr.write( 'sorting ... ' )

    rows.sort( key = lambda row: row.start )
    rows.sort( key = lambda row: row.seqname )

    stderr.write( 'done.\n' )

    stderr.write( 'writing ... ' )

    print_gtf( stdout, rows )

    stderr.write( 'done.\n' )

if __name__ == '__main__':
    main()



#!/usr/bin/python


#
#   gtf_transcripts
#   ---------------
#   Make a bed file giving the extents of all the transcripts
#   given in a gtf file.
#
#



import re
from collections import namedtuple, defaultdict
from copy        import copy
from sys         import stdout, stdin, stderr
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



def get_transcripts( rows ):

    # organize by transcript

    stderr.write( 'getting transcripts ... ' )

    transcripts = defaultdict(list)
    for row in rows:
        if 'transcript_id' in row.attributes:
            transcripts[row.attributes['transcript_id']].append( row )

    stderr.write( 'done (%d transcripts).\n' % len(transcripts) )

    # get introns

    stderr.write( 'getting introns ... ' )

    bed_rows = []
    for (transcript_id,rows) in transcripts.iteritems():
        if len(set( (row.seqname,row.strand) for row in rows )) != 1:
            stderr.write( 'Malformed transcript "%s". Skipping.' % transcript_id)

        starts = [ row.start for row in rows ]
        ends   = [ row.end   for row in rows ]

        starts.sort()
        ends.sort()

        bed_rows.append( bed_row( rows[0].seqname, starts[0]-1, ends[-1],
                                  transcript_id, 0, rows[0].strand ) )


    return bed_rows


def print_bed( fout, rows ):

    bed_str = '{seqname}\t{start}\t{end}\t{name}\t{score}\t{strand}\n'

    for row in rows:
        fout.write( bed_str.format(
            seqname = row.seqname,
            start   = row.start,
            end     = row.end,
            name    = row.name,
            score   = row.score,
            strand  = row.strand ) )


def main():

    rows = parse_gtf(stdin)
    transcripts = get_transcripts(rows)
    print_bed(stdout,transcripts)



if __name__ == '__main__':
    main()

#!/usr/bin/env python


#
#   gtf_annotated
#   ---------------
#   Make a bed file giving disjoint intervals of all the regions covered by
#   annotations.
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

    stderr.write( 'done (%d regions).\n' % len(rows) )

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



def get_annotated( rows, stranded=True ):

    # organize by transcript

    stderr.write( 'sorting by position ... ' )

    if stranded:
        rows.sort( key = lambda row: (row.seqname,row.strand,row.start) )
    else:
        rows.sort( key = lambda row: (row.seqname,row.start) )

    stderr.write( 'done.\n' )



    stderr.write( 'getting annotated regions ... ' )

    bed_rows = []

    last_seq = None
    i = None
    j = None

    for row in rows:
        if stranded:
            seq = (row.seqname,row.strand)
        else:
            seq = row.seqname

        if last_seq != seq:
            last_seq = seq
            i = row.start
            j = row.end
            continue

        if j >= row.start:
            j = row.end
        else:
            bed_rows.append( bed_row( row.seqname, i-1, j, '.', 0, row.strand ))
            i = row.start
            j = row.end



    if i and j and rows:
        row = rows[-1]
        bed_rows.append( bed_row( row.seqname, i-1, j, '.', 0, row.strand ))

    stderr.write( 'done (%d regions).' % len(bed_rows) )

    return bed_rows



def print_gtf( fout, rows ):
    gtf_str = '{seqname}\t.\ttranscript\t{start}\t{end}\t{score}\t{strand}\t.\tgene_id "{name}"; transcript_id "{name}"\n'

    for row in rows:
        fout.write( gtf_str.format(
            seqname = row.seqname,
            start   = row.start + 1,
            end     = row.end,
            name    = row.name,
            score   = row.score,
            strand  = row.strand ) )




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
    annotated   = get_annotated(transcripts)
    print_bed(stdout,annotated)
    #print_gtf( stdout, transcripts )



if __name__ == '__main__':
    main()

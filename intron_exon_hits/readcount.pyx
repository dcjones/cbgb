

from sys         import (argv,stdout,stderr)
from collections import defaultdict
from warnings    import warn
from functools   import partial
import re
import csv

cdef extern from "stdio.h":
    int snprintf(char *str, size_t size, char *format, ...)


cdef extern from "samtools/bam.h":
    ctypedef struct bam1_t:
        pass

    ctypedef struct bam_header_t:
        pass

    ctypedef struct bam_index_t:
        pass

    ctypedef void*    bamFile

    ctypedef bam1_t* const_bam1_ptr "const bam1_t*"
    ctypedef char* const_char_ptr "const char*"
    ctypedef int (*bam_fetch_f)(const_bam1_ptr b, void *data)

    int bam1_strand( bam1_t* )
    void bam_close( bamFile )
    bamFile bam_open( char* fn, char* mode )
    bam_index_t *bam_index_load(char *fn)
    bam_header_t *bam_header_read(bamFile fp)


    void bam_header_destroy( bam_header_t* )
    int bam_parse_region(bam_header_t *header, const_char_ptr str, int *ref_id, int *begin, int *end)
    int bam_fetch(bamFile fp, bam_index_t *idx, int tid, int beg, int end, void *data, bam_fetch_f func)

cdef struct count_state:
    unsigned char strand
    size_t count




cdef int tally_read( const_bam1_ptr b, void* a ):
    if bam1_strand(b) == (<count_state*>a).strand:
        (<count_state*>a).count += 1
    return 0



cdef size_t _readcount(
        bamFile read_f, bam_header_t* read_header,
        bam_index_t* read_index,
        char* seq, size_t start, size_t end, int strand ):

    cdef count_state s
    s.strand = strand
    s.count = 0

    cdef int ref_id,start_,end_
    cdef char region[256]
    snprintf( region, 256, "%s:%d-%d", seq, start, end )

    bam_parse_region( read_header, <const_char_ptr>region, &ref_id, &start_, &end_ )
    bam_fetch( read_f, read_index, ref_id, start_, end_, <void*>&s, &tally_read )

    return s.count



cdef class Bam:
    cdef bamFile       read_file
    cdef bam_index_t*  read_index
    cdef bam_header_t* read_header

    def __cinit__( self ):
        self.read_file   = NULL
        self.read_index  = NULL
        self.read_header = NULL

    def __init__( self, filename, index_filename ):
        self.read_file   = bam_open( filename, "rb" )
        if self.read_file == NULL:
            raise Exception('Can\'t open BAM file "%s"' % filename )

        self.read_index  = bam_index_load( index_filename )
        if self.read_index == NULL:
            raise Exception('Can\'t open BAM index "%s"' % index_filename )

        self.read_header = bam_header_read(self.read_file)
        if self.read_header == NULL:
            raise Exception('Can\'t parse BAM header.')


    def __dealloc__( self ):
        bam_header_destroy(self.read_header)
        # bam_index_destroy(idx);*/ /* XXX: why does this cause a crash? */
        bam_close(self.read_file)

    def readcount( self, seq, start, end, strand ):
        return _readcount( self.read_file, self.read_header, self.read_index,
                           seq, start, end, 0 if strand == '+' else 1 )



def gtf_exons_introns_cds( f ):
    '''
    Given a gtf file, create a dictionary of each exon associated with a
    transcript.
    '''

    exons = defaultdict(list)
    cds   = defaultdict(list)

    attribute_pat = re.compile( r'\s*(\w+)\s+(([\.\w]+)|"([\.\w]+)")\s*;' )

    for row in csv.reader(f,delimiter='\t'):
        if len(row) < 8:
            warn( 'Only {n} fields found.'.format( n=len(row) ) )
            continue

        (seqname, source, feature, start, end, score, strand, fname) = row[:8]

        (attributes_str, comments) = ('','')
        if len(row) > 8:
            attributes_str = row[8]
        if len(row) > 9:
            comments = row[9]

        attributes = {}
        for mat in attribute_pat.finditer( attributes_str ):
            attributes[mat.group(1)] = mat.group(2)
        if not 'transcript_id' in attributes:
            warn( 'Missing transcript_id from {seqname}.'.format( seqname=seqname ) )
            continue


        if feature == 'exon':
            exons[(seqname,strand,attributes['transcript_id'])].append( \
                          (int(start),int(end)) )
        elif feature == 'CDS':
            cds[(seqname,strand,attributes['transcript_id'])].append( \
                          (int(start),int(end)) )


    for i in exons.itervalues():
        i.sort()

    for i in cds.itervalues():
        i.sort()


    introns = introns_from_exons(exons)

    return (exons,introns,cds)


def introns_from_exons(exons):
    introns = {}
    for k in exons.iterkeys():
        (seqname,strand,transcript_id) = k
        (exon_starts,exon_ends) = zip(*exons[k])
        introns[k] = zip(exon_ends[:-1],exon_starts[1:])

    return introns


def chrom_num( seqname ):
    ''' Return the sequences chrom number, if it has one. '''
    mat =  re.search(r'(\d+)',seqname)
    if mat:
        return int(mat.group(1))
    else:
        return 99999 # arbitrary large number so unnumbered chromosomes are ordered at the end


def fst(x):
    return x[0]

def chrom_num_fst(x):
    return chrom_num(x[0])


def start_pos( dct, key ):
    if not dct[key]:
        return 0
    return dct[key][0]


def reads_per_nt( Bam B, seq, strand, intervals ):


    length = 0
    reads  = 0

    for (start,end) in intervals:
        reads  += B.readcount( seq, start, end, strand )
        length += end-start # XXX: are these inclusive ??

    if length == 0:
        return None
    else:
        return (reads,length)


if __name__ == '__main__':
    if len(argv) < 4:
        stderr.write( 'Usage: readcount annotations.gtf reads.bam reads_index.bai\n' )
        exit(1)

    B = Bam( argv[2], argv[3] )

    (exons,introns,cds) = gtf_exons_introns_cds( open(argv[1]) )

    print len(exons)

    keys = [ key for key in exons.keys() if len(exons[key]) > 0 ]
    keys.sort( key=partial( start_pos, cds ) )
    keys.sort( key=fst )           # seqname + strand

    for (i,key) in enumerate(keys):

        (seq,strand,transcript) = key

        pred = exons[keys[i-1]] if i   > 0 and keys[i-1][:2] == key[:2] else None
        succ = exons[keys[i+1]] if i+1 < len(keys) and keys[i+1][:2] == key[:2] else None

        if pred and pred[-1][1] >= exons[key][0][0]:
            continue
        if succ and succ[0][0] <= exons[key][-1][1]:
            continue

        e = reads_per_nt( B, seq, strand, exons[key] )
        i = reads_per_nt( B, seq, strand, introns[key] )
        if e and i:
            stdout.write( '%s\t%d\t%d\t%d\t%d\n' %
                    (transcript,e[0],i[0],e[1],i[1]) )


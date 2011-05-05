

'''
A function to build numpy arrays of read coverage over a given interval.
'''


import numpy as np
from sys                import stdout, stderr, stdin
from Bio.Seq            import Seq, reverse_complement
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA


cdef extern from 'stdlib.h':
    void* malloc( size_t )
    void free( void* )



cdef extern from 'stdint.h':
    ctypedef unsigned long long uint64_t
    ctypedef unsigned int  uint32_t
    ctypedef signed int    int32_t
    ctypedef unsigned char uint8_t


cdef extern from 'samtools/sam.h':

    ctypedef void* bamFile

    ctypedef struct bam_header_t:
        int32_t n_targets
        char** target_name
        uint32_t* target_len
        pass

    ctypedef struct _dummy_:
        bamFile bam
        pass

    ctypedef struct samfile_t:
        _dummy_ x
        bam_header_t* header
        pass

    ctypedef struct bam_index_t:
        pass

    struct __bam_iter_t:
        pass

    ctypedef __bam_iter_t* bam_iter_t


    ctypedef struct bam1_core_t:
        int32_t tid
        int32_t pos
        int32_t n_cigar
        int32_t l_qseq

        pass

    ctypedef struct bam1_t:
        bam1_core_t core
        pass

    samfile_t* samopen( char* fn, char* mode, void* aux )
    bam_index_t* bam_index_load( char* fn )
    void bam_index_destroy( bam_index_t* idx )
    void samclose( samfile_t* )

    int bam_parse_region( bam_header_t *header, \
                         char *s, int *ref_id, int *begin, int *end )

    void bam_iter_destroy( bam_iter_t )

    bam_iter_t bam_iter_query( bam_index_t *idx, int tid, int beg, int end )

    bam1_t* bam_init1()
    void bam_destroy1( bam1_t* )

    int bam_iter_read( bamFile fp, bam_iter_t iter, bam1_t *b )

    int bam1_strand( bam1_t* )
    void* bam1_seq( bam1_t* )
    char  bam1_seqi( void* s, int i )

    uint32_t* bam1_cigar( bam1_t* )

    uint32_t bam_calend( bam1_core_t *c, uint32_t *cigar )

    uint32_t BAM_CIGAR_MASK
    uint32_t BAM_CIGAR_SHIFT
    #uint8_t BAM_CMATCH




cdef extern from 'samtools/bam.h':
    int32_t bam_get_tid(bam_header_t *header, char *seq_name)


cdef extern from "bam_init_header_hash.h":
    void bam_init_header_hash(bam_header_t *header)


# To force the incomplete type into a complete type.
# This is a bit risky, as samtools may change their code.
cdef struct complete_bam_iter_t__:
    int from_first
    int tid, beg, end, n_off, i, finished
    uint64_t curr_off
    void*    off


# Convert 4bit encodings of nucleotides to characters.
cdef char bam1_basei( bam1_t* b, int i ):
    cdef int k
    k = bam1_seqi( bam1_seq(b), i )
    if k == 1:
        return b'A'
    elif k == 2:
        return b'C'
    elif k == 4:
        return b'G'
    elif k == 8:
        return b'T'
    else:
        return b'N'



# Cigar operations.
BAM_CMATCH     = 0
BAM_CINS       = 1
BAM_CDEL       = 2
BAM_CREF_SKIP  = 3
BAM_CSOFT_CLIP = 4
BAM_CHARD_CLIP = 5
BAM_CPAD       = 6




cdef class Read:
    '''
    A python class wrapping bam1_t.
    '''

    cdef bam1_t* b
    cdef object seq
    cdef object cigar

    def __cinit__( self ):
        self.b = bam_init1()
        self.cigar = None
        self.seq   = None

    def __dealloc__( self ):
        bam_destroy1( self.b )

    property strand:
        def __get__( self ):
            return int(bam1_strand(self.b))

    property pos:
        def __get__( self ):
            return self.b.core.pos

    property seq:
        def __get__( self ):
            if self.seq is not None:
                return self.seq

            cdef int32_t n = self.b.core.l_qseq
            cdef char* c_seq = <char*>malloc((n+1)*sizeof(char))
            for i in range(n):
                c_seq[i] = bam1_basei( self.b, i )
            c_seq[n] = 0

            self.seq = Seq( c_seq, IUPACAmbiguousDNA() )
            free(c_seq)

            return self.seq

    property cigar:
        def __get__( self ):
            if self.cigar is not None:
                return self.cigar

            cdef uint32_t* cigar = bam1_cigar(self.b)
            cdef int i
            cdef uint8_t  c_op
            cdef uint32_t c_len

            self.cigar = []

            for i in range(self.b.core.n_cigar):
                c_op  = cigar[i] & BAM_CIGAR_MASK
                c_len = cigar[i] >> BAM_CIGAR_SHIFT

                self.cigar.append( ( c_op, c_len) )

            return self.cigar





    # TODO: access to other fields



cdef class Bam


cdef class BamIter:
    cdef bam_iter_t it
    cdef Bam bam
    cdef int strand

    def __cinit__( self, Bam bam, \
                   chrom=None, start=None, end=None, strand=None ):

        cdef int bam_tid, bam_start, bam_end

        if   strand == '+' or strand == 0: self.strand = 0
        elif strand == '-' or strand == 1: self.strand = 1
        else: self.strand = -1

        if chrom is not None:
            if start is None: start = 0
            if end   is None: end   = bam.seqlens[chrom]

        self.bam = bam
        self.it  = NULL

        if chrom is not None and start is not None and end is not None:
            region = '%s:%d-%d' % (chrom,start,end)
            bam_parse_region( bam.reads_f.header, region,
                              &bam_tid, &bam_start, &bam_end )
            if bam_tid == -1:
                return

            self.it = bam_iter_query( bam.reads_index, bam_tid, bam_start, bam_end )
        else:
            self.it = <bam_iter_t>malloc( sizeof(complete_bam_iter_t__) )
            (<complete_bam_iter_t__*>self.it).finished   = 0
            (<complete_bam_iter_t__*>self.it).off        = NULL
            (<complete_bam_iter_t__*>self.it).from_first = 1


    def __dealloc( self ):
        bam_iter_destroy(self.it)

    def __iter__( self ):
        return self

    def __next__( self ):
        if self.it == NULL:
            return StopIteration

        cdef int m
        cdef Read read = Read()

        while True:
            m = bam_iter_read( self.bam.reads_f.x.bam, self.it, read.b )
            if m <= 0:
                raise StopIteration

            if self.strand == -1 or bam1_strand(read.b) == self.strand:
                break

        return read



cdef class Bam:
    cdef samfile_t* reads_f
    cdef bam_index_t* reads_index
    cdef object seqlens

    def __cinit__( self, fn ):
        self.reads_f = samopen( fn, 'rb', NULL )
        if( self.reads_f == NULL ):
            raise IOError( 'Can\'t open file %s\n' % fn )
            exit(1)

        self.reads_index = bam_index_load( fn )
        if( self.reads_index == NULL ):
            raise IOError( 'Can\'t open bam index %s.bai\n' % fn )
            exit(1)

        bam_init_header_hash( self.reads_f.header );

        self.seqlens = dict()
        cdef int i = 0
        for i in range(self.reads_f.header.n_targets):
            self.seqlens[ self.reads_f.header.target_name[i] ] = \
                            self.reads_f.header.target_len[i]

    property seqlens:
        def __get__( self ):
            return self.seqlens

    def __dealloc__( self ):
        bam_index_destroy( self.reads_index )
        samclose( self.reads_f )


    def reads( self, chrom = None, start = None, end = None, strand = None ):
        '''
        Iterate through all the reads, or if (chrom, start, end, strand) are
        specified, all the reads that overlap that particular region.
        '''

        return BamIter( self, chrom, start, end, strand )


    def __iter__( self ):
        return self.reads()


    def counts( self, chrom, start, end, strand=None ):

        if strand == '+': strand = 0
        if strand == '-': strand = 1

        cdef int bam_tid, bam_start, bam_end
        bam_start = start
        bam_end   = end

        bam_tid   = bam_get_tid( self.reads_f.header, chrom )
        if bam_tid == -1:
            return None

        cdef bam_iter_t it = bam_iter_query( self.reads_index, bam_tid, bam_start, bam_end )

        xs = np.zeros( end-start+1, dtype=np.float )

        if it == NULL: return xs

        cdef bam1_t* read = bam_init1()

        while bam_iter_read( self.reads_f.x.bam, it, read ) >= 0:
            if strand is not None and strand != bam1_strand(read): continue

            if bam1_strand(read) == 0:
                i = read.core.pos
                if start <= i <= end: xs[i-start] += 1
            else:
                i = bam_calend( &read.core, bam1_cigar(read) ) - 1
                if start <= i <= end: xs[i-start] += 1

        bam_destroy1(read)
        bam_iter_destroy(it)

        if strand == 1: return xs[::-1]
        else:          return xs


    def count(self, chrom, start, end, strand = None):

        if strand == '+': strand = 0
        if strand == '-': strand = 1

        cdef int bam_tid, bam_start, bam_end
        bam_start = start
        bam_end   = end

        bam_tid   = bam_get_tid( self.reads_f.header, chrom )
        if bam_tid == -1:
            return 0

        cdef bam_iter_t it = bam_iter_query( self.reads_index, bam_tid, bam_start, bam_end )

        if it == NULL: return 0

        cdef bam1_t* read = bam_init1()

        x = 0

        while bam_iter_read( self.reads_f.x.bam, it, read ) >= 0:
            if strand is not None and strand != bam1_strand(read): continue

            if bam1_strand(read) == 0:
                i = read.core.pos
                if start <= i <= end: x += 1
            else:
                i = bam_calend( &read.core, bam1_cigar(read) ) - 1
                if start <= i <= end: x += 1

        bam_destroy1(read)
        bam_iter_destroy(it)

        return x


    def coverage( self, chrom, start, end, strand=None ):

        if strand == '+': strand = 0
        if strand == '-': strand = 1

        cdef int bam_tid, bam_start, bam_end
        region = '%s:%d-%d' % (chrom,start,end)
        bam_parse_region( self.reads_f.header, region,
                          &bam_tid, &bam_start, &bam_end )
        if bam_tid == -1:
            return None

        cdef bam_iter_t it = bam_iter_query( self.reads_index, bam_tid, bam_start, bam_end )

        cdef bam1_t* read = bam_init1()
        cdef size_t i, j

        cdef uint32_t* cigar
        cdef int32_t pos
        cdef uint8_t op
        cdef uint32_t clen

        xs = np.zeros( end-start+1, dtype=np.float )

        while bam_iter_read( self.reads_f.x.bam, it, read ) >= 0:
            if strand is not None and strand != bam1_strand(read): continue

            cigar = bam1_cigar(read)
            pos   = read.core.pos

            for i in range(read.core.n_cigar):
                if pos > end: break
                op  = cigar[i] & BAM_CIGAR_MASK
                clen = cigar[i] >> BAM_CIGAR_SHIFT

                if op == BAM_CMATCH:
                    for j in range(clen):
                        if start <= pos <= end: xs[pos - start] += 1
                        pos += 1
                else:
                    pos += clen


        bam_destroy1(read)
        bam_iter_destroy(it)

        return xs




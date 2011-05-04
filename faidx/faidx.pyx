

'''
Accessing indexed fasta (faidx) files, using samtools.


Daniel Jones <dcjones@cs.washington.edu>
Nov. 10, 2010
'''


from sys                import stderr
from Bio.Seq            import Seq, reverse_complement
from Bio.Alphabet.IUPAC import IUPACAmbiguousDNA


cdef extern from 'stdlib.h':
    void free( void* )



cdef extern from 'samtools/faidx.h':
    ctypedef struct faidx_t:
        pass

    faidx_t* fai_load( char* fn )
    void fai_destroy( faidx_t* )
    char* faidx_fetch_seq( faidx_t* fai, char* name, int start, int end, int* seq_len )




cdef class Faidx:
    cdef faidx_t* fai

    def __cinit__( self, fn ):
        self.fai = fai_load( fn )
        if self.fai == NULL:
            raise IOError( 'Can\'t open faidx file %s.' % fn )


    def __dealloc__( self ):
        if self.fai != NULL:
            fai_destroy( self.fai )


    def fetch( self, seqname, start, end, strand = '+' ):

        cdef int len_out
        cdef char* c_seq

        c_seq = faidx_fetch_seq( self.fai, seqname, start, end, &len_out )

        if c_seq == NULL:
            raise Exception('Invalid region: %s(%s):%d-%d' % \
                                (seqname, strand, start, end) )

        seq = Seq( c_seq, IUPACAmbiguousDNA() )
        free(<void*>c_seq)

        if strand != '+' and strand != 0:
            seq = reverse_complement( seq )

        return seq






from libc.stdlib cimport *
from libc.stdio cimport *


cdef extern from "math.h":
    int isfinite(double x)


cdef extern from "common.h":
    ctypedef long pos_t

    ctypedef enum strand_t:
        strand_pos
        strand_neg
        strand_na


cdef extern from "str_map.h":

    ctypedef struct str_map_pair:
        char*         key
        size_t        keylen
        void*         value
        str_map_pair* next

    ctypedef struct str_map:
        str_map_pair** A
        size_t n
        size_t m
        size_t max_m


cdef extern from "gtf_parse.h":
    ctypedef struct str_t:
        char*  s
        size_t n
        size_t size

    ctypedef struct gtf_row_t:
        str_t*   seqname
        str_t*   source
        str_t*   feature
        pos_t    start
        pos_t    end
        double   score
        strand_t strand
        int      frame
        str_map* attributes

    gtf_row_t* gtf_row_alloc()
    void       gtf_row_free(gtf_row_t*)


    ctypedef struct gtf_file_t:
        pass

    gtf_file_t* gtf_file_alloc(FILE*)
    void        gtf_file_free(gtf_file_t*)
    bint        gtf_next(gtf_file_t*, gtf_row_t*)



import sys



cdef class gtf_row:
    '''
    A representation of a GTF row.
    '''

    cdef public str seqname
    cdef public str source
    cdef public str feature
    cdef public int start
    cdef public int end
    cdef public double score
    cdef strand_t strand
    cdef public int  frame
    cdef public dict attributes

    property strand:
        def __get__(self):
            if   self.strand == strand_pos: return '+'
            elif self.strand == strand_neg: return '-'
            else: return '.'

        def __set__(self, x):
            if x == '+' or x == 0: self.strand = strand_pos
            elif x == '-' or x > 1: self.strand = strand_neg
            else: self.strand = strand_na


    def __cinit__(self):
        self.attributes = dict()
        pass


    def __repr__(self):
        return self.print_line()


    def clear(self):
        self.attributes.clear()


    cdef void set(self, gtf_row_t* row):
        self.seqname = row.seqname.s
        self.source  = row.source.s
        self.feature = row.feature.s
        self.start   = row.start
        self.end     = row.end
        self.score   = row.score
        self.frame   = row.frame

        self.attributes.clear()

        cdef str k, v

        cdef size_t i
        cdef str_map_pair* u
        for i in range(row.attributes.n):
            u = row.attributes.A[i]
            while u != NULL:
                if u.value != NULL:
                    k = <str>u.key[:u.keylen]
                    v = <str>(<str_t*>u.value).s[:(<str_t*>u.value).n]
                    self.attributes[k] = v
                u = u.next


    def print_line(self):
        s = '{seqname}\t{source}\t{feature}\t{start}\t{end}\t' \
                '{score}\t{strand}\t{frame}\t{attributes}\n'

        attr_str   = ' '.join(['{k} "{v}";'.format(k = k, v = v)
                                  for (k, v) in self.attributes.iteritems()])
        score_str  = ('%0.4e' % self.score) if isfinite(self.score) else '.'
        if self.strand == strand_pos:   strand_str = '+'
        elif self.strand == strand_neg: strand_str = '-'
        else:                          strand_str = '.'
        frame_str  = ('%d' % self.frame) if self.frame >= 0 else '.'

        return s.format(seqname    = self.seqname,
                        source     = self.source,
                        feature    = self.feature,
                        start      = self.start,
                        end        = self.end,
                        score      = score_str,
                        strand     = strand_str,
                        frame      = frame_str,
                        attributes = attr_str)





cdef class gtf_file:
    '''
    A generator for GTF rows.
    '''


    cdef FILE* file
    cdef gtf_file_t* f
    cdef gtf_row_t* row

    def __cinit__(self, fn):
        cdef char* c_fn = fn

        self.file = fopen(c_fn, 'r')
        if self.file == NULL:
            raise IOError('Can\'t open file %r' % fn)

        self.f   = gtf_file_alloc(self.file)
        self.row = gtf_row_alloc()


    def __dealloc__( self ):
        gtf_row_free(self.row)
        gtf_file_free(self.f)
        fclose(self.file)


    def __next__( self ):
        if gtf_next(self.f, self.row) <= 0:
            raise StopIteration

        r = gtf_row()
        r.set(self.row)
        return r

    def __iter__(self):
        return self



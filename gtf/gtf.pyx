
from libc.stdlib cimport *
from libc.stdio cimport *

cdef extern from 'string.h':
    char* strtok_r( char* s, char* delim, char** saveptr )

cdef extern from 'ctype.h':
    int isspace( char )

import sys



cdef size_t max_line = 4096

ctypedef long pos

cdef class gtf_row:
    '''
    A representation of a GTF row.
    '''

    cdef public str seqname
    cdef public str source
    cdef public str feature
    cdef public int start
    cdef public int end
    cdef public str score
    cdef int strand
    cdef public int  frame
    cdef public dict attributes

    property strand:
        def __get__( self ):
            if self.strand == 0: return '+'
            elif self.strand > 0: return '-'
            else: return '.'

        def __set__( self, x ):
            if x == '+' or x == 0: self.strand = 0
            elif x == '-' or x > 1: self.strand = 1
            else: self.strand = -1


    def __cinit__( self ):
        self.attributes = dict()
        pass


    def __repr__( self ):
        cdef char* buf = <char*>malloc( max_line )
        self.write( buf )
        pybuf = buf
        free( <void*>buf )

        return pybuf


    def clear( self ):
        self.attributes.clear()


    cdef int parse( self, char* line ):
        if line == NULL or line[0] == b'\0' or line[0] == b'\n': return False

        self.clear()

        cdef int k = 0
        cdef char* saveptr
        cdef char *s, *c

        cdef char *key_start, *key_end, *val_start, *val_end
        cdef str key, val

        s = strtok_r( line, '\t\n',  &saveptr )
        while s:
            if k == 0:
                self.seqname = s
            elif k == 1:
                self.source = s
            elif k == 2:
                self.feature = s
            elif k == 3:
                self.start = <pos>atoi(s)
            elif k == 4:
                self.end   = <pos>atoi(s)
            elif k == 5:
                self.score = s
            elif k == 6:
                if s[0] == b'+':
                    self.strand = 0
                elif s[0] == b'-':
                    self.strand = 1
                else:
                    self.strand = -1
            elif k == 7:
                if s[0] == b'.':
                    self.frame = -1
                else:
                    self.frame = atoi(s)

            # attribute parsing
            else:

                key_start = key_end = val_start = val_end = NULL
                c =  s
                while c[0] != b'\0':
                    if key_start == NULL:
                        if not isspace(c[0]): key_start = c
                    elif key_end == NULL:
                        if isspace(c[0]): key_end = c
                    elif val_start == NULL:
                        if not isspace(c[0]) and c[0] != b'\"':
                            val_start = c
                    elif val_end == NULL:
                        if isspace(c[0]) or c[0] == b'\"':
                            val_end = c
                    else:
                        break

                    c += 1

                if key_start == NULL and \
                     key_end == NULL and \
                   val_start == NULL and \
                     val_end == NULL:
                         break

                if key_start == NULL or \
                   key_end == NULL or \
                   val_start == NULL or \
                   val_end == NULL:
                       sys.stderr.write( 'Malformed GTF attribute: %s\n' % s )

                key = str(key_start[ : key_end - key_start ])
                val = str(val_start[ : val_end - val_start ])
                self.attributes[key] = val
            k += 1
            if k < 8: s = strtok_r( NULL, '\t\n', &saveptr )
            else:     s = strtok_r( NULL, ';\n', &saveptr )

        return True


    cdef int parse_bed( self, char* line ):
        if line == NULL or line[0] == b'\0' or line[0] == b'\n': return False

        self.clear()

        cdef char* saveptr
        cdef char* s
        cdef int i = 0

        s = strtok_r( line, '\t\n', &saveptr )
        while s:
            if i == 0: self.seqname = s
            elif i == 1: self.start = <pos>atoi(s)
            elif i == 2: self.end = <pos>atoi(s)
            elif i == 3: self.feature = s
            elif i == 4: self.score = s
            elif i == 5:
                if s[0] == b'+':
                    self.strand = 0
                elif s[0] == b'-':
                    self.strand = 1
                else:
                    self.strand = -1
            else: break

            s = strtok_r( NULL, '\t\n', &saveptr )
            i += 1

        return True

    cdef bint write( self, char* buf ):
        pass # TODO



cdef class gtf_file:
    '''
    A generator for GTF rows.
    '''


    cdef FILE* f
    cdef char* buf

    def __cinit__( self, fn ):
        cdef char* c_fn = fn

        self.f = fopen( c_fn, 'r' )
        if self.f == NULL:
            raise IOError( 'Can\'t open file %r' % fn )

        self.buf = <char*>malloc( max_line )


    def __dealloc__( self ):
        fclose(self.f)
        free( <void*>self.buf )



    def __next__( self ):
        if fgets( self.buf, max_line, self.f ) == NULL: raise StopIteration

        cdef gtf_row row = gtf_row()
        row.parse( self.buf )

        return row


    def __iter__( self ):
        return self


cdef class bed_file:
    '''
    Generate GTF rows from a bed file.
    '''

    cdef FILE* f
    cdef char* buf


    def __cinit__( self, fn ):

        cdef char* c_fn = fn

        self.f = fopen( c_fn, 'r' )
        if self.f == NULL:
            raise IOError( 'Can\'t open file %r' % fn )

        self.buf = <char*>malloc( max_line )


    def __dealloc__( self ):
        fclose(self.f)
        free( <void*>self.buf )



    def __next__( self ):
        if fgets( self.buf, max_line, self.f ) == NULL: raise StopIteration

        cdef gtf_row row = gtf_row()

        row.parse_bed( self.buf )

        return row


    def __iter__( self ):
        return self




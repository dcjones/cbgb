
from libc.stdlib cimport malloc, free
import sys


cdef extern from "stdlib.h":
    void srandom(unsigned int seed)


cdef extern from 'ushuffle.h':
    void c_shuffle  "shuffle"  (char *s, char *t, int l, int k)




def seed(long s):
    srandom(s)


def shuffle(bytes s, k = 2):
    '''
    ushuffle(seq, k = 2) -> seq

    Generate a uniform random permutation of a biological sequenc (e.g. DNA,
    RNA, AA) than preserves the exect k-let count.
    '''

    l = len(s)

    cdef bytes out
    cdef char* c_out = <char*>malloc((l + 1) * sizeof(char))
    c_out[l] = b'\0'

    c_shuffle(s, c_out, l, k)

    out = c_out
    free(c_out)

    return out




#include "common.h"
#include <stdio.h>


void* malloc_or_die(size_t n)
{
    void* ptr = malloc(n);
    if (ptr == NULL) {
        fprintf(stderr, "Falied to allocate %zu bytes. Out of memory.\n", n);
        return NULL;
    }

    return ptr;
}




void* realloc_or_die(void* ptr, size_t n)
{
    ptr = realloc(ptr, n);
    if (ptr == NULL) {
        fprintf(stderr,"Falied to (re)allocate %zu bytes. Out of memory.\n", n);
        return NULL;
    }

    return ptr;
}


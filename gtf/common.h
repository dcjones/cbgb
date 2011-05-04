
#ifndef ISOLATOR_COMMON_H
#define ISOLATOR_COMMON_H

#ifdef __cplusplus
extern "C" {
#endif

#include <stdlib.h>

/* genomic position */
typedef long pos_t;

/* sequence identifier */
typedef int seqid_t;

typedef enum {
    strand_pos = 0,
    strand_neg = 1,
    strand_na  = 2
}strand_t;


void* malloc_or_die(size_t);
void* realloc_or_die(void*, size_t);


#ifdef __cplusplus
}
#endif

#endif




#ifndef BAMHASH_HASH
#define BAMHASH_HASH


#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <samtools/sam.h>

struct hashed_value
{
    char*    seq;
    uint32_t count;

    struct hashed_value* next;
};


/* Hash table structure. */
struct table
{
    struct hashed_value** A; /* table proper */
    size_t n;                /* table size */
    size_t m;                /* hashed items */
    size_t max_m;            /* max hashed items before rehash */
    size_t min_m;            /* min hashed items before rehash */
    size_t read_len;         /* fixed read length */
    size_t read_bytes;       /* number of bytes used to encode read sequence */
};


void table_create( struct table* T, size_t read_len );
void table_destroy( struct table* T );

void     table_inc( struct table*, char* seq );
uint64_t table_get( struct table*, char* seq );

void sort_by_count( struct table* T,
                    struct hashed_value*** S );


#endif



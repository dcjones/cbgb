

#ifndef BAMHASH_HASH
#define BAMHASH_HASH


#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <samtools/sam.h>

struct hashed_value
{
    uint8_t* seq;
    int32_t tid;
    int32_t pos;

    uint32_t pos_count;
    uint32_t neg_count;

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

void table_inc( struct table*, bam1_t* read );



void sort_by_count( struct table* T,
                    struct hashed_value*** S_pos,
                    struct hashed_value*** S_neg );



#endif



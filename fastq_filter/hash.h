

#ifndef BAMHASH_HASH
#define BAMHASH_HASH


#include <stdint.h>
#include <stdbool.h>
#include <stdlib.h>
#include <samtools/sam.h>

struct hashed_value
{
    char* read_id;
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
};


void table_create( struct table* T );
void table_destroy( struct table* T );

void     table_add( struct table*, char* read_id );
bool     table_member( struct table*, char* read_id );


#endif



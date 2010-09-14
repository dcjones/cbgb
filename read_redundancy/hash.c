

#include "hash.h"
#include "superfasthash.h"


#define INITIAL_TABLE_SIZE 128
#define MAX_LOAD 0.75
#define MIN_LOAD 0.05 /* make sure this is less than MAX_LOAD/2 */


/* Test string equality. (Slightly faster than memcmp, since we don't care about
 * ordering. Note: n is in nibbles! */
bool memeq( const void* x, const void* y, size_t n )
{
    uint8_t* a = (uint8_t*)x;
    uint8_t* b = (uint8_t*)y;

    if( a == b ) return true;

    while( n >= 2*sizeof(uint32_t) ) {
        if( *((uint32_t*)a) != *((uint32_t*)b) ) return false;
        a += sizeof(uint32_t);
        b += sizeof(uint32_t);
        n -= 2*sizeof(uint32_t);
    }

    while( n ) {
        if( (*a & 0xf) != (*b & 0xf) )   return false;
        if( --n == 0 ) break;
        if( (*a & 0xf0) != (*b & 0xf0) ) return false;

        a++;
        b++;
        n--;
    }

    return true;
}


void rehash( struct table* T, size_t new_n );



/* Create an empty hash table. */
void table_create( struct table* T, size_t read_len )
{
    T->A = malloc(sizeof(struct hashed_value*)*INITIAL_TABLE_SIZE);
    memset( T->A, 0, sizeof(struct hashed_value*)*INITIAL_TABLE_SIZE );
    T->n = INITIAL_TABLE_SIZE;
    T->m = 0;
    T->max_m = T->n * MAX_LOAD;
    T->min_m = T->n * MIN_LOAD;
    T->read_len = read_len;
    T->read_bytes = read_len * sizeof(char);
}



/* Remove all elements in the table. */
void table_clear( struct table* T )
{
    struct hashed_value* u;
    size_t i;
    for( i = 0; i < T->n; i++ ) {
        while( T->A[i] ){
            u = T->A[i]->next;
            free(T->A[i]->seq);
            free(T->A[i]);
            T->A[i] = u;
        }
    }
    T->m = 0;
}



/* Free all memory associated with a table. */
void table_destroy( struct table* T )
{
    table_clear(T);
    free(T->A);
}



void table_inc( struct table* T, char* seq )
{
    if( T->m == T->max_m ) rehash( T, T->n*2 );

    uint32_t h = hash(seq, T->read_bytes) % T->n;

    struct hashed_value* u = T->A[h];

    while(u) {
        if( memeq( u->seq, seq, T->read_bytes ) ) {
            u->count++;
            return;
        }

        u = u->next;
    }

    u = malloc(sizeof(struct hashed_value));
    u->seq = strdup(seq);

    u->count = 1;

    u->next = T->A[h];
    T->A[h] = u;

    T->m++;
}

uint64_t table_get( struct table* T, char* seq )
{
    uint32_t h = hash(seq, T->read_bytes) % T->n;

    struct hashed_value* u = T->A[h];
    while(u) {
        if( memeq( u->seq, seq, T->read_bytes ) ) {
            return u->count;
        }

        u = u->next;
    }

    return 0;
}



/* Insert existing entries without copying sequences. Used for rehashing. */
bool table_insert_without_copy( struct table* T, struct hashed_value* V )
{
    if( T->m == T->max_m ) rehash( T, T->n*2 );

    uint32_t h = hash(V->seq,T->read_bytes) % T->n;

    V->next = T->A[h];
    T->A[h] = V;

    T->m++;

    return true;
}


/* Rezise the table T to new_n. */
void rehash( struct table* T, size_t new_n )
{
    fprintf( stderr, "\t(rehashing)\n" );
    struct table U;
    U.n = new_n;
    U.m = 0;
    U.max_m = U.n*MAX_LOAD;
    U.min_m = U.n*MIN_LOAD;
    U.A = malloc( sizeof(struct hashed_value*) * U.n );
    U.read_len   = T->read_len;
    U.read_bytes = T->read_bytes;
    memset( U.A, 0, sizeof(struct hashed_value*) * U.n );


    struct hashed_value *j,*k;
    size_t i;
    for( i = 0; i < T->n; i++ ) {
        j = T->A[i];
        while( j ) {
            k = j->next;
            table_insert_without_copy( &U, j );
            j = k;
        }
        T->A[i] = NULL;
    }

    free(T->A);
    T->A = U.A;
    T->n = U.n;
    T->max_m = T->n*MAX_LOAD;
    T->min_m = T->n*MIN_LOAD;
}

int comp_hashed_value( const void* x, const void* y )
{
    struct hashed_value* const * a = x;
    struct hashed_value* const * b = y;

    if( (*a)->count > (*b)->count ) return -1;
    if( (*a)->count < (*b)->count ) return 1;
    return 0;
}

void sort_by_count( struct table* T,
                    struct hashed_value*** S_ )
{
    struct hashed_value** S = malloc( sizeof(struct hashed_value*) * T->m );
    memset( S, 0, sizeof(struct hashed_value*) * T->m );

    struct hashed_value* j;
    size_t i,k;
    for( i=0, k=0; i < T->n; i++ ) {
        j = T->A[i];
        while( j ) {
            S[k] = j;
            k++;
            j = j->next;
        }
    }

    qsort( S, T->m, sizeof(struct hashed_value*), comp_hashed_value );
    *S_ = S;
}


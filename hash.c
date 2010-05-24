

#include "hash.h"


#define INITIAL_TABLE_SIZE 128
#define MAX_LOAD 0.75
#define MIN_LOAD 0.05 /* make sure this is less than MAX_LOAD/2 */



/* This is Jenkin's hash. The implementation is stolen from the Linux kernel,
 * modified only slightly.
 */

#define __jhash_mix(a, b, c) \
{ \
  a -= b; a -= c; a ^= (c>>13); \
  b -= c; b -= a; b ^= (a<<8); \
  c -= a; c -= b; c ^= (b>>13); \
  a -= b; a -= c; a ^= (c>>12);  \
  b -= c; b -= a; b ^= (a<<16); \
  c -= a; c -= b; c ^= (b>>5); \
  a -= b; a -= c; a ^= (c>>3);  \
  b -= c; b -= a; b ^= (a<<10); \
  c -= a; c -= b; c ^= (b>>15); \
}

/* The golden ration: an arbitrary value */
#define JHASH_GOLDEN_RATIO	0x9e3779b9

/* The most generic version, hashes an arbitrary sequence
 * of bytes.  No alignment or length assumptions are made about
 * the input key.
 */
uint32_t hash( const void *key, uint32_t length )
{
	uint32_t a, b, c, len;
	const uint8_t *k = key;

	len = length;
	a = b = JHASH_GOLDEN_RATIO;

    c = 0;

	while (len >= 12) {
		a += (k[0] +((uint32_t)k[1]<<8) +((uint32_t)k[2]<<16) +((uint32_t)k[3]<<24));
		b += (k[4] +((uint32_t)k[5]<<8) +((uint32_t)k[6]<<16) +((uint32_t)k[7]<<24));
		c += (k[8] +((uint32_t)k[9]<<8) +((uint32_t)k[10]<<16)+((uint32_t)k[11]<<24));

		__jhash_mix(a,b,c);

		k += 12;
		len -= 12;
	}

	c += length;
	switch (len) {
	case 11: c += ((uint32_t)k[10]<<24);
	case 10: c += ((uint32_t)k[9]<<16);
	case 9 : c += ((uint32_t)k[8]<<8);
	case 8 : b += ((uint32_t)k[7]<<24);
	case 7 : b += ((uint32_t)k[6]<<16);
	case 6 : b += ((uint32_t)k[5]<<8);
	case 5 : b += k[4];
	case 4 : a += ((uint32_t)k[3]<<24);
	case 3 : a += ((uint32_t)k[2]<<16);
	case 2 : a += ((uint32_t)k[1]<<8);
	case 1 : a += k[0];
	};

	__jhash_mix(a,b,c);

	return c;
}




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
    T->read_bytes = ((read_len-1)/2)+1;
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



void table_inc( struct table* T, bam1_t* read )
{
    if( T->m == T->max_m ) rehash( T, T->n*2 );

    uint32_t h = hash(bam1_seq(read), T->read_bytes) % T->n;

    struct hashed_value* u = T->A[h];

    while(u) {
        if( memeq( u->seq, bam1_seq(read), T->read_len ) ) {

            if( u->tid != read->core.tid || u->pos != read->core.pos ) {
                fputs( "Warning: inconsistant mapping found.\n", stderr );
            }

            if( bam1_strand(read) == 1 ) u->pos_count++;
            else                         u->neg_count++;
            return;
        }

        u = u->next;
    }

    u = malloc(sizeof(struct hashed_value));
    u->seq = malloc(T->read_bytes);
    memcpy( u->seq, bam1_seq(read), T->read_bytes );
    u->tid    = read->core.tid;
    u->pos    = read->core.pos;

    if( bam1_strand(read) == 1 ) { u->pos_count = 1; u->neg_count = 0; }
    else                         { u->pos_count = 0; u->neg_count = 1; }

    u->next = T->A[h];
    T->A[h] = u;

    T->m++;
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
    struct table U;
    U.n = new_n;
    U.m = 0;
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



int comp_hashed_value_pos( const void* x, const void* y )
{
    struct hashed_value* const * a = x;
    struct hashed_value* const * b = y;

    if( (*a)->pos_count < (*b)->pos_count ) return 1;
    if( (*a)->pos_count > (*b)->pos_count ) return -1;
    return 0;
}

int comp_hashed_value_neg( const void* x, const void* y )
{
    struct hashed_value* const * a = x;
    struct hashed_value* const * b = y;

    if( (*a)->neg_count < (*b)->neg_count ) return 1;
    if( (*a)->neg_count > (*b)->neg_count ) return -1;
    return 0;
}

void sort_by_count( struct table* T,
                    struct hashed_value*** S_pos_,
                    struct hashed_value*** S_neg_ )
{
    struct hashed_value** S_pos = malloc( sizeof(struct hashed_value*) * T->m );
    memset( S_pos, 0, sizeof(struct hashed_value*) * T->m );

    struct hashed_value** S_neg = malloc( sizeof(struct hashed_value*) * T->m );
    memset( S_neg, 0, sizeof(struct hashed_value*) * T->m );


    /* read off values from the table */
    struct hashed_value* j;
    size_t i,k;
    for( i=0, k=0; i < T->n; i++ ) {
        j = T->A[i];
        while( j ) {
            S_pos[k] = j;
            S_neg[k] = j;
            k++;
            j = j->next;
        }
    }


    /* sort */
    qsort( S_pos, T->m, sizeof(struct hashed_value*), comp_hashed_value_pos );
    qsort( S_neg, T->m, sizeof(struct hashed_value*), comp_hashed_value_neg );

    *S_pos_ = S_pos;
    *S_neg_ = S_neg;
}

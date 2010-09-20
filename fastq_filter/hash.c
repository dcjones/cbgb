

#include "hash.h"
#include "superfasthash.h"


#define INITIAL_TABLE_SIZE 128
#define MAX_LOAD 0.75
#define MIN_LOAD 0.05 /* make sure this is less than MAX_LOAD/2 */


void rehash( struct table* T, size_t new_n );

/* Create an empty hash table. */
void table_create( struct table* T )
{
    T->A = malloc(sizeof(struct hashed_value*)*INITIAL_TABLE_SIZE);
    memset( T->A, 0, sizeof(struct hashed_value*)*INITIAL_TABLE_SIZE );
    T->n = INITIAL_TABLE_SIZE;
    T->m = 0;
    T->max_m = T->n * MAX_LOAD;
    T->min_m = T->n * MIN_LOAD;
}



/* Remove all elements in the table. */
void table_clear( struct table* T )
{
    struct hashed_value* u;
    size_t i;
    for( i = 0; i < T->n; i++ ) {
        while( T->A[i] ){
            u = T->A[i]->next;
            free(T->A[i]->read_id);
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



void table_add( struct table* T, char* read_id )
{
    if( T->m == T->max_m ) rehash( T, T->n*2 );

    size_t read_id_len = strlen(read_id);

    uint32_t h = hash( read_id, read_id_len ) % T->n;

    struct hashed_value* u = T->A[h];

    while(u) {
        if( strcmp( u->read_id, read_id ) == 0 ) {
            return;
        }

        u = u->next;
    }

    u = malloc(sizeof(struct hashed_value));
    u->read_id = strdup(read_id);

    u->next = T->A[h];
    T->A[h] = u;

    T->m++;
}

bool table_member( struct table* T, char* read_id )
{
    size_t read_id_len = strlen(read_id);
    uint32_t h = hash( read_id, read_id_len ) % T->n;

    struct hashed_value* u = T->A[h];
    while(u) {
        if( strcmp( u->read_id, read_id ) == 0 ) {
            return true;
        }

        u = u->next;
    }

    return false;
}



/* Insert existing entries without copying sequences. Used for rehashing. */
bool table_insert_without_copy( struct table* T, struct hashed_value* V )
{
    if( T->m == T->max_m ) rehash( T, T->n*2 );

    size_t read_id_len = strlen(V->read_id);
    uint32_t h = hash( V->read_id, read_id_len ) % T->n;

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


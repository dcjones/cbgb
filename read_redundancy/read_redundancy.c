
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <zlib.h>
#include <samtools/sam.h>
#include "hash.h"
#include "dibase_decode.h"

const size_t MAX_LINE_WIDTH=4096;
bool ignore_N = true;
bool decode_colorspace = false;

typedef union {
    gzFile     rawf;
    samfile_t* samf;
} READ_FILE;


typedef size_t (*readget)( READ_FILE* f, char* read );

size_t csfasta_getread( READ_FILE* f, char* read );
size_t fastq_getread  ( READ_FILE* f, char* read );
size_t sam_getread    ( READ_FILE* f, char* read );


struct table* hash_reads( READ_FILE* f, readget getread, size_t* num_ignored )
{
    fprintf( stderr, "hashing reads ...\n" );

    char* read = malloc(MAX_LINE_WIDTH*sizeof(char*));
    size_t n, m;
    *num_ignored = 0;

    /* get one read to guess the read length */
    if( (m = getread( f, read )) == 0 ) return NULL;

    /* create table */
    struct table* T = malloc(sizeof(struct table));
    table_create( T, m );

    n = 0;

    /* hash */
    do {
        if( read[0] == '\0' ) {
            (*num_ignored)++;
            continue;
        }
        n++;
        table_inc( T, read );
        if( n % 100000 == 0 ) fprintf( stderr, "\t%zu reads (%zu unique, %zu ignored).\n",
                                       n, T->m, *num_ignored );
    } while( getread( f, read ) );


    free(read);

    fprintf( stderr, "done. (%zu reads hashed, %zu are unique, %zu were ignored)\n",
                      n, T->m, *num_ignored );

    return T;
}


void usage()
{
    fprintf( stderr,
             "Usage: read_redundancy [OPTIONS] input_file\n\n"
             "-C, -Q, -S, -B        input is csfasta, fastq, sam, or bam, respectively.\n"
             "                      exactly one must be specified!\n" 
             "-n                    count reads containing 'N' or '.' characters\n"
             "-d                    decode colorspace reads\n\n"
             );
}


void file_not_found( const char* fn )
{
    fprintf( stderr, "Could not open file '%s' for reading.\n", fn );
    exit(1);
}

void not_implemented( const char* func )
{
    fprintf( stderr, "Function '%s' not yet implemented.\n", func );
    exit(1);
}



int main( int argc, char* argv[] )
{
    size_t num_ignored = 0;
    int C_opt, Q_opt, S_opt, B_opt;
    C_opt = Q_opt = S_opt = B_opt = 0;

    int c;

    const char* optstring = "CQSBnd";
    do {
        c = getopt( argc, argv, optstring );
        switch( c ) {
            case 'C': C_opt = 1; break;
            case 'Q': Q_opt = 1; break;
            case 'S': S_opt = 1; break;
            case 'B': B_opt = 1; break;
            case 'n': ignore_N = false; break;
            case 'd': decode_colorspace = true; break;
        }
    } while( c != -1 );

    if( optind >= argc ) {
        usage();
        fprintf( stderr, "Too few arguments.\n" );
        exit(1);
    }

    if( C_opt + Q_opt + S_opt + B_opt != 1 ) {
        usage();
        fprintf( stderr, "Exactly one input option must be specified.\n" );
        exit(1);
    }

    const char* fn = argv[optind];
    READ_FILE f;
    struct table* T = NULL;

    if( C_opt ) {
        if( (f.rawf = gzopen( fn, "r" )) == NULL ) file_not_found( fn );
        T = hash_reads( &f, csfasta_getread, &num_ignored );
        gzclose( f.rawf );
    }
    else if( Q_opt ) {
        if( (f.rawf = gzopen( argv[optind], "r" )) == NULL ) file_not_found( fn );
        T = hash_reads( &f, fastq_getread, &num_ignored );
        gzclose( f.rawf );
    }
    else if( S_opt ) {
        if( (f.samf = samopen( argv[optind], "r", NULL )) == NULL ) file_not_found( fn );
        T = hash_reads( &f, sam_getread, &num_ignored );
        samclose(f.samf);
    }
    else if( B_opt ) {
        if( (f.samf = samopen( argv[optind], "rb", NULL )) == NULL ) file_not_found( fn );
        T = hash_reads( &f, sam_getread, &num_ignored );
        samclose(f.samf);
    }

    if( T == NULL ) {
        fprintf( stderr, "hashreads inexplicably failed.\n" );
        exit(1);
    }


    struct hashed_value** S;

    fprintf( stderr, "sorting ... " );
    sort_by_count( T, &S );
    fprintf( stderr, "done.\n" );

    fprintf( stderr, "printing ... " );
    size_t i;
    for( i = 0; i < T->m; i++ ) {
        fprintf( stdout, "%s\t%d\n", S[i]->seq, S[i]->count );
    }
    fprintf( stdout, "ignored\t%zu\n", num_ignored );

    fprintf( stderr, "done.\n" );

    fprintf( stderr, "dismantling table ... " );
    free(S);
    table_destroy( T );
    fprintf( stderr, "done.\n" );
    return 0;
}


/* make one pass through the string to determine the length and whether it
 * contains the given character. */
bool strchrlen( const char* s, char c, size_t* len ) 
{
    bool  c_found = false;
    const char* s0 = s;

    while( *s ) if( *s++ == c ) c_found = true;

    *len = s - s0;
    return c_found;
}


size_t csfasta_getread( READ_FILE* f, char* read )
{
    /* get read name */
    if( !gzgets( f->rawf, read, MAX_LINE_WIDTH ) ) return 0;

    /* get read sequence */
    if( !gzgets( f->rawf, read, MAX_LINE_WIDTH ) ) return 0;

    size_t n = 0;
    bool N_found = strchrlen( read, '.', &n );

    if( N_found ) read[0] = '\0';
    read[ n-1 ] = '\0'; /* trim newline */

    if( decode_colorspace ) dibase_decode( read, read );

    return n-1;
};




size_t fastq_getread  ( READ_FILE* f, char* read )
{
    not_implemented( __func__ );
    return 0;
};

size_t sam_getread    ( READ_FILE* f, char* read )
{
    not_implemented( __func__ );
    return 0;
};




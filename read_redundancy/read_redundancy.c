
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <samtools/sam.h>
#include "hash.h"

const size_t MAX_LINE_WIDTH=4096;

typedef union {
    FILE*      rawf;
    samfile_t* samf;
} READ_FILE;


typedef bool (*readget)( READ_FILE* f, char* read );

bool csfasta_getread( READ_FILE* f, char* read );
bool fastq_getread  ( READ_FILE* f, char* read );
bool sam_getread    ( READ_FILE* f, char* read );


struct table* hash_reads( READ_FILE* f, readget getread  )
{
    fprintf( stderr, "hashing reads ...\n" );

    char* read = malloc(MAX_LINE_WIDTH*sizeof(char*));
    size_t n, m;

    /* get one read to guess the read length */
    if( !getread( f, read ) ) return NULL;
    m = strlen(read);

    /* create table */
    struct table* T = malloc(sizeof(struct table));
    table_create( T, m );

    n = 0;

    /* hash */
    do {
        n++;
        table_inc( T, read );
        if( n % 100000 == 0 ) fprintf( stderr, "\t%d reads.\n", n );
    } while( getread( f, read ) );


    free(read);

    fprintf( stderr, "done. (%d reads hashed, %d are unique)\n", n, T->m );

    return T;
}


void usage()
{
    fprintf( stderr,
             "Usage: read_redundancy [OPTIONS] input_file\n\n"
             "-C, -Q, -S, -B        input is csfasta, fastq, sam, or bam, respectively.\n"
             "                      exactly one must be specified!\n\n" 
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
    int C_opt, Q_opt, S_opt, B_opt;
    C_opt = Q_opt = S_opt = B_opt = 0;

    int c;

    const char* optstring = "CQSB";
    do {
        c = getopt( argc, argv, optstring );
        switch( c ) {
            case 'C': C_opt = 1; break;
            case 'Q': Q_opt = 1; break;
            case 'S': S_opt = 1; break;
            case 'B': B_opt = 1; break;
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
    struct table* T;

    if( C_opt ) {
        if( (f.rawf = fopen( fn, "r" )) == NULL ) file_not_found( fn );
        T = hash_reads( &f, csfasta_getread );
        fclose( f.rawf );
    }
    else if( Q_opt ) {
        if( (f.rawf = fopen( argv[optind], "r" )) == NULL ) file_not_found( fn );
        T = hash_reads( &f, fastq_getread );
        fclose( f.rawf );
    }
    else if( S_opt ) {
        if( (f.samf = samopen( argv[optind], "r", NULL )) == NULL ) file_not_found( fn );
        T = hash_reads( &f, sam_getread );
        samclose(f.samf);
    }
    else if( B_opt ) {
        if( (f.samf = samopen( argv[optind], "rb", NULL )) == NULL ) file_not_found( fn );
        T = hash_reads( &f, sam_getread );
        samclose(f.samf);
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

    fprintf( stderr, "done.\n" );

    table_destroy( T );
    return 0;
}




bool csfasta_getread( READ_FILE* f, char* read )
{
    /* get read name */
    if( !fgets( read, MAX_LINE_WIDTH, f->rawf ) ) return false;

    /* get read sequence */
    if( !fgets( read, MAX_LINE_WIDTH, f->rawf ) ) return false;

    /* trim newline */
    read[ strlen(read)-1 ] = '\0';

    return true;
};




bool fastq_getread  ( READ_FILE* f, char* read )
{
    not_implemented( __func__ );
    return false;
};

bool sam_getread    ( READ_FILE* f, char* read )
{
    not_implemented( __func__ );
    return false;
};




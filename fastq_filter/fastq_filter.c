/*
 *      fastq_filter
 *      ------------------
 *      Filter a fastq file by the read name using the mapped or unmapped reads
 *      in a SAM/BAM file.
 *
 *      September 2010 / Daniel Jones <dcjones@cs.washington.edu>
 *
 */




#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
#include <unistd.h>
#include <samtools/sam.h>
#include "hash.h"



void usage()
{
    fprintf( stderr, "Usage: fastq_filter [OPTIONS] filter.sam/bam [in.fastq] > out.fastq \n\n"
                     "Filters the reads in in.fastq, printing only those that are\n"
                     "aligned in filter.sam/bam, or only those that are unaligned\n"
                     "if the -v switch is used.\n\n"
                     "Options:\n"
                     "-v         invert, i.e. keep only unaligned reads from in.fastq\n"
                     "-S         filter is in SAM format\n"
                     "-b         filter is in BAM format (default)\n\n" );
}



void double_encode( int n, char* read_seq, char* qual_seq, int* quals )
{
    /* convert read sequence */
    char* s = read_seq+2;
    int i = 0;

    while( *s ) {
        switch( *s ) {
            case '0' : read_seq[i++] = 'A'; break;
            case '1' : read_seq[i++] = 'C'; break;
            case '2' : read_seq[i++] = 'G'; break;
            case '3' : read_seq[i++] = 'T'; break;
            case '.' : read_seq[i++] = 'N'; break;
            default: break;
        };
        s++;
    }
    read_seq[i] = '\0';

    /* convert quality sequence */
    for( i = 1; i < n; i++ ) {
        if( quals[i] == -1 ) quals[i] = 0;
        qual_seq[i-1] = (char)(quals[i]+33);
    }
    qual_seq[i-1] = '\0';
}

char* fgets_noncomment( char* buf, size_t buf_size, FILE* f )
{
    char* r;
    while( (r = fgets( buf, buf_size, f)) && r[0] == '#' );

    return r;
}



struct table* hash_read_ids( const char* filter_fn, bool input_bam, bool invert )
{
    samfile_t* filter_f = samopen( filter_fn, input_bam ? "rb" : "r", NULL );
    if( filter_f == NULL ) {
        fprintf( stderr, "Can't open SAM/BAM file '%s'.\n", filter_fn );
        exit(1);
    }

    size_t n = 0;

    struct table* T = malloc(sizeof(struct table));
    table_create( T );

    bam1_t *b = bam_init1();
    bool unmapped;

    while( samread( filter_f, b ) ) {
        n++;
        unmapped = (b->core.flag & BAM_FUNMAP) > 0;
        if( unmapped == invert ) table_add( T, bam1_qname(b) );
        if( n % 100000 == 0 ) fprintf( stderr, "\t%zd reads proccessed\n", n );
    }

    bam_destroy1(b);
    samclose( filter_f );

    return T;
}


int main( int argc, char* argv[] )
{
    bool input_bam = true;
    bool invert    = false;
    bool using_stdin = false;

    const char* optstring = "vbS";
    int opt;

    do {
        opt = getopt( argc, argv, optstring );
        switch( opt ) {
            case 'v':
                invert = true;
                break;
            case 'b':
                input_bam = true;
                break;
            case 'S':
                input_bam = false;
                break;
        }
    } while( opt != -1 );


    if( argc - optind < 1 ) {
        usage();
        fprintf( stderr, "Too few arguments.\n" );
        exit(1);
    }

    const char* filter_fn = argv[optind++];

    FILE* input_f;
    if( optind >= argc ) {
        input_f = stdin;
        using_stdin = true;
    }
    else {
        const char* input_fn  = argv[optind++];
        input_f  = fopen( input_fn, "r" );
        if( input_f == NULL ) {
            fprintf( stderr, "Can't open fastq file '%s'.\n", input_fn );
            exit(1);
        }
        using_stdin = false;
    }



    samfile_t* filter_f = samopen( filter_fn, input_bam ? "rb" : "r", NULL );
    if( filter_f == NULL ) {
        fprintf( stderr, "Can't open SAM/BAM file '%s'.\n", filter_fn );
        exit(1);
    }


    fprintf( stderr, "hashinging read id's ... " );
    struct table* T = hash_read_ids( filter_fn, input_bam, invert );
    fprintf( stderr, "done. (%zd hashed)\n", T->m );


    const size_t buf_size = 4096;

    char* buf = malloc( 4*buf_size );


    char* read_name = buf + 0*buf_size;
    char* read_seq  = buf + 1*buf_size;
    char* qual_name = buf + 2*buf_size;
    char* qual_seq  = buf + 3*buf_size;

    size_t n = 0;

    fprintf( stderr, "filtering ... " );

    /* iterate over groups of four lines */
    while( fgets_noncomment( read_name, buf_size, input_f ) &&
           fgets_noncomment( read_seq,  buf_size, input_f ) &&
           fgets_noncomment( qual_name, buf_size, input_f ) &&
           fgets_noncomment( qual_seq,  buf_size, input_f ) )
    {
        n++;        
        read_name[strlen(read_name)] = '\0';
        if( !table_member( T, read_name+1 ) ) continue;

        printf( "%s\n%s%s%s", read_name, read_seq, qual_name, qual_seq );
    }

    fprintf( stderr, "done. (%zd reads processed)\n", n );


    if( !using_stdin ) fclose( input_f );
    free(buf);
    table_destroy(T);
    free(T);

    return 0;
}





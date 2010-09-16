/*
 *           bam_binned_coverage
 *           -------------------
 *           A program to print number of reads falling within
 *           bins of size k across the echo chromosome.
 *
 */

#include <stdio.h>
#include <samtools/sam.h>
#include <samtools/bam.h>

uint32_t k = 100000;


void usage()
{
    fprintf( stderr,
             "Usage: bam_binned_coverage [OPTIONS] in1.bam [in2.bam ...]\n\n"
             "Options:\n"
             "-k K      size of bins (default: %u)\n\n", k );
}


typedef struct chrom_count_
{
    char* seqname;
    int tid;
    uint32_t n;
    uint32_t* count;
    struct chrom_count_* next;
} chrom_count;


chrom_count* init_counts( bam_header_t* header )
{
    chrom_count* C = NULL;
    chrom_count* B;
    int32_t i;
    for( i = 0; i < header->n_targets; i++ ) {
        B = C;
        C = malloc( sizeof(chrom_count) );
        C->seqname = strdup(header->target_name[i]);
        C->tid     = i;
        C->n       = (header->target_len[i]/k)+1;
        C->count   = malloc( sizeof(uint32_t)*C->n );
        memset( C->count, 0, sizeof(uint32_t)*C->n );
        C->next = B;
    }

    return C;
}

void set_tids( chrom_count* C, bam_header_t* header )
{
    chrom_count* B;
    int i;

    for( i = 0; i < header->n_targets; i++ ) {
        for( B = C; B != NULL; B = B->next ) {
            if( strcmp( B->seqname, header->target_name[i] ) == 0 ) {
                B->tid = i;
                break;
            }
        }
    }
}


void destroy_counts( chrom_count** C )
{
    chrom_count* B;
    while( *C ) {
        B = (*C)->next;
        free( (*C)->seqname );
        free( (*C)->count );
        free( *C );
        *C = B;
    }
}

void print_counts( chrom_count* C )
{
    uint32_t i;
    while( C ) {
        for( i = 0; i < C->n; i++ ) {
            printf( "%s\t%u\t%u\n",  C->seqname, i*k, C->count[i] );
        }

        C = C->next;
    }
}



int tally( const bam1_t *b, void* data )
{
    chrom_count* C = (chrom_count*)data;
    if( b->core.tid != C->tid ) return 0;
    C->count[ b->core.pos / k ]++;

    return 0;
}



int main( int argc, char* argv[] )
{
    const char* optstring = "k:";
    int c;

    do {
        c = getopt( argc, argv, optstring );
        switch( c ) {
            case 'k':
                k = atoi(optarg);
                break;
        }
    } while( c != -1 );

    if( optind >= argc ) {
        usage();
        fprintf( stderr, "Too few arguments.\n" );
        exit(1);
    }

    samfile_t* f; 
    bam_index_t* idx;

    chrom_count* C = NULL;
    chrom_count* B;

    for( ; optind < argc; optind++ ) {
        f = samopen( argv[optind], "rb", NULL );
        if( f == NULL ) {
            fprintf( stderr, "Can't open file '%s'.\n", argv[optind] );
            exit(1);
        }


        if( C == NULL ) {
            if( f->header == NULL ) {
                fprintf( stderr, "Error: bam file '%s' is missing a header.\n", argv[optind] );
                exit(1);
            }

            C = init_counts( f->header );
        }

        set_tids( C, f->header );

        idx = bam_index_load( argv[optind] );

        if( idx == NULL ) {
            bam_index_build( argv[optind] );
            idx = bam_index_load( argv[optind] );

            if( idx == NULL ) {
                fprintf( stderr, "Error: can't load index for bam file '%s'.\n", argv[optind] );
                exit(1);
            }
        }

        for( B = C; B != NULL; B = B->next ) {
            bam_fetch( f->x.bam, idx, B->tid, 0, B->n*k, (void*)B, tally );
        }

        bam_index_destroy(idx);
        samclose(f);
    }

    print_counts( C );

    return 0;
}


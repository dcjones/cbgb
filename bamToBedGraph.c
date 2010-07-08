
/*
 *                  bamToBedGraph
 *                  -------------
 *                  Create a coverage histogram from a BAM file.
 *
 *                  (Note: this is what genomeCoverageBed from BEDTools does,
 *                  but it does so incorrectly: it does not take reads with gaps
 *                  into account, and lacks to option to consider strand.)
 *
 *                  July 2010  /  Daniel Jones <dcjones@cs.washington.edu>
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <getopt.h>
#include <samtools/sam.h>



/* I. Run length encoding of an array, to covert the coverage array to bedGraph. */
void print_rle( uint32_t* A, size_t n, const char* chrom )
{
    size_t i,j;
    uint32_t k;

    i = 0;

    while( i < n ) {
        j = i;
        k = A[i];
        while( j < n && A[j] == k ) j++;

        printf( "%s\t%d\t%d\t%d\n", chrom, i, j, k );

        i = j;
    }
}


void print_usage()
{
    fprintf( stderr,
            "Usage: bamToBedGraph [options] in.bam\n"
            "Covert a BAM file to a Bed Graph file.\n\n"
            "Options:\n"
            "  -s, --strand=STRAND        Where STRAND is '+' or '-'. Count only reads aligned\n"
            "                             to this strand. (By default reads on both strands\n"
            "                             are counted.)\n"
           );
}


void count_read( bam1_t* b, uint32_t* A, size_t n )
{
    uint32_t* cigar = bam1_cigar(b);

    int32_t pos = b->core.pos;

    size_t i, j;
    uint8_t  op;
    uint32_t len;

    /* TODO: must I do something else to handle strand properly?? */

    for( i = 0; i < b->core.n_cigar; i++ ) {
        op  = cigar[i] & BAM_CIGAR_MASK;
        len = cigar[i] >> BAM_CIGAR_SHIFT;

        switch( op ) {
            case BAM_CMATCH:
                for( j = 0; j < len; j++ ) A[pos++]++;
                break;
            case BAM_CDEL:
            case BAM_CREF_SKIP:
                pos += len;
                break;

            case BAM_CSOFT_CLIP:
            case BAM_CHARD_CLIP:
            /* note: insertions are not plotted, because there is no reasonable
             * way to do so. */
            case BAM_CINS:
            default:
                break;
        }
    }
}


int main( int argc, char* argv[] )
{
    char strand = '*';


    /* 1. Parse Args */

    static struct option long_options[] = {
        {"strand", 1, 0, 0},
        {0,0,0,0} };

    int c, option_index;
    do {
        c = getopt_long( argc, argv, "s:", long_options, &option_index );

        if( c == 's' || (c == 0 && option_index == 0) ) {
            strand = optarg[0];
            if( strand != '-' && strand != '+' && strand != '*' ) {
                fprintf( stderr, "Error: Invalid strand: %c\n", strand );
                exit(1);
            }
        }
        else if( c == '?' ) {
            print_usage();
            exit(1);
        }

    } while( c != -1 );

    if( optind >= argc ) {
        print_usage();
        exit(1);
    }

    char* bam_fn   = argv[optind];


    samfile_t* bam_f = samopen( bam_fn, "rb", NULL );

    if( bam_f == NULL ) {
        fprintf( stderr, "Error: Can't open file '%s'.\n", bam_fn );
        exit(1);
    }

    bam1_t* b = bam_init1();

    int32_t curr_tid = -1;
    uint32_t* A = NULL;
    size_t    n = 0;

    while( samread( bam_f, b ) > 0 ) {
        if( strand == '+' && bam1_strand(b) == 1 ) continue;
        if( strand == '-' && bam1_strand(b) == 0 ) continue;

        if( b->core.tid != curr_tid ) {
            if( curr_tid != -1 ) {
                print_rle( A, n, bam_f->header->target_name[curr_tid] );
                free(A);
            }

            curr_tid = b->core.tid;

            n = bam_f->header->target_len[curr_tid];
            A = calloc( n, sizeof(uint32_t) );
            memset( A, 0, n*sizeof(uint32_t) );

        }

        count_read( b, A, n );
    }

    if( curr_tid != -1 ) {
        print_rle( A, n, bam_f->header->target_name[curr_tid] );
        free(A);
    }

    bam_destroy1(b);

    samclose( bam_f );

    return 0;
}



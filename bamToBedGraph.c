
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


/* 
 *  I. Reading/storing tab seperated chromosome length file.
 *  (Like those used by UCSC and BEDTools.)
 */

typedef struct _chrom_size
{
    struct _chrom_size* next;
    char*               chrom;
    uint32_t            size;
} chrom_size;


void chrom_sizes_push( chrom_size** S, const char* chrom, uint32_t size )
{
    chrom_size* T = malloc(sizeof(chrom_size));
    T->chrom  = strdup(chrom);
    T->size   = size;
    T->next   = *S;

    *S = T;
}

void chrom_sizes_clear( chrom_size** S )
{
    chrom_size* T;

    while( *S ) {
        T = (*S)->next;
        free((*S)->chrom);
        free(*S);
        *S = T;
    }
}

uint32_t chrom_sizes_get( chrom_size* S, const char* chrom )
{
    /* simply linear search: the number of chroms is typicall small */
    while( S && strcmp( S->chrom, chrom ) ) S = S->next;

    if( S ) return S->size;
    else    return 0;
}


chrom_size* read_chrom_sizes( const char* fn )
{
    fprintf( stderr, "Reading chrom sizes..." );

    FILE* f = fopen( fn, "r" );
    if( f == NULL ) {
        fprintf( stderr, "Can't open file '%s'\n", fn );
        exit(1);
    }

    size_t buf_size = 16384;
    char* buf = malloc( sizeof(char) * buf_size );

    size_t n = 0;

    char*    chrom = NULL;
    uint32_t size;

    chrom_size* cs = NULL;

    while( fgets( buf, buf_size, f ) ) {
        if( sscanf( buf, "%s\t%d", chrom, &size ) != 2 ) {
            if( chrom ) {
                free(chrom);
                chrom = NULL;
            }
            continue;
        }

        chrom_sizes_push( &cs, chrom, size );
        chrom = NULL;
        n++;
    }

    free(buf);
    fclose(f);

    fprintf( stderr, "Done. (%d read.)\n", n );
    return cs;
}


/* II. Run length encoding of an array, to covert the coverage array to bedGraph. */
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
            "Usage: bamToBedGraph [options] in.bam chrom.sizes\n"
            "Covert a BAM file to a Bed Graph file\n\n."
            "options:\n"
            "  -s, --strand=STRAND    Either '+', '-', or '*', where '*' is both strands.\n" 
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

    if( argc - optind < 2 ) {
        fprintf( stderr, "Error: Too few arguments.\n" );
        print_usage();
        exit(1);
    }

    char* bam_fn   = argv[optind];
    char* chrom_fn = argv[optind+1];


    chrom_size* cs = read_chrom_sizes( chrom_fn ); 


    samfile_t* bam_f = samopen( bam_fn, "rb", NULL );

    if( bam_f == NULL ) {
        fprintf( "Error: Can't open file '%s'.\n", bam_fn );
        exit(1);
    }

    bam1_t* b = bam_init1();

    int32_t curr_tid = -1;
    uint32_t* A = NULL;
    size_t    n = 0;

    /* TODO: HOLY FUCK: I don't need to read in a chrom.sizes file. */

    while( samread( bam_f, b ) > 0 ) {
        if( b->core.tid != curr_tid ) {
            if( curr_tid != -1 ) {
                print_rle( A, n, bam_f->header->targen_name[curr_tid] );
                free(A);
            }

            n = bam_f->header->target_len[tid];
            A = calloc( n, sizeof(uint32_t) );
            memset( A, 0, n*sizeof(uint32_t) );

            curr_tid = tid;
        }

        count_read( b, A, n );
    }

    if( curr_tid != -1 ) {
        print_rle( A, n, bam_f->header->targen_name[curr_tid] );
        free(A);
    }

    bam_destroy1(b);

    samclose( bam_f );
    chrom_sizes_clear(&cs);

    return 0;
}





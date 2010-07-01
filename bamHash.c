
#include "hash.h"
#include <samtools/sam.h>


void usage()
{
    fputs( "Usage: bamHash <reads.bam> [k]\n", stderr );
}


void print_hashed_value( bam_header_t* header, struct table* T, struct hashed_value* v, int strand )
{ 
    size_t count = 0;
    if( strand == 0 )       count = v->pos_count;
    else if( strand == 1 )  count = v->neg_count;
    else if( strand == -1 ) count = v->pos_count + v->neg_count;

    if( count == 0 ) return;

    char* seq = malloc((sizeof(char)+1) * T->read_len);
    int i;
    for( i = 0; i < T->read_len; i++ ) {
        seq[i] = bam_nt16_rev_table[bam1_seqi( v->seq, i )];
    }
    seq[i] = '\0';


    if( strand >= 0 ) {
        fprintf( stdout, "%d\t%s\t%c\n",
                count,
                seq,
                strand == 0 ? '+' : '-' );
    } 
    else {
        fprintf( stdout, "%d\t%s\n", count, seq );
    }

    
    free(seq);
}
void print_top_reads( bam_header_t* header, struct table* T,
                      struct hashed_value** S )
{
    size_t i;
    for( i = 0; i < T->m; i++ ) {
        print_hashed_value( header, T, S[i], -1 );
    }
}

void print_top_reads_stranded( bam_header_t* header, struct table* T,
                               struct hashed_value** S_pos, struct hashed_value** S_neg )
{
    size_t n = T->m;

    size_t i=0,j=0;
    while( i<n || j<n ) {
        if( j>=n ) {
            print_hashed_value( header, T, S_pos[i++], 0 );
        }
        else if( i>=n ) {
            print_hashed_value( header, T, S_neg[j++], 1 );
        }
        else {
            if( S_pos[i]->pos_count > S_neg[j]->neg_count )
                print_hashed_value( header, T, S_pos[i++], 0 );
            else
                print_hashed_value( header, T, S_neg[j++], 1 );
        }
    }
}




int main( int argc, char* argv[] )
{
    if( argc < 2 ) {
        usage();
        exit(1);
    }

    size_t n = 0;
    struct table T;

    samfile_t* fp = NULL;
    bam1_t* read  = bam_init1();

    fp = samopen( argv[1], "rb", NULL );
    if( fp == NULL ) {
        fprintf( stderr, "Can't open file '%s'.\n", argv[1] );
        exit(1);
    }

    if( argc > 2 ) {
        table_create( &T, atoi(argv[2]) );
    }
    else {
        /* assume all reads are of the same length, so we can read one to determine
         * the length */
        samread(fp,read);
        table_create( &T, read->core.l_qseq );
        table_inc( &T, read );
        n = 1;
    }




    fprintf( stderr, "Reading...\n" );
    while( samread(fp,read) > 0 ) {
        table_inc( &T, read );
        n++;

        if( n % 1000000 == 0 ) {
            fprintf( stderr, "\t%d reads (%d unique)...\n", n, T.m );
        }
    }
    fprintf( stderr, "%d reads (%d unique). Finished.\n", n, T.m );


    struct hashed_value** S;

    fprintf( stderr, "Sorting...\n" );
    sort_by_count( &T, &S );

    print_top_reads( fp->header, &T, S );

    free(S);

    table_destroy(&T);

    bam_destroy1(read);
    samclose(fp);

    return 0;
}






#include "hash.h"
#include <samtools/sam.h>


void usage()
{
    fprintf( stderr,
            "Usage: bamHash <reads.bam> [k] [d]\n"
            "Hashes the first k nucleotides of each read and prints counts\n"
            "in the form:\n"
            "   seq1   count1\n"
            "   seq2   count2\n"
            "If d is specified, produce counts for the sequence and all sequence\n"
            "within Hamming distance d\n" );
}

const char* NT = "ATCGN";

char* bam_seq_to_seq( bam1_t* read )
{
    char* seq = malloc(sizeof(char)*(1+read->core.l_qseq));
    int i;
    for( i = 0; i < read->core.l_qseq; i++ ) {
        seq[i] = bam_nt16_rev_table[bam1_seqi( bam1_seq(read), i )];
    }
    seq[i] = '\0';

    return seq;
}


void bam_seq_to_seq_n( bam1_t* read, char* out, size_t n )
{
    if( read->core.l_qseq < n ) n = read->core.l_qseq;

    int i;
    for( i = 0; i < n; i++ ) {
        out[i] = bam_nt16_rev_table[bam1_seqi( bam1_seq(read), i )];
    }

    out[i] = '\0';
}



void print_reads( bam_header_t* header, struct table* T,
                      struct hashed_value** S )
{
    size_t i;
    for( i = 0; i < T->m; i++ ) {
        fprintf( stdout, "%s\t%llu\n", S[i]->seq, S[i]->count );
    }
}



uint64_t subst_count( struct table* T, char* seq, size_t d )
{
    if( d == 0 ) return 0;

    char c;
    uint64_t n = 0;
    size_t i, j;
    for( i = 0; i < T->read_len; i++ ) {
        c = seq[i];

        for( j = 0; j < 5; j++ ) {
            if( c == NT[j] ) continue;
            seq[i] = NT[j];
            n += table_get( T, seq );
            n += subst_count( T, seq, d-1 );
        }

        seq[i] = c;
    }

    return n;
}


void add_cluster_counts( struct table* T, size_t d )
{
    size_t i;
    struct hashed_value *j;
    for( i = 0; i < T->n; i++ ) {
        j = T->A[i];
        while( j ) {
            j->count += subst_count(T,j->seq,d);
            j = j->next;
        }
        if( i > 0 && i % 100000 == 0 ) {
            fprintf( stderr, "\t%d reads...\n", i );
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
    size_t d = 0;
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
        char* tmp = bam_seq_to_seq(read);
        table_inc( &T, tmp );
        free(tmp);
        n = 1;
    }

    if( argc > 3 ) d = atoi(argv[3]);


    char* seq = malloc(sizeof(char)*(1+T.read_len));

    fprintf( stderr, "Reading...\n" );
    while( samread(fp,read) > 0 ) {
        bam_seq_to_seq_n( read, seq, T.read_len );
        table_inc( &T, seq );
        n++;

        if( n % 1000000 == 0 ) {
            fprintf( stderr, "\t%d reads (%d unique)...\n", n, T.m );
        }
    }
    fprintf( stderr, "%d reads (%d unique). Finished.\n", n, T.m );


    struct hashed_value** S;

    fprintf( stderr, "Clustering...\n" );

    if( d > 0 ) {
        add_cluster_counts( &T, d );
    }

    fprintf( stderr, "Sorting...\n" );

    sort_by_count( &T, &S );

    print_reads( fp->header, &T, S );

    free(S);

    table_destroy(&T);

    bam_destroy1(read);
    samclose(fp);

    return 0;
}





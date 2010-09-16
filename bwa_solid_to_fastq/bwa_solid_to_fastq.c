/*
 *      bwa_solid_to_fastq
 *      ------------------
 *      A C rewrite of the solid2fastq.pl script that comes with BWA.
 *      (Because I want this to finish in my lifetime.)
 *
 *      September 2010 / Daniel Jones <dcjones@cs.washington.edu>
 *
 */




#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>



void usage()
{
    fprintf( stderr, "Usage: bwa_solid_to_fastq in.csfasta in.qual\n\n"
                     "Writes to stdout a fastq file suitable for use with BWA.\n\n" );
}

FILE* safe_fopen( const char* fn, const char* mode ) {
    FILE* f = fopen( fn, mode );
    if( f == NULL ) {
        fprintf( stderr, "Error: Can't open file '%s'.\n", fn );
        exit(1);
    }

    return f;
}

/* convert the quals string into an array of doubles, return the number of
 * qualities extracted. */
int get_quals( char* in, int* out )
{
    int i = 0, j, l = 0;
    char c;
    while( in[i] ) {
        while( isspace(in[i]) ) i++;
        if( !isdigit(in[i]) ) break;
        j = i;
        while( isdigit(in[j]) ) j++;

        c = in[j];
        in[j] = '\0';
        out[l++] = atoi( in + i );
        in[j] = c;

        i = j;
    }

    return l;
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


int main( int argc, char* argv[] )
{
    if( argc < 3 ) {
        usage();
        exit(1);
    }

    FILE* reads_in = safe_fopen( argv[1], "r" );
    FILE* quals_in = safe_fopen( argv[2], "r" );

    const size_t buf_size = 4096;

    char* read_name = malloc( buf_size*sizeof(char) );
    char* read_seq  = malloc( buf_size*sizeof(char) );
    char* qual_name = malloc( buf_size*sizeof(char) );
    char* qual_seq  = malloc( buf_size*sizeof(char) );

    int* quals = malloc( buf_size*sizeof(int) );

    int n;

    /* skip initial comments */
    while( fgets_noncomment( read_name, buf_size, reads_in ) &&
           fgets_noncomment( read_seq,  buf_size, reads_in ) &&
           fgets_noncomment( qual_name, buf_size, quals_in ) &&
           fgets_noncomment( qual_seq,  buf_size, quals_in ) )
    {
        if( strcmp( read_name, qual_name ) != 0 ) {
            fprintf( stderr, "Error: Mismatching read/quality pair.\n" );
            exit(1);
        }

        n = get_quals( qual_seq, quals );

        read_name[strlen(read_name)-1] = '\0';
        double_encode( n, read_seq, qual_seq, quals );

        printf( "@%s/1\n%s\n+\n%s\n", read_name+1, read_seq, qual_seq );
    }


    fclose( reads_in );
    fclose( quals_in );

    free(read_name);
    free(read_seq);
    free(qual_name);
    free(qual_seq);


    return 0;
}





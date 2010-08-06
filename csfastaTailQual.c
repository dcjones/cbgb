/*
 *
 * Sample the mean base qualities of poly-0 tails found in reads.
 *
 */

#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>

void usage()
{

}


/* return the offset at which the poly-0 tail begins, or -1 if it has none. */
int poly_tail( const char* seq )
{
    int n = strlen(seq)-1; /* -1 to skip newline */
    int i = n-1;
    while( i >= 0 && (seq[i] == '0' || seq[i] == '.') ) i--;

    if( i == n-1 ) return -1;
    else           return  i-1;  /* -1 to ignore than last non-0 color as well */
}

/* convert the quals string into an array of doubles, return the number of
 * qualities extracted. */
int get_quals( char* in, double* out )
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
        out[l++] = atof( in + i );
        in[j] = c;

        i = j;
    }

    return l;
}



double mean_qual( double* quals, size_t i, size_t j )
{
    double n = 0.0;
    double mu = 0.0;
    while( i <= j ) {
        mu += quals[i];
        n++;
        i++;
    }

    return mu/n;
}


int main( int argc, char* argv[] )
{
    if( argc < 3 ) {
        usage();
        exit(1);
    }

    FILE* reads_in = fopen( argv[1], "r" );
    FILE* quals_in = fopen( argv[2], "r" );

    const size_t buf_size = 4096;

    char *read_name, *read_seq, *qual_name, *qual_seq;

    read_name = malloc( buf_size*sizeof(char) );
    read_seq  = malloc( buf_size*sizeof(char) );
    qual_name = malloc( buf_size*sizeof(char) );
    qual_seq  = malloc( buf_size*sizeof(char) );
    double* quals = malloc( buf_size*sizeof(double) );

    int k;
    int n;

    while( fgets( read_name, buf_size, reads_in ) &&
           fgets( read_seq,  buf_size, reads_in ) &&
           fgets( qual_name, buf_size, quals_in ) &&
           fgets( qual_seq,  buf_size, quals_in ) )
    {
        if( strcmp( read_name, qual_name ) != 0 ) {
            fprintf( stderr, "Mismatching read, quality pair.\n" );
            exit(1);
        }

        k = poly_tail( read_seq );

        if( k < 5 ) continue;

        n = get_quals( qual_seq, quals );
        double mu = mean_qual( quals, k, n-1 );

        printf( "%0.5f\n", mu );
    }


    fclose(reads_in);
    fclose(quals_in);

    return 0;
}


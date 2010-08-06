/*
 *                     csfastaUntail
 *                     -------------
 *                     For all reads with trailing 0's, chop this tail off and
 *                     write to a file.
 *
 *                     July 2010  /  Daniel Jones <dcjones@cs.washington.edu>
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>


const int    MIN_LEN  = 6;
const double MIN_QUAL = 12.0;


void usage()
{
    fprintf( stderr, "Usage: csfastaUntail in.csfasta in.qual out.csfasta out.qual\n" );
}


FILE* safe_fopen( const char* fn, const char* mode ) {
    FILE* f = fopen( fn, mode );
    if( f == NULL ) {
        fprintf( stderr, "Can't open file '%s'.\n", fn );
        exit(1);
    }

    return f;
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


/* return the offset at which the poly-0 tail begins, or -1 if it has none. */
int poly_tail( const char* seq )
{
    int n = strlen(seq)-1; /* -1 to skip newline */
    int i = n-1;
    while( i >= 0 && (seq[i] == '0' || seq[i] == '.') ) i--;

    if( i == n-1 ) return -1;
    else           return  i-1;  /* -1 to ignore than last non-0 color as well */
}


/* print k space seperated quality values from the sequence 'seq'. */
void printn_quals( FILE* out, const char* name, const double* quals, int k )
{
    fputs( name, out );

    int i;

    if( k ) {

        fprintf( out, "%d", (int)(quals[0]) );
        i = 1;
    
        while( i < k ) {
            fprintf( out, " %d", (int)(quals[i]) );
            i++;
        }
    }

    fputc( '\n', out );
}


int main( int argc, char* argv[] )
{
    if( argc < 5 ) {
        usage();
        exit(1);
    }

    const size_t buf_size = 4096;

    FILE *reads_in, *quals_in, *reads_out, *quals_out;
    char *read_name, *read_seq, *qual_name, *qual_seq;

    reads_in = safe_fopen( argv[1], "r" );
    quals_in = safe_fopen( argv[2], "r" );
    reads_out = safe_fopen( argv[3], "w" );
    quals_out = safe_fopen( argv[4], "w" );


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
        if( k < 0 ) continue;

        n = get_quals( qual_seq, quals );
        if( mean_qual( quals, k, n-1 ) < MIN_QUAL ) continue;

        if( n - k < MIN_LEN ) continue;

        read_seq[k+1] = '\0';
        fprintf( reads_out, "%s%s\n", read_name, read_seq );

        printn_quals( quals_out, qual_name, quals, k );
    }


    free(read_name);
    free(read_seq);
    free(qual_name);
    free(qual_seq);

    free(quals);

    fclose(reads_in);
    fclose(quals_in);
    fclose(reads_out);
    fclose(quals_out);

    return 0;
}


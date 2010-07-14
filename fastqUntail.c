/*
 *                     fastaqUntail
 *                     -------------
 *                     For all reads with trailing or leading A's, chop this tail off and
 *                     write to a file.
 *
 *                     July 2010  /  Daniel Jones <dcjones@cs.washington.edu>
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <unistd.h>

void usage()
{
    fprintf( stderr,
             "Usage: fastaqUntail [OPTIONS] in.fastq\n"
             "Trim particular nucleotides of either end of a read.\n"
             "Options:\n"
             "-5            trim 5' end\n"
             "-3            trim 3' end\n"
             "-n nt         trim the given nucleotides (default 'ATN')\n"
             );
}


const char* all_nt = "ATN";


/* is a character in a string */
bool is_in( char c, char* S )
{
    while( *S ) {
        if( c == *S ) return true;
        S++;
    }

    return false;
}


int main( int argc, char* argv[] )
{
    const char* optstring = "53n:";
    char* nt = NULL;
    int c;

    bool trim_5 = false, trim_3 = false;

    do {
        c = getopt( argc, argv, optstring );
        switch(c) {
            case '5':
                trim_5 = true;
                break;
            case '3':
                trim_3 = true;
                break;
            case 'n':
                nt = strdup(optarg);
                break;
        }

    } while( c != -1 );

    if( nt == NULL ) nt = strdup(all_nt);

    if( optind >= argc ) {
        usage();
        exit(1);
    }

    FILE* f = fopen( argv[optind], "r" );
    if( f == NULL ) {
        fprintf( stderr, "Can't open file '%s'.\n", argv[optind] );
        exit(1);
    }

    const size_t max_line_len = 4096;
    char* id1  = malloc( max_line_len*sizeof(char) );
    char* seq  = malloc( max_line_len*sizeof(char) );
    char* id2  = malloc( max_line_len*sizeof(char) );
    char* qual = malloc( max_line_len*sizeof(char) );

    int i,j,n;

    while( fgets( id1, max_line_len, f ) &&
           fgets( seq, max_line_len, f ) &&
           fgets( id2, max_line_len, f ) &&
           fgets( qual, max_line_len, f ) )
    {
        if( !( id2[1] == '\n' || strcmp( id1+1, id2+1 ) == 0 ) ) {
            fprintf( stderr, "Mismatched read/quality ids:\nRead ID: %s\nQual ID%s\n",
                    id1, id2 );
            exit(1);
        }
        
        n = strlen(seq)-1;

        i = 0;   /* 5' trim */
        j = n-1; /* 3' trim */

        if( trim_5 ) {
            while( i < j && is_in( seq[i], nt ) ) i++;
        }

        if( trim_3 ) {
            while( j > i && is_in( seq[j], nt ) ) j--;
        }

        /* supress untrimmed reads */
        if( i == 0 && j == n-1 ) continue;

        seq[j+1] = '\0';
        qual[j+1] = '\0';
        printf( "%s%s\n%s%s\n", id1, seq+i, id2, qual+i );
    }


    fclose(f);
    free(id1);
    free(seq);
    free(id2);
    free(qual);
    free(nt);
    return 0;
}







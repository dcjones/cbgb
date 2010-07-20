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
#include <regex.h>

void usage()
{
    fprintf( stderr,
             "Usage: fastaqUntail [OPTIONS] in.fastq\n"
             "Trim poly-A/T/C/G stretches from either end of a read.\n"
             "Options:\n"
             "-m            minimum prefix/suffix length (default: 4)\n"
             "-5            trim 5' end\n"
             "-3            trim 3' end\n"
             "-n nt         trim the given nucleotides (default: 'AT')\n"
             );
}


const char* all_nt = "AT";


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
    const char* optstring = "53n:m:";
    char* nt = NULL;
    int c;
    int m = 4;

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
            case 'm':
                m = atoi(optarg);
                break;
        }

    } while( c != -1 );

    if( nt == NULL ) nt = strdup(all_nt);

    if( optind >= argc ) {
        usage();
        exit(1);
    }


    /* build prefix and suffix regular expressions */
    char pat[512];

    int offset=1;

    offset += sprintf( pat+1, "([%cN]{%d,}", nt[0], m );

    char* s = nt+1;
    while( *s ) {
        offset += sprintf( pat+offset, "|[%cN]{%d,}", *s, m );
        s++;
    }


    pat[offset++] = ')';
    pat[offset]   = '$';
    pat[offset+1] = '\0';

    /*fprintf( stderr, "Compiling regex: '%s'\n", pat+1 );*/

    regex_t suf_re;
    if( regcomp( &suf_re, pat+1, REG_ICASE | REG_EXTENDED ) ) {
        fprintf( stderr, "Regular expression failed to compile. Sorry, please report this.\n" );
    }

    pat[0] = '^';
    pat[offset] = '\0';

    /*fprintf( stderr, "Compiling regex: '%s'\n", pat );*/

    regex_t pre_re;
    if( regcomp( &pre_re, pat, REG_ICASE | REG_EXTENDED ) ) {
        fprintf( stderr, "Regular expression failed to compile. Sorry, please report this.\n" );
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

    regmatch_t mat;

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
        seq[n] = '\0';

        i = 0;   /* 5' trim */
        j = n-1; /* 3' trim */

        if( trim_5 && regexec( &pre_re, seq, 1, &mat, 0 ) == 0 ) {
            i = mat.rm_eo;
        }

        if( trim_3 && regexec( &suf_re, seq+i, 1, &mat, 0 ) == 0 ) {
            j = mat.rm_so-1;
        }

        
        /* supress untrimmed reads */
        if( i == 0 && j == n-1 ) continue;

        seq[j+1] = '\0';
        qual[j+1] = '\0';
        printf( "%s%s\n%s%s\n", id1, seq+i, id2, qual+i );
    }


    regfree(&suf_re);
    regfree(&pre_re);
    fclose(f);
    free(id1);
    free(seq);
    free(id2);
    free(qual);
    free(nt);
    return 0;
}







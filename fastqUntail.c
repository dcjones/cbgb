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
#include <unistd.h>

void usage()
{
    fprintf( stderr,
             "Usage: fastaqUntail [OPTIONS] in.fastq\n\n"
             "Options:\n"
             "-5            trim 5' end\n"
             "-3            trim 3' end\n"
             "-n nt         trim the given nucleotides (default 'AT')\n"
             );
}


const char* all_nt = "ATCGN";

int main( int argc, char* argv[] )
{
    const char* optstring = "53n:";
    char* nt = all_nt;
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


    if( nt != all_nt ) free(nt);
    return 0;
}







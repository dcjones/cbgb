
/* Replace ensembl chromosome names with UCSC chromosome names, very quickly. */

#include <pcre.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>


int main( int arg, char* argv[] )
{
    pcre* seq_pat;
    const char* error;
    int erroffset;
    seq_pat = pcre_compile( "^SEQ \\S+ (\\S+)", 0, &error, &erroffset, NULL );
    if( seq_pat == NULL ) {
        fprintf( stderr, "Can't compile RE pattern.\n" );
        exit(1);
    }

    pcre_extra* seq_pat_extra = pcre_study( seq_pat, 0, &error );

    int mat;
    int ovec[6];
    size_t maxline = 4096;
    char* line = malloc( maxline*sizeof(char) ); 

    /* for converting integers */
    char c;
    char *endptr;
    long chrom;

    while( fgets( line, maxline, stdin ) ) {
        if( line[0] != 'S' ) {
            fputs( line, stdout );
            continue;
        }

        mat = pcre_exec( seq_pat,        /* code */
                         seq_pat_extra,  /* extra */ 
                         line,           /* subject */
                         strlen(line),   /* length */
                         0,              /* startoffset */
                         0,              /* options */
                         ovec,           /* ovector */
                         6 );            /* ovector size */ 


        if( mat == 2 ) {
            fprintf( stderr, "MATCH: %d-%d\n", ovec[2], ovec[3] );

            if( strncmp( "MT", line+ovec[2], ovec[3]-ovec[2] ) == 0 ) {
                fwrite( line, sizeof(char), ovec[2], stdout );
                fprintf( stdout, "chrM%s", line+ovec[3] );
            }
            else {
                c = line[ovec[3]];
                line[ovec[3]] = '\0';
                chrom = strtol( line+ovec[2], &endptr, 10 );
                if( *endptr == '\0' ) {
                    line[ovec[3]] = c;
                    fwrite( line, sizeof(char), ovec[2], stdout );
                    fprintf( stdout, "chr%ld%s", chrom, line+ovec[3] );
                }
                else {
                    line[ovec[3]] = c;
                    fputs( line, stdout );
                }
            }

        }
        else if( mat == PCRE_ERROR_NOMATCH ) fputs( line, stdout );
        else {
            fprintf( stderr, "Error in pcre_exec: %d\n", mat );
            exit(1);
        }
    }

    free(line);
    pcre_free( seq_pat );
    pcre_free( seq_pat_extra );
    return 0;
}


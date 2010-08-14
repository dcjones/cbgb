/*
 * Compute base pair probabilities for RNA secondary structure, according to
 * ViennaRNA.
 *
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <stdbool.h>
#include <ctype.h>

#include <ViennaRNA/fold.h>
#include <ViennaRNA/part_func.h>
#include <ViennaRNA/fold_vars.h>


bool is_nt( char c )
{
    switch(c) {
        case 'A': case 'T': case 'C': case 'G': case 'U': return true;
        default: return false;
    }
}


char* read_fa( FILE* f )
{
    const int bufsize = 2048;
    char* buf = malloc(bufsize*sizeof(char));
    char* seq;

    int seqlen = 0;
    char* c;

    /* get sequence length */
    while( fgets(buf,bufsize,f) ) {
        if( buf[0] == '>' ) continue;
        for( c = buf; *c; c++ ) {
            *c = toupper(*c);
            if( is_nt(*c) ) seqlen++;
        }
    }

    rewind(f);

    /* get sequence */
    seq = malloc((seqlen+1)*sizeof(char));
    int i = 0;

    while( fgets(buf,bufsize,f) ) {
        if( buf[0] == '>' ) continue;
        for( c = buf; *c; c++ ) {
            *c = toupper(*c);
            if( is_nt(*c) ) seq[i++] = *c;
        }
    }
    seq[i] = '\0';

    free(buf);

    return seq;
}


void print_pr( int len )
{
    int i,j;
    for( i = 0; i < len; i++ ) {
        for( j = 0; j < len; j++ ) {
            printf( "%e", pr[iindx[i]-j] );
            j == len-1 ? fputc('\n',stdout) : fputc('\t',stdout);
        }
    }
}

double* marginals( int len )
{
    double* M = malloc(sizeof(double)*len);

    int i,j;
    for( i = 1; i <= len; i++ ) {
        M[i-1] = 0.0;
        for( j = i+1; j <= len; j++ ) {
            M[i-1] += pr[iindx[i]-j];
        }
        for( j = 1; j < i; j++ ) {
            M[i-1] += pr[iindx[j]-i];
        }
    }

    return M;
}


int main( int argc, char* argv[] )
{
    if( argc < 2 ) {
        fprintf( stderr, "Usage: RNApart in.fa\n" );
        exit(1);
    }

    FILE* f = fopen(argv[1],"r");
    char* seq = read_fa(f);
    fclose(f);

    int n = strlen(seq);

    char* struc = malloc(sizeof(char)*(n+1));
    fprintf( stderr, "Folding...\n" );
    double min_en = fold(seq,struc);
    free(struc);
    free_arrays();

    do_backtrack = 1;

    /* set pf_scale the same way RNAfold does */
    double sfact=1.07;
    double kT = (temperature+273.15)*1.98717/1000.; /* in Kcal */
    pf_scale = exp(-(sfact*min_en)/kT/(double)n);
    if (n>2000) fprintf(stderr, "scaling factor %f\n", pf_scale);


    init_pf_fold(n);

    fprintf( stderr, "Computing pair probabilities...\n" );
    pf_fold(seq,NULL);

    fprintf( stderr, "Computing marginals...\n" );
    double* M = marginals(n);

    int i;
    for( i = 0; i < n; i++ ) {
        printf( "%e\n", M[i] );
    }

    free(M);
    free(seq);
    free_pf_arrays();

    return 0;
}


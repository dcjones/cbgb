/*
 *                      resample
 *                      --------
 *                      Efficiently print random lines (or clumps of lines),
 *                      sampling uniformly with replacement. Useful for
 *                      bootstraping, etc.
 *
 *                      Daniel Jones
 *                      dcjones@cs.washington.edu
 *                      April/2010
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>


/* count the number of lines in the file */
size_t count_lines( FILE* f )
{
    long pos = ftell(f);
    rewind(f);

    char *p, *buf = malloc(BUFSIZ);

    size_t cnt;
    size_t n = 0;
    while( (cnt = fread( buf, sizeof(char), BUFSIZ, f )) ) {
        p = buf;
        while( (p = strchr( p, '\n' )) ) { p++; n++; }
    }

    return n;

    free(buf);
    fseek(f,pos,SEEK_SET);
}


/* jobustly generate a random integers in {1,...,n-1} */
size_t uniform(size_t n)
{
    size_t r; 
    do{
        r = random();
    } while( r >= n*(RAND_MAX/n) );

    return r % n;
}


/* compare two unsigned integers */
int cmpul( const void* a, const void* b ) {
    if(      *(size_t*)a < *(size_t*)b ) return -1;
    else if( *(size_t*)a > *(size_t*)b ) return 1;
    else                                 return 0;
}



/* read chunks of k lines from a file */
char* fgetsk( char* s, int size, FILE* stream, size_t k )
{
    size_t cnt;
    s[0] = '\0';
    while( k-- ) {
        cnt = strlen(s);
        size -= cnt;
        s    += cnt;
        if( !fgets( s, size, stream ) ) return NULL;
    }

    return s;
}


int main( int argc, char* argv[] )
{
    if( argc < 3 ) {
        fprintf( stderr,
                "Usage: resample file.in m [k]\n\n"
                "Print m random lines from the specified file, with replacement,\n"
                "and in the order in which they appear in the file.\n\n"
                "If 'k' is specified, consider groups of k lines.\n\n"
                "Works in O(m) space and O(m+n) time, where n is the size of the file.\n"
                );
        exit(EXIT_FAILURE);
    }


    srandom((unsigned int)time(NULL));


    /* 1. Count lines */
    FILE*  f = fopen(argv[1],"r");

    if( !f ){
        fprintf( stderr, "Can't open file '%s'\n", argv[1] );
        exit(EXIT_FAILURE);
    }


    size_t k = 1;
    if( argc > 3 ) k = strtoul(argv[3],NULL,10);
    size_t n = count_lines(f);
    size_t m = strtoul(argv[2],NULL,10);
    rewind(f);

    if( n % k ) {
        fprintf( stderr,
                 "Warning: number of lines n=%d is not divisible by k=%d\n",
                 n, k );
    }
                 


    /* 2. Generate a sorted stream of random numbers. */
    size_t* lines = malloc(sizeof(size_t)*m);

    size_t i,j;
    for( i = 0; i < m; i++ ) lines[i] = uniform(n/k);

    qsort( lines, m, sizeof(size_t), cmpul );



    /* 3. Scan through the input file outputing line as enumerated. */
    const size_t buf_size = 4096*k;
    char* buf = malloc(buf_size);

    i = 0;
    j = 0;

    while( i < m && fgetsk(buf,buf_size,f,k) ) {
        while( lines[i] == j ) {
            fwrite(buf,sizeof(char),strlen(buf),stdout);
            i++;
        }
        j++;
    }



    free(lines);
    free(buf);

    return 0;
}




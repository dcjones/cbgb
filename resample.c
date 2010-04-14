

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


/* Robustly generate a random integers in {1,...,n-1} */
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


int main( int argc, char* argv[] )
{
    srandom((unsigned int)time(NULL));


    /* 1. Count lines */
    FILE*  f = fopen(argv[1],"r");
    size_t n = count_lines(f);
    size_t m = strtoul(argv[2],NULL,10);
    rewind(f);



    /* 2. Generate a sorted stream of random numbers. */
    size_t* lines = malloc(sizeof(size_t)*m);

    size_t i,j;
    for( i = 0; i < m; i++ ) lines[i] = uniform(n);

    qsort( lines, m, sizeof(size_t), cmpul );



    /* 3. Scan through the input file outputing line as enumerated. */
    char* buf = malloc(BUFSIZ);

    i = 0;
    j = 0;

    while( i < m && fgets(buf,BUFSIZ,f) ) {
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




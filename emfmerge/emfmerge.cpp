


#include "emfparse.hpp"

#include <cstdlib>
#include <cstdio>
#include <cstring>


void usage()
{
    fprintf( stderr, "mafmerge species.txt ref.fasta < input.emf\n" );
}


vector<string> parse_species( const char* fn )
{
    vector<string> species;
    FILE* f = fopen( fn, "r" );
    if( f == NULL ) {
        fprintf( stderr, "Can't open file %s.\n", fn );
        exit(1);
    }

    size_t max_line = 2048;
    char* buf = new char[max_line];

    while( fgets( buf, max_line, f ) != NULL ) {
        if( buf[0] != '\0' && buf[0] != '\n' ) {
            buf[strlen(buf)-1] = '\0';
            species.push_back( string(buf) );
        }
    }

    fclose(f);
    delete[] buf;

    return species;
}


FILE* fopen_or_die( const char* fn, const char* mode )
{
    FILE* f = fopen( fn, mode );
    if( f == NULL ) {
        fprintf( stderr, "Can't open file %s\n", fn );
        exit(1);
    }

    return f;
}



int main( int argc, char* argv[] )
{
    if( argc < 3 ) {
        usage();
        exit(1);
    }

    vector<string> species = parse_species( argv[1] );

    emfparse* aln;

    if( argc >= 3 ) {
        FILE* f = fopen_or_die( argv[3], "r" );
        aln = new emfparse( species, argv[2], f );
        fclose(f);
    }
    else {
        aln = new emfparse( species, argv[2], stdin );
    }

    aln->print_fasta( stdout );


    free(aln);
    return 0;
}







#include "emfparse.hpp"

#include <pcrecpp.h>

#include <string>

#include <cstdlib>
#include <cstdio>
#include <cstring>

using namespace std;


bool isnt( char c ) {
    switch( c ) {
        case 'A':
        case 'T':
        case 'C':
        case 'G':
        case 'N':
        case 'a':
        case 't':
        case 'c':
        case 'g':
        case 'n':
            return true;
        default:
            return false;
    }
}

char complement( char c ) {
    switch( c ) {
        case 'A': return 'T';
        case 'C': return 'G';
        case 'T': return 'A';
        case 'G': return 'C';
        case 'a': return 't';
        case 'c': return 'g';
        case 't': return 'a';
        case 'g': return 'c';
        default: return c;
    }
}

void* safe_realloc( void* ptr, size_t size )
{
    void* p;
    if( (p = realloc( ptr, size )) == NULL ) {
        fprintf( stderr, "Out of memory!\n" );
        exit(1);
    }

    return p;
}



emfparse::emfparse( vector<string> species, const char* ref_fn, FILE* emf_f )
{
    parse_ref( ref_fn, &ref_seqname, &ref_seq );

    // assume the first species in the list is the reference
    msa[string(species[0])] = ref_seq;
    ref_name = species[0];

    n = strlen(ref_seq);

    // allocate space for the rest of the species
    fprintf( stderr, "allocating memory ... " );
    char* seq;
    vector<string>::const_iterator i = species.begin();
    i++; // skip reference
    for( ; i != species.end(); i++ ) {
        seq = (char*)calloc( n, sizeof(char) );
        memset( seq, '-', n );
        msa[*i] = seq;
    }
    fprintf( stderr, "done.\n" );

    parse( emf_f );
}

emfparse::~emfparse()
{
    map<string,char*>::iterator i;
    for( i = msa.begin(); i != msa.end(); i++ ) {
        free(i->second);
    }
}


void emfparse::parse_ref( const char* ref_fn, char** seqname_out, char** seq_out )
{
    fprintf( stderr, "reading reference from %s ... ", ref_fn );
    FILE* f = fopen( ref_fn, "r" );
    if( f == NULL ) {
        fprintf( stderr, "Can't open file %s\n", ref_fn );
        exit(1);
    }

    char* buf = new char[BUFSIZ];
    char* c;

    enum fa_state { fa_state_eol, fa_state_seqname, fa_state_seq };

    fa_state s = fa_state_eol;

    string seqname;
    size_t n = 1024;
    size_t j = 0;
    char* seq = (char*)calloc( n, sizeof(char) );

    while( fgets( buf, BUFSIZ, f ) != NULL ) {
        for( c = buf; *c; c++ ) {
            if( *c == '\n' ) {
                s = fa_state_eol;
                continue;
            }

            switch( s ) {
                case fa_state_eol:
                    if( *c == '>' ) {
                        if( !seqname.empty() ) goto parse_ref_done;
                        s = fa_state_seqname;
                        break;
                    }
                    s = fa_state_seq;
                    // fall through

                case fa_state_seq:
                    if( isnt(*c) ) {
                        if( j+1 >= n ) {
                            n *= 2;
                            seq = (char*)safe_realloc( seq, n*sizeof(char) );
                        }
                        seq[j++] = *c;
                    }
                    break;

                case fa_state_seqname:
                    seqname += *c;
                    break;
            }
        }
    }
parse_ref_done:

    seq = (char*)realloc( seq, (j+1)*sizeof(char) );
    seq[j] = '\0';

    *seqname_out = strdup(seqname.c_str());
    *seq_out     = seq;

    fprintf( stderr, "done. (%zunt)\n", strlen(seq) );
}




struct pos
{
    int col, row;
    size_t start, end;
    int strand;
};




void emfparse::parse( FILE* f )
{
    fprintf( stderr, "parsing EMF file ... \n" );

    size_t maxline = 10000;
    char* line = new char[maxline];
    size_t linelen;


    size_t blocknum = 0;
 

    // regular expression for SEQ lines
    pcrecpp::RE seq_re( "^SEQ\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s+(\\S+)\\s++(\\S+)" );


    // parser states
    enum emf_state {
        emf_state_header, // reading block header
        emf_state_data,   // reading alignment
    };

    emf_state s = emf_state_header;

    pos tmp_pos;
    std::vector<pos> ref_pos;
    std::vector<pos>::iterator i;


    string seq_name, seq_chrom, seq_start, seq_end, seq_strand;

    // keep track of the order of the columns
    int seq_col_num = 0;
    std::vector<char*> seqs;


    // keep track of the current row and column, relative to the start of the block
    int col;


    while( fgets( line, maxline, f ) != NULL ) {

        linelen = strlen(line);

        if( line[0] == '#' ) continue;                     // skip comments
        if( line[0] == '\0' || line[0] == '\n' ) continue; // skip blank lines

        if( s == emf_state_header ) {
            if( linelen >= 3 && strncmp( line, "SEQ", 3 ) == 0 ) {
                if( !seq_re.PartialMatch( line, &seq_name, &seq_chrom,
                                             &seq_start,&seq_end, &seq_strand ) ) {
                    fprintf( stderr, "Malformed SEQ line: %s", line ); 
                    exit(1);
                }

                if( seq_name == ref_name ) {

                    seqs.push_back(msa[seq_name]);

                    if( seq_chrom == ref_seqname ) {
                        tmp_pos.col    = seq_col_num;
                        tmp_pos.row    = 0;
                        tmp_pos.start  = atoi(seq_start.c_str());
                        tmp_pos.end    = atoi(seq_end.c_str());
                        tmp_pos.strand = atoi(seq_strand.c_str());
                        ref_pos.push_back(tmp_pos);
                    }
                }
                else if( msa.find(seq_name) != msa.end() ) {
                    seqs.push_back(msa[seq_name]);
                }
                else {
                    seqs.push_back(NULL);
                }

                seq_col_num++;
            }
            if( linelen >= 4 && strncmp( line, "DATA", 4 ) == 0 ) {
                s = emf_state_data;
            }
        }


        else if( s == emf_state_data ) {
            if( linelen >= 2 && strncmp( line, "//", 2 ) == 0 ) {
                blocknum++;
                s = emf_state_header;
                seqs.clear();
                ref_pos.clear();
                seq_col_num = 0;
                if( blocknum % 100 == 0 ) fprintf( stderr, "\t%zd blocks read ...\n", blocknum );
                continue;
            }

            if( ref_pos.size() > 1 ) {
                fprintf( stderr, "(%d)\n", ref_pos.size() );
            }

            for( i = ref_pos.begin(); i != ref_pos.end(); i++ ) {
                if( line[i->col] == '-' ) continue;

                for( col = 0; !isspace(line[col]); col++ ) {
                    if( seqs[col] == NULL ) continue;

                    // positive strand
                    else if( i->strand == 1 ) {
                        if( i->start + i->row >= n ) {
                            fprintf( stderr, "Index out of bounds: %zd / %zd\n",
                                     i->start + i->row - 1, n );
                            exit(1);
                        }

                        if( seqs[col] == ref_seq ) {
                            if( tolower(seqs[col][i->start + i->row - 1]) != tolower(line[col]) ) {
                                fprintf( stderr, "Ref mismatch (+): %c <-> %c\n", seqs[col][i->start + i->row - 1], line[col] );
                            }
                        }
                        else {
                            seqs[col][i->start + i->row - 1] = line[col];
                        }
                    }

                    // negative strand
                    else {
                        if( i->end - i->row - 1 >= n ) {
                            fprintf( stderr, "Index out of bounds: %zd / %zd\n",
                                     i->end - i->row - 1, n );
                            exit(1);
                        }

                        if( seqs[col] == ref_seq ) {
                            if( tolower(seqs[col][i->end - i->row - 1]) != tolower(complement(line[col])) ) {
                                fprintf( stderr, "Ref mismatch (+): %c <-> %c\n", seqs[col][i->end - i->row - 1], complement(line[col]) );
                            }
                        }
                        else {
                            seqs[col][i->end- i->row - 1] = complement(line[col]);
                        }
                    }
                }

                i->row++;
            }
        }
    }

    delete[] line;

    fprintf( stderr, "done.\n" );
}


void emfparse::print_fasta( FILE* f )
{
    fprintf( stderr, "printing fasata ... " );

    /* print reference first */
    char* tmp;
    asprintf( &tmp, "%s.%s", ref_name.c_str(), ref_seqname );
    print_fasta_sub( f, tmp, msa[ref_name] );
    free( tmp );

    map<string,char*>::iterator i;
    for( i = msa.begin(); i != msa.end(); i++ ) {
        if( i->first == ref_name ) continue;
        print_fasta_sub( f, i->first.c_str(), i->second );
    }

    fprintf( stderr, "done.\n" );
}

void emfparse::print_fasta_sub( FILE* f, const char* seqname, const char* seq )
{
    const int linewidth = 79;
    char* fmt;
    asprintf( &fmt, "%%.%ds\n", linewidth );

    fprintf( f, ">%s\n", seqname );
    size_t i = 0;

    while( i < n ) {
        fprintf( f, fmt, seq+i );
        i += linewidth;
    }
        
    free(fmt);
}



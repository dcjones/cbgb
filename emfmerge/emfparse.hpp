

#include <vector>
#include <string>
#include <map>

using namespace std;

class emfparse
{
    public:
        emfparse( vector<string> species, const char* ref_fn, FILE* emf_f );
        ~emfparse();

        void print_fasta( FILE* f );

    private:
        void parse_ref( const char* ref_fn, char** seqname, char** seq );
        void parse( FILE* f );

        void print_fasta_sub( FILE* f, const char* seqname, const char* seq );

        size_t n;
        char* ref_seqname;
        string ref_name;
        map<string,char*> msa;
};


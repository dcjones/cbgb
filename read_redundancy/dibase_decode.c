

const char decode_dict[5][5] =
  /*            0    1    2    3    .       */
  /* A */  { { 'A', 'C', 'G', 'T', 'N' },
  /* C */    { 'C', 'A', 'T', 'G', 'N' },
  /* G */    { 'G', 'T', 'A', 'C', 'N' },
  /* T */    { 'T', 'G', 'C', 'A', 'N' },
  /* N */    { 'N', 'N', 'N', 'N', 'N' } };




/* convert characters to number for indexing into arrays */
int to_num( char c )
{
    switch( c ) {
        case 'A': case 'a': case '0': return 0;
        case 'C': case 'c': case '1': return 1;
        case 'G': case 'g': case '2': return 2;
        case 'T': case 't': case '3': return 3;
        case 'N': case 'n': case '.': default: return 4;
    }
}



void dibase_decode( char* seq_in, char* seq_out )
{
    /* copy the primer base */
    seq_out[0] = seq_in[0];

    int i,j;
    int k = 1;

    while( seq_in[k] ) {
        i = to_num( seq_out[k-1] );
        j = to_num( seq_in[k] );
        seq_out[k] = decode_dict[i][j];
        k++;
    }
    seq_out[k] = '\0';
}



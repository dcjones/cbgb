#!/usr/bin/env python


from sys import stdin, stdout, stderr
import re

seq_re = re.compile( r'^SEQ' )
int_re = re.compile( r'\d+' )

for line in stdin:
    if seq_re.match( line ):
        row = re.split( '\s+', line )
        if row[2] == 'MT':
            row[2] = 'chrM'
        elif int_re.match(row[2]):
            row[2] = 'chr' + row[2]
        line = ' '.join( row ) + '\n'

    stdout.write( line )




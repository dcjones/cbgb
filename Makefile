
all : fastqUntail

fastqUntail : fastqUntail.c
	gcc -O3 -g -Wall -D_FILE_OFFSET_BITS=64 -o fastqUntail fastqUntail.c

clean : 
	rm -f fastqUntail

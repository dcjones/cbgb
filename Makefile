
all : fastqUntail

fastqUntail : fastqUntail.c
	gcc -O3 -g -Wall -o fastqUntail fastqUntail.c

clean : 
	rm -f fastqUntail

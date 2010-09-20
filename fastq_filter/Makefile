
CFLAGS=-g -Wall -O2
INC=
LIBS=-lbam -lz

all: fastq_filter

fastq_filter : fastq_filter.o hash.o superfasthash.o
	gcc $(CFLAGS) -o fastq_filter fastq_filter.o hash.o superfasthash.o ${LIBS}

.c.o :
	gcc ${CFLAGS} ${INC} -c $<

clean:
	rm -f *.o fastq_filter






CFLAGS=-g -Wall -O0

all : bamToBedGraph

bamToBedGraph : bamToBedGraph.c
	gcc $(CFLAGS) -o bamToBedGraph bamToBedGraph.c -lbam -lz

clean :
	rm -f bamToBedGraph


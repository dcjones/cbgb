


CFLAGS=-g -Wall -O3

all : bamToBedGraph

bamToBedGraph : bamToBedGraph.c
	gcc $(CFLAGS) -o bamToBedGraph bamToBedGraph.c -lbam -lz

clean :
	rm -f bamToBedGraph


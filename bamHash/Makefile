
CFLAGS=-Wall -g -O3
LIBS=-lbam -lz
OBJ=hash.o bamHash.o

all : bamHash summarize

.c.o :
	gcc ${CFLAGS} ${INC} -c $<

bamHash :  ${OBJ}
	gcc ${CFLAGS} -o $@ ${OBJ} ${LIBS}

summarize : summarize.o
	gcc ${CFLAGS} -o $@ $^

clean:
	rm -rf ${OBJ} bamHash summarize




CFLAGS=-Wall -g -O3
LIBS=-lbam -lz
OBJ=hash.o bamHash.o

all : bamHash

.c.o :
	gcc ${CFLAGS} ${INC} -c $<

bamHash :  ${OBJ}
	gcc ${CFLAGS} -o bamHash ${OBJ} ${LIBS}

clean:
	rm -rf ${OBJ} bamHash



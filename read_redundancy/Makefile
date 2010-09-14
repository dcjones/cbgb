
CFLAGS=-g -Wall -O3
INC=-I${HOME}/include
LIB=${HOME}/lib/libbam.a -lz

OBJ=read_redundancy.o hash.o superfasthash.o

all : read_redundancy

clean :
	rm -f ${OBJ}

.c.o :
	gcc ${CFLAGS} ${INC} -c $<

read_redundancy : ${OBJ}
	gcc ${CFLAGS} ${INC} -o read_redundancy ${OBJ} ${LIB}




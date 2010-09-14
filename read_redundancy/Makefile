
CFLAGS=-g -Wall
INC=
LIB=-lbam -lz

OBJ=read_redundancy.o hash.o

all : read_redundancy

clean :
	rm -f ${OBJ}

.c.o :
	gcc ${CFLAGS} ${INC} -c $<

read_redundancy : ${OBJ}
	gcc ${CFLAGS} ${INC} -o read_redundancy ${OBJ} ${LIB}






all : readcount


CFLAGS=$(shell python-config --cflags )
INCLUDES=$(shell python-config --includes )
LIBS=$(shell python-config --libs ) -lbam -lz


readcount.c : readcount.pyx
	cython --embed readcount.pyx

readcount : readcount.c
	gcc ${CFLAGS} ${INCLUDES} ${LIBS} -o readcount readcount.c /usr/local/lib/libbam.a


clean:
	rm -f readcount.c readcount


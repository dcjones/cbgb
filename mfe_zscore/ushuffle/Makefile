
CFLAGS=-O3 -Wall
PY_INCLUDE=/usr/include/python2.6/
PY_INSTALL=/usr/local/lib/python2.6/dist-packages/

all: ushuffle ushuffle.o

ushuffle.o: ushuffle.c ushuffle.h
	gcc $(CFLAGS) -c ushuffle.c

ushuffle: main.c ushuffle.o 
	gcc $(CFLAGS) main.c ushuffle.o -o ushuffle

ushuffle.so: ushuffle.o ushufflemodule.c
	gcc $(CFLAGS) -I$(PY_INCLUDE) -shared ushuffle.o ushufflemodule.c -o ushuffle.so

install: all
	install ushuffle /usr/local/bin
	install ushuffle.so $(PY_INSTALL)


clean:
	rm -rf *.o ushuffle




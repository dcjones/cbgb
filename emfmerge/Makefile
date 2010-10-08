
CPPFLAGS=-Wall -g -D_FILE_OFFSET_BITS=64 -I${HOME}/include -O3

all : emfmerge


%.o : %.cpp
	g++ -c $(CPPFLAGS) $<


emfmerge : emfmerge.o emfparse.o
	g++ -o $@ $^ ${HOME}/lib/libpcrecpp.a ${HOME}/lib/libpcre.a

clean :
	rm -f emfmerge *.o


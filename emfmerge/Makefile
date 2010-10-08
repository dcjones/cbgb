
CPPFLAGS=-Wall -g -D_FILE_OFFSET_BITS=64

all : emfmerge


%.o : %.cpp
	g++ -c $(CPPFLAGS) $<


emfmerge : emfmerge.o emfparse.o
	g++ -o $@ $^ -lpcrecpp

clean :
	rm -f emfmerge *.o


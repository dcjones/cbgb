
CPPFLAGS=-Wall -g -D_FILE_OFFSET_BITS=64 -I${HOME}/include -O3

all : emfmerge emf_ensembl_to_ucsc


%.o : %.cpp
	g++ -c $(CPPFLAGS) $<


emfmerge : emfmerge.o emfparse.o
	g++ -o $@ $^ ${HOME}/lib/libpcrecpp.a ${HOME}/lib/libpcre.a


emf_ensembl_to_ucsc : emf_ensembl_to_ucsc.c
	gcc $(CPPFLAGS) -o $@ $< -lpcre


clean :
	rm -f emfmerge *.o


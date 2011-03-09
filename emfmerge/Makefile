
CPPFLAGS=-Wall -g -D_FILE_OFFSET_BITS=64 -I${HOME}/include -O3

pcre_libs=-lpcre
pcre_cpp_libs=-lpcre -lpcrecpp

all : emfmerge emf_ensembl_to_ucsc


%.o : %.cpp
	g++ -c $(CPPFLAGS) $<


emfmerge : emfmerge.o emfparse.o
	g++ -o $@ $^ $(pcre_cpp_libs)


emf_ensembl_to_ucsc : emf_ensembl_to_ucsc.c
	gcc $(CPPFLAGS) -o $@ $< $(pcre_libs)


clean :
	rm -f emfmerge *.o emf_ensembl_to_ucsc


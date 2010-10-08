
CPPFLAGS=-Wall -g -D_FILE_OFFSET_BITS=64

all : emfmerge emf_ensembl_to_ucsc


%.o : %.cpp
	g++ -c $(CPPFLAGS) $<


emfmerge : emfmerge.o emfparse.o
	g++ -o $@ $^ -lpcrecpp


emf_ensembl_to_ucsc : emf_ensembl_to_ucsc.c
	gcc $(CPPFLAGS) -o $@ $< -lpcre


clean :
	rm -f emfmerge *.o


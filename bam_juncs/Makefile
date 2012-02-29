

samtools_obj = $(subst .c,.o, $(shell ls samtools/*.c))
obj = bam-juncs.o \


all : bam-juncs bam-gsnap-juncs

CFLAGS=-D_USE_KNETFILE -D_FILE_OFFSET_BITS=64 -g -Wall -O2

%.o : %.c
	gcc $(CFLAGS) -c $< -o $@


bam-juncs : bam-juncs.o $(samtools_obj)
	gcc -o $@ $^ -lz

bam-gsnap-juncs : bam-gsnap-juncs.o $(samtools_obj)
	gcc -o $@ $^ -lz

clean :
	rm -f *.o samtools/*.o bam-juncs bam-gsnap-juncs


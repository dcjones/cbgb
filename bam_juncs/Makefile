

obj = bam-juncs.o \
	  $(subst .c,.o, $(shell ls samtools/*.c))


all : bam-juncs

CFLAGS=-D_USE_KNETFILE -D_FILE_OFFSET_BITS=64 -g -Wall -O2

%.o : %.c
	gcc $(CFLAGS) -c $< -o $@


bam-juncs : $(obj)
	gcc -o $@ $^ -lz


clean :
	rm -f *.o samtools/*.o bam-juncs


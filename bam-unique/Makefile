

obj = bam-unique.o hash.o common.o \
	  $(subst .c,.o, $(shell ls samtools/*.c))



all : bam-unique

CFLAGS=-D_USE_KNETFILE -D_FILE_OFFSET_BITS=64 -g -Wall -O2

%.o : %.c
	gcc $(CFLAGS) -c $< -o $@


bam-unique : $(obj)
	gcc -o $@ $^ -lz


clean :
	rm -f *.o samtools/*.o bam-unique


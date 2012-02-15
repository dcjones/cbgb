

obj = bam-summarize.o \
	  $(subst .c,.o, $(shell ls samtools/*.c)) \
	  $(subst .c,.o, $(shell ls hat-trie/*.c))



all : bam-summarize

CFLAGS=-D_USE_KNETFILE -D_FILE_OFFSET_BITS=64 -g -Wall -O2

%.o : %.c
	gcc $(CFLAGS) -c $< -o $@


bam-summarize : $(obj)
	gcc -o $@ $^ -lz


clean :
	rm -f *.o samtools/*.o bam-summarize


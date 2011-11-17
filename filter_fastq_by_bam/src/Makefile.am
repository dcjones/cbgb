
SUBDIRS = hat-trie samtools

bin_PROGRAMS = ffbb

ffbb_SOURCES = ffbb.c parse.h parse.c common.h common.c
ffbb_LDADD = samtools/libbam.la hat-trie/libhat-trie.la -lz



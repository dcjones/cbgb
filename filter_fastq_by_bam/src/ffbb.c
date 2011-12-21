
#include "parse.h"
#include "hat-trie/hat-trie.h"
#include "samtools/sam.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>


void print_usage(FILE* fout)
{
    fprintf(fout,
           "Usage: ffbb output_prefix alignments1.bam[,alignments2.bam,...] reads2.fastq [reads2.fastq]\n");
}


int main(int argc, char* argv[])
{
    if (argc < 4) {
        print_usage(stderr);
        exit(EXIT_FAILURE);
    }

    char* prefix = argv[1];

    char* bam_fns = argv[2];

    char* reads1_fn = argv[3];
    char* reads2_fn = NULL;
    if (argc > 4) reads2_fn = argv[4];




    /* hash ids of mapped reads */

    printf("hashing ids ...\n");
    hattrie_t* ids = hattrie_create();
    samfile_t* sf;
    bam1_t* read = bam_init1();

    size_t count = 0;
    char* id;
    char* bam_fn = strtok(bam_fns, ",");
    while (bam_fn) {
        sf = samopen(bam_fn, "rb", NULL);
        if (sf == NULL) {
            fprintf(stderr, "Can't open BAM file %s.\n", bam_fn);
            continue;
        }

        while (samread(sf, read) > 0) {
            if (read->core.flag & BAM_FUNMAP) continue;

            id = bam1_qname(read);
            hattrie_get(ids, id, strlen(id));

            if (++count % 100000 == 0) {
                printf("\t%zu reads, %zu ids\n", count, hattrie_size(ids));
            }
        }

        samclose(sf);
        bam_fn = strtok(NULL, ",");
    }

    bam_destroy1(read);
    printf("done. (%zu reads, %zu ids)\n", count, hattrie_size(ids));



    /* filter fastq files */
    FILE* f1 = fopen(reads1_fn, "r");
    if (f1 == NULL) {
        fprintf(stderr, "Can't open FASTQ file %s.\n", reads1_fn);
        exit(EXIT_FAILURE);
    }
    fastq_t* fq1 = fastq_open(f1);
    seq_t* read1 = fastq_alloc_seq();

    char fn[256];

    if (reads2_fn != NULL) snprintf(fn, sizeof(fn), "%s_1.fastq", prefix);
    else                   snprintf(fn, sizeof(fn), "%s.fastq", prefix);

    FILE* fout1 = fopen(fn, "w");
    if (fout1 == NULL) {
        fprintf(stderr, "Can't open %s for writing.\n", fn);
        exit(EXIT_FAILURE);
    }


    FILE* f2 = NULL;
    fastq_t* fq2 = NULL;
    seq_t* read2 = NULL;
    FILE* fout2 = NULL;

    if (reads2_fn != NULL) {
        f2 = fopen(reads2_fn, "r");
        if (f2 == NULL) {
            fprintf(stderr, "Can't open FASTQ file %s.\n", reads2_fn);
            exit(EXIT_FAILURE);
        }
        fq2 = fastq_open(f2);
        read2 = fastq_alloc_seq();

        snprintf(fn, sizeof(fn), "%s_2.fastq", prefix);
        fout2 = fopen(fn, "w");
        if (fout2 == NULL) {
            fprintf(stderr, "Can't open %s for writing.\n", fn);
            exit(EXIT_FAILURE);
        }

        printf("writing unmapped reads to [%s] and [%s] ...\n", reads1_fn, reads2_fn);
    }
    else {
        printf("writing unmapped reads to [%s] ...\n", reads1_fn);
    }


    count = 0;

    char* s;
    char* t;
    size_t idlen;

    if (reads2_fn) {
        while (fastq_next(fq1, read1) && fastq_next(fq2, read2)) {
            s = strchr(read1->id1.s, '/');
            t = strchr(read1->id1.s, ' ');
            if (s && t) idlen = (s < t ? s : t) - read1->id1.s;
            else if (s) idlen = s - read1->id1.s;
            else if (t) idlen = t - read1->id1.s;
            else        idlen = read1->id1.n;


            if (hattrie_tryget(ids, read1->id1.s, idlen) == 0) {
                if (++count % 100000 == 0) printf("\t%zu reads.\n", count);
                fastq_print(fout1, read1);
                fastq_print(fout2, read2);
            }
        }
    }
    else {
        while (fastq_next(fq1, read1)) {
            if (hattrie_tryget(ids, read1->id1.s, read1->id1.n) == 0) {
                if (++count % 100000 == 0) printf("\t%zu reads.\n", count);
                fastq_print(fout1, read1);
            }
        }
    }


    printf("done. (%zu reads)\n", count);


    fclose(fout1);
    fastq_close(fq1);
    fastq_free_seq(read1);

    if (reads2_fn != NULL) {
        fclose(fout2);
        fastq_close(fq2);
        fastq_free_seq(read2);
    }

    hattrie_free(ids);

    return 0;
}


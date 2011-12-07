
#include "hash.h"
#include "samtools/sam.h"
#include "samtools/bam.h"
#include <stdio.h>



void usage()
{
    fprintf(stderr, "Usage: bam-unique in.bam\n");
}



hash_table* hash_ids(const char* fn)
{
    fprintf(stderr, "hashing ... \n");

    hash_table* T = create_hash_table();

    samfile_t* f = samopen(fn, "rb", NULL);
    if (f == NULL) {
        fprintf(stderr, "can't open bam file %s\n", fn);
        exit(1);
    }

    bam1_t* b = bam_init1();

    uint32_t n = 0;

    char* qname = NULL;
    size_t qname_size = 0;

    while (samread(f, b) >= 0) {
        if (++n % 1000000 == 0) {
            fprintf(stderr, "\t%d reads\n", n);
        }

        if (qname_size < b->core.l_qname + 3) {
            qname_size = b->core.l_qname + 3;
            qname = realloc(qname, qname_size);
        }

        memcpy(qname, bam1_qname(b), b->core.l_qname);

        if (b->core.flag & BAM_FREAD2) {
            qname[b->core.l_qname]     = '/';
            qname[b->core.l_qname + 1] = '2';
            qname[b->core.l_qname + 2] = '\0';
        }
        else {
            qname[b->core.l_qname]     = '/';
            qname[b->core.l_qname + 1] = '1';
            qname[b->core.l_qname + 2] = '\0';
        }


        inc_hash_table(T, qname, b->core.l_qname + 2);
    }

    free(qname);

    bam_destroy1(b);
    samclose(f);

    fprintf(stderr, "done.\n");
    return T;
}



void filter_by_id(const char* fn, hash_table* T)
{
    fprintf(stderr, "filtering ... \n");

    samfile_t* fin = samopen(fn, "rb", NULL);
    if (fin == NULL) {
        fprintf(stderr, "can't open bam file %s\n", fn);
        exit(1);
    }

    samfile_t* fout = samopen("-", "w", (void*)fin->header);
    if (fout == NULL) {
        fprintf(stderr, "can't open stdout, for some reason.\n");
        exit(1);
    }

    fputs(fin->header->text, stdout);

    bam1_t* b = bam_init1();
    uint32_t n = 0;

    char* qname = NULL;
    size_t qname_size = 0;

    while (samread(fin, b) >= 0) {
        if (++n % 1000000 == 0) {
            fprintf(stderr, "\t%d reads\n", n);
        }


        if (qname_size < b->core.l_qname + 3) {
            qname_size = b->core.l_qname + 3;
            qname = realloc(qname, qname_size);
        }

        memcpy(qname, bam1_qname(b), b->core.l_qname);

        if (b->core.flag & BAM_FREAD2) {
            qname[b->core.l_qname]     = '/';
            qname[b->core.l_qname + 1] = '2';
            qname[b->core.l_qname + 2] = '\0';
        }
        else {
            qname[b->core.l_qname]     = '/';
            qname[b->core.l_qname + 1] = '1';
            qname[b->core.l_qname + 2] = '\0';
        }

        if (get_hash_table(T, qname, b->core.l_qname + 2) == 1) {
            samwrite(fout, b);
        }
    }

    free(qname);

    bam_destroy1(b);
    samclose(fout);
    samclose(fin);

    fprintf(stderr, "done.\n");
}



int main(int argc, char* argv[])
{
    if (argc < 2) {
        usage();
        exit(1);
    }

    hash_table* T = hash_ids(argv[1]);
    filter_by_id(argv[1], T);

    destroy_hash_table(T);

    return 0;
}




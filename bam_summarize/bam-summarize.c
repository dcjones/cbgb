

#include "samtools/sam.h"
#include "hat-trie/hat-trie.h"
#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>


/* any gap at or above this length is considered a splice junction */
const int min_splice_length = 50;


/*
 * This is a ugly hodge-podge script that generate a bunch of information from a
 * BAM file.
 * 
 *   1. Count the number of reads with perfect alignments.
 *   2. Output the number of reads with multiple alignments.
 *   3. Output junctions.
 *   4. Output SNPs.
 *
 */

typedef struct
{
    /* how many alignments exist for the read */
    uint32_t aln_count;

    uint32_t unspliced_perfect_cnt;
    uint32_t spliced_perfect_cnt;

    uint32_t spliced_cnt;
    uint32_t gapped_cnt;

} read_stat_t;



int main(int argc, char* argv[])
{
    if (argc < 2) {
        fprintf(stderr, "Usage: bam-summarize reads.bam\n");
        exit(EXIT_FAILURE);
    }

    samfile_t* f = samopen(argv[1], "rb", NULL);
    if (f == NULL) {
        fprintf(stderr, "can't open bam file %s\n", argv[1]);
        exit(1);
    }

    bam1_t* b = bam_init1();


    hattrie_t* T = hattrie_create();

    char* qname = NULL;
    size_t qname_size = 0;

    size_t j, n = 0;
    uint32_t* cigar;
    uint32_t cigar_op, cigar_len;

    read_stat_t** val;

    while (samread(f, b) >= 0) {
        if (++n % 1000000 == 0) {
            fprintf(stderr, "\t%zu alignments\n", n);
        }

        bool perfect = true;
        bool spliced = false;
        bool gapped  = false;

        cigar = bam1_cigar(b);
        for (j = 0; j < b->core.n_cigar; ++j) {
            cigar_op  = cigar[j] & BAM_CIGAR_MASK;
            cigar_len = cigar[j] >> BAM_CIGAR_SHIFT;

            if (cigar_op == BAM_CREF_SKIP) {
                if (cigar_len < min_splice_length) gapped = true;
                else                               spliced = true;
            }
            else if (cigar_op != BAM_CMATCH)  perfect = false;

            if (cigar_op == BAM_CSOFT_CLIP || cigar_op == BAM_CHARD_CLIP) break;
        }

        /* Skip any clipped alignments. We don't want your kind! */
        if (cigar_op == BAM_CSOFT_CLIP || cigar_op == BAM_CHARD_CLIP) continue;

        /* Hack the read to include mate information. */
        if (b->core.flag & BAM_FPAIRED) {
            if (qname_size < b->core.l_qname + 3) {
                qname_size = b->core.l_qname + 3;
                qname = realloc(qname, qname_size);
            }
            memcpy(qname, bam1_qname(b), b->core.l_qname);

            if (b->core.flag & BAM_FREAD1) {
                qname[b->core.l_qname]     = '/';
                qname[b->core.l_qname + 1] = '2';
                qname[b->core.l_qname + 2] = '\0';
            }
            else {
                qname[b->core.l_qname]     = '/';
                qname[b->core.l_qname + 1] = '1';
                qname[b->core.l_qname + 2] = '\0';
            }

            val = (read_stat_t**) hattrie_get(T, qname, b->core.l_qname + 2);
        }
        else {
            val = (read_stat_t**) hattrie_get(T, bam1_qname(b), b->core.l_qname);
        }


        if (*val == NULL) {
            *val = malloc(sizeof(read_stat_t));
            memset(*val, 0, sizeof(read_stat_t));
        }

        (*val)->aln_count++;
        if (perfect) {
            if (spliced) (*val)->spliced_perfect_cnt++;
            else         (*val)->unspliced_perfect_cnt++;
        }

        if (spliced) (*val)->spliced_cnt++;
        if (gapped) (*val)->gapped_cnt++;
    }

    printf("alignment_count\t%zu\n", n);
    printf("read_count\t%zu\n", hattrie_size(T));


    /* print stats from the table */

    uint32_t multi_count = 0;

    uint32_t unspliced_perfect_cnt = 0;
    uint32_t spliced_perfect_cnt = 0;

    uint32_t spliced_cnt = 0;
    uint32_t gapped_cnt = 0;

    /* excluding multireads */
    uint32_t unique_unspliced_perfect_cnt = 0;
    uint32_t unique_spliced_perfect_cnt = 0;

    uint32_t unique_spliced_cnt = 0;
    uint32_t unique_gapped_cnt = 0;



    hattrie_iter_t* i;
    for (i = hattrie_iter_begin(T);
         !hattrie_iter_finished(i);
         hattrie_iter_next(i))
    {
        val = (read_stat_t**) hattrie_iter_val(i);

        if ((*val)->aln_count == 1) {
            unique_unspliced_perfect_cnt += (*val)->unspliced_perfect_cnt;
            unique_spliced_perfect_cnt   += (*val)->spliced_perfect_cnt;

            unique_spliced_cnt += (*val)->spliced_cnt;
            unique_gapped_cnt  += (*val)->gapped_cnt;
        }
        else multi_count++;

        unspliced_perfect_cnt += (*val)->unspliced_perfect_cnt;
        spliced_perfect_cnt   += (*val)->spliced_perfect_cnt;

        spliced_cnt += (*val)->spliced_cnt;
        gapped_cnt  += (*val)->gapped_cnt;
    }

    hattrie_iter_free(i);


    printf("multi_count\t%u\n", multi_count);
    printf("unspliced_perfect_cnt\t%u\n", unspliced_perfect_cnt);
    printf("spliced_perfect_cnt\t%u\n", spliced_perfect_cnt);
    printf("spliced_cnt\t%u\n", spliced_cnt);
    printf("gapped_cnt\t%u\n", gapped_cnt);

    printf("unique_unspliced_perfect_cnt\t%u\n", unique_unspliced_perfect_cnt);
    printf("unique_spliced_perfect_cnt\t%u\n", unique_spliced_perfect_cnt);
    printf("unique_spliced_cnt\t%u\n", unique_spliced_cnt);
    printf("unique_gapped_cnt\t%u\n", unique_gapped_cnt);


    /* free the table */
    for (i = hattrie_iter_begin(T);
         !hattrie_iter_finished(i);
         hattrie_iter_next(i))
    {
        free(* (read_stat_t**) hattrie_iter_val(i));
    }

    hattrie_iter_free(i);
    hattrie_free(T);
    free(qname);

    bam_destroy1(b);

    return 0;
}





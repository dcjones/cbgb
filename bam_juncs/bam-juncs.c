

#include "samtools/sam.h"
#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>


/* any gap at or above this length is considered a splice junction */
const int min_splice_length = 50;




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

    size_t j, n = 0;
    uint32_t* cigar;
    uint32_t cigar_op, cigar_len;
    int32_t off;

    while (samread(f, b) >= 0) {
        if (++n % 1000000 == 0) {
            fprintf(stderr, "\t%zu alignments\n", n);
        }

        cigar = bam1_cigar(b);
        for (j = 0, off = 0; j < b->core.n_cigar; ++j) {
            cigar_op  = cigar[j] & BAM_CIGAR_MASK;
            cigar_len = cigar[j] >> BAM_CIGAR_SHIFT;

            if (cigar_op == BAM_CREF_SKIP) {
                printf("%s\t%d\t%d\n",
                        f->header->target_name[b->core.tid],
                        b->core.pos + off,
                        b->core.pos + off + cigar_len);
            }

            off += cigar_len;
        }
    }
    bam_destroy1(b);

    return 0;
}






#include "samtools/sam.h"
#include <stdio.h>
#include <stdbool.h>
#include <stdint.h>


int main(int argc, char* argv[])
{
    if (argc < 2) {
        fprintf(stderr, "Usage: bam-gsnap-juncs reads.bam\n");
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
    uint8_t* xs_tag;
    uint8_t strand;
    const char* seqname;

    while (samread(f, b) >= 0) {
        if (++n % 1000000 == 0) {
            fprintf(stderr, "\t%zu alignments\n", n);
        }

        /* Look for an XS (splice strand) field. */
        if ((xs_tag = bam_aux_get(b, "XS")) == 0) continue;
        strand = xs_tag[1];
        seqname = f->header->target_name[b->core.tid];

        cigar = bam1_cigar(b);
        for (j = 0, off = 0; j < b->core.n_cigar; ++j) {
            cigar_op  = cigar[j] & BAM_CIGAR_MASK;
            cigar_len = cigar[j] >> BAM_CIGAR_SHIFT;

            /* Print introns */
            if (cigar_op == BAM_CREF_SKIP) {
                if (strand == '+') {
                    printf(">_ %s:%d-%d\n",
                        seqname,
                        b->core.pos + off,
                        b->core.pos + off + cigar_len);
                }
                else {
                    printf(">_ %s:%d-%d\n",
                        seqname,
                        b->core.pos + off + cigar_len,
                        b->core.pos + off);
                }
            }

            /* Print splice sites. */
#if 0
            if (cigar_op == BAM_CREF_SKIP) {
                if (strand == '+') {
                    printf(">_ %s:%d..%d donor\n",
                        seqname,
                        b->core.pos + off,
                        b->core.pos + off + 1);
                    printf(">_ %s:%d..%d acceptor\n",
                        seqname,
                        b->core.pos + off + cigar_len - 1,
                        b->core.pos + off + cigar_len);
                }
                else {
                    printf(">_ %s:%d..%d donor\n",
                        seqname,
                        b->core.pos + off + 1,
                        b->core.pos + off);
                    printf(">_ %s:%d..%d acceptor\n",
                        seqname,
                        b->core.pos + off + cigar_len,
                        b->core.pos + off + cigar_len - 1);
                }
            }
#endif

            off += cigar_len;
        }
    }
    bam_destroy1(b);

    return 0;
}

#ifndef P_UTILS_H
#define P_UTILS_H

#include "htslib/htslib/sam.h"
#include "htslib/htslib/khash.h"

KHASH_MAP_INIT_INT(occ, int)

int read_bam(void *data, bam1_t *b);
void count_pos(int n, const bam_pileup1_t *plp, khint_t key, khash_t(occ) *hash);

typedef struct aux_t {

  samFile *in;
  bam_hdr_t *h;
  int min_mapQ;

} aux_t;

#endif

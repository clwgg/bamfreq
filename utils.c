
#include "utils.h"

int read_bam(void *data, bam1_t *b)
{

  aux_t *aux = (aux_t*)data;
  int ret;

  while (1) {
    ret = sam_read1(aux->in, aux->h, b);
    if (ret < 0) break;
    if ( (int)b->core.qual < aux->min_mapQ ) continue;
    if ( b->core.flag & (BAM_FUNMAP | BAM_FDUP) ) continue;
    break;
  }
  return ret;

}

void count_pos(int n, const bam_pileup1_t *plp, khint_t key, khash_t(occ) *hash)
{

  int i;
  for (i = 0; i < n; ++i) {

    const bam_pileup1_t *p = plp + i;
    uint8_t *seq = bam_get_seq(p->b);
    uint8_t nuc = bam_seqi(seq, p->qpos);

    if (p->is_del) nuc = 32;

    int ret;
    key = kh_put(occ, hash, nuc, &ret);
    if (ret) kh_value(hash, key) = 0;
    kh_value(hash, key)++;
  }

}

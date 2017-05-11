#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <string.h>

#include "utils.h"

static int usage(char **argv)
{
  printf("\nPrint base frequencies by position.\n");
  printf("\nUsage: %s [options] file.bam\n\n", argv[0]);
  printf("Options:\n");

  printf("\t-h\tonly print heterozygous positions\n");
  printf("\t-f\tmin fraction of read support\n");
  printf("\t-q\tmin map quality\n");
  printf("\t-c\tmin coverage of position to be reported\n");
  printf("\t-m\tmax coverage of position to be reported\n\n");
  printf("\t-r\tonly print raw counts without bases, coverage or position\n\n");

  printf("Default output:\nSeq\tPos\tCov\tBase1\tCount1\t...\tBaseN\tCountN\n\n");
  printf("Postitions are 1-based\n\n");

  return 1;
}

char t_seq(int i)
{
  char c;
  switch (i) {
    case 1:  c = 'A'; break;
    case 2:  c = 'C'; break;
    case 4:  c = 'G'; break;
    case 8:  c = 'T'; break;
    case 15: c = 'N'; break;
    case 32: c = 'D'; break;
  }
  return c;
}

void p_count(khash_t(occ) *hash)
{

  khint_t key;
  int i = 0;

  for (key = 0; key != kh_end(hash); ++key) {
    if (kh_exist(hash, key)) {
      int count = hash->vals[key];
      if (i > 0) fprintf(stdout, "\t");
      fprintf(stdout, "%d", count);
      ++i;
    }
  }
  fprintf(stdout, "\n");

}

void p_pos(int n, khash_t(occ) *hash, bam_hdr_t *h, int tid, int pos)
{

  pos = pos + 1; // make 1-based for output
  char *seq = h->target_name[tid];

  fprintf(stdout, "%s\t%d\t%d", seq, pos, n);

  khint_t key;

  for (key = 0; key != kh_end(hash); ++key) {
    if (kh_exist(hash, key)) {
      int base = hash->keys[key];
      int count = hash->vals[key];
      fprintf(stdout, "\t%c\t%d", t_seq(base), count);
    }
  }
  fprintf(stdout, "\n");

}

int main(int argc, char **argv)
{

  float minfrac = 0.0;
  int het = 0;
  int minc = 0;
  int maxc = 0;
  int raw = 0;
  int elem;

  bam_plp_t iter;
  const bam_pileup1_t *plp;
  int tid, pos, n;

  aux_t data;
  data.min_mapQ = 0;

  while (( elem = getopt(argc, argv, "m:f:q:hc:r") ) >= 0) {
    switch(elem) {
      case 'f': minfrac = atof(optarg); break;
      case 'q': data.min_mapQ = atoi(optarg); break;
      case 'h': het = 1; break;
      case 'c': minc = atoi(optarg); break;
      case 'm': maxc = atoi(optarg); break;
      case 'r': raw = 1; break;
    }
  }

  if (argc - optind != 1) {
    return usage(argv);
  }

  data.in = sam_open(argv[optind], "r");
  if (!data.in) return usage(argv);

  data.h = sam_hdr_read(data.in);

  iter = bam_plp_init(read_bam, (void*)&data);

  while ((plp = bam_plp_auto(iter, &tid, &pos, &n)) != 0) {

    if (minc && n < minc) continue;
    if (maxc && n > maxc) continue;

    khash_t(occ) *hash;
    khint_t key = 0;
    hash = kh_init(occ);

    count_pos( n, plp, key, hash );

    if (minfrac) {
      for (key = 0; key != kh_end(hash); ++key) {
        if (kh_exist(hash, key) && kh_val(hash, key) < (int)(minfrac * n)) {
          kh_del(occ, hash, key);
        }
      }
    }

    if (het && hash->size < 2) {
      kh_destroy(occ, hash);
    }
    else if (raw) {
      p_count(hash);
      kh_destroy(occ, hash);
    }
    else {
      p_pos(n, hash, data.h, tid, pos);
      kh_destroy(occ, hash);
    }
  }

  bam_plp_destroy(iter);
  bam_hdr_destroy(data.h);
  sam_close(data.in);

  return 0;
}


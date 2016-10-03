/*
 * Obtain pairing information from XF bam field.
 *
 * Copyright (c) 2016 Graham Gower <graham.gower@gmail.com>
 *
 * Permission to use, copy, modify, and distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>

#include <htslib/sam.h>

#include "kseq.h"
KSTREAM_INIT(int, read, 16384);

typedef struct {
	char *bam_fn;
	char *bed_fn;
	int min_mapq;
	int min_baseq;
	int hp_len;
	FILE *metrics_fp;
} opt_t;

typedef struct {
	opt_t *opt;
	samFile *bam_fp;
	bam_hdr_t *bam_hdr;
	hts_idx_t *bam_idx;
	hts_itr_t *bam_iter;
} bam_aux_t;


static int nt2int[256] = {['A']=0, ['C']=1, ['G']=2, ['T']=3};
static char cmap[] = {['A']='T', ['C']='G', ['G']='C', ['T']='A', ['N']='N', ['n']='N',
				['a']='t', ['c']='g', ['g']='c', ['t']='a'};

/*
 * 0<<2|0	0	AA	AA	0
 * 0<<2|1	1	AC	CA	4
 * 0<<2|2	2	AG	GA	8
 * 0<<2|3	3	AT	TA	12
 * 1<<2|0	4	CA	AC	1
 * 1<<2|1	5	CC	CC	5
 * 1<<2|2	6	CG	GC	9
 * 1<<2|3	7	CT	TC	13
 * 2<<2|0	8	GA	AG	2
 * 2<<2|1	9	GC	CG	6
 * 2<<2|2	10	GG	GG	10
 * 2<<2|3	11	GT	TG	14
 * 3<<2|0	12	TA	AT	3
 * 3<<2|1	13	TC	CT	7
 * 3<<2|2	14	TG	GT	11
 * 3<<2|3	15	TT	TT	15
 */
static int nt16rev[16] = { 0, 4, 8, 12, 1, 5, 9, 13,
			2, 6, 10, 14, 3, 7, 11, 15 };

/*
 * Get next alignment.
 */
static int
next_aln(void *data, bam1_t *b)
{
	bam_aux_t *bat = data;
	int ret;

	while (1) {
		ret = sam_itr_next(bat->bam_fp, bat->bam_iter, b);
		if (ret < 0) {
			// iterator exhausted
			break;
		}
	

		if (b->core.flag & (BAM_FUNMAP | BAM_FQCFAIL | BAM_FDUP) ||
				b->core.qual < bat->opt->min_mapq)
			// skip these
			continue;

		break;
	}

	return ret;
}

int
xf2ssqq(char *xf, char **s1, char **s2, char **q1, char **q2)
{
	char *p = xf;
	int len;

	while (*p != 0 && *p != '|')
		p++;

	if (*p == 0)
		return -1;

	len = p-xf;

	*s1 = xf;
	*s2 = xf+len+1;
	*q1 = xf+2*(len+1);
	*q2 = xf+3*(len+1);

	return len;
}

#define WIN_SZ (128*1024)
#define MAX_DP 1024
#define CTX_SZ 30

int
wcache_select(uint8_t *win, uint16_t *win_dp, int wi)
{
	int dp;
	if (wi < 0 || wi >= WIN_SZ)
		return -1;
	dp = win_dp[wi];
	if (dp == 0)
		return -1;
	if (dp == 1)
		return win[wi*MAX_DP];
	return win[wi*MAX_DP +(rand()%dp)];
}

int
scan_pairs(opt_t *opt)
{
	bam_aux_t bat;
	bam_plp_t plpiter;
	const bam_pileup1_t *plp;
	bam1_t *b;
	int tid, pos, n;
	int x, y;
	int i;
	int ret;

	int win_tid, win_pos;
	uint16_t win_dp[WIN_SZ];
	uint8_t *win;

	win = calloc(WIN_SZ*MAX_DP, sizeof(uint8_t));
	if (win == NULL) {
		perror("calloc");
		ret = -1;
		goto err0;
	}

	/*
	 * 5' and 3' contextual pairing info.
	 * ctx5p: The first CTX_SZ slots are upstream of the most 5' position
	 *        and the latter CTX_SZ slots are downstream of the most 5'.
	 * ctx3p: The first CTX_SZ slots are upstream of the most 3' position
	 *        and the latter CTX_SZ slots are downstream of the most 3'.
	 */
	uint64_t ctx5p[2*CTX_SZ][16] = {{0,},}, ctx3p[2*CTX_SZ][16] = {{0,},};

	bat.opt = opt;

	bat.bam_fp = sam_open(opt->bam_fn, "r");
	if (bat.bam_fp == NULL) {
		fprintf(stderr, "bam_open: %s: %s\n", opt->bam_fn, strerror(errno));
		ret = -2;
		goto err1;
	}

	bat.bam_hdr = sam_hdr_read(bat.bam_fp);
	if (bat.bam_hdr == NULL) {
		fprintf(stderr, "%s: couldn't read header\n", opt->bam_fn);
		ret = -3;
		goto err2;
	}

	bat.bam_idx = sam_index_load(bat.bam_fp, opt->bam_fn);
	if (bat.bam_idx == NULL) {
		fprintf(stderr, "Couldn't load index for %s\n", opt->bam_fn);
		ret = -4;
		goto err3;
	}

	b = bam_init1();

	win_tid = 0;
	win_pos = 0;

	while (win_tid < bat.bam_hdr->n_targets) {

		bat.bam_iter = sam_itr_queryi(bat.bam_idx, win_tid, win_pos, win_pos+WIN_SZ);
		if (bat.bam_iter == NULL)
			goto next_window;

		plpiter = bam_plp_init(next_aln, &bat);

		memset(win_dp, 0, sizeof(win_dp));

		/*
		 * Iterate through the pileup to cache nucleotide pairing info
		 * for a window of WIN_SZ nucleotide pairs.
		 */
		while ((plp = bam_plp_auto(plpiter, &tid, &pos, &n)) != 0) {

			if (pos < win_pos || pos >= win_pos+WIN_SZ)
				continue;

			int wi = pos - win_pos;
			int dp = 0;

			for (i=0; i<n; i++) {
				const bam_pileup1_t *p = plp + i;
				uint8_t *xf_aux;
				char *xf;
				char *s1, *s2, *q1, *q2;
				int ci, cj;

				if (p->is_del || p->is_refskip)
					continue;
				
				if (bam_get_qual(p->b)[p->qpos] < opt->min_baseq)
					continue;

				xf_aux = bam_aux_get(p->b, "XF");
				if (xf_aux == NULL)
					continue;
				xf = bam_aux2Z(xf_aux);
				if (xf == NULL)
					continue;

				if (xf2ssqq(xf, &s1, &s2, &q1, &q2) < 0)
					continue;

				if (bam_is_rev(p->b)) {
					ci = cmap[(int)s2[p->b->core.l_qseq - p->qpos-1]];
					cj = s1[p->b->core.l_qseq - p->qpos-1];
					//fprintf(stderr, "[-]%s  %d:%d:%d  %c/%c\n", bam_get_qname(p->b), pos, p->qpos, i, ci, cj);
				} else {
					ci = s1[p->qpos];
					cj = cmap[(int)s2[p->qpos]];
					//fprintf(stderr, "[+]%s  %d:%d:%d  %c/%c\n", bam_get_qname(p->b), pos, p->qpos, i, ci, cj);
				}

				if (ci == 'N' || cj == 'N')
					continue;

				win[wi*MAX_DP + dp] = nt2int[(int)ci]<<2 | nt2int[(int)cj];
				if (++dp == MAX_DP)
					break;
			}
			win_dp[wi] = dp;
		}


		bam_plp_destroy(plpiter);

		// reset the iterator
		sam_itr_destroy(bat.bam_iter);
		bat.bam_iter = sam_itr_queryi(bat.bam_idx, win_tid, win_pos, win_pos+WIN_SZ);
		if (bat.bam_iter == NULL)
			goto next_window;

		/*
		 * Iterate back through the window, accumulating
		 * pairing information for CTX_SZ bases up/downstream of both
		 * ends of each read.  Nucleotide pairing info from outside of
		 * the reads is obtained from win[WIN_SZ], cached above.
		 */
		while (sam_itr_next(bat.bam_fp, bat.bam_iter, b) >= 0) {
			uint8_t *xf_aux;
			char *xf;
			char *s1, *s2, *q1, *q2;
			int len;
			int hairpin;
			int ci, cj;

			//if (b->core.pos < win_pos || b->core.pos >= win_pos+WIN_SZ)
			//	continue;

			if (b->core.l_qseq < 2*CTX_SZ)
				continue;

			if (b->core.flag & (BAM_FUNMAP|BAM_FSECONDARY|BAM_FQCFAIL|BAM_FDUP|BAM_FSUPPLEMENTARY))
				continue;

			xf_aux = bam_aux_get(b, "XF");
			if (xf_aux == NULL)
				continue;
			xf = bam_aux2Z(xf_aux);
			if (xf == NULL)
				continue;

			if ((len = xf2ssqq(xf, &s1, &s2, &q1, &q2)) < 0)
				continue;

			if (b->core.l_qseq+opt->hp_len < len)
				hairpin = 1;
			else
				hairpin = 0;

			if (bam_is_rev(b)) {
				// within the read
				for (x=0; x<CTX_SZ; x++) {
					// 5'
					ci = cmap[(int)s2[b->core.l_qseq-x-1]];
					cj = s1[b->core.l_qseq-x-1];
					if (ci != 'N' && cj != 'N')
						ctx5p[CTX_SZ+x][(nt2int[(int)ci]<<2)|nt2int[(int)cj]]++;

					if (!hairpin)
						continue;

					// 3'
					ci = cmap[(int)s2[x]];
					cj = s1[x];
					if (ci != 'N' && cj != 'N')
						ctx3p[CTX_SZ-x-1][(nt2int[(int)ci]<<2)|nt2int[(int)cj]]++;
				}

				// upstream and downstream context of the read
				int off5p = b->core.pos - win_pos + b->core.l_qseq;
				int off3p = b->core.pos - win_pos;
				//fprintf(stderr, "[-]%s  %d, %d, %d:%d:%d\n", bam_get_qname(b), win_pos, b->core.pos, off5p, b->core.l_qseq, off3p);
				for (x=0; x<CTX_SZ; x++) {
					// upstream
					y = wcache_select(win, win_dp, off5p+CTX_SZ-x-1);
					if (y > 0)
						ctx5p[x][nt16rev[y]]++;

					if (!hairpin)
						continue;

					// downstream
					y = wcache_select(win, win_dp, off3p-x-1);
					if (y > 0)
						ctx3p[CTX_SZ+x][nt16rev[y]]++;
				}
			} else {
				// within the read
				for (x=0; x<CTX_SZ; x++) {
					// 5'
					ci = s1[x];
					cj = cmap[(int)s2[x]];
					if (ci != 'N' && cj != 'N')
						ctx5p[CTX_SZ+x][(nt2int[(int)ci]<<2)|nt2int[(int)cj]]++;

					if (!hairpin)
						continue;

					// 3'
					ci = s1[b->core.l_qseq-x-1];
					cj = cmap[(int)s2[b->core.l_qseq-x-1]];
					if (ci != 'N' && cj != 'N')
						ctx3p[CTX_SZ-x-1][(nt2int[(int)ci]<<2)|nt2int[(int)cj]]++;
				}

				// upstream and downstream context of the read
				int off5p = b->core.pos - win_pos;
				int off3p = b->core.pos+b->core.l_qseq - win_pos;
				//fprintf(stderr, "[+]%s  %d, %d, %d:%d:%d\n", bam_get_qname(b), win_pos, b->core.pos, off5p, b->core.l_qseq, off3p);
				for (x=0; x<CTX_SZ; x++) {
					// upstream
					y = wcache_select(win, win_dp, off5p-CTX_SZ+x);
					if (y > 0)
						ctx5p[x][y]++;

					if (!hairpin)
						continue;

					// downstream
					y = wcache_select(win, win_dp, off3p+x);
					if (y > 0)
						ctx3p[CTX_SZ+x][y]++;
				}
			}
		}

		sam_itr_destroy(bat.bam_iter);
next_window:
		win_pos += WIN_SZ;
		if (win_pos > bat.bam_hdr->target_len[win_tid]) {
			win_pos = 0;
			win_tid++;
		}
	}

	fprintf(opt->metrics_fp, "#CTX5p\tpos");
	for (x=0; x<4; x++) {
		for (y=0; y<4; y++)
			fprintf(opt->metrics_fp, "\t%c/%c", "ACGT"[x], "ACGT"[y]);
	}
	fprintf(opt->metrics_fp, "\n");

	for (x=0; x<2*CTX_SZ; x++) {
		fprintf(opt->metrics_fp, "CTX5p\t%d", x-CTX_SZ);
		for (y=0; y<16; y++) {
			fprintf(opt->metrics_fp, "\t%jd", (uintmax_t)ctx5p[x][y]);
		}
		fprintf(opt->metrics_fp, "\n");
	}


	fprintf(opt->metrics_fp, "\n");
	fprintf(opt->metrics_fp, "#CTX3p\tpos");
	for (x=0; x<4; x++) {
		for (y=0; y<4; y++)
			fprintf(opt->metrics_fp, "\t%c/%c", "ACGT"[x], "ACGT"[y]);
	}
	fprintf(opt->metrics_fp, "\n");

	for (x=0; x<2*CTX_SZ; x++) {
		fprintf(opt->metrics_fp, "CTX3p\t%d", x-CTX_SZ);
		for (y=0; y<16; y++) {
			fprintf(opt->metrics_fp, "\t%jd", (uintmax_t)ctx3p[x][y]);
		}
		fprintf(opt->metrics_fp, "\n");
	}

	ret = 0;

	bam_destroy1(b);
	hts_idx_destroy(bat.bam_idx);
err3:
	bam_hdr_destroy(bat.bam_hdr);
err2:
	sam_close(bat.bam_fp);
err1:
	free(win);
err0:
	return ret;
}

void
usage(char *argv0)
{
	fprintf(stderr, "usage: %s [-h HPLEN] in.bam\n", argv0);
	exit(1);
}

int
main(int argc, char **argv)
{
	opt_t opt;
	int c;

	memset(&opt, 0, sizeof(opt_t));
	opt.min_mapq = 25;
	opt.min_baseq = 10;
	opt.hp_len = 29;
	opt.metrics_fp = stdout;

	while ((c = getopt(argc, argv, "h:")) != -1) {
		switch (c) {
			case 'h':
				opt.hp_len = strtoul(optarg, NULL, 0);
				if (opt.hp_len < 1 || opt.hp_len > 150) {
					fprintf(stderr, "Error: invalid hairpin length: `%s'\n", optarg);
					usage(argv[0]);
				}
				break;
			default:
				usage(argv[0]);
		}
	}

	if (argc-optind != 1) {
		usage(argv[0]);
	}

	opt.bam_fn = argv[optind];

	srand(31415);

	if (scan_pairs(&opt) < 0) {
		return 1;
	}

	return 0;
}

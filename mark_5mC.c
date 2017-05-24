/*
 * Print per-site methylation counts.
 *
 * Copyright (c) 2016,2017 Graham Gower <graham.gower@gmail.com>
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
#include <ctype.h>

#include <htslib/sam.h>
#include <htslib/faidx.h>

#include "kseq.h"
KSTREAM_INIT(int, read, 16384);

#include "fold.h"

typedef struct {
	char *bam_fn;
	char *bed_fn;
	char *fasta_fn;
	char *hairpin;
	char *rhairpin;
	size_t hlen;
	int min_mapq;
	int min_baseq;
	int min_5;
	int min_3;
} opt_t;

typedef struct _bed_entry {
	uint32_t tid;
	uint32_t start;
	uint32_t end;
	struct _bed_entry *next;
} bed_t;

typedef struct {
	opt_t *opt;
	samFile *bam_fp;
	bam_hdr_t *bam_hdr;
	hts_idx_t *bam_idx;
	hts_itr_t *bam_iter;
	bed_t *bed_head;
	char *ref;
	faidx_t *fai;
} bam_aux_t;

// defined in the SAM spec
#define PHRED_SCALE 33

static int nt2int[256] = {['A']=0, ['C']=1, ['G']=2, ['T']=3};
static char cmap[] = {['A']='T', ['C']='G', ['G']='C', ['T']='A', ['N']='N', ['n']='N',
				['a']='t', ['c']='g', ['g']='c', ['t']='a'};

/*
 * Load .bed loci into a singly linked list. Takes the chromosome name from
 * the first column, matching it to one from the bam header
 */
bed_t *
load_bed(const char *bed_fn, bam_hdr_t *bam_hdr)
{
	int fd;
	int dret;
	int lineno = 0;
	kstream_t *ks;
	kstring_t *str;
	bed_t *head, *tail, *tmp;

	head = tail = tmp = NULL;

	fd = open(bed_fn, O_RDONLY);
	if (fd == -1) {
		fprintf(stderr, "open: %s: %s\n", bed_fn, strerror(errno));
		goto err0;
	}

	str = calloc(1, sizeof(kstring_t));
	if (str == NULL) {
		perror("load_bed: calloc");
		goto err1;
	}

	ks = ks_init(fd);

	while (ks_getuntil(ks, KS_SEP_SPACE, str, &dret) >= 0) {
		int tid, start, end;

		lineno++;

		tid = bam_name2id(bam_hdr, str->s);
		if (tid < 0) {
			fprintf(stderr, "%s: line %d: couldn't find RNAME `%s' in sam header\n",
					bed_fn, lineno, str->s);
			goto err2;
		}

		if (dret == '\n') {
			fprintf(stderr, "%s: line %d: expected 3 columns, only got 1\n",
					bed_fn, lineno);
			goto err2;
		}

		if (ks_getuntil(ks, KS_SEP_SPACE, str, &dret) < 0)
			break;

		if (dret == '\n') {
			fprintf(stderr, "%s: line %d: expected 3 columns, only got 2\n",
					bed_fn, lineno);
			goto err2;
		}

		start = atoi(str->s);
		if (ks_getuntil(ks, KS_SEP_SPACE, str, &dret) < 0)
			break;
		end = atoi(str->s);

		// got 3 columns
		tmp = malloc(sizeof(bed_t));
		if (tmp == NULL) {
			perror("load_bed: malloc");
			goto err2;
		}

		tmp->tid = tid;
		tmp->start = start;
		tmp->end = end;
		tmp->next = NULL;

		if (tail == NULL) {
			head = tail = tmp;
		} else {
			tail->next = tmp;
			tail = tmp;
		}

		if (dret != '\n') {
			// consume the rest of the line
			while ((dret = ks_getc(ks)) > 0 && dret != '\n')
				;
		}
		if (dret < 0)
			break;
	}

	goto done;
err2:
	for (tmp=head; tmp!=NULL; ) {
		bed_t *x = tmp;
		tmp = tmp->next;
		free(x);
	}
	head = NULL;
done:
	ks_destroy(ks);
	free(str->s);
	free(str);
err1:
	close(fd);
err0:
	return head;
}


/*
 * Get next alignment.
 */
static int
next_aln(void *data, bam1_t *b)
{
	bam_aux_t *bat = data;
	int ret;

	uint8_t *xf_aux;
	char *xf;
	char *s1, *s2, *q1, *q2;
	int len;

	while (1) {
		if (bat->opt->bed_fn) {
			if (bat->bed_head == NULL) {
				ret = -1;
				break;
			}
			if (bat->bam_iter == NULL) {
				bed_t *x = bat->bed_head;
				bat->bed_head = x->next;
				bat->bam_iter = sam_itr_queryi(bat->bam_idx, x->tid, x->start, x->end);
				free(x);
				if (bat->bam_iter == NULL) {
					ret = -2;
					break;
				}
			}
			ret = sam_itr_next(bat->bam_fp, bat->bam_iter, b);
			if (ret < 0) {
				// iterator exhausted
				hts_itr_destroy(bat->bam_iter);
				bat->bam_iter = NULL;
				continue;
			}
		} else {
			ret = sam_read1(bat->bam_fp, bat->bam_hdr, b);
			if (ret < 0)
				break;
		}

		if (b->core.flag & (BAM_FUNMAP | BAM_FQCFAIL | BAM_FDUP) ||
				b->core.qual < bat->opt->min_mapq)
			// skip these
			continue;

		xf_aux = bam_aux_get(b, "XF");
		if (xf_aux == NULL)
			continue;
		xf = bam_aux2Z(xf_aux);
		if (xf == NULL)
			continue;

		len = xf2ssqq(xf, &s1, &s2, &q1, &q2);
		if (len < 0)
			continue;

		if (correct_s1s2(s1, q1, len, s2, q2, len,
				bat->opt->hairpin,
				bat->opt->rhairpin,
				bat->opt->hlen,
				PHRED_SCALE, PHRED_SCALE) == -1)
			continue;

		break;
	}

	return ret;
}

/*
 * Get base in the reference.
 */
char
get_refbase(bam_aux_t *bat, int tid, int pos)
{
	static int ref_tid = -1;
	static int ref_len = -1;

	if (ref_tid != tid) {
		if (bat->ref)
			free(bat->ref);
		bat->ref = NULL;
		ref_len = faidx_seq_len(bat->fai, bat->bam_hdr->target_name[tid]);
		if (ref_len == -1) {
			fprintf(stderr, "%s has region `%s', which is not in %s\n",
					bat->opt->bam_fn,
					bat->bam_hdr->target_name[tid],
					bat->opt->fasta_fn);
			return -1;
		}
		bat->ref = faidx_fetch_seq(bat->fai, bat->bam_hdr->target_name[tid], 0, ref_len, &ref_len);
		if (bat->ref == NULL)
			return -1;
		ref_tid = tid;
	}

	if (pos < 0 || pos >= ref_len)
		return 'N';

	return toupper(bat->ref[pos]);
}

int
print_m(bam_aux_t *bat, int tid, int pos, char refbase, int depth, int C, int mC)
{
	char strand, ctx1, ctx2;

	if (refbase == 'C') {
		strand = '+';
		ctx1 = get_refbase(bat, tid, pos+1);
		ctx2 = get_refbase(bat, tid, pos+2);
		if (ctx1 == -1 || ctx2 == -1)
			return -1;
	} else {
		strand = '-';
		ctx1 = get_refbase(bat, tid, pos-1);
		ctx2 = get_refbase(bat, tid, pos-2);
		if (ctx1 == -1 || ctx2 == -1)
			return -1;
		ctx1 = cmap[(int)ctx1];
		ctx2 = cmap[(int)ctx2];
	}

	printf("%s\t%d\t%d\t%c\t%d\t%d\t%d\t%c%c%c\n",
			bat->bam_hdr->target_name[tid],
			pos,
			pos+1,
			strand,
			depth,
			C,
			mC,
			'C',
			ctx1,
			ctx2);
	return 0;
}

int
mark_5mC(opt_t *opt)
{
	bam_aux_t bat;
	bam_plp_t plpiter;
	const bam_pileup1_t *plp;
	int tid, pos, n;
	int i;
	int ret;

	memset(&bat, 0, sizeof(bam_aux_t));
	bat.opt = opt;

	bat.bam_fp = sam_open(opt->bam_fn, "r");
	if (bat.bam_fp == NULL) {
		fprintf(stderr, "bam_open: %s: %s\n", opt->bam_fn, strerror(errno));
		ret = -1;
		goto err0;
	}

	bat.bam_hdr = sam_hdr_read(bat.bam_fp);
	if (bat.bam_hdr == NULL) {
		fprintf(stderr, "%s: couldn't read header\n", opt->bam_fn);
		ret = -2;
		goto err1;
	}

	bat.bed_head = NULL;
	if (opt->bed_fn) {
		bat.bam_idx = sam_index_load(bat.bam_fp, opt->bam_fn);
		if (bat.bam_idx == NULL) {
			fprintf(stderr, "Couldn't load index for %s\n", opt->bam_fn);
			ret = -3;
			goto err2;
		}
		bat.bed_head = load_bed(opt->bed_fn, bat.bam_hdr);
		if (bat.bed_head == NULL) {
			ret = -4;
			goto err3;
		}
	}

	bat.fai = fai_load(opt->fasta_fn);
	if (bat.fai == NULL) {
		ret = -5;
		goto err4;
	}

	bat.bam_iter = NULL;
	plpiter = bam_plp_init(next_aln, &bat);

	printf("chrom\tpos-0\tpos-1\tstrand\tdepth\tC\tmC\tcontext\n");

	while ((plp = bam_plp_auto(plpiter, &tid, &pos, &n)) != 0) {

		if (bat.bam_iter && (pos < bat.bam_iter->beg || pos > bat.bam_iter->end))
			continue;

		int ntpair[16] = {0,};
		int depth = 0;
		char refbase = 0;

		for (i=0; i<n; i++) {
			const bam_pileup1_t *p = plp + i;
			uint8_t *xf_aux;
			char *xf;
			char *s1, *s2, *q1, *q2;
			int ci, cj;
			int qi, qj;

			if (p->is_del || p->is_refskip)
				continue;

			if (p->qpos < opt->min_5 || p->b->core.l_qseq-p->qpos < opt->min_3)
				continue;

			if (bam_get_qual(p->b)[p->qpos] < opt->min_baseq)
				continue;

			xf_aux = bam_aux_get(p->b, "XF");
			if (xf_aux == NULL) {
				fprintf(stderr, "%s:%d: missing auxiliary field XF:Z\n", bam_get_qname(p->b), pos);
				ret = -6;
				goto err5;
			}
			xf = bam_aux2Z(xf_aux);
			if (xf == NULL) {
				fprintf(stderr, "%s:%d: invalid auxiliary field XF, not XF:Z\n", bam_get_qname(p->b), pos);
				ret = -7;
				goto err5;
			}

			if (xf2ssqq(xf, &s1, &s2, &q1, &q2) < 0) {
				fprintf(stderr, "%s:%d: invalid auxiliary field XF:Z, not created by foldreads\n", bam_get_qname(p->b), pos);
				ret = -8;
				goto err5;
			}

			int hclip = 0;
			if (bam_is_rev(p->b)) {
				if (p->b->core.n_cigar > 1) {
					if (bam_cigar_op(bam_get_cigar(p->b)[p->b->core.n_cigar-1]) == BAM_CHARD_CLIP)
						hclip = bam_cigar_oplen(bam_get_cigar(p->b)[p->b->core.n_cigar-1]);
				}
				ci = cmap[(int)s2[p->b->core.l_qseq - p->qpos-1 +hclip]];
				cj = s1[p->b->core.l_qseq - p->qpos-1 +hclip];
				qi = q2[p->b->core.l_qseq - p->qpos-1 +hclip];
				qj = q1[p->b->core.l_qseq - p->qpos-1 +hclip];
			} else {
				if (p->b->core.n_cigar > 1) {
					if (bam_cigar_op(bam_get_cigar(p->b)[0]) == BAM_CHARD_CLIP)
						hclip = bam_cigar_oplen(bam_get_cigar(p->b)[0]);
				}
				ci = s1[p->qpos+hclip];
				cj = cmap[(int)s2[p->qpos+hclip]];
				qi = q1[p->qpos+hclip];
				qj = q2[p->qpos+hclip];
			}

			/*
			if (pos == 1524196) {
				fprintf(stderr, "[%c] %s  %d:%d S=%c:Q=%d, S=%c/%c Q=%d/%d %s[%d]\n",
						"+-"[bam_is_rev(p->b)],
						bam_get_qname(p->b),
						pos,
						p->qpos,
						"=ACMGRSVTWYHKDBN"[bam_seqi(bam_get_seq(p->b),p->qpos)],
						bam_get_qual(p->b)[p->qpos],
						ci, cj,
						qi, qj,
						hclip>0?"HCLIP":"", hclip);
			}*/

			if (qi < PHRED_SCALE+opt->min_baseq || qj < PHRED_SCALE+opt->min_baseq)
				continue;

			int pair = nt2int[(int)ci]<<2 | nt2int[(int)cj];
			ntpair[pair]++;
			depth++;
		}

		if (!depth)
			continue;

		refbase = get_refbase(&bat, tid, pos);
		if (refbase == -1) {
			ret = -9;
			goto err5;
		}

		if (refbase == 'C') {
			int f_C = ntpair[nt2int['T']<<2 | nt2int['G']];
			int f_mC = ntpair[nt2int['C']<<2 | nt2int['G']];
			ret = print_m(&bat, tid, pos, refbase, depth, f_C, f_mC);
		} else if (refbase == 'G') {
			int r_C = ntpair[nt2int['G']<<2 | nt2int['T']];
			int r_mC = ntpair[nt2int['G']<<2 | nt2int['C']];
			ret = print_m(&bat, tid, pos, refbase, depth, r_C, r_mC);
		}

		if (ret == -1) {
			ret = -9;
			goto err5;
		}
	}


	ret = 0;
err5:
	if (bat.ref)
		free(bat.ref);

	if (opt->fasta_fn)
		fai_destroy(bat.fai);
	bam_plp_destroy(plpiter);
err4:
	if (bat.bed_head) {
		bed_t *x, *xx;
		for (x=bat.bed_head; x!=NULL; ) {
			xx = x;
			x = x->next;
			free(xx);
		}
	}
err3:
	if (bat.bam_idx)
		hts_idx_destroy(bat.bam_idx);
err2:
	bam_hdr_destroy(bat.bam_hdr);
err1:
	sam_close(bat.bam_fp);
err0:
	return ret;
}

void
usage(char *argv0, opt_t *opt)
{
	fprintf(stderr, "mark_5mC v9\n\n");
	fprintf(stderr, " Print the methylation status of cytosines in CpG/CHG/CHH contexts,\n"
			" where context is determined from the reference sequence.\n\n");

	fprintf(stderr, "usage: %s [...] in.bam ref.fasta\n", argv0);
	fprintf(stderr, " -p SEQ            The hairpin SEQuence\n");
	fprintf(stderr, " -b REGIONS.BED    Count methylation levels for specified regions\n");
	fprintf(stderr, " -M MAPQ           Minimum mapping quality for a read to be counted [%d]\n", opt->min_mapq);
	fprintf(stderr, " -B BASEQ          Minimum base quality for a base to be counted [%d]\n", opt->min_baseq);
	fprintf(stderr, " -5 N              Only count bases at least N bp from the 5' end of a read [%d]\n", opt->min_5);
	fprintf(stderr, " -3 M              Only count bases at least M bp from the 3' end of a read [%d]\n", opt->min_3);
	exit(1);
}

int
main(int argc, char **argv)
{
	opt_t opt;
	int c;
	int i;
	int ret;

	memset(&opt, 0, sizeof(opt_t));
	opt.min_mapq = 25;
	opt.min_baseq = 10;
	opt.min_5 = 0;
	opt.min_3 = 0;

	while ((c = getopt(argc, argv, "p:b:M:B:5:3:")) != -1) {
		switch (c) {
			case 'p':
				opt.hairpin = optarg;
				break;
			case 'b':
				opt.bed_fn = optarg;
				break;
			case 'M':
				opt.min_mapq = strtoul(optarg, NULL, 0);
				if (opt.min_mapq < 0 || opt.min_mapq > 100) {
					fprintf(stderr, "-M ``%s'' is invalid\n", optarg);
					usage(argv[0], &opt);
				}
				break;
			case 'B':
				opt.min_baseq = strtoul(optarg, NULL, 0);
				if (opt.min_baseq < 0 || opt.min_baseq > 100) {
					fprintf(stderr, "-B ``%s'' is invalid\n", optarg);
					usage(argv[0], &opt);
				}
				break;
			case '5':
				opt.min_5 = strtoul(optarg, NULL, 0);
				if (opt.min_5 < 0 || opt.min_5 > 10000) {
					fprintf(stderr, "-5 ``%s'' is invalid\n", optarg);
					usage(argv[0], &opt);
				}
				break;
			case '3':
				opt.min_3 = strtoul(optarg, NULL, 0);
				if (opt.min_3 < 0 || opt.min_3 > 10000) {
					fprintf(stderr, "-3 ``%s'' is invalid\n", optarg);
					usage(argv[0], &opt);
				}
				break;
			default:
				usage(argv[0], &opt);
		}
	}

	if (opt.hairpin == NULL) {
		fprintf(stderr, "Error: must specify a hairpin sequence.\n");
		usage(argv[0], &opt);
	}

	opt.hlen = strlen(opt.hairpin);
	if (opt.hlen < 5 || opt.hlen > 1000) {
		fprintf(stderr, "Error: hairpin too %s (len=%zd).\n",
				opt.hlen<5?"short":"long", opt.hlen);
		usage(argv[0], &opt);
	}

	for (i=0; i<opt.hlen; i++)
		opt.hairpin[i] = toupper(opt.hairpin[i]);
	opt.rhairpin = strdup(opt.hairpin);
	revcomp(opt.rhairpin, opt.hlen);

	if (argc-optind != 2) {
		usage(argv[0], &opt);
	}

	opt.bam_fn = argv[optind];
	opt.fasta_fn = argv[optind+1];

	ret = (mark_5mC(&opt) != 0);

	free(opt.rhairpin);

	return ret;
}

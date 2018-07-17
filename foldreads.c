/*
 * Fold r1/r2 sequences back together at the hairpin.
 *
 * Copyright (c) 2016-2018 Graham Gower <graham.gower@gmail.com>
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
#include <ctype.h>
#include <errno.h>
#include <getopt.h>
#include <stdint.h>
#include <math.h>
#include <assert.h>
#include <zlib.h>
#include "khash.h"

#include "kseq.h"
KSEQ_INIT(gzFile, gzread);

#include "fold.h"

#define min(a,b) ((a)<(b)?(a):(b))
#define max(a,b) ((a)>(b)?(a):(b))

int fit_lognorm(const uint64_t *hist, int len, uint64_t area, double *mu, double *sigma);

struct adapter {
	char *s;
	size_t l;
	double *pv;
};

struct hairpin {
	struct adapter *fwd, *rev;
};

typedef struct {
	char *fn1, *fn2; // Input filenames.
	char *metrics_fn;
	FILE *fos; // File pointer for outputting sequences.
	FILE *f_unmatched_r1, *f_unmatched_r2; // File pointers for unfolded reads
	int phred_scale_in, phred_scale_out;
	size_t adapter_matchlen;

	int n_hairpins;
	struct hairpin *hplist;
	struct adapter *a1, *a2;
} opt_t;

typedef struct {
	uint64_t total_reads;
	uint64_t folded;
	uint64_t hairpin_complete;
	uint64_t hairpin_dislocated;

	// fragment length histogram
	uint64_t *fl_hist;
	size_t fl_hist_len;

	// hairpin counts
	uint64_t *hp_hist;
} metrics_t;

/*
 * Fold the input sequences, s1 and s2.
 * Folded sequence is placed in s_out, folded quality scores placed in q_out.
 * Returns sequence length if folding was successful, 0 otherwise.
 */
int
fold(const opt_t *opt, metrics_t *metrics,
		const char *s1, const char *q1, size_t len1,
		const char *s2, const char *q2, size_t len2,
		char **_s_out, char **_q_out, int *hindex)
{
	struct hairpin *h = NULL;
	char *s_out, *q_out;
	double mm;
	int i;

	// hairpin indices
	int h1, h2;

	if (len1 != len2) {
		fprintf(stderr, "Error: input r1/r2 reads have different lengths.\n");
		exit(1);
	}

	s_out = malloc(len1+1);
	q_out = malloc(len1+1);

	if (s_out == NULL || q_out == NULL) {
		perror("malloc");
		exit(-2);
	}


	*hindex = -1;

	for (i=0; i<opt->n_hairpins; i++) {
		h = opt->hplist + i;

		// look for hairpin
		find_adapters(s1, q1, len1, s2, q2, len2,
				h->fwd->pv, h->fwd->l,
				h->rev->pv, h->rev->l, &h1, &h2);

		if (h1 != h2) {
			// these reads are untrustworthy
			if (h1+opt->adapter_matchlen < len1 && h2+opt->adapter_matchlen < len2) {
				// polymerase slippage? chimera?
				metrics->hairpin_dislocated++;
			}
			goto discard_reads;
		}

		if (h1+opt->adapter_matchlen < len1) {
			// found the hairpin
			*hindex = i;
			if (h1+h->fwd->l <= len1) // && h2+h->rev->l <= len2)
				metrics->hairpin_complete++;
			break;
		}
	}

	if (*hindex == -1) {
		/*
		 * No hairpin sequences found, so look for trailing adapters.
		 * If found, this means the molecule lacks a hairpin.
		 */
		int a1, a2;
		find_adapters(s1, q1, len1, s2, q2, len2,
				opt->a1->pv, opt->a1->l,
				opt->a2->pv, opt->a2->l, &a1, &a2);
		if (a1+opt->adapter_matchlen < len1 || a2+opt->adapter_matchlen < len2)
			goto discard_reads;
	}

	if (*hindex != -1 && h1+h->fwd->l < len1) {// && h2+h->rev->l < len2) {
		// short molecule, there are valid bases after the hairpin
		const char *s3 = s1 + h1 + h->fwd->l;
		const char *q3 = q1 + h1 + h->fwd->l;
		const char *s4 = s2 + h2 + h->fwd->l;
		const char *q4 = q2 + h2 + h->fwd->l;
		int len3 = 2*h1+h->fwd->l > len1 ? len1 - h1 - h->fwd->l : h1;
		int len4 = 2*h2+h->rev->l > len2 ? len2 - h2 - h->rev->l : h2;

		len1 = h1;
		len2 = h2;

		match4(s1, q1, len1,
				s2, q2, len2,
				s3, q3, len3,
				s4, q4, len4,
				s_out, q_out);
	} else {
		// long molecule, just match up to the hairpin
		if (*hindex != -1) {
			len1 = h1;
			len2 = h2;
		}
		match2(s1, q1, len1,
				s2, q2, len2,
				s_out, q_out,
				1);
	}

	size_t len = min(len1, len2);

	s_out[len] = '\0';
	q_out[len] = '\0';

	mm = posterior_error(q_out, len);

	if (mm < maxdiff(len, AVG_ERR, MAXDIFF_THRES)) {
		metrics->folded++;
		*_s_out = s_out;
		*_q_out = q_out;
		return len;
	} else {
discard_reads:
		free(s_out);
		free(q_out);
		*_s_out = *_q_out = NULL;
		*hindex = -1;
		return 0;
	}
}

/*
 * Convert phred scaled Q score to ASCII and print to file.
 */
void
fput_qual(FILE *fp, int phred_scale_out, const char *q, size_t len)
{
	int i;
	for (i=0; i<len; i++)
		fputc(q[i]+phred_scale_out, fp);
}

/*
 */
int
foldreads_pe(const opt_t *opt, metrics_t *metrics)
{
	gzFile fp1, fp2;
	kseq_t *seq1, *seq2;
	int len1, len2;
	int ret, fraglen;
	char *s_out, *q_out;
	int hindex;

	metrics->hp_hist = calloc(opt->n_hairpins, sizeof(*metrics->hp_hist));
	if (metrics->hp_hist == NULL) {
		perror("calloc");
		ret = -1;
		goto err0;
	}

	fp1 = gzopen(opt->fn1, "r");
	if (fp1 == NULL) {
		fprintf(stderr, "%s: %s\n", opt->fn1, strerror(errno));
		ret = -2;
		goto err0;
	}

	fp2 = gzopen(opt->fn2, "r");
	if (fp2 == NULL) {
		fprintf(stderr, "%s: %s\n", opt->fn2, strerror(errno));
		ret = -3;
		goto err1;
	}

	seq1 = kseq_init(fp1);
	seq2 = kseq_init(fp2);

	for (;;) {
		len1 = kseq_read(seq1);
		len2 = kseq_read(seq2);
		if (len1 < 0 || len2 < 0)
			break;

		if ((seq1->name.s[seq1->name.l-2] == '/' && seq1->name.s[seq1->name.l-1] == '1' &&
			seq2->name.s[seq2->name.l-2] == '/' && seq2->name.s[seq2->name.l-1] == '2')
		  || (seq1->name.s[seq1->name.l-2] == '.' && seq1->name.s[seq1->name.l-1] == '1' &&
			seq2->name.s[seq2->name.l-2] == '.' && seq2->name.s[seq2->name.l-1] == '2')) {
			seq1->name.l -= 2;
			seq1->name.s[seq1->name.l] = 0;
			seq2->name.l -= 2;
			seq2->name.s[seq2->name.l] = 0;
		}

		if (strcmp(seq1->name.s, seq2->name.s)) {
			fprintf(stderr, "R1 and R2 sequence names don't match:\n%s\n%s\n",
					seq1->name.s, seq2->name.s);
			ret = -4;
			goto err3;
		}

		metrics->total_reads++;

		//if (len1 < opt->hlen || len2 < opt->hlen)
		//	continue;

		if (seq1->qual.l == 0) {
			fprintf(stderr, "%s: qual scores required.\n", opt->fn1);
			ret = -5;
			goto err3;
		}

		if (seq2->qual.l == 0) {
			fprintf(stderr, "%s: qual scores required.\n", opt->fn2);
			ret = -6;
			goto err3;
		}

		if (metrics->fl_hist == NULL) {
			metrics->fl_hist_len = len1+1;
			metrics->fl_hist = calloc(metrics->fl_hist_len,
						sizeof(*metrics->fl_hist));
			if (metrics->fl_hist == NULL) {
				perror("malloc");
				ret = -7;
				goto err3;
			}
		} else if (metrics->fl_hist_len != len1+1) {
			fprintf(stderr, "Reads have different length to previous reads: %s\n",
					seq1->name.s);
			ret = -8;
			goto err3;
		}

		clean_quals(seq1->seq.s, seq1->qual.s, seq1->seq.l, opt->phred_scale_in);
		clean_quals(seq2->seq.s, seq2->qual.s, seq2->seq.l, opt->phred_scale_in);

		fraglen = fold(opt, metrics, seq1->seq.s, seq1->qual.s, seq1->seq.l,
			seq2->seq.s, seq2->qual.s, seq2->seq.l, &s_out, &q_out, &hindex);
		if (fraglen > 0) {
			fprintf(opt->fos, "@%s", seq1->name.s);

			fprintf(opt->fos, "\tr1:Z:%s", seq1->seq.s);
			fprintf(opt->fos, "\tq1:Z:");
			fput_qual(opt->fos, opt->phred_scale_out, seq1->qual.s, seq1->qual.l);

			fprintf(opt->fos, "\tr2:Z:%s", seq2->seq.s);
			fprintf(opt->fos, "\tq2:Z:");
			fput_qual(opt->fos, opt->phred_scale_out, seq2->qual.s, seq2->qual.l);

			if (hindex != -1) {
				struct hairpin *h = opt->hplist + hindex;
				fprintf(opt->fos, " hp:Z:%s", h->fwd->s);
				metrics->hp_hist[hindex]++;
			}

			fprintf(opt->fos, "\n%s\n+\n", s_out);
			fput_qual(opt->fos, opt->phred_scale_out, q_out, strlen(s_out));

			fputc('\n', opt->fos);

			free(s_out);
			free(q_out);

			fraglen = min(fraglen, len1);
			metrics->fl_hist[fraglen]++;

		} else if (opt->f_unmatched_r1) {

			if (seq1->comment.l)
				fprintf(opt->f_unmatched_r1, "@%s %s\n%s\n+\n", seq1->name.s, seq1->comment.s, seq1->seq.s);
			else
				fprintf(opt->f_unmatched_r1, "@%s\n%s\n+\n", seq1->name.s, seq1->seq.s);
			fput_qual(opt->f_unmatched_r1, opt->phred_scale_out, seq1->qual.s, seq1->qual.l);
			fputc('\n', opt->f_unmatched_r1);

			if (seq2->comment.l)
				fprintf(opt->f_unmatched_r2, "@%s %s\n%s\n+\n", seq2->name.s, seq2->comment.s, seq2->seq.s);
			else
				fprintf(opt->f_unmatched_r2, "@%s\n%s\n+\n", seq2->name.s, seq2->seq.s);
			fput_qual(opt->f_unmatched_r2, opt->phred_scale_out, seq2->qual.s, seq2->qual.l);
			fputc('\n', opt->f_unmatched_r2);
		}
	}

	if (len1 != -1) {
		fprintf(stderr, "%s: unexpected end of file.\n", opt->fn1);
		ret = -9;
		goto err3;
	}

	if (len2 != -1) {
		fprintf(stderr, "%s: unexpected end of file.\n", opt->fn2);
		ret = -10;
		goto err3;
	}

	ret = 0;

err3:
	kseq_destroy(seq2);
	kseq_destroy(seq1);
//err2:
	gzclose(fp2);
err1:
	gzclose(fp1);
err0:
	return ret;
}

void
print_metrics(const opt_t *opt, const metrics_t *metrics)
{
	FILE *fp = stderr;
	int i;

	if (opt->metrics_fn) {
		fp = fopen(opt->metrics_fn, "w");
		if (fp == NULL) {
			fprintf(stderr, "Error: %s: %s\n", opt->metrics_fn, strerror(errno));
			fp = stderr;
		}
	}

	fprintf(fp, "# foldreads version\n");
	fprintf(fp, "VV\t%s\n\n", FOLDREADS_VERSION);

	fprintf(fp, "# Number of read pairs processed.\n");
	fprintf(fp, "NR\t%jd\n\n", (uintmax_t)metrics->total_reads);

	fprintf(fp, "# Number of read pairs successfully folded.\n");
	fprintf(fp, "NF\t%jd\n\n", (uintmax_t)metrics->folded);

	fprintf(fp, "# Number of read pairs with complete hairpin.\n");
	fprintf(fp, "NH\t%jd\n\n", (uintmax_t)metrics->hairpin_complete);

	fprintf(fp, "# Number of read pairs with dislocated hairpin.\n");
	fprintf(fp, "ND\t%jd\n\n", (uintmax_t)metrics->hairpin_dislocated);

	fprintf(fp, "# Number of hairpins identified for each hairpin.\n");
	fprintf(fp, "#HP\tseq\tnum\n");
	for (i=0; i<opt->n_hairpins; i++) {
		struct hairpin *h = opt->hplist+i;
		fprintf(fp, "HP\t%s\t%jd\n", h->fwd->s, (uintmax_t)metrics->hp_hist[i]);
	}
	fprintf(fp, "\n");


	int histlen = metrics->fl_hist_len -1 -opt->adapter_matchlen;

	double mu, sigma;
	fit_lognorm(metrics->fl_hist, histlen, metrics->folded, &mu, &sigma);
	fprintf(fp, "# Best fit LogNormal distribution for fragment lengths.\n");
	fprintf(fp, "#LN\tmu\tsigma\n");
	fprintf(fp, "LN\t%.3lf\t%.3lf\n\n", mu, sigma);

	fprintf(fp, "# Fragment length counts.\n");
	fprintf(fp, "# The last line is for fragments greater than or equal to that length.\n");
	fprintf(fp, "#FL\tfraglen\tnseqs\n");
	for (i=1; i<histlen; i++) {
		fprintf(fp, "FL\t%d\t%jd\n",
				i, (uintmax_t)metrics->fl_hist[i]);
	}

	uintmax_t longmol = 0;
	for (; i<metrics->fl_hist_len; i++)
		longmol += metrics->fl_hist[i];

	fprintf(fp, "FL\t%d\t%jd\n",
			histlen,
			longmol);
	//fprintf(fp, "\n");

	free(metrics->fl_hist);
	free(metrics->hp_hist);

	if (fp != stderr)
		fclose(fp);
}

struct adapter *
adapter_init(const char *s)
{
	struct adapter *a;
	int i;

	a = calloc(sizeof(*a), 1);
	if (a == NULL) {
		perror("adapter_init:calloc");
		goto err0;
	}

	a->l = strlen(s);

	a->s = malloc(a->l+1);
	if (a->s == NULL) {
		perror("adapter_init:malloc");
		goto err1;
	}

	for (i=0; i<a->l; i++)
		a->s[i] = toupper(s[i]);

	a->s[a->l] = '\0';

	str2pvec(a->s, a->l, &a->pv);
	if (a->pv == NULL)
		goto err2;

	return a;
//err3:
//	free(a->pv);
err2:
	free(a->s);
err1:
	free(a);
err0:
	return NULL;
}

void
adapter_free(struct adapter *a)
{
	if (a == NULL)
		return;
	free(a->pv);
	free(a->s);
	free(a);
}

int
hairpin_init(struct hairpin *h, const char *s)
{
	char *s_rev;
	int ret;

	h->fwd = adapter_init(s);
	if (h->fwd == NULL) {
		ret = -1;
		goto err0;
	}

	s_rev = strdup(s);
	revcomp(s_rev, h->fwd->l);

	h->rev = adapter_init(s_rev);
	free(s_rev);
	if (h->fwd == NULL) {
		ret = -2;
		goto err1;
	}

	return 0;
//err2:
//	adapter_free(h->rev);
err1:
	adapter_free(h->fwd);
err0:
	return ret;
}

void
hairpin_free(struct hairpin *h)
{
	adapter_free(h->rev);
	adapter_free(h->fwd);
}

void
usage(char *argv0)
{
	fprintf(stderr, "foldreads v%s\n", FOLDREADS_VERSION);
	fprintf(stderr, "usage: %s [...] -1 IN1.FQ -2 IN2.FQ\n", argv0);
	fprintf(stderr, " -o OUT.FQ         Fastq output file [stdout]\n");
	fprintf(stderr, " -m FILE           Metrics output file [stderr]\n");
	fprintf(stderr, " -u PREFIX         Filename prefix for unfolded reads []\n");
	fprintf(stderr, " -p SEQ            The hairpin SEQuence [ACGCCGGCGGCAAGTGAAGCCGCCGGCGT]\n");
	fprintf(stderr, " -1 IN1.FQ[.GZ]    R1 fastq input file\n");
	fprintf(stderr, " -2 IN2.FQ[.GZ]    R2 fastq input file\n");
	fprintf(stderr, " -T SEQ            Adapter SEQuence trailing R1 (p7) []\n");
	fprintf(stderr, " -B SEQ            Adapter SEQuence trailing R2 (p5) []\n");
	exit(-1);
}

int
main(int argc, char **argv)
{
	opt_t opt;
	metrics_t metrics;
	struct hairpin *htmp, *h = NULL;
	int c, ret;
	int i;
	char *unmatched_pfx = NULL;
	char *fos_fn = NULL;
	char *a1, *a2;

	// Y-adapters
	// p7, might read into this at end of r1
	a1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG";
	// p5, might read into this at end of r2
	a2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT";

	memset(&opt, '\0', sizeof(opt_t));
	memset(&metrics, '\0', sizeof(metrics_t));
	opt.fos = stdout;
	opt.phred_scale_in = 33;
	opt.phred_scale_out = 33;
	opt.adapter_matchlen = 9;

	while ((c = getopt(argc, argv, "o:m:p:u:1:2:T:B:")) != -1) {
		switch (c) {
			case 'o':
				fos_fn = optarg;
				break;
			case 'm':
				opt.metrics_fn = optarg;
				break;
			case 'p':
				if (strlen(optarg) < opt.adapter_matchlen) {
					fprintf(stderr, "Hairpin sequence `%s' too short for reliable identification\n", optarg);
					exit(1);
				}
				opt.n_hairpins++;
				htmp = realloc(h, sizeof(*h)*opt.n_hairpins);
				if (htmp == NULL) {
					free(h);
					perror("realloc");
					exit(1);
				}
				h = htmp;
				if (hairpin_init(h + (opt.n_hairpins-1), optarg) < 0)
					exit(1);
				break;
			case 'u':
				unmatched_pfx = optarg;
				break;
			case '1':
				opt.fn1 = optarg;
				break;
			case '2':
				opt.fn2 = optarg;
				break;
			case 'T':
				a1 = optarg;
				break;
			case 'B':
				a2 = optarg;
				break;
			default:
				usage(argv[0]);
		}
	}

	if (opt.fn1 == NULL || opt.fn2 == NULL) {
		fprintf(stderr, "Error: must specify input files.\n");
		usage(argv[0]);
	}

	if (h == NULL) {
		opt.n_hairpins = 3;
		h = malloc(sizeof(*h)*opt.n_hairpins);
		if (h == NULL) {
			perror("malloc");
			exit(1);
		}

		// default hairpin sequence, central T is biotinylated
		if (hairpin_init(h, "ACGCCGGCGGCAAGTGAAGCCGCCGGCGT") < 0)
			exit(1);

		// polymerase skipped a base after biotin
		if (hairpin_init(h+1, "ACGCCGGCGGCAAGTAAGCCGCCGGCGT") < 0)
			exit(1);

		// polymerase skipped two bases after biotin
		if (hairpin_init(h+2, "ACGCCGGCGGCAAGTAGCCGCCGGCGT") < 0)
			exit(1);
	}

	opt.hplist = h;

	opt.a1 = adapter_init(a1);
	if (opt.a1 == NULL)
		exit(1);
	opt.a2 = adapter_init(a2);
	if (opt.a2 == NULL)
		exit(1);

	if (fos_fn) {
		opt.fos = fopen(fos_fn, "w");
		if (opt.fos == NULL) {
			fprintf(stderr, "Error: %s: %s\n", fos_fn, strerror(errno));
			usage(argv[0]);
		}
	}

	if (unmatched_pfx) {
		char *fn;
		fn = malloc(strlen(unmatched_pfx)+10);
		sprintf(fn, "%s_r1.fq", unmatched_pfx);
		opt.f_unmatched_r1 = fopen(fn, "w");
		if (opt.f_unmatched_r1 == NULL) {
			fprintf(stderr, "Error: %s: %s\n", fn, strerror(errno));
			usage(argv[0]);
		}
		sprintf(fn, "%s_r2.fq", unmatched_pfx);
		opt.f_unmatched_r2 = fopen(fn, "w");
		if (opt.f_unmatched_r2 == NULL) {
			fprintf(stderr, "Error: %s: %s\n", fn, strerror(errno));
			usage(argv[0]);
		}
		free(fn);
	}

	ret = foldreads_pe(&opt, &metrics);

	if (unmatched_pfx) {
		fclose(opt.f_unmatched_r1);
		fclose(opt.f_unmatched_r2);
	}

	if (fos_fn)
		fclose(opt.fos);

	adapter_free(opt.a1);
	adapter_free(opt.a2);

	print_metrics(&opt, &metrics);

	for (i=0; i<opt.n_hairpins; i++) {
		hairpin_free(h+i);
	}
	free(h);

	return ret;
}

/*
 * Fold r1/r2 sequences back together at the hairpin.
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

typedef struct {
	char *hairpin, *rhairpin; // Hairpin sequence, reverse complement.
	size_t hlen; // Hairpin length.
	char *fn1, *fn2; // Input filenames.
	char *metrics_fn;
	FILE *fos; // File pointer for outputting sequences.
	FILE *f_unmatched_r1, *f_unmatched_r2; // File pointers for unfolded reads
	int phred_scale_in, phred_scale_out;
	size_t adapter_matchlen;
} opt_t;

typedef struct {
	uint64_t total_reads;
	uint64_t folded;
	uint64_t hairpin_complete;
	uint64_t hairpin_dislocated;
} metrics_t;

/*
 * Fold the input sequences, s1 and s2.
 * Folded sequence is placed in s_out, folded quality scores placed in q_out.
 * Returns 1 if folding was successful, 0 otherwise.
 */
int
fold(const opt_t *opt, metrics_t *metrics,
		const char *s1, const char *q1, size_t len1,
		const char *s2, const char *q2, size_t len2,
		char **_s_out, char **_q_out)
{
	char *s_out, *q_out;
	int mm;

	// hairpin indices
	int h1, h2;

	if (len1 != len2) {
		fprintf(stderr, "Error: input r1/r2 reads have different lengths.\n");
		exit(1);
		/*
		*_s_out = *_q_out = NULL;
		return 0;
		*/
	}

	s_out = malloc(len1+1);
	q_out = malloc(len1+1);

	if (s_out == NULL || q_out == NULL) {
		perror("malloc");
		exit(-2);
	}

	find_hp_adapter(s1, len1, s2, len2,
			opt->hairpin, opt->rhairpin, opt->hlen,
			&h1, &h2);

	if (h1 && h2) {
		if (h1 != h2 && h1 < len1-opt->adapter_matchlen && h2 < len2-opt->adapter_matchlen) {
			// polymerase slippage? chimera?
			metrics->hairpin_dislocated++;
			//fprintf(stderr, "PCR SLIP: %d\n", h1-h2);
			goto discard_reads;
		}
	}
	if (h1 && !h2)
		h2 = h1;
	if (h2 && !h1)
		h1 = h2;

	h1 = h2 = min(h1, h2);

	if (h1+opt->hlen <= len1)// && h2+opt->hlen <= len2)
		metrics->hairpin_complete++;

	if (h1+opt->hlen < len1) {// && h2+opt->hlen < len2) {
		// short molecule, there are valid bases after the hairpin
		const char *s3 = s1 + h1 + opt->hlen;
		const char *q3 = q1 + h1 + opt->hlen;
		const char *s4 = s2 + h2 + opt->hlen;
		const char *q4 = q2 + h2 + opt->hlen;
		int len3 = 2*h1+opt->hlen > len1 ? len1 - h1 - opt->hlen : h1;
		int len4 = 2*h2+opt->hlen > len2 ? len2 - h2 - opt->hlen : h2;

		len1 = h1;
		len2 = h2;

		mm = match4(s1, q1, len1,
				s2, q2, len2,
				s3, q3, len3,
				s4, q4, len4,
				s_out, q_out,
				opt->phred_scale_in,
				opt->phred_scale_out);
	} else {
		// long molecule, just match up to the hairpin
		len1 = h1;
		len2 = h2;
		mm = match2(s1, q1, len1,
				s2, q2, len2,
				s_out, q_out,
				1,
				opt->phred_scale_in,
				opt->phred_scale_out);
	}

	s_out[min(len1,len2)] = '\0';
	q_out[min(len1,len2)] = '\0';

	if (mm <= maxdiff(len1+len2, AVG_ERR, MAXDIFF_THRES)) {
		metrics->folded++;
		*_s_out = s_out;
		*_q_out = q_out;
		return 1;
	} else {
discard_reads:
		free(s_out);
		free(q_out);
		*_s_out = *_q_out = NULL;
		return 0;
	}
}

/*
 * Write out a single fastq entry.
 */
void
seq_write(const char *name, const char *comment, const char *seq, const char *qual, FILE *fs)
{
	if (comment)
		fprintf(fs, "@%s %s\n%s\n+\n%s\n", name, comment, seq, qual);
	else
		fprintf(fs, "@%s\n%s\n+\n%s\n", name, seq, qual);
}

/*
 */
int
foldreads_pe(const opt_t *opt, metrics_t *metrics)
{
	gzFile fp1, fp2;
	kseq_t *seq1, *seq2;
	int len1, len2;
	int ret, f;
	char *s_out, *q_out;

	fp1 = gzopen(opt->fn1, "r");
	if (fp1 == NULL) {
		fprintf(stderr, "%s: %s\n", opt->fn1, strerror(errno));
		ret = 1;
		goto err0;
	}

	fp2 = gzopen(opt->fn2, "r");
	if (fp2 == NULL) {
		fprintf(stderr, "%s: %s\n", opt->fn2, strerror(errno));
		ret = 1;
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
			ret = 3;
			goto err3;
		}

		metrics->total_reads++;

		if (len1 < opt->hlen || len2 < opt->hlen)
			continue;

		if (seq1->qual.l == 0) {
			fprintf(stderr, "%s: qual scores required.\n", opt->fn1);
			ret = 4;
			goto err3;
		}

		if (seq2->qual.l == 0) {
			fprintf(stderr, "%s: qual scores required.\n", opt->fn2);
			ret = 4;
			goto err3;
		}

		f = fold(opt, metrics, seq1->seq.s, seq1->qual.s, seq1->seq.l,
			seq2->seq.s, seq2->qual.s, seq2->seq.l, &s_out, &q_out);
		if (f == 1) {
			fprintf(opt->fos,
					"@%s XF:Z:%s|%s|%s|%s\n%s\n+\n%s\n",
					seq1->name.s,
					seq1->seq.s, seq2->seq.s,
					seq1->qual.s, seq2->qual.s,
					s_out, q_out);
			free(s_out);
			free(q_out);
		} else if (opt->f_unmatched_r1) {
			char *comment = seq1->comment.l==0 ? NULL : seq1->comment.s;
			seq_write(seq1->name.s, comment, seq1->seq.s, seq1->qual.s, opt->f_unmatched_r1);
			comment = seq2->comment.l==0 ? NULL : seq2->comment.s;
			seq_write(seq2->name.s, comment, seq2->seq.s, seq2->qual.s, opt->f_unmatched_r2);
		}
	}

	if (len1 != -1) {
		fprintf(stderr, "%s: unexpected end of file.\n", opt->fn1);
		ret = 5;
		goto err3;
	}

	if (len2 != -1) {
		fprintf(stderr, "%s: unexpected end of file.\n", opt->fn2);
		ret = 5;
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

	if (opt->metrics_fn) {
		fp = fopen(opt->metrics_fn, "w");
		if (fp == NULL) {
			fprintf(stderr, "Error: %s: %s\n", opt->metrics_fn, strerror(errno));
			fp = stderr;
		}
	}

	fprintf(fp, "Total read pairs: %jd\n", (uintmax_t)metrics->total_reads);
	fprintf(fp, "Number of successfully folded read pairs: %jd (%lf)\n",
			(uintmax_t)metrics->folded,
			(double)metrics->folded/metrics->total_reads);
	fprintf(fp, "Number of read pairs with complete hairpin: %jd\n",
			(uintmax_t)metrics->hairpin_complete);
	fprintf(fp, "Number of read pairs with dislocated hairpins: %jd\n",
			(uintmax_t)metrics->hairpin_dislocated);

	if (fp != stderr)
		fclose(fp);
}

void
usage(char *argv0)
{
	fprintf(stderr, "foldreads v7\n");
	fprintf(stderr, "usage: %s [...] -p SEQ -1 IN1.FQ -2 IN2.FQ\n", argv0);
	fprintf(stderr, " -o OUT.FQ         Fastq output file [stdout]\n");
	fprintf(stderr, " -m FILE           Metrics output file [stderr]\n");
	fprintf(stderr, " -u PREFIX         Filename prefix for unfolded reads []\n");
	fprintf(stderr, " -p SEQ            The hairpin SEQuence\n");
	fprintf(stderr, " -1 IN1.FQ[.GZ]    R1 fastq input file\n");
	fprintf(stderr, " -2 IN2.FQ[.GZ]    R2 fastq input file\n");
	exit(-1);
}

int
main(int argc, char **argv)
{
	opt_t opt;
	metrics_t metrics;
	int c, ret;
	int i;
	char *unmatched_pfx = NULL;
	char *fos_fn = NULL;

	memset(&opt, '\0', sizeof(opt_t));
	memset(&metrics, '\0', sizeof(metrics_t));
	opt.fos = stdout;
	opt.phred_scale_in = 33;
	opt.phred_scale_out = 33;
	opt.adapter_matchlen = 9;

	while ((c = getopt(argc, argv, "o:m:p:u:1:2:")) != -1) {
		switch (c) {
			case 'o':
				fos_fn = optarg;
				break;
			case 'm':
				opt.metrics_fn = optarg;
				break;
			case 'p':
				opt.hairpin = optarg;
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
			default:
				usage(argv[0]);
		}
	}

	if (opt.fn1 == NULL || opt.fn2 == NULL) {
		fprintf(stderr, "Error: must specify input files.\n");
		usage(argv[0]);
	}

	if (opt.hairpin == NULL) {
		fprintf(stderr, "Error: must specify a hairpin sequence.\n");
		usage(argv[0]);
	}

	opt.hlen = strlen(opt.hairpin);

	if (opt.hlen < 5 || opt.hlen > 1000) {
		fprintf(stderr, "Error: hairpin too %s (len=%zd).\n",
				opt.hlen<5?"short":"long", opt.hlen);
		usage(argv[0]);
	}

	for (i=0; i<opt.hlen; i++)
		opt.hairpin[i] = toupper(opt.hairpin[i]);

	opt.rhairpin = strdup(opt.hairpin);
	revcomp(opt.rhairpin, opt.hlen);

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

	free(opt.rhairpin);

	if (unmatched_pfx) {
		fclose(opt.f_unmatched_r1);
		fclose(opt.f_unmatched_r2);
	}

	if (fos_fn)
		fclose(opt.fos);

	print_metrics(&opt, &metrics);

	return ret;
}

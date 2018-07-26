/*
 * Create a quality score profile from a fastq file.
 *
 * Copyright (c) 2018 Graham Gower <graham.gower@gmail.com>
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
#include <stdint.h>
#include <errno.h>
#include <zlib.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread);

#define MAX_READS 100000
#define PHRED_SCALE 33

int
qualprofile(char *fn)
{
	int i, j;
	int ret;
	int len;
	int qlen = 0;
	uint64_t nseqs;
	gzFile fp;
	kseq_t *seq;

	double *mu = NULL; // MVN mean
	double *Sigma = NULL; // MVN covariance matrix

	fp = gzopen(fn, "r");
	if (fp == NULL) {
		fprintf(stderr, "%s: %s\n", fn, strerror(errno));
		ret = -1;
		goto err0;
	}

	seq = kseq_init(fp);
	nseqs = 0;

	while (nseqs<MAX_READS && (len = kseq_read(seq)) >= 0) {
		if (seq->qual.l == 0) {
			fprintf(stderr, "%s: qual scores required.\n", fn);
			ret = -2;
			goto err4;
		}

		if (nseqs++ == 0) {
			qlen = seq->qual.l;

			mu = calloc(qlen, sizeof(double));
			if (mu == NULL) {
				perror("calloc");
				ret = -3;
				goto err4;
			}
			Sigma = calloc(qlen*qlen, sizeof(double));
			if (Sigma == NULL) {
				perror("calloc");
				ret = -4;
				goto err4;
			}
		}
		
		if (qlen != seq->qual.l) {
			fprintf(stderr, "%s: read length changed at %s.\n", fn, seq->name.s);
			ret = -5;
			goto err4;
		}

		/*
		 * Qual scores ought to be correlated along a sequence,
		 * so model each read as having MVN distributed qual scores.
		 *
		 * NOTE: we really estimate mu and Sigma for a truncated MVN,
		 * and would like to recover the pre-trucated parameters in
		 * order that simulations (subsequently truncated also)
		 * reflect the correct empirical distribution.  The effect
		 * of not taking this into account is that our parameter
		 * estimates are biased low.  Possible solutions:
		 * https://stats.stackexchange.com/questions/58081/pre-truncation-moments-for-truncated-multivariate-normal
		 *
		 * Cumulative mean and covariance adapted from
		 * Knuth TAOCP vol 2, 3rd edition, p 232.
		 */
		for (i=0; i<qlen; i++) {
			double qi = seq->qual.s[i];
			mu[i] += (qi - mu[i]) / nseqs;
			for (j=0; j<=i; j++) {
				double qj = seq->qual.s[j];
				Sigma[i*qlen+j] += (qj - mu[j]) * (qi - mu[i]);
			}
		}
	}

	if (len < -1) {
		fprintf(stderr, "%s: unexpected end of file.\n", fn);
		ret = -6;
		goto err4;
	}

	if (nseqs == 0) {
		fprintf(stderr, "%s: no sequences found.\n", fn);
		ret = -7;
		goto err4;
	}

	printf("# Average phred scaled quality score for each position in the reads.\n");
	printf("MEAN");
	for (i=0; i<qlen; i++)
		printf(" %f", mu[i]-PHRED_SCALE);
	printf("\n\n");

	printf("# Covariance matrix (lower triangle) of the quality scores.\n");
	for (i=0; i<qlen; i++) {
		printf("COV");
		for (j=0; j<=i; j++)
			printf(" %f", Sigma[i*qlen+j] / (nseqs-1));
		printf("\n");
	}

	ret = 0;
err4:
	if (Sigma)
		free(Sigma);
//err3:
	if (mu)
		free(mu);
//err2:
	kseq_destroy(seq);
//err1:
	gzclose(fp);
err0:
	return ret;
}

int
main(int argc, char **argv)
{
	if (argc != 2) {
		fprintf(stderr, "usage: %s file.fq", argv[0]);
		return -1;
	}

	if (qualprofile(argv[1]) < 0)
		return -2;

	return 0;
}

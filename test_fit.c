/*
 * Test fit_lognorm.c with randomly drawn read lengths.
 *
 * Copyright (c) 2017 Graham Gower <graham.gower@gmail.com>
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
#include <stdint.h>
#include <sys/types.h>
#include <math.h>
#include "kmath.h"

int fit_lognorm(const uint64_t *hist, int len, uint64_t area, double *mu, double *sigma);

static krand_t *kr;

uint
rlognorm(double mu, double sigma)
{
	double x = kr_normal(kr);
	return round(exp(mu + x*sigma));
}

void
sample(uint64_t *hist, int len, size_t n,
		double mu, double sigma)
{
	int i;
	for (i=0; i<n; i++) {
		uint x = rlognorm(mu, sigma);
		if (x > len)
			x = len;
		hist[x]++;
	}
}

#define HIST_SIZE 141
#define NSAMPLES 100
#define SAMPLE_SIZE 10000

int
main()
{
	uint64_t hist[HIST_SIZE+1];
	double mu, sigma;
	double est_mu, est_sigma;
	double ssqe_mu, ssqe_sigma;
	uint64_t tsum;
	int i;

	kr = kr_srand(31415);

	for (mu=4.0; mu<5.55; mu+=0.1) {
		for (sigma=0.2; sigma<0.55; sigma+=0.1) {
			ssqe_mu = ssqe_sigma = 0;
			tsum = 0;
			for (i=0; i<NSAMPLES; i++) {
				memset(hist, 0, sizeof(hist));
				sample(hist, HIST_SIZE, SAMPLE_SIZE, mu, sigma);
				fit_lognorm(hist, HIST_SIZE, SAMPLE_SIZE, &est_mu, &est_sigma);
				ssqe_mu += (mu-est_mu)*(mu-est_mu);
				ssqe_sigma += (sigma-est_sigma)*(sigma-est_sigma);
				tsum += hist[HIST_SIZE];
			}
			printf("mu=%.3lf (msqe=%.3lg), sigma=%.3lf (msqe=%.3lg), tail=%.3lg\n",
					mu, ssqe_mu/NSAMPLES,
					sigma, ssqe_sigma/NSAMPLES,
					(double)tsum/NSAMPLES/SAMPLE_SIZE);
		}
	}

	free(kr);
}

/*
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
#include <sys/types.h>
#include <stdint.h>
#include <math.h>
#include <float.h>
#include "kmath.h"

// doesn't seem to work very well
#define USE_TAIL 0

typedef struct {
	double *p; // normalised histogram array
	double p_tail; // proportion that is right censored (longer than `len')
	size_t len; // length of p
	size_t min, max; // minimum and maximum indices for valid data
} hist_t;

/*
 * PDF of LogNormal distribution.
 */
static double
pdf_lognorm(int x, double mu, double sigma)
{
	if (x == 0)
		return 0;
	double a = erf((log(x+0.5) - mu) / (M_SQRT2*sigma));
	double b = erf((log(x-0.5) - mu) / (M_SQRT2*sigma));
	return 0.5*(a-b);
}

#if USE_TAIL
/*
 * 1-CDF of LogNormal distribution.
 */
static double
sf_lognorm(int x, double mu, double sigma)
{
	if (x == 0)
		return 1;
	double a = erf((log(x-0.5) - mu) / (M_SQRT2*sigma));
	return 0.5*(1.0-a);
}
#endif

/*
 * Sum of Squared Error between data and LogNormal.
 */
static double
ssqe_lognorm(int n, double *x, void *data)
{
	hist_t *hist = data;
	double mu = x[0];
	double sigma = x[1];
	double ssqe = 0;
	int i = 0;

	if (mu < 0 || sigma < 0)
		return DBL_MAX;

	for (i=hist->min; i<hist->max; i++) {
		double p_i = pdf_lognorm(i, mu, sigma);
		double diff = p_i - hist->p[i];
		ssqe += diff*diff;
		if (isnan(ssqe))
			return DBL_MAX;
	}

#if USE_TAIL
	// tail
	double p_tail = sf_lognorm(hist->len, mu, sigma);
	double diff = p_tail - hist->p_tail;
	ssqe += diff*diff;
	if (isnan(ssqe))
		return DBL_MAX;
#endif

	return ssqe;
}

/*
 * Fit histogram to a LogNormal.  The histogram may be right censored due
 * to the read length, so we fit by minimising the squared error for the
 * region of the histogram that was observed.
 */
int
fit_lognorm(const uint64_t *hist, int len, uint64_t area, double *mu, double *sigma)
{
	hist_t hist_normed;
	int ret;
	uint64_t tail_area;
	int i;

	*mu = *sigma = 0;
	hist_normed.min = -1;

	if (area == 0) {
		// no reads
		ret = -1;
		goto err0;
	}

	tail_area = area;
	for (i=0; i<len; i++) {
		tail_area -= i*hist[i];
		if (hist[i] > 0) {
			if (hist_normed.min == -1)
				hist_normed.min = i;
			hist_normed.max = i;
		}
	}

	if (tail_area < 0) {
		fprintf(stderr, "negative tail_area!\n");
		ret = -2;
		goto err0;
	}

	hist_normed.len = len;
	hist_normed.p = malloc(len * sizeof(*hist_normed.p));
	if (hist_normed.p == NULL) {
		ret = -3;
		goto err0;
	}
	// normalise hist by area
	hist_normed.p_tail = (double)tail_area / area;
	for (i=0; i<len; i++)
		hist_normed.p[i] = (double)hist[i]/area;

	double x[] = {5.0, 0.5}; // initial values for mu and sigma
	kmin_hj(ssqe_lognorm, 2, x, &hist_normed,
			KMIN_RADIUS, KMIN_EPS, KMIN_MAXCALL);

	*mu = x[0];
	*sigma = x[1];
	ret = 0;

	free(hist_normed.p);
err0:
	return ret;
}

/*
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
#include <string.h>
#include <stdlib.h>
#include <sys/types.h>
#include <ctype.h>
#include <assert.h>
#include <math.h>
#include "fold.h"

static char cmap[] = {['A']='T', ['C']='G', ['G']='C', ['T']='A', ['N']='N', ['n']='N',
				['a']='t', ['c']='g', ['g']='c', ['t']='a'};

#define min(a,b) ((a)<(b)?(a):(b))
#define max(a,b) ((a)>(b)?(a):(b))

static uint _n2i[] = {['A']=0, ['C']=1, ['G']=2, ['T']=3, ['N']=4};
#define n2i(nt) _n2i[(uint)nt]


/*
 * Inverse poisson CDF, stolen from bwa: bwtaln.c
 */
static int
bwa_cal_maxdiff(int l, double err, double thres)
{
	double elambda = exp(-l * err);
	double sum, y = 1.0;
	int k, x = 1;
	for (k = 1, sum = elambda; k < MAXLEN; ++k) {
		y *= l * err;
		x *= k;
		sum += elambda * y / x;
		if (1.0 - sum < thres) return k;
	}
	return 2;
}

int
maxdiff(int l, double err, double thres)
{
	static int maxdiff[MAXLEN];
	static int init = 0;

	if (!init) {
		int i;
		for (i=0; i<MAXLEN; i++)
			maxdiff[i] = bwa_cal_maxdiff(i, err, thres);
		init = 1;
	}

	return l<MAXLEN ? maxdiff[l] : maxdiff[MAXLEN-1];
}

/*
 * Reverse complement of s.
 */
void
revcomp(char *s, size_t len)
{
	int i, j;
	int tmp;

	for (i=0, j=len-1; i<len/2; i++, j--) {
		tmp = cmap[(int)s[i]];
		s[i] = cmap[(int)s[j]];
		s[j] = tmp;
	}

	if (len % 2 == 1)
		s[i] = cmap[(int)s[i]];
}

/*
 * Reverse s (no complement).
 */
void
reverse(char *s, size_t len)
{
	int i, j;
	int tmp;

	for (i=0, j=len-1; i<len/2; i++, j--) {
		tmp = s[i];
		s[i] = s[j];
		s[j] = tmp;
	}
}

/*
 * Covert PHRED scaled Q value to a probability.
 */
static inline double
q2p(int q)
{
	static int q2p_cache_inited = 0;
	static double q2p_cache[PHRED_MAX];

	if (q2p_cache_inited)
		return q2p_cache[q];

	int x;

	// probability of error for random nucleotide
	q2p_cache[0] = q2p_cache[1] = 0.75;

	for (x=2; x<PHRED_MAX; x++)
		q2p_cache[x] = pow(10, -x/10.0);

	q2p_cache_inited = 1;
	return q2p_cache[q];
}

static inline int
p2q(double p)
{
	if (p == 0)
		return Q_MAX;
	int q = -round(10.0*log10(p));
	if (q > Q_MAX)
		q = Q_MAX;
	return q;
}

/* 
 * Compare len bytes of s1 and s2, allowing mismatches.
 * Return summed mismatch probabilities, or -1 if Ns of the same length do
 * better.
 */
static double
mmcmp(const char *s, const char *q, const double *pv, size_t len)
{
	int i;
	double mm;

	// must do better than a string of Ns
	mm = (q2p(2)/3) * len;

	for (i=0; i<len; i++) {
		int x = i<<(PHRED_SHIFT+3) | q[i]<<(3) | n2i(s[i]);
		mm -= pv[x];
		if (mm < 0)
			return -1;
	}

	return mm;
}

void
str2pvec(const char *s, size_t len, double **pv)
{
	int i, x;
	int q, n;

	*pv = calloc(len*PHRED_MAX*8, sizeof(double));
	if (*pv == NULL) {
		perror("calloc");
		return;
	}

	for (i=0; i<len; i++) {
		for (q=0; q<PHRED_MAX; q++) {
			for (n=0; n<8; n++) {
				x = i<<(PHRED_SHIFT+3) | q<<(3) | n;

				if (s[i] == 'N')
					(*pv)[x] = 0.75;
				else if (n != n2i(s[i]))
					(*pv)[x] = 1-q2p(q);
				else
					(*pv)[x] = q2p(q)/3;
			}
		}
	}
}

/*
 * Search s1 and s2 for adapter sequences.
 */
void
find_adapters(const char *s1, const char *q1, size_t len1,
		const char *s2, const char *q2, size_t len2,
		const double *pv1, size_t pv1_len,
		const double *pv2, size_t pv2_len,
		int *h1, int *h2)
{
	int i;
	double mm1, mm2;
	double mm1_max = -1, mm2_max = -1;

	*h1 = len1;
	*h2 = len2;

	for (i=0; i<len1; i++) {
		if ((mm1=mmcmp(s1+i, q1+i, pv1, min(len1-i, pv1_len))) > 0) {
			if (mm1 > mm1_max) {
				mm1_max = mm1;
				*h1 = i;
			}
		}
	}

	for (i=0; i<len2; i++) {
		if ((mm2=mmcmp(s2+i, q2+i, pv2, min(len2-i, pv2_len))) > 0) {
			if (mm2 > mm2_max) {
				mm2_max = mm2;
				*h2 = i;
			}
		}
	}
}

/*
 * Expected number of wrong bases.
 */
double
posterior_error(const char *qvec, int len)
{
	int i;
	double sum_p = 0;
	for (i=0; i<len; i++)
		sum_p += q2p(qvec[i]);
	return sum_p;
}

/*
 * Match a single base pair using posterior base probability as
 * described in Renaud et al. 2014, with minor modifications to
 * allow bisulfite converted base pairs.
 */
static void
match1bp_slow(char c1, char c2, char q1, char q2,
		char *c_out, char *q_out,
		int allow_bs)
{
	double p[4], p1[4], p2[4];
	uint q[4], q_max;
	double p_sum;
	int i, i_max;
	int n_maxs;

	/*
	 * Pr(nt | N_obs=c1, Q_obs=q1)
	 *   Probability of nt given one stranded observation.
	 */
	p1[n2i('A')] = c1=='A'? 1-q2p(q1) : q2p(q1)/3;
	p1[n2i('C')] = c1=='C'? 1-q2p(q1) : q2p(q1)/3;
	p1[n2i('G')] = c1=='G'? 1-q2p(q1) : q2p(q1)/3;
	p1[n2i('T')] = c1=='T'? 1-q2p(q1) : q2p(q1)/3;

	p2[n2i('A')] = c2=='A'? 1-q2p(q2) : q2p(q2)/3;
	p2[n2i('C')] = c2=='C'? 1-q2p(q2) : q2p(q2)/3;
	p2[n2i('G')] = c2=='G'? 1-q2p(q2) : q2p(q2)/3;
	p2[n2i('T')] = c2=='T'? 1-q2p(q2) : q2p(q2)/3;

	if (allow_bs) {
		/*
		 * Pr(nt | N_obs={n1,n2}, Q_obs={q1,q2}, bisulfite treated)
		 *   Probability of nt given two stranded observation, and
		 *   we can have differences caused by bisulfite treatement.
		 */
		p[n2i('A')] = p1[n2i('A')] * p2[n2i('A')];
		p[n2i('T')] = p1[n2i('T')] * p2[n2i('T')];
		p[n2i('C')] = 0.5*(p1[n2i('C')]+p1[n2i('T')]) * p2[n2i('C')];
		p[n2i('G')] = p1[n2i('G')] * 0.5*(p2[n2i('G')]+p2[n2i('A')]);
	} else {
		/*
		 * Pr(nt | N_obs={n1,n2}, Q_obs={q1,q2})
		 *   Probability of nt given two stranded observation.
		 */
		p[n2i('A')] = p1[n2i('A')] * p2[n2i('A')];
		p[n2i('T')] = p1[n2i('T')] * p2[n2i('T')];
		p[n2i('C')] = p1[n2i('C')] * p2[n2i('C')];
		p[n2i('G')] = p1[n2i('G')] * p2[n2i('G')];
	}

	p_sum = p[n2i('A')] + p[n2i('C')] + p[n2i('G')] + p[n2i('T')];
	assert(p_sum != 0);

	q[n2i('A')] = p2q((p[n2i('C')] + p[n2i('G')] + p[n2i('T')])/p_sum);
	q[n2i('C')] = p2q((p[n2i('A')] + p[n2i('G')] + p[n2i('T')])/p_sum);
	q[n2i('G')] = p2q((p[n2i('A')] + p[n2i('C')] + p[n2i('T')])/p_sum);
	q[n2i('T')] = p2q((p[n2i('A')] + p[n2i('C')] + p[n2i('G')])/p_sum);

	q_max = 0, i_max = -1;
	n_maxs = 0;
	for (i=0; i<4; i++) {
		if (q[i] > q_max) {
			q_max = q[i];
			i_max = i;
			n_maxs = 1;
		} else if (q[i] == q_max)
			n_maxs++;
	}

	assert(n_maxs != 0);

	if (n_maxs > 1) {
		*c_out = 'N';
		*q_out = 0;
	} else {
		*c_out = "ACGT"[i_max];
		*q_out = q[i_max];
	}

	if (allow_bs) {
		// check if methylated
		if (*c_out == 'C') {
			if (p1[n2i('C')] > p1[n2i('T')])
				*c_out = 'c';
		} else if (*c_out == 'G') {
			if (p2[n2i('G')] > p2[n2i('A')])
				*c_out = 'g';
		}
	}
}

static char _nt_MAP_n1n2q1q2[PHRED_MAX*PHRED_MAX*4*4];
static char _nt_MAP_n1n2q1q2bs[PHRED_MAX*PHRED_MAX*4*4];
static char _Q_MAP_n1n2q1q2[PHRED_MAX*PHRED_MAX*4*4];
static char _Q_MAP_n1n2q1q2bs[PHRED_MAX*PHRED_MAX*4*4];

/*
 * Cache the most probable nucleotide, and its Q value, for all possible
 * two stranded observations.
 */
static void
init_match1bp_cache()
{
	int n1, n2, q1, q2;
	int x;
	char c, q;

	for (q1=0; q1<PHRED_MAX; q1++) {
		for (q2=0; q2<PHRED_MAX; q2++) {
			for (n1=0; n1<4; n1++) {
				for (n2=0; n2<4; n2++) {
					x = q1<<(PHRED_SHIFT+2+2) | q2<<(2+2) | n1<<(2) | n2;

					// no bs
					match1bp_slow("ACGT"[n1], "ACGT"[n2], q1, q2, &c, &q, 0);
					_nt_MAP_n1n2q1q2[x] = c;
					_Q_MAP_n1n2q1q2[x] = q;

					// bs
					match1bp_slow("ACGT"[n1], "ACGT"[n2], q1, q2, &c, &q, 1);
					_nt_MAP_n1n2q1q2bs[x] = c;
					_Q_MAP_n1n2q1q2bs[x] = q;
				}
			}
		}
	}
}

/*
 * Match a single base pair.
 */
static void
match1bp(char c1, char c2, char q1, char q2,
		char *c_out, char *q_out,
		int allow_bs)
{
	int x;
	static int cached_inited = 0;

	if (!cached_inited) {
		init_match1bp_cache();
		cached_inited = 1;
	}

	c1 = toupper(c1);
	c2 = toupper(c2);

	if (c1 == 'N' || c2 == 'N') {
		if (allow_bs) {
			// not much can be said with confidence from a
			// single stranded observation

			switch (c1) {
				case 'A':
					// have an A
					*c_out = c1;
					*q_out = q1;
					break;
				case 'C':
					// must be 5mC
					*c_out = 'c';
					*q_out = q1;
					break;
				case 'G':
					// might have C or 5mC on other strand
					*c_out = 'N';
					break;
				case 'T':
					// might have a T or C
					*c_out = 'N';
					break;
			}
			switch (c2) {
				case 'A':
					// might be an A or a G
					*c_out = 'N';
					break;
				case 'C':
					// might be C or 5mC
					*c_out = 'N';
					break;
				case 'G':
					// must have 5mC on other strand
					*c_out = 'g';
					*q_out = q2;
					break;
				case 'T':
					// have a T
					*c_out = c2;
					*q_out = q2;
					break;
			}

			if (c1 == 'N' && c2 == 'N')
				*c_out = 'N';

		} else {
			*c_out = c1 == 'N' ? c2 : c1;
			*q_out = c1 == 'N' ? q2 : q1;
		}

		if (*c_out == 'N')
			*q_out = 0;

	} else {
		x = q1<<(PHRED_SHIFT+2+2) | q2<<(2+2) | n2i(c1)<<(2) | n2i(c2);
		if (allow_bs) {
			*c_out = _nt_MAP_n1n2q1q2bs[x];
			*q_out = _Q_MAP_n1n2q1q2bs[x];
		} else {
			*c_out = _nt_MAP_n1n2q1q2[x];
			*q_out = _Q_MAP_n1n2q1q2[x];
		}
	}
}

void
match2(const char *s1, const char *q1, size_t len1,
	const char *s2, const char *q2, size_t len2,
	char *s_out, char *q_out,
	int allow_bs)
{
	int i;
	int len = min(len1, len2);

	for (i=0; i<len; i++) {
		match1bp(s1[i], s2[i], q1[i], q2[i],
				s_out+i, q_out+i,
				allow_bs);
	}
}

/*
 * R1 -> --s1---==hairpin==--s3---
 *       --s4---==hairpin==--s2--- <- R2
 *
 * We first match the top strand to the bottom strand, correcting errors
 * induced by the sequencing platform (match s1 to s4 and s2 to s3).
 * Sequences upstream and downstream of the hairpin are then matched, where
 * discordance is mostly from bisulfite conversion but might also be
 * from errors during sequencing or polymerase copying.
 */
void
match4(const char *_s1, const char *_q1, size_t len1,
	const char *_s2, const char *_q2, size_t len2,
	const char *_s3, const char *_q3, size_t len3,
	const char *_s4, const char *_q4, size_t len4,
	char *s_out, char *q_out)
{
	char *mem, *s1, *s2, *s3, *s4, *q1, *q2, *q3, *q4;

	assert(len1 >= len4);
	assert(len2 >= len3);

	mem = malloc(2*(len1+len2+len3+len4));
	if (mem == NULL) {
		perror("malloc");
		exit(1);
	}

	s1 = mem;
	s2 = s1+len1;
	s3 = s2+len2;
	s4 = s3+len3;
	q1 = s4+len4;
	q2 = q1+len1;
	q3 = q2+len2;
	q4 = q3+len3;

	memcpy(s1, _s1, len1);
	memcpy(s2, _s2, len2);
	memcpy(s3, _s3, len3);
	memcpy(s4, _s4, len4);
	memcpy(q1, _q1, len1);
	memcpy(q2, _q2, len2);
	memcpy(q3, _q3, len3);
	memcpy(q4, _q4, len4);

	revcomp(s4, len4);
	reverse(q4, len4);
	match2(s1+len1-len4, q1+len1-len4, len4,
		s4, q4, len4,
		s1+len1-len4, q1+len1-len4,
		0);

	revcomp(s3, len3);
	reverse(q3, len3);
	match2(s2+len2-len3, q2+len2-len3, len3,
		s3, q3, len3,
		s2+len2-len3, q2+len2-len3,
		0);

	match2(s1, q1, len1,
			s2, q2, len2,
			s_out, q_out,
			1);

	/*{
		int i;
		for (i=0; i<len1; i++)
			putchar(s1[i]);
		printf("\n");
		for (i=0; i<len2; i++)
			putchar(s2[i]);
		printf("\n+\n");
		for (i=0; i<len1; i++)
			putchar(q1[i]);
		printf("\n");
		for (i=0; i<len2; i++)
			putchar(q2[i]);
		printf("\n+\n");
		for (i=0; i<min(len1, len2); i++)
			putchar(s_out[i]);
		printf("\n");
		for (i=0; i<min(len1, len2); i++)
			putchar(q_out[i]);
		printf("\n\n", mm);
	}*/

	free(mem);
}

/*
 * Like match4(), but correct the sequence and qualities in place
 * instead of copying.  And there is no matching of s1 to s2 here.
 */
void
correct_s1s2(char *s1, char *q1, size_t len1,
		char *s2, char *q2, size_t len2,
		size_t hplen, size_t hppos)
{
	if (hppos+hplen < len1) {
		// short molecule, there are valid bases after the hairpin
		char *s3 = s1 + hppos + hplen;
		char *q3 = q1 + hppos + hplen;
		char *s4 = s2 + hppos + hplen;
		char *q4 = q2 + hppos + hplen;

		int len3 = 2*hppos+hplen > len1 ? len1 - hppos - hplen : hppos;
		int len4 = 2*hppos+hplen > len2 ? len2 - hppos - hplen : hppos;
		int len1 = hppos;
		int len2 = hppos;

		assert(len1 >= len4);
		assert(len2 >= len3);

		revcomp(s4, len4);
		reverse(q4, len4);
		match2(s1+len1-len4, q1+len1-len4, len4,
			s4, q4, len4,
			s1+len1-len4, q1+len1-len4,
			0);

		revcomp(s3, len3);
		reverse(q3, len3);
		match2(s2+len2-len3, q2+len2-len3, len3,
			s3, q3, len3,
			s2+len2-len3, q2+len2-len3,
			0);
	}
}

/*
 * Convert ASCII value to phred scaled Q score.
 */
void
clean_quals(const char *s, char *q, size_t len, int phred_scale_in)
{
	int i;

	for (i=0; i<len; i++) {
		char *qi = q+i;

		/*
		if (s[i] == 'N') {
			*qi = 0;
			continue;
		}*/

		*qi -= phred_scale_in;
		if (*qi < 2)
			*qi = 2;
		if (*qi >= PHRED_MAX)
			*qi = PHRED_MAX-1;
	}
}

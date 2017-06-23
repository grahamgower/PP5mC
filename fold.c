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
 * Compare len bytes of s1 and s2, allowing mismatches.
 * Return 0 for a match, -1 otherwise.
 */
static int
mmcmp(const char *s1, const char *s2, size_t len)
{
	int i, mm;

	mm = maxdiff(len, AVG_ERR, MAXDIFF_THRES);

	for (i=0; i<len; i++) {
		mm -= (s1[i] != s2[i]);
		if (mm < 0)
			return -1;
	}

	return 0;
}

/*
 * Search s1 and s2 for hairpin sequences.
 */
void
find_hp_adapter(const char *s1, size_t len1,
		const char *s2, size_t len2,
		const char *hairpin, const char *rhairpin,
		size_t hlen,
		int *h1, int *h2)
{
	int i;

	*h1 = 0;
	*h2 = 0;

	for (i=0; i<len1; i++) {
		if (!mmcmp(s1+i, hairpin, min(len1-i, hlen))) {
			*h1 = i;
			if (!mmcmp(s2+i, rhairpin, min(len2-i, hlen)))
				*h2 = i;
			break;
		}
	}

	if (!*h2) {
		for (i=0; i<len2; i++) {
			if (!mmcmp(s2+i, rhairpin, min(len2-i, hlen))) {
				*h2 = i;
				break;
			}
		}
	}
}

static inline double
q2p(int q)
{
	if (q <= 1)
		return 0.75;
	if (q > 40)
		q = 40;
	return pow(10.0, -q/10.0);
}

static inline int
p2q(double p)
{
	if (p == 0)
		return 40;
	int q = -round(10.0*log10(p));
	if (q > 40)
		q = 40;
	return q;
}

static uint _n2i[] = {['A']=0, ['C']=1, ['G']=2, ['T']=3};
#define n2i(nt) _n2i[(uint)nt]

/*
 * Match a single base pair using posterior base probability as
 * described in Renaud et al. 2014, with minor modifications to
 * allow bisulfite converted base pairs.
 */
static int
match1bp(char c1, char c2, char q1, char q2,
		char *c_out, char *q_out,
		int allow_bs,
		int phred_scale_in, int phred_scale_out)
{
	double p[4], p1[4], p2[4];
	uint q[4], q_max;
	double p_sum;
	int ret;
	int i, i_max;
	int n_maxs;

	c1 = toupper(c1);
	c2 = toupper(c2);
	q1 = q1 -phred_scale_in;
	q2 = q2 -phred_scale_in;

	if (c1 == 'N' || c2 == 'N') {
		*c_out = c1 == 'N' ? c2 : c1;
		*q_out = (c1 == 'N' ? q2 : q1);
		if (*c_out == 'N')
			ret = -1;
		else
			ret = 1;
		goto done;
	}

	p1[n2i('A')] = c1=='A'? 1-q2p(q1) : q2p(q1)/3;
	p1[n2i('C')] = c1=='C'? 1-q2p(q1) : q2p(q1)/3;
	p1[n2i('G')] = c1=='G'? 1-q2p(q1) : q2p(q1)/3;
	p1[n2i('T')] = c1=='T'? 1-q2p(q1) : q2p(q1)/3;

	p2[n2i('A')] = c2=='A'? 1-q2p(q2) : q2p(q2)/3;
	p2[n2i('C')] = c2=='C'? 1-q2p(q2) : q2p(q2)/3;
	p2[n2i('G')] = c2=='G'? 1-q2p(q2) : q2p(q2)/3;
	p2[n2i('T')] = c2=='T'? 1-q2p(q2) : q2p(q2)/3;

	p[n2i('A')] = p1[n2i('A')] * p2[n2i('A')];
	p[n2i('T')] = p1[n2i('T')] * p2[n2i('T')];

	if (allow_bs) {
		p[n2i('C')] = 0.5*(p1[n2i('C')]+p1[n2i('T')]) * p2[n2i('C')];
		p[n2i('G')] = p1[n2i('G')] * 0.5*(p2[n2i('G')]+p2[n2i('A')]);
	} else {
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
		ret = -1;
		goto done;
	}

	*c_out = "ACGT"[i_max];
	*q_out = q[i_max];

	if (allow_bs) {
		// check if methylated
		if (i_max == n2i('C')) {
			if (p1[n2i('C')] > p1[n2i('T')])
				*c_out = 'c';
		} else if (i_max == n2i('G')) {
			if (p2[n2i('G')] > p2[n2i('A')])
				*c_out = 'g';
		}
	}

	if (c1 == c2) {
		ret = 0;
	} else if (allow_bs && ((c1 == 'T' && c2 == 'C') || (c1 == 'G' && c2 == 'A'))) {
		// bisulfite converted.
		ret = 0;
	} else {
		ret = 1;
	}
done:
	*q_out += phred_scale_out;
	return ret;
}

/*
 * Match a single base pair, naiively.
 */
static int
match1bp_naiive(char c1, char c2, char q1, char q2,
		char *c_out, char *q_out,
		int allow_bs,
		int phred_scale_in, int phred_scale_out)
{
	int ret;

	c1 = toupper(c1);
	c2 = toupper(c2);
	q1 = q1 -phred_scale_in;
	q2 = q2 -phred_scale_in;

	if (c1 == 'N' || c2 == 'N') {
		*c_out = c1 == 'N' ? c2 : c1;
		*q_out = (c1 == 'N' ? q2 : q1);
		if (*c_out == 'N')
			ret = -1;
		else
			ret = 1;
	} else if (c1 == c2) {
		if (allow_bs && (c1 == 'C' || c1 == 'G'))
			// Unconverted C<->G.
			c1 = tolower(c1);
		*c_out = c1;
		*q_out = min(q1+q2, 40);
		ret = 0;
	} else if (allow_bs && ((c1 == 'T' && c2 == 'C') || (c1 == 'G' && c2 == 'A'))) {
		// bisulfite converted.
		*c_out = c2 == 'C' ? 'C' : 'G';
		*q_out = min(q1+q2, 40);
		ret = 0;
	} else {
		if (q1 > q2) {
			*c_out = c1;
			*q_out = q1-q2;
		} else if (q2 > q1) {
			*c_out = c2;
			*q_out = q2-q1;
		} else {
			*c_out = 'N';
			*q_out = 0;
		}
		ret = 1;
	}

	*q_out += phred_scale_out;
	return ret;
}

int
match2(const char *s1, const char *q1, size_t len1,
	const char *s2, const char *q2, size_t len2,
	char *s_out, char *q_out,
	int allow_bs,
	int phred_scale_in, int phred_scale_out)
{
	int i;
	int mm = 0;
	int n_count = 0;
	int len = min(len1, len2);

	for (i=0; i<len; i++) {
		int m = match1bp_naiive(s1[i], s2[i], q1[i], q2[i],
				s_out+i, q_out+i,
				allow_bs,
				phred_scale_in, phred_scale_out);
		switch (m) {
			case -1:
				n_count++;
				break;
			case 0:
				break;
			case 1:
				mm++;
				break;
		}
	}

	// add length mismatch
	//mm += abs(len1-len2);

	if (3*n_count > len)
		// I don't trust this pair of reads
		mm += n_count;

	return mm;
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
int
match4(const char *_s1, const char *_q1, size_t len1,
	const char *_s2, const char *_q2, size_t len2,
	const char *_s3, const char *_q3, size_t len3,
	const char *_s4, const char *_q4, size_t len4,
	char *s_out, char *q_out,
	int phred_scale_in, int phred_scale_out)
{
	int mm;
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
		0, phred_scale_in, phred_scale_out);

	revcomp(s3, len3);
	reverse(q3, len3);
	match2(s2+len2-len3, q2+len2-len3, len3,
		s3, q3, len3,
		s2+len2-len3, q2+len2-len3,
		0, phred_scale_in, phred_scale_out);

	mm = match2(s1, q1, len1,
			s2, q2, len2,
			s_out, q_out,
			1, phred_scale_in, phred_scale_out);

	free(mem);

	return mm;
}

/*
 * Like match4(), but correct the sequence and qualities in place
 * instead of copying.  And there is no matching of s1 to s2 here.
 */
int
correct_s1s2(char *s1, char *q1, size_t len1,
		char *s2, char *q2, size_t len2,
		const char *hairpin, const char *rhairpin, size_t hlen,
		int phred_scale_in, int phred_scale_out)
{
	int h1, h2;

	find_hp_adapter(s1, len1, s2, len2,
			hairpin, rhairpin, hlen,
			&h1, &h2);

	if (h1 != 0 && h2 != 0 && h1 != h2)
		// dislocated hairpins
		return -1;

	if (h1+hlen < len1 && h2+hlen < len2) {
		// short molecule, there are valid bases after the hairpin
		char *s3 = s1 + h1 + hlen;
		char *q3 = q1 + h1 + hlen;
		char *s4 = s2 + h2 + hlen;
		char *q4 = q2 + h2 + hlen;

		int len3 = 2*h1+hlen > len1 ? len1 - h1 - hlen : h1;
		int len4 = 2*h2+hlen > len2 ? len2 - h2 - hlen : h2;
		int len1 = h1;
		int len2 = h2;

		assert(len1 >= len4);
		assert(len2 >= len3);

		revcomp(s4, len4);
		reverse(q4, len4);
		match2(s1+len1-len4, q1+len1-len4, len4,
			s4, q4, len4,
			s1+len1-len4, q1+len1-len4,
			0, phred_scale_in, phred_scale_out);

		revcomp(s3, len3);
		reverse(q3, len3);
		match2(s2+len2-len3, q2+len2-len3, len3,
			s3, q3, len3,
			s2+len2-len3, q2+len2-len3,
			0, phred_scale_in, phred_scale_out);
	}

	return 0;
}

/*
 * Parse the XF:Z sam/bam field, containing pipe delimited sequence
 * and quality scores, i.e. XF:Z:s1|s2|q1|q2.
 */
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

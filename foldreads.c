/*
 * Fold r1/r2 sequences back together at the hairpin.
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
#include <ctype.h>
#include <errno.h>
#include <getopt.h>
#include <stdint.h>
#include <math.h>
#include <zlib.h>
#include "kseq.h"
#include "khash.h"

KSEQ_INIT(gzFile, gzread);
KHASH_SET_INIT_STR(str);

#define min(a,b) ((a)<(b)?(a):(b))
#define max(a,b) ((a)>(b)?(a):(b))

typedef struct {
	char *hairpin, *rhairpin; // Hairpin sequence, reverse complement.
	size_t hlen; // Hairpin length.
	char *fn1, *fn2; // Input filenames.
	FILE *fos; // File pointer for outputting sequences.
	FILE *fom; // File pointer for outputting metrics.
	FILE *f_unmatched_r1, *f_unmatched_r2; // File pointers for unfolded reads
	int phred_scale;
} opt_t;

#define N_POS 30
typedef uint64_t pairs_t[N_POS][16];

// count things
typedef struct {
	uint64_t total_reads;
	uint64_t folded_pairs;
	uint64_t hairpin_complete;
	uint64_t unique_hairpin_complete;
	uint64_t ntcomp[4]; // nucleotide composition

	uint64_t match_total;
	uint64_t damage_total;
	uint64_t mm_total;

	uint64_t mm_hairpin[16]; // X->Y mismatches in the hairpin

	khash_t(str) *seqmap;
	uint64_t clones;

	struct {
		struct {
			// unmethylated, methylated, hemimethylated
			uint64_t u, m, h;
		} cpg, chg, chh; // CpG, CHG, CHH contexts
	} p, m; // plus, minus strand

	uint64_t hemi_cpg, hemi_chg;

#define DS_DIST_FROM_END 20
	// X<->Y pairings which are putatively in the double stranded region
	// of the fragment - i.e. at least DS_DIST_FROM_END from an end.
	uint64_t ds_pairs[16];

	// X<->Y pairs for each position within a fragment,
	// for up to N_POS bases from an end.
	pairs_t pairs_l, pairs_r;

	uint64_t **pairs;
} metrics_t;

static int nt2int[] = {['A']=0, ['C']=1, ['G']=2, ['T']=3, ['a']=0, ['c']=1, ['g']=2, ['t']=3};
static char cmap[] = {['A']='T', ['C']='G', ['G']='C', ['T']='A', ['N']='N', ['n']='N',
				['a']='t', ['c']='g', ['g']='c', ['t']='a'};

/*
 * Inverse poisson CDF, stolen from bwa: bwtaln.c
 */
#define AVG_ERR 0.02
#define MAXDIFF_THRES 0.01
int
bwa_cal_maxdiff(int l, double err, double thres)
{
	double elambda = exp(-l * err);
	double sum, y = 1.0;
	int k, x = 1;
	for (k = 1, sum = elambda; k < 1000; ++k) {
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
	static int maxdiff[1000];
	static int init = 0;

	if (!init) {
		int i;
		for (i=0; i<1000; i++)
			maxdiff[i] = bwa_cal_maxdiff(i, err, thres);
		init = 1;
	}

	return l<1000 ? maxdiff[l] : maxdiff[999];
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
 * Compare s1 and s2 to the hairpin h1 and its reverse complement h2.
 * Return 0 for a match, -1 otherwise.
 */
int
hpcmp(const char *s1, const char *s2, const char *h1, const char *h2, size_t len)
{
	int i, mm;

	mm = maxdiff(len*2, AVG_ERR, MAXDIFF_THRES);

	for (i=0; i<len; i++) {
		mm -= (s1[i] != h1[i]) + (s2[i] != h2[i]);
		if (mm < 0)
			return -1;
	}

	return 0;
}

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
	int i, j;
	int hairpin = 0;

	if (len1 != len2) {
		fprintf(stderr, "Error: input r1/r2 reads have different lengths.\n");
		exit(1);
		/*
		*_s_out = *_q_out = NULL;
		return 0;
		*/
	}

	// metrics
	uint32_t ntcomp[4] = {0,};
	uint32_t match_total = 0;
	uint32_t damage_total = 0;
	uint32_t mm_total = 0;
	uint32_t ds_pairs[16] = {0,};
	uint32_t pairs[len1][16];
	memset(pairs, 0, sizeof(pairs));

	s_out = malloc(len1+1);
	q_out = malloc(len1+1);

	if (s_out == NULL || q_out == NULL) {
		perror("malloc");
		exit(-2);
	}

	for (i=0; i<len1; i++) {

		// Look for the hairpin in both reads.
		if (!hpcmp(s1+i, s2+i, opt->hairpin, opt->rhairpin, min(len1-i, opt->hlen))) {
			if (len1-i >= opt->hlen) {
				// full hairpin
				hairpin = 2;
				metrics->hairpin_complete++;
			} else
				// incomplete hairpin
				hairpin = 1;
			break;
		}

		char c1 = toupper(s1[i]);
		char c2 = toupper(s2[i]);

		if (c1 == 'N' || c2 == 'N') {
			// Use the other read, if we can.
			s_out[i] = c2 == 'N' ? c1 : c2;
			q_out[i] = c2 == 'N' ? q1[i] : q2[i];
			continue;
		}

		int nt_idx = nt2int[(int)c1] << 2 | nt2int[(int)cmap[(int)c2]];
		pairs[i][nt_idx]++;

		if (i > DS_DIST_FROM_END && i < (len1-DS_DIST_FROM_END))
			// Assume this position is not derived from a
			// single strand overhang.
			ds_pairs[nt_idx]++;

		if (c1 == c2) {
			if (c2 == 'C' || c2 == 'G')
				// Unconverted C<->G.
				c2 = tolower(c2);
			s_out[i] = c2;
			q_out[i] = min((int)q1[i]+(int)q2[i]-2*opt->phred_scale, 40) +33;
			match_total++;
			ntcomp[nt2int[(int)c2]]++;
		} else if ((c1 == 'T' && c2 == 'C') || (c1 == 'G' && c2 == 'A')) {
			// Putatively damaged, or bisulfite converted.
			s_out[i] = c2 == 'C' ? 'C' : 'G';
			q_out[i] = min((int)q1[i]+(int)q2[i]-2*opt->phred_scale, 40) +33;
			damage_total++;
			ntcomp[nt2int[(int)s_out[i]]]++;
		} else {
			// Mismatch, take highest quality base.
			if (q1[i] > q2[i]) {
				s_out[i] = c1;
				q_out[i] = max((int)q1[i]-(int)q2[i], 0) +33;
			} else if (q1[i] < q2[i]) {
				s_out[i] = c2;
				q_out[i] = max((int)q2[i]-(int)q1[i], 0) +33;
			} else {
				s_out[i] = 'N';
				q_out[i] = '!';
			}
			mm_total++;
		}
	}

	int last_idx = i;

	if (mm_total > maxdiff(len1, AVG_ERR, MAXDIFF_THRES)) {
		// Sequences are not complementary.
		free(s_out);
		free(q_out);
		*_s_out = *_q_out = NULL;
		return 0;
	}

	if (i == 0) {
		// Hairpin at the start of both reads.
		free(s_out);
		free(q_out);
		*_s_out = *_q_out = NULL;
		return 0;
	}

	{
		// Check for duplicate sequences in the seqmap.
		khint_t k;
		int ret;
		char c = s_out[i]; // save first character of the hairpin
		s_out[i] = '\0';

		k = kh_put(str, metrics->seqmap, s_out, &ret);
		switch (ret) {
			case -1:
				fprintf(stderr, "Out of memory\n");
				exit(1);
			case 0:
				// duplicate
				metrics->clones++;
				goto no_metrics;
			case 1:
			case 2:
				// unique
				kh_key(metrics->seqmap, k) = strdup(s_out);
				break;
		}
		s_out[i] = c;
	}

	/*
	 * Record metrics.
	 */

	// For full hairpins only.
	if (hairpin == 2) {
		// Hairpin mismatches
		int end = last_idx+opt->hlen;
		for (i=last_idx, j=end-1; i<end; i++, j--) {
			char c1 = s1[i];
			char c2 = cmap[(int)s2[j]];
			char h = opt->hairpin[i-last_idx];
			int nt_idx;
			if (c1 != h && c1 != 'N') {
				// r1 has incorrect base
				nt_idx = nt2int[(int)h] << 2 | nt2int[(int)c1];
				metrics->mm_hairpin[nt_idx]++;
			}
			if (c2 != h && c2 != 'N') {
				// r2 has incorrect base
				nt_idx = nt2int[(int)cmap[(int)h]] << 2 | nt2int[(int)s2[j]];
				metrics->mm_hairpin[nt_idx]++;
			}
			static int nterr[256] = {['A']=1, ['C']=1, ['G']=1, ['T']=1, ['N']=1};
			if (!nterr[(int)c1])
				fprintf(stderr, "%s\n", s1);
		}

		if (last_idx >= 2*N_POS) {
			// Folded watson/crick pairs.
			for (i=0; i<N_POS; i++) {
				for (j=0; j<16; j++) {
					metrics->pairs_l[i][j] += pairs[i][j];
					metrics->pairs_r[i][j] += pairs[last_idx-i-1][j];
				}
			}
			for (j=0; j<16; j++)
				metrics->ds_pairs[j] += ds_pairs[j];

		}

		// Contextual methylation metrics.
		for (i=0, j=last_idx-1; i<last_idx-1; i++, j--) {

			// (+) strand
			switch (s_out[i]) {
				case 'C':
					// CpG?
					if (s_out[i+1] == 'G' || s_out[i+1] == 'g')
						metrics->p.cpg.u++;
					else {
						// CHG?
						if (i+2 < last_idx && s_out[i+1] != 'N') {
							if (s_out[i+2] == 'G' || s_out[i+2] == 'g')
								metrics->p.chg.u++;
							else {
								// CHH?
								if (s_out[i+2] != 'N') {
									metrics->p.chh.u++;
								}
							}
						}
					}
					break;
				case 'c':
					// cpG?
					if (s_out[i+1] == 'G')
						metrics->p.cpg.h++;
					else if (s_out[i+1] == 'g')
						metrics->p.cpg.m++;
					else {
						// cHG?
						if (i+2 < last_idx && s_out[i+1] != 'N') {
							if (s_out[i+2] == 'G')
								metrics->p.chg.h++;
							else if (s_out[i+2] == 'g')
								metrics->p.chg.m++;
							else {
								// cHH?
								if (s_out[i+2] != 'N') {
									metrics->p.chh.m++;
								}
							}
						}
					}
					break;
			}

			// (-) strand
			switch (s_out[j]) {
				case 'G':
					// GpC?
					if (s_out[j-1] == 'C' || s_out[j-1] == 'c')
						metrics->m.cpg.u++;
					else {
						// GHC?
						if (j-2 > 0 && s_out[j-1] != 'N') {
							if (s_out[j-2] == 'C' || s_out[j-2] == 'c')
								metrics->m.chg.u++;
							else {
								// HHC?
								if (s_out[j-2] != 'N') {
									metrics->m.chh.u++;
								}
							}
						}
					}
					break;
				case 'g':
					// Gpc?
					if (s_out[j-1] == 'C')
						metrics->m.cpg.h++;
					else if (s_out[j-1] == 'c')
						metrics->m.cpg.m++;
					else {
						// GHc?
						if (j-2 > 0 && s_out[j-1] != 'N') {
							if (s_out[j-2] == 'C')
								metrics->m.chg.h++;
							else if (s_out[j-2] == 'c')
								metrics->m.chg.m++;
							else {
								// HHc?
								if (s_out[j-2] != 'N') {
									metrics->m.chh.m++;
								}
							}
						}
					}
					break;
			}
		}

		metrics->unique_hairpin_complete++;
	}

	metrics->match_total += match_total;
	metrics->damage_total += damage_total;
	metrics->mm_total += mm_total;

	for (i=0; i<4; i++)
		metrics->ntcomp[i] += ntcomp[i];

no_metrics:
	s_out[last_idx] = '\0';
	q_out[last_idx] = '\0';
	*_s_out = s_out;
	*_q_out = q_out;

	return 1;
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
		
		metrics->total_reads++;

		if (len1 < opt->hlen || len2 < opt->hlen)
			continue;

		if (seq1->qual.l == 0) {
			fprintf(stderr, "%s: qual scores required.\n", opt->fn1);
			ret = 2;
			goto err2;
		}

		if (seq2->qual.l == 0) {
			fprintf(stderr, "%s: qual scores required.\n", opt->fn2);
			ret = 2;
			goto err2;
		}

		f = fold(opt, metrics, seq1->seq.s, seq1->qual.s, seq1->seq.l,
			seq2->seq.s, seq2->qual.s, seq2->seq.l, &s_out, &q_out);
		if (f) {
			char *comment = seq1->comment.l==0 ? NULL : seq1->comment.s;
			seq_write(seq1->name.s, comment, s_out, q_out, opt->fos);
			free(s_out);
			free(q_out);
			metrics->folded_pairs++;
		} else if (opt->f_unmatched_r1) {
			char *comment = seq1->comment.l==0 ? NULL : seq1->comment.s;
			seq_write(seq1->name.s, comment, seq1->seq.s, seq1->qual.s, opt->f_unmatched_r1);
			comment = seq2->comment.l==0 ? NULL : seq2->comment.s;
			seq_write(seq2->name.s, comment, seq2->seq.s, seq2->qual.s, opt->f_unmatched_r2);
		}
	}

	if (len1 != -1) {
		fprintf(stderr, "%s: unexpected end of file.\n", opt->fn1);
		ret = 3;
		goto err2;
	}

	if (len2 != -1) {
		fprintf(stderr, "%s: unexpected end of file.\n", opt->fn2);
		ret = 3;
		goto err2;
	}

	ret = 0;

err2:
	kseq_destroy(seq2);
	kseq_destroy(seq1);
	gzclose(fp2);
err1:
	gzclose(fp1);
err0:
	return ret;
}

void
print_metrics(const opt_t *opt, const metrics_t *metrics)
{
	FILE *fp = opt->fom;
	const pairs_t *pairs;
	uint32_t hp_comp[4] = {0,};
	uint64_t total_bases, total_bases2, total_ds_pairs;
	int i, j, k;

	// Total nucleotide count, excluding mismatches.
	total_bases = 0;
	for(i=0; i<4; i++) {
		total_bases += metrics->ntcomp[i];
	}

	// Total nucleotide count, including mismatches.
	total_bases2 = metrics->match_total + metrics->damage_total + metrics->mm_total;

	fprintf(fp, "Total read pairs: %jd\n", (uintmax_t)metrics->total_reads);
	fprintf(fp, "Number of successfully folded read pairs: %jd\n", (uintmax_t)metrics->folded_pairs);
	fprintf(fp, "Number of read pairs with complete hairpin: %jd\n",
			(uintmax_t)metrics->hairpin_complete);
	fprintf(fp, "Clonality of folded pairs: %lf\n\n", (double)metrics->clones/metrics->folded_pairs);

	fprintf(fp, "Base pair +/- concordance: %lf\n",
			(double)metrics->match_total/total_bases2);
	fprintf(fp, "Base pair +/- discordance (C<->T or G<->A): %lf\n",
			(double)metrics->damage_total/total_bases2);
	fprintf(fp, "Base pair +/- discordance (other): %lf\n\n",
			(double)metrics->mm_total/total_bases2);

	fprintf(fp, "Base frequency: A=%lf, C=%lf, G=%lf, T=%lf, GC=%lf\n",
			(double)metrics->ntcomp[0]/total_bases,
			(double)metrics->ntcomp[1]/total_bases,
			(double)metrics->ntcomp[2]/total_bases,
			(double)metrics->ntcomp[3]/total_bases,
			(double)(metrics->ntcomp[1]+metrics->ntcomp[2])/total_bases);
	fprintf(fp, "Methylation frequency (+): CpG=%lf, CHG=%lf, CHH=%lf\n",
			(double)(metrics->p.cpg.m+metrics->p.cpg.h) / (metrics->p.cpg.m+metrics->p.cpg.h+metrics->p.cpg.u),
			(double)(metrics->p.chg.m+metrics->p.chg.h) / (metrics->p.chg.m+metrics->p.chg.h+metrics->p.chg.u),
			(double)(metrics->p.chh.m) / (metrics->p.chg.m+metrics->p.chg.u));
	fprintf(fp, "Methylation frequency (-): CpG=%lf, CHG=%lf, CHH=%lf\n",
			(double)(metrics->m.cpg.m+metrics->m.cpg.h) / (metrics->m.cpg.m+metrics->m.cpg.h+metrics->m.cpg.u),
			(double)(metrics->m.chg.m+metrics->m.chg.h) / (metrics->m.chg.m+metrics->m.chg.h+metrics->m.chg.u),
			(double)(metrics->m.chh.m) / (metrics->m.chg.m+metrics->m.chg.u));
	fprintf(fp, "Hemimethylation frequency (+/- discordance): CpG=%lf, CHG=%lf\n",
			(double)(metrics->p.cpg.h+metrics->m.cpg.h)/(metrics->p.cpg.m+metrics->p.cpg.h+metrics->p.cpg.u + metrics->m.cpg.m+metrics->m.cpg.h+metrics->m.cpg.u),
			(double)(metrics->p.chg.h+metrics->m.chg.h)/(metrics->p.chg.m+metrics->p.chg.h+metrics->p.chg.u + metrics->m.chg.m+metrics->m.chg.h+metrics->m.chg.u));

	// Nucleotide composition of the hairpin.
	for (i=0; i<opt->hlen; i++)
		hp_comp[nt2int[(int)opt->hairpin[i]]]++;

	fprintf(fp, "\nHairpin sequencing mismatches (Hairpin->Observed shown as Row->Column):\n");
	fprintf(fp, " \tA\t\tC\t\tG\t\tT\t\tTotal\n");
	double col[4] = {0,};
	for (i=0; i<4; i++) {
		double row = 0;
		fprintf(fp, "%c", "ACGT"[i]);
		for (j=0; j<4; j++) {
			double mm = (double)metrics->mm_hairpin[i<<2|j]/(2.0 * hp_comp[i] * metrics->unique_hairpin_complete);
			col[j] += mm;
			row += mm;
			fprintf(fp, "\t%lf", mm);
		}
		fprintf(fp, "\t%lf\n", row);

	}
	fprintf(fp, "Total\t%lf\t%lf\t%lf\t%lf\t%lf\n",
			col[0], col[1], col[2], col[3],
			col[0]+col[1]+col[2]+col[3]);

	total_ds_pairs = 0;
	for (i=0; i<16; i++)
		total_ds_pairs += metrics->ds_pairs[i];

	fprintf(fp, "\nWatson<->Crick [>%dbp from fragment end]\n", DS_DIST_FROM_END);
	for (i=0; i<4; i++) {
		for (j=0; j<4; j++) {
			fprintf(fp, "%c<->%c: %lf\n", "ACGT"[i], "ACGT"[j],
				(double)metrics->ds_pairs[i<<2|j]/total_ds_pairs);
		}
	}


	pairs = &metrics->pairs_l;
	for (;;) {
		uint64_t bases[N_POS] = {0,};
		for (i=0; i<4; i++) {
			for (j=0; j<4; j++) {
				for (k=0; k<N_POS; k++) {
					bases[k] += (*pairs)[k][i<<2|j];
				}
			}
		}

		fprintf(fp, "\nWatson<->Crick [%c']", pairs == &metrics->pairs_l ? '5' : '3');
		for (k=0; k<N_POS; k++)
			fprintf(fp, "\t%d", k);
		fprintf(fp, "\n");

		for (i=0; i<4; i++) {
			for (j=0; j<4; j++) {
				fprintf(fp, "%c<->%c", "ACGT"[i], "ACGT"[j]);
				for (k=0; k<N_POS; k++) {
					fprintf(fp, "\t%lf", (double)(*pairs)[k][i<<2|j] / bases[k]);
				}
				fprintf(fp, "\n");
			}
		}

		if (pairs == &metrics->pairs_r)
			break;
		else
			pairs = &metrics->pairs_r;
	}
}

void
usage(char *argv0)
{
	fprintf(stderr, "foldreads v3\n");
	fprintf(stderr, "usage: %s [-o OUT.fq] [-m FILE] [-u PFX] -p SEQ -1 IN1.fq -2 IN2.fq\n", argv0);
	fprintf(stderr, " -o OUT.fq         Fastq output file [stdout].\n");
	fprintf(stderr, " -m FILE           Metrics output file [stderr].\n");
	fprintf(stderr, " -u PFX            Filename prefix for unfolded reads []");
	fprintf(stderr, " -p SEQ            The hairpin SEQuence.\n");
	fprintf(stderr, " -1 IN1.fq[.gz]    R1 fastq input file.\n");
	fprintf(stderr, " -2 IN2.fq[.gz]    R2 fastq input file.\n");
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
	char *fos_fn = NULL, *fom_fn = NULL;

	memset(&opt, '\0', sizeof(opt_t));
	memset(&metrics, '\0', sizeof(metrics_t));
	opt.fos = stdout;
	opt.fom = stderr;
	opt.phred_scale = 33; // Input phred scale; we always output phred+33.

	while ((c = getopt(argc, argv, "o:m:p:u:1:2:")) != -1) {
		switch (c) {
			case 'o':
				fos_fn = optarg;
				break;
			case 'm':
				fom_fn = optarg;
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

	metrics.seqmap = kh_init(str);
	ret = foldreads_pe(&opt, &metrics);

	khint_t k;
	for (k=kh_begin(metrics.seqmap); k!=kh_end(metrics.seqmap); k++) {
		if (kh_exist(metrics.seqmap, k))
			free((char *)kh_key(metrics.seqmap, k));
	}
	kh_destroy(str, metrics.seqmap);

	free(opt.rhairpin);

	if (unmatched_pfx) {
		fclose(opt.f_unmatched_r1);
		fclose(opt.f_unmatched_r2);
	}

	if (fos_fn)
		fclose(opt.fos);

	if (fom_fn) {
		opt.fom = fopen(fom_fn, "w");
		if (opt.fom == NULL) {
			fprintf(stderr, "Error: %s: %s\n", fom_fn, strerror(errno));
			usage(argv[0]);
		}
	}

	print_metrics(&opt, &metrics);

	if (fom_fn)
		fclose(opt.fom);


	return ret;
}

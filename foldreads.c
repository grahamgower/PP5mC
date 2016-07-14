/*
 * Fold input sequences at the hairpin.
 *
 * Copyright (c) 2015 Graham Gower <graham.gower@gmail.com>
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
#include <zlib.h>
#include "kseq.h"

KSEQ_INIT(gzFile, gzread);

#define min(a,b) ((a)<(b)?(a):(b))
#define max(a,b) ((a)>(b)?(a):(b))

typedef struct {
	int bisulfite; // Input data is from a bisulfite treated library.
	char *hairpin, *rhairpin; // Hairpin sequence, reverse complement.
	size_t hlen; // Hairpin length.
	char *fn1, *fn2; // Input filenames.
	FILE *fos; // File pointer for outputting sequences.
	FILE *fom; // File pointer for outputting metrics.
} opt_t;

#define N_POS 20
typedef uint64_t pairs_t[N_POS][16];

// count things
typedef struct {
	uint64_t total_reads;
	uint64_t hairpin_missing;
	uint64_t palindrome_missing;
	uint64_t ntcomp[4]; // nucleotide composition

	uint64_t match_total;
	uint64_t damage_total;
	uint64_t mm_total;

	uint64_t mm_hairpin[16]; // X->Y mismatches in the hairpin

#define DS_DIST_FROM_END 15
	// X<->Y pairings which are putatively in the double stranded region
	// of the fragment - i.e. at least DS_DIST_FROM_END from an end.
	uint64_t ds_pairs[16];

	// X<->Y pairs for each position within a fragment,
	// for up to N_POS bases from an end.
	pairs_t pairs_l, pairs_r;
} metrics_t;

static int nt2int[] = {['A']=0, ['C']=1, ['G']=2, ['T']=3, ['a']=0, ['c']=1, ['g']=2, ['t']=3};

/*
 * Reverse complement of s.
 */
void
revcomp(char *s, size_t len)
{
	static char cmap[] = {['A'] = 'T', ['C'] = 'G', ['G'] = 'C', ['T'] = 'A', ['N'] = 'N'};
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
 * Reverse string s.
 */
void
reverse(char *s, size_t len)
{
	int i, j;
	char tmp;

	for (i=0, j=len-1; i<len/2; i++, j--) {
		tmp = s[i];
		s[i] = s[j];
		s[j] = tmp;
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
 * Fold the input sequences, s1 and s2.
 * Folded sequence is placed in s_out, folded quality scores placed in q_out.
 * Returns 1 if folding was successful, 0 otherwise.
 *
 * TODO: try gapped alignment if the simple case fails.
 */
int
fold(const opt_t *opt, metrics_t *metrics,
		const char *s1, const char *q1, size_t len1,
		const char *s2, const char *q2, size_t len2,
		char **_s_out, char **_q_out)
{
	static char cmap[] = {['A']='T', ['C']='G', ['G']='C', ['T']='A',
				['N']='N', ['n']='N',
				['a']='t', ['c']='g', ['g']='c', ['t']='a'};
	char *s_out, *q_out;
	int i, j;

	if (len1 != len2) {
		// TODO: may need gapped alignment
		*_s_out = *_q_out = NULL;
		metrics->palindrome_missing++;
		return 0;
	}

	// metrics
	uint32_t ntcomp[4] = {0,};
	uint32_t match_total = 0;
	uint32_t damage_total = 0;
	uint32_t mm_total = 0;
	uint32_t ds_pairs[16] = {0,};
	pairs_t pairs_l = {{0,},};
	pairs_t pairs_r = {{0,},};

	s_out = malloc(len1+1);
	q_out = malloc(len1+1);

	if (s_out == NULL || q_out == NULL) {
		perror("malloc");
		exit(-2);
	}

	for (i=0, j=len2-1; i<len1 && j>=0; i++, j--) {
		char c1 = s1[i];
		char c2 = cmap[(int)s2[j]];

		if (c1 == 'N' || c2 == 'N') {
			s_out[i] = c2 == 'N' ? c1 : c2;
			q_out[i] = c2 == 'N' ? q1[i] : q2[j];
			continue;
		}

		int nt_idx = nt2int[(int)c1] << 2 | nt2int[(int)s2[j]];

		if (i < N_POS)
			pairs_l[i][nt_idx]++;
		if (j < N_POS)
			pairs_r[j][nt_idx]++;

		if (i > DS_DIST_FROM_END && j > DS_DIST_FROM_END)
			// Assume this position is not derived from a
			// single strand overhang.
			ds_pairs[nt_idx]++;

		if (c1 == c2) {
			if (opt->bisulfite && (c2 == 'C' || c2 == 'G'))
				// Unconverted C<->G.
				c2 = tolower(c2);
			s_out[i] = c2;
			// TODO: combine quality scores
			q_out[i] = max(q1[i], q2[j]);
			match_total++;
			ntcomp[nt2int[(int)c2]]++;
		} else if ((c1 == 'T' && c2 == 'C') || (c1 == 'G' && c2 == 'A')) {
			// Putatively damaged, or bisulfite converted.
			s_out[i] = c2 == 'C' ? 'C' : 'G';
			q_out[i] = c2 == 'C' ? q2[j] : q1[i];
			damage_total++;
			ntcomp[nt2int[(int)s_out[i]]]++;
		} else {
			// Mismatch, probably sequencing error.

			// Take the highest quality base.
			// TODO: combine quality scores
			if (q1[i] > q2[j]) {
				s_out[i] = c1;
				q_out[i] = q1[i];
			} else {
				s_out[i] = c2;
				q_out[i] = q2[j];
			}
			mm_total++;
		}
	}

	if (mm_total > match_total) {
		// This is not a palindromic sequence.
		// TODO: may need gapped alignment
		free(s_out);
		free(q_out);
		*_s_out = *_q_out = NULL;
		metrics->palindrome_missing++;
		return 0;
	}

	s_out[i] = '\0';
	q_out[i] = '\0';

	metrics->match_total += match_total;
	metrics->damage_total += damage_total;
	metrics->mm_total += mm_total;

	for (i=0; i<4; i++)
		metrics->ntcomp[i] += ntcomp[i];

	for (j=0; j<16; j++) {
		for (i=0; i<N_POS; i++) {
			metrics->pairs_l[i][j] += pairs_l[i][j];
			metrics->pairs_r[i][j] += pairs_r[i][j];
		}
		metrics->ds_pairs[j] += ds_pairs[j];
	}

	*_s_out = s_out;
	*_q_out = q_out;
	
	return 1;
}

/*
 * Find the hairpin in sequence s.
 * Returns the index where the hairpin is found, or -1 otherwise.
 * TODO: look harder, don't assume hairpin is always in the middle.
 */
int
find_hairpin_se(const opt_t *opt, metrics_t *metrics, const char *s, size_t slen)
{
	int len;
	int i;
	int h_mm = 0, h_mm_i = -1;

	if ((slen - opt->hlen) % 2 != 0) {
		// TODO: look harder
		metrics->hairpin_missing++;
		return -1;
	}

	// Folded fragment length.
	len = (slen - opt->hlen) / 2;

	for (i=0; i<opt->hlen; i++) {
		if (s[i+len] != opt->hairpin[i]) {
			// Allow one mismatch only in the hairpin.
			if (++h_mm > 1) {
				// Missing hairpin.
				metrics->hairpin_missing++;
				return -1;
			}
			h_mm_i = i;
		}
	}

	if (h_mm) {
		int idx = nt2int[(int)opt->hairpin[h_mm_i]] << 2 | nt2int[(int)s[len+h_mm_i]];
		metrics->mm_hairpin[idx]++;
	}

	return len;
}

int
strlcmp_mm(const char *s1, const char *s2, size_t len)
{
	int i, mm=0;
	for (i=0; i<len; i++) {
		if (s1[i] != s2[i]) {
			if (mm++ == 0)
				return -1;
		}
	}

	return 0;
}

/*
 * Search for a hairpin in sequence s, starting from position start.
 * Mismatches are not allowed.  If partial is 1, accept partial matches to the
 * hairpin at the end of the sequence.
 * Returns the index where the hairpin is found, or -1 otherwise.
 */
int
find_hairpin(char *hairpin, size_t hlen, const char *s, size_t slen, int start, int partial)
{
	int i;
	int end = partial ? slen : slen-hlen;

	for (i=start; i<end; i++) {
		if (strlcmp_mm(s+i, hairpin, min(hlen, slen-i)) == 0)
			return i;
	}

	return -1;
}

/*
 * Paired end reads which remained uncollapsed after read merging
 * by e.g. AdapterRemoval.  We expect the two reads to be palindromic.
 */
int
foldreads_pe(const opt_t *opt, metrics_t *metrics)
{
	gzFile fp1, fp2;
	kseq_t *seq1, *seq2;
	int len1, len2;
	int ret, f;
	char *s_out, *q_out;
	int hp_pos;

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

		hp_pos = find_hairpin(opt->hairpin, opt->hlen, seq1->seq.s, seq1->seq.l, seq1->seq.l-opt->hlen, 1);
		if (hp_pos > 0) {
			seq1->seq.l = hp_pos;
			seq1->qual.l = hp_pos;
		}

		hp_pos = find_hairpin(opt->rhairpin, opt->hlen, seq2->seq.s, seq2->seq.l, seq2->seq.l-opt->hlen, 1);
		if (hp_pos > 0) {
			seq2->seq.l = hp_pos;
			seq2->qual.l = hp_pos;
		}

		revcomp(seq2->seq.s, seq2->seq.l);
		reverse(seq2->qual.s, seq2->qual.l);

		f = fold(opt, metrics, seq1->seq.s, seq1->qual.s, seq1->seq.l,
			seq2->seq.s, seq2->qual.s, seq2->seq.l, &s_out, &q_out);
		if (f) {
			char *comment = seq1->comment.l==0 ? NULL : seq1->comment.s;
			seq_write(seq1->name.s, comment, s_out, q_out, opt->fos);
			free(s_out);
			free(q_out);
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


/*
 * Single ended reads are assumed to be paired end reads collapsed
 * by e.g. AdapterRemoval.  Hence we can expect the terminal regions
 * to be palindromic.
 */
int
foldreads_se(const opt_t *opt, metrics_t *metrics)
{
	gzFile fp1;
	kseq_t *seq;
	int len;
	int ret, f;
	char *s_out, *q_out;
	int hp_pos; // hairpin position
	int s2_start;

	fp1 = gzopen(opt->fn1, "r");
	if (fp1 == NULL) {
		fprintf(stderr, "%s: %s\n", opt->fn1, strerror(errno));
		ret = 1;
		goto err0;
	}

	seq = kseq_init(fp1);

	while ((len = kseq_read(seq)) >= 0) {
		metrics->total_reads++;

		if (len < opt->hlen)
			continue;

		if (seq->qual.l == 0) {
			fprintf(stderr, "%s: qual scores required.\n", opt->fn1);
			ret = 2;
			goto err1;
		}

		hp_pos = find_hairpin_se(opt, metrics, seq->seq.s, seq->seq.l);
		if (hp_pos < 0)
			continue;

		s2_start = hp_pos + opt->hlen;

		f = fold(opt, metrics, seq->seq.s, seq->qual.s, hp_pos,
			seq->seq.s+s2_start, seq->qual.s+s2_start,
			seq->seq.l-s2_start, &s_out, &q_out);
		if (f) {
			char *comment = seq->comment.l==0 ? NULL : seq->comment.s;
			seq_write(seq->name.s, comment, s_out, q_out, opt->fos);
			free(s_out);
			free(q_out);
		}
	}

	if (len != -1) {
		fprintf(stderr, "%s: unexpected end of file.\n", opt->fn1);
		ret = 3;
		goto err1;
	}

	ret = 0;

err1:
	kseq_destroy(seq);
	gzclose(fp1);
err0:
	return ret;
}

void
print_metrics(const opt_t *opt, const metrics_t *metrics)
{
	FILE *fp = opt->fom;
	const pairs_t *pairs;
	uint32_t hp_comp[4];
	uint64_t n_hairpins = metrics->total_reads - metrics->hairpin_missing;
	uint64_t total_bases, total_bases2, total_ds_pairs;
	int i, j, k;

	// Total nucleotide count, excluding mismatches.
	total_bases = 0;
	for(i=0; i<4; i++) {
		total_bases += metrics->ntcomp[i];
	}

	// Total nucleotide count, including mismatches.
	total_bases2 = metrics->match_total + metrics->damage_total + metrics->mm_total;

	fprintf(fp, "Total reads: %jd\n", (uintmax_t)metrics->total_reads);
	fprintf(fp, "Fraction of reads with missing hairpin: %lf\n",
			(double)metrics->hairpin_missing/metrics->total_reads);
	fprintf(fp, "Fraction of reads with missing palindrome: %lf\n",
			(double)metrics->palindrome_missing/metrics->total_reads);
	fprintf(fp, "Matches=%lf, Damaged=%lf, Mismatches=%lf\n",
			(double)metrics->match_total/total_bases2,
			(double)metrics->damage_total/total_bases2,
			(double)metrics->mm_total/total_bases2);
	fprintf(fp, "Base composition (folded sequences): A=%lf, C=%lf, G=%lf, T=%lf, GC=%lf\n",
			(double)metrics->ntcomp[0]/total_bases,
			(double)metrics->ntcomp[1]/total_bases,
			(double)metrics->ntcomp[2]/total_bases,
			(double)metrics->ntcomp[3]/total_bases,
			(double)(metrics->ntcomp[1]+metrics->ntcomp[2])/total_bases);

	// Nucleotide composition of the hairpin.
	for(i=0; i<4; i++)
		hp_comp[i] = 0;
	for (i=0; i<opt->hlen; i++)
		hp_comp[nt2int[(int)opt->hairpin[i]]]++;

	for (i=0; i<4; i++) {
		uint64_t mm = 0;
		for (j=0; j<4; j++)
			mm += metrics->mm_hairpin[i<<2|j];
		fprintf(fp, "Hairpin mismatch: %c->X = %lg\n",
					"ACGT"[i],
					(double)mm / (hp_comp[i] * n_hairpins));
	}

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

		fprintf(fp, "\nWatson<->Crick [%c']", pairs == &metrics->pairs_l ? '5' : '3');
		for (k=0; k<N_POS; k++)
			fprintf(fp, "\t%d", k);
		fprintf(fp, "\n");

		for (i=0; i<4; i++) {
			for (j=0; j<4; j++) {
				fprintf(fp, "%c<->%c", "ACGT"[i], "ACGT"[j]);
				for (k=0; k<N_POS; k++) {
					fprintf(fp, "\t%lf", (double)(*pairs)[k][i<<2|j] / (metrics->total_reads - metrics->hairpin_missing - metrics->palindrome_missing));
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
	fprintf(stderr, "usage: %s [-b] [-o OUT.fq] [-m FILE] -p SEQ -1 IN1.fq [-2 IN2.fq]\n", argv0);
	fprintf(stderr, " -b                Mark putative methylated cytosines with lowercase 'c',\n"
			"                   or lowercase 'g' for methylation on the other strand.\n");
	fprintf(stderr, " -o OUT.fq         Fastq output file [stdout].\n");
	fprintf(stderr, " -m FILE           Metrics output file [stderr].\n");
	fprintf(stderr, " -p SEQ            The hairpin SEQuence.\n");
	fprintf(stderr, " -1 IN1.fq[.gz]    Fastq input file (paired or singled ended reads).\n");
	fprintf(stderr, " -2 IN2.fq[.gz]    Secondary fastq input file (for paired reads only).\n");
	exit(-1);
}

int
main(int argc, char **argv)
{
	opt_t opt;
	metrics_t metrics;
	int c, ret;
	int i;
	char *fos_fn = NULL, *fom_fn = NULL;

	memset(&opt, '\0', sizeof(opt_t));
	memset(&metrics, '\0', sizeof(metrics_t));
	opt.fos = stdout;
	opt.fom = stderr;

	while ((c = getopt(argc, argv, "bo:m:p:1:2:")) != -1) {
		switch (c) {
			case 'b':
				opt.bisulfite = 1;
				break;
			case 'o':
				fos_fn = optarg;
				break;
			case 'm':
				fom_fn = optarg;
				break;
			case 'p':
				opt.hairpin = optarg;
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

	if (opt.fn1 == NULL) {
		fprintf(stderr, "Error: must specify at least one input file.\n");
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

	if (opt.fn2)
		ret = foldreads_pe(&opt, &metrics);
	else
		ret = foldreads_se(&opt, &metrics);

	free(opt.rhairpin);

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

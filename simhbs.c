/*
 * Simulate various types of reads obtained from hairpin bs-seq libraries.
 *
 * Copyright (c) 2017,2018 Graham Gower <graham.gower@gmail.com>
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
#include <errno.h>
#include <stdint.h>
#include <sys/types.h>
#include <fcntl.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include "fold.h"
#include "kmath.h"
#include "kseq.h"
KSEQ_INIT(int, read);

#define min(a,b) ((a)<(b)?(a):(b))
#define max(a,b) ((a)>(b)?(a):(b))


typedef struct {
	char *refseq;
	uint64_t reflen;

	size_t n_seqs;
	size_t readlen;
	int phred_scale_out;

	char *oprefix;

#define MOL_BSSEQ	0
#define MOL_HAIRPIN	1
#define MOL_PAL_STAR	2
#define MOL_PAL_MODHP	3
	int mol_type;

	char *hairpin, *rhairpin; // hairpin and reverse complement
	char *a1, *a2; // Y adapter sequences
	size_t hlen, alen; // length of hairpin and Y adapters

	double mu, sigma; // mean and std. of log(read length)

	double err_rate; // sequencing error rate
	double bs_rate; // bisulfite conversion rate

	// random state
	krand_t *state;
} opt_t;

static char cmap[] = {['A']='T', ['C']='G', ['G']='C', ['T']='A', ['N']='N', ['n']='N',
				['a']='t', ['c']='g', ['g']='c', ['t']='a'};

double
rlognorm(krand_t *state, double mu, double sigma)
{
	double x = kr_normal(state);
	return exp(mu + x*sigma);
}

int
randbit(opt_t *opt)
{
	static krint64_t x;
	static int bit = 0;

	if (bit == 0) {
		x = kr_rand(opt->state);
		bit = 63;
	} else {
		bit--;
	}

	return (x >> bit) & 0x1;
}

char
randnt(opt_t *opt)
{
	return "ACGT"[(randbit(opt)<<1) + randbit(opt)];
}


// add sequencing errors
void
overlay_err(opt_t *opt, char *seq, size_t n)
{
	int i;
	for (i=0; i<n; i++) {
		double e = kr_drand(opt->state);
		if (e < opt->err_rate)
			seq[i] = randnt(opt);
	}
}

// copy, with bisulfite conversion outside CpG context
void
bs_copy(opt_t *opt, char *dst, char *src, size_t n, int rev, int noctx)
{
	int i;

	if (!rev) {
		// forward strand
		for (i=0; i<n; i++) {
			if (src[i] == 'C' && (noctx || src[i+1] != 'G')
			    && kr_drand(opt->state) < opt->bs_rate) {
				dst[i] = 'T';
			} else {
				dst[i] = src[i];
			}
		}
	} else {
		// reverse
		for (i=0; i<n; i++) {
			if (src[i] == 'G' && (noctx || src[i-1] != 'C')
			    && kr_drand(opt->state) < opt->bs_rate) {
				dst[n-i-1] = 'T';
			} else {
				dst[n-i-1] = cmap[(int)src[i]];
			}
		}
	}
}

char *
sim_hp_bs_mol(opt_t *opt, char *seq, size_t mol_len, size_t *full_len)
{
	int rev;
	char *fullseq;

	*full_len = 2*mol_len + opt->hlen;
	fullseq = malloc(*full_len +1);
	if (fullseq == NULL) {
		perror("malloc");
		fullseq = NULL;
		goto err0;
	}

	rev = randbit(opt);
	bs_copy(opt, fullseq, seq, mol_len, rev, 0);
	memcpy(fullseq+mol_len, opt->hairpin, opt->hlen);
	bs_copy(opt, fullseq+mol_len+opt->hlen, seq, mol_len, !rev, 0);
	fullseq[*full_len] = '\0';

err0:
	return fullseq;
}

char *
sim_bs_mol(opt_t *opt, char *seq, size_t mol_len)
{
	char *fullseq;

	fullseq = malloc(mol_len+1);
	if (fullseq == NULL) {
		perror("malloc");
		fullseq = NULL;
		goto err0;
	}

	bs_copy(opt, fullseq, seq, mol_len, randbit(opt), 0);
	fullseq[mol_len] = '\0';

err0:
	return fullseq;
}

/*
 * Simulate interrupted palindromes.
 *
 * ====s1====---loop---====s2=====---y-adapter---
 * s1 and s2 are reverse complements, like for proper hairpin molecules,
 * but the loop is not a proper hairpin sequence.
 * The loop sequence might be hairpin derived (star==0),
 * or it might be derived from the original molecule (star==1).
 */
char *
sim_bs_pal_mol(opt_t *opt, char *seq, size_t mol_len, size_t *full_len, int star)
{
	int rev;
	char *fullseq;
	char *loopseq;
	size_t pal_len;
	size_t loop_len;

#define MIN_COMPLEMENT 4
#define MIN_LOOP 1
	if (mol_len < MIN_COMPLEMENT+MIN_LOOP) {
		*full_len = mol_len;
		return sim_bs_mol(opt, seq, mol_len);
	}

	if (star) {
		/* Take loop from the originating molecule.
		 * I.e. make palindrome look like Star et al. 2014
		 * doi://10.1371/journal.pone.0089676
		 */
		size_t a, b; // loop start, loop end
		do {
			a = kr_rand(opt->state) % (mol_len-2*MIN_COMPLEMENT)
				+ MIN_COMPLEMENT;
			b = kr_rand(opt->state) % (mol_len-2*MIN_COMPLEMENT)
				+ MIN_COMPLEMENT;
			pal_len = min(a,b);
			loop_len = max(a,b) - pal_len;
		} while (loop_len < MIN_LOOP);
		//fprintf(stderr, "mol_len=%zd, pal_len=%zd, loop_len=%zd\n", mol_len, pal_len, loop_len);
		loopseq = seq + pal_len;
	} else {
		/* Take loop from hairpin sequence.
		 * I.e. polymerase jumped across the hairpin stem, leaving a
		 * partial/modified hairpin between the palindromic segments.
		 *
		 * Original:        /---hairpin---\
		 * ====s1====------/               \
		 * ====s2====----------hairpin-----/
		 *
		 * Polymerase jumped, so loop is short:
		 * ====s1====---------\
		 * ====s2====---------/
		 *
		 * Actually, the empirical cases look a bit weird, and the
		 * polymerase is probably slipping around making shit up, so
		 * the loop is often not recognisable as coming from the
		 * hairpin sequence, other than being GC and methylation rich.
		 * Here, I simplify the loop sequence, by just taking a random
		 * subsequence of the hairpin.
		 */
		pal_len = mol_len;
		do {
			// 10-15 bp
			//loop_len = round(rlognorm(opt->state, 2.5, 0.1));
			// 6-20 bp
			loop_len = round(rlognorm(opt->state, 2.4, 0.25));
		} while (loop_len < MIN_LOOP || loop_len >= opt->hlen);
		loopseq = opt->hairpin + kr_rand(opt->state) % (opt->hlen-loop_len);
	}

	*full_len = 2*pal_len + loop_len;
	fullseq = malloc(*full_len +1);
	if (fullseq == NULL) {
		perror("malloc");
		fullseq = NULL;
		goto err0;
	}

	rev = randbit(opt);
	bs_copy(opt, fullseq, seq, pal_len, rev, 0);
	if (star) {
		// loop comes from original molecule, which will be BS-treated
		bs_copy(opt, fullseq+pal_len, loopseq, loop_len, rev, 0);

		// Star et al. palindromes are Bst filled, so methylation
		// signal is lost on the bottom strand (bases in all contexts
		// are BS converted, including CpGs).
		bs_copy(opt, fullseq+pal_len+loop_len, seq, pal_len, !rev, 1);
	} else {
		// loop comes from hairpin, which is methylated
		memcpy(fullseq+pal_len, loopseq, loop_len);

		bs_copy(opt, fullseq+pal_len+loop_len, seq, pal_len, !rev, 0);
	}

	fullseq[*full_len] = '\0';
err0:
	return fullseq;
}

void
sim_quals(opt_t *opt, char *q, size_t len)
{
	int i;
	int e = opt->err_rate==0 ? 40 : -10 * log10(opt->err_rate);

	if (e > 40)
		e = 40;

	for (i=0; i<len; i++)
		q[i] = e + opt->phred_scale_out;
}


int
simhbs(opt_t *opt)
{
	int ret;
	int n;
	int i;
	uint64_t pos;
	size_t mol_len;

	char *s1, *s2, *q1, *q2;

	char *tmpfn;
	FILE *ofp1, *ofp2;

	s1 = malloc(opt->readlen + 1);
	if (s1 == NULL) {
		perror("malloc");
		ret = -1;
		goto err0;
	}

	s2 = malloc(opt->readlen + 1);
	if (s2 == NULL) {
		perror("malloc");
		ret = -2;
		goto err1;
	}

	q1 = malloc(opt->readlen + 1);
	if (q1 == NULL) {
		perror("malloc");
		ret = -3;
		goto err2;
	}

	q2 = malloc(opt->readlen + 1);
	if (q2 == NULL) {
		perror("malloc");
		ret = -4;
		goto err3;
	}

	s1[opt->readlen] = '\0';
	s2[opt->readlen] = '\0';
	q1[opt->readlen] = '\0';
	q2[opt->readlen] = '\0';

	tmpfn = malloc(strlen(opt->oprefix)+64);
	if (tmpfn == NULL) {
		perror("malloc");
		ret = -5;
		goto err4;
	}
	sprintf(tmpfn, "%s.r1.fq", opt->oprefix);
	ofp1 = fopen(tmpfn, "w");
	if (ofp1 == NULL) {
		free(tmpfn);
		ret = -6;
		goto err4;
	}
	sprintf(tmpfn, "%s.r2.fq", opt->oprefix);
	ofp2 = fopen(tmpfn, "w");
	if (ofp2 == NULL) {
		fclose(ofp1);
		free(tmpfn);
		ret = -7;
		goto err4;
	}
	free(tmpfn);

	for (n=0; n<opt->n_seqs; n++) {
		char *s;
		size_t len;

		do {
			pos = kr_rand(opt->state) % opt->reflen;
			mol_len = round(rlognorm(opt->state, opt->mu, opt->sigma));
		} while (mol_len == 0 || pos < 1 || pos+mol_len >= opt->reflen);

		switch (opt->mol_type) {
			default:
			case MOL_BSSEQ:
				s = sim_bs_mol(opt, opt->refseq+pos,
						mol_len);
				len = mol_len;
				break;
			case MOL_HAIRPIN:
				s = sim_hp_bs_mol(opt, opt->refseq+pos,
						mol_len, &len);
				break;
			case MOL_PAL_STAR:
				s = sim_bs_pal_mol(opt, opt->refseq+pos,
						mol_len, &len, 1);
				break;
			case MOL_PAL_MODHP:
				s = sim_bs_pal_mol(opt, opt->refseq+pos,
						mol_len, &len, 0);
				break;
		}

		if (s == NULL) {
			ret = -8;
			goto err5;
		}

		memcpy(s1, s, min(opt->readlen, len));
		for (i=0; i<min(opt->readlen, len); i++)
			s2[i] = cmap[(int)s[len-i-1]];

		free(s);

		if (opt->readlen > len) {
			memcpy(s1+len, opt->a1, min(opt->readlen-len,opt->alen));
			memcpy(s2+len, opt->a2, min(opt->readlen-len,opt->alen));
			if (opt->readlen > len+opt->alen) {
				for (i=0; i<opt->readlen-len-opt->alen; i++) {
					s1[len+opt->alen+i] = randnt(opt);
					s2[len+opt->alen+i] = randnt(opt);
				}
			}
		}

		overlay_err(opt, s1, opt->readlen);
		overlay_err(opt, s2, opt->readlen);

		sim_quals(opt, q1, opt->readlen);
		sim_quals(opt, q2, opt->readlen);

		fprintf(ofp1, "@simhbs:%d LEN=%zd\n%s\n+\n%s\n", n, mol_len, s1, q1);
		fprintf(ofp2, "@simhbs:%d LEN=%zd\n%s\n+\n%s\n", n, mol_len, s2, q2);
	}

	ret = 0;
err5:
	fclose(ofp2);
	fclose(ofp1);
err4:
	free(q2);
err3:
	free(q1);
err2:
	free(s2);
err1:
	free(s1);
err0:
	return ret;
}

int
load_ref(opt_t *opt, char *fn)
{
	int ret;
	int i;
	int fd;
	int len;
	kseq_t *ks;

	uint64_t reflen;
	char *refseq;
	void *tmp;

	fd = open(fn, O_RDONLY);
	if (fd == -1) {
		fprintf(stderr, "%s: %s\n", fn, strerror(errno));
		ret = -1;
		goto err0;
	}

	ks = kseq_init(fd);

	reflen = 0;
	refseq = NULL;

	for (;;) {
		len = kseq_read(ks);
		if (len < 0)
			break;

		tmp = realloc(refseq, sizeof(*refseq)*(reflen+len));
		if (tmp == NULL) {
			perror("realloc");
			ret = -3;
			goto err1;
		}
		refseq = tmp;

		//memcpy(refseq+reflen, ks->seq.s, len);
		for (i=0; i<len; i++)
			refseq[reflen+i] = toupper(ks->seq.s[i]);

		reflen += len;
	}

	if (len != -1) {
		fprintf(stderr, "%s: unexpected end of file.\n", fn);
		ret = -5;
		goto err1;
	}

	opt->refseq = refseq;
	opt->reflen = reflen;
	ret = 0;
err1:
	if (ret < 0) {
		free(refseq);
	}
	kseq_destroy(ks);
	close(fd);
err0:
	return ret;
}

void
usage(char *argv0, opt_t *opt)
{
	fprintf(stderr, "simhbs v1\n");
	fprintf(stderr, "usage: %s [...] ref.fa\n", argv0);
	fprintf(stderr, " -b FLOAT       Bisulfite conversion rate [%.3f]\n", opt->bs_rate);
	fprintf(stderr, " -e FLOAT       Sequencing error rate [%.3f]\n", opt->err_rate);
	fprintf(stderr, " -h             Simulate hairpin data [%s]\n",
					((char *[]){"no", "yes"})[opt->mol_type==MOL_HAIRPIN]);
	fprintf(stderr, " -p             Simulate palindromes (as in Star et al. 2014) [%s]\n",
					((char *[]){"no", "yes"})[opt->mol_type==MOL_PAL_STAR]);
	fprintf(stderr, " -P             Simulate palindromes (modified hairpin) [%s]\n",
					((char *[]){"no", "yes"})[opt->mol_type==MOL_PAL_MODHP]);
	fprintf(stderr, " -n INT         Number of sequences to simulate [%zd]\n", opt->n_seqs);
	fprintf(stderr, " -o STR         Output prefix for fastq files [%s]\n", opt->oprefix);
	fprintf(stderr, " -r INT         Read length [%zd]\n", opt->readlen);
	fprintf(stderr, " -s FLOAT       std.dev. of LogNormal molecule length [%.3f]\n", opt->sigma);
	fprintf(stderr, " -u FLOAT       mean of LogNormal molecule length [%.3f]\n", opt->mu);
	exit(1);
}

int
main(int argc, char **argv)
{
	opt_t opt;
	int c;
	int ret;
	char *ref_fn;
	int a1len, a2len;

	memset(&opt, '\0', sizeof(opt_t));
	opt.oprefix = "simhbs";
	opt.n_seqs = 100;
	opt.readlen = 150;
	opt.phred_scale_out = 33;
	opt.mol_type = MOL_BSSEQ;

	opt.hairpin = "ACGCCGGCGGCAAGTGAAGCCGCCGGCGT";
	opt.a1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACNNNNNNATCTCGTATGCCGTCTTCTGCTTG";
	opt.a2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT";

	// sane aDNA defaults
	opt.mu = 4.0;
	opt.sigma = 0.25;

	opt.err_rate = 0.001;
	opt.bs_rate = 0.98;

	// seed the random state
	opt.state = kr_srand(time(NULL));

	while ((c = getopt(argc, argv, "b:e:hpPn:o:u:s:")) != -1) {
		switch (c) {
			case 'b':
				opt.bs_rate = atof(optarg);
				if (opt.bs_rate < 0 || opt.bs_rate > 1) {
					fprintf(stderr, "b=`%s' out of range\n", optarg);
					ret = -1;
					goto err0;
				}
				break;
			case 'e':
				opt.err_rate = atof(optarg);
				if (opt.err_rate < 0 || opt.err_rate > 1) {
					fprintf(stderr, "e=`%s' out of range\n", optarg);
					ret = -1;
					goto err0;
				}
				break;
			case 'h':
				if (opt.mol_type != MOL_BSSEQ) {
					fprintf(stderr, "-h, -p, and -P are mutually exclusive\n");
					ret = -2;
					goto err0;
				}
				opt.mol_type = MOL_HAIRPIN;
				break;
			case 'p':
				if (opt.mol_type != MOL_BSSEQ) {
					fprintf(stderr, "-h, -p, and -P are mutually exclusive\n");
					ret = -2;
					goto err0;
				}
				opt.mol_type = MOL_PAL_STAR;
				break;
			case 'P':
				if (opt.mol_type != MOL_BSSEQ) {
					fprintf(stderr, "-h, -p, and -P are mutually exclusive\n");
					ret = -2;
					goto err0;
				}
				opt.mol_type = MOL_PAL_MODHP;
				break;
			case 'n':
				opt.n_seqs = strtoul(optarg, NULL, 0);
				if (opt.n_seqs < 0) {
					fprintf(stderr, "n=`%s' out of range\n", optarg);
					ret = -1;
					goto err0;
				}
				break;
			case 'o':
				opt.oprefix = optarg;
				break;
			case 'r':
				opt.readlen = strtoul(optarg, NULL, 0);
				if (opt.readlen < 10 || opt.readlen > 100000) {
					fprintf(stderr, "r=`%s' out of range\n", optarg);
					ret = -1;
					goto err0;
				}
				break;
			case 's':
				opt.sigma = atof(optarg);
				if (opt.sigma < 0) {
					fprintf(stderr, "s=`%s' out of range\n", optarg);
					ret = -1;
					goto err0;
				}
				break;
			case 'u':
				opt.mu = atof(optarg);
				if (opt.mu < 0) {
					fprintf(stderr, "u=`%s' out of range\n", optarg);
					ret = -1;
					goto err0;
				}
				break;
			default:
				usage(argv[0], &opt);
				break;
		}
	}

	if (argc-optind != 1) {
		usage(argv[0], &opt);
	}

	ref_fn = argv[optind];

	opt.hlen = strlen(opt.hairpin);
	if (opt.hlen < 5 || opt.hlen > 100) {
		fprintf(stderr, "Error: hairpin is too %s (len=%zd).\n",
				opt.hlen<5?"short":"long", opt.hlen);
		usage(argv[0], &opt);
	}

	opt.rhairpin = strdup(opt.hairpin);
	revcomp(opt.rhairpin, opt.hlen);

	a1len = strlen(opt.a1);
	if (a1len < 5 || a1len > 100) {
		fprintf(stderr, "Error: adapter1 is too %s (len=%d).\n",
				a1len<5?"short":"long", a1len);
		usage(argv[0], &opt);
	}

	a2len = strlen(opt.a2);
	if (a2len < 5 || a2len > 100) {
		fprintf(stderr, "Error: adapter2 is too %s (len=%d).\n",
				a2len<5?"short":"long", a2len);
		usage(argv[0], &opt);
	}

	opt.alen = min(a1len, a2len);



	if (load_ref(&opt, ref_fn) < 0) {
		ret = -1;
		goto err0;
	}

	ret = simhbs(&opt);
	free(opt.rhairpin);
	free(opt.refseq);
	free(opt.state);
err0:
	return ret;
}

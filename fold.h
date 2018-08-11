#ifndef _FOLD_H
#define _FOLD_H

#define FOLDREADS_VERSION "12"

// max Q value for output
#define Q_MAX 40

// max Q value for internal caches
#define PHRED_MAX 64
//#define PHRED_SHIFT (log(PHRED_MAX)/log(2.0))
#define PHRED_SHIFT 6

#define AVG_ERR 0.01
#define MAXDIFF_THRES 0.04
#define MAXLEN 1000

int
maxdiff(int l, double err, double thres);

void
revcomp(char *s, size_t len);

void
reverse(char *s, size_t len);

void
str2pvec(const char *s, size_t len, double **pv);

void
find_adapters(const char *s1, const char *q1, size_t len1,
		const char *s2, const char *q2, size_t len2,
		const double *pv1, size_t pv1_len,
		const double *pv2, size_t pv2_len,
		int *h1, int *h2);

double
posterior_error(const char *qvec, int len);

void
match2(const char *s1, const char *q1, size_t len1,
	const char *s2, const char *q2, size_t len2,
	char *s_out, char *q_out,
	int allow_bs);

void
match4(const char *_s1, const char *_q1, size_t len1,
	const char *_s2, const char *_q2, size_t len2,
	const char *_s3, const char *_q3, size_t len3,
	const char *_s4, const char *_q4, size_t len4,
	char *s_out, char *q_out);

void
correct_s1s2(char *s1, char *q1, size_t len1,
		char *s2, char *q2, size_t len2,
		size_t hplen, size_t hppos);

void
clean_quals(const char *s, char *q, size_t len, int phred_scale_in);

#endif // _FOLD_H

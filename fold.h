#ifndef _FOLD_H
#define _FOLD_H

#define AVG_ERR 0.02
#define MAXDIFF_THRES 0.01
#define MAXLEN 1000

int
maxdiff(int l, double err, double thres);

void
revcomp(char *s, size_t len);

void
reverse(char *s, size_t len);

void
find_hp_adapter(const char *s1, size_t len1,
		const char *s2, size_t len2,
		const char *hairpin, const char *rhairpin,
		size_t hlen,
		int *h1, int *h2);

int
match2(const char *s1, const char *q1, size_t len1,
	const char *s2, const char *q2, size_t len2,
	char *s_out, char *q_out,
	int allow_bs,
	int phred_scale_in, int phred_scale_out);

int
match4(const char *_s1, const char *_q1, size_t len1,
	const char *_s2, const char *_q2, size_t len2,
	const char *_s3, const char *_q3, size_t len3,
	const char *_s4, const char *_q4, size_t len4,
	char *s_out, char *q_out,
	int phred_scale_in, int phred_scale_out);

int
correct_s1s2(char *s1, char *q1, size_t len1,
		char *s2, char *q2, size_t len2,
		const char *hairpin, const char *rhairpin, size_t hlen,
		int phred_scale_in, int phred_scale_out);

int
xf2ssqq(char *xf, char **s1, char **s2, char **q1, char **q2);

#endif // _FOLD_H

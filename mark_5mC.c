#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>

#include <htslib/sam.h>

#include "kseq.h"
KSTREAM_INIT(int, read, 16384);

typedef struct {
	char *bam_fn;
	char *bed_fn;
	int min_mapq;
	int min_baseq;
} opt_t;

typedef struct _bed_entry {
	uint32_t tid;
	uint32_t start;
	uint32_t end;
	struct _bed_entry *next;
} bed_t;

typedef struct {
	opt_t *opt;
	samFile *bam_fp;
	bam_hdr_t *bam_hdr;
	hts_idx_t *bam_idx;
	hts_itr_t *bam_iter;
	bed_t *bed_head;
} bam_aux_t;


static char cmap[] = {['A']='T', ['C']='G', ['G']='C', ['T']='A', ['N']='N', ['n']='N',
				['a']='t', ['c']='g', ['g']='c', ['t']='a'};

/*
 * Load .bed loci into a singly linked list. Takes the chromosome name from
 * the first column, matching it to one from the bam header
 */
bed_t *
load_bed(const char *bed_fn, bam_hdr_t *bam_hdr)
{
	int fd;
	int dret;
	int lineno = 0;
	kstream_t *ks;
	kstring_t *str;
	bed_t *head, *tail, *tmp;

	head = tail = tmp = NULL;

	fd = open(bed_fn, O_RDONLY);
	if (fd == -1) {
		fprintf(stderr, "open: %s: %s\n", bed_fn, strerror(errno));
		goto err0;
	}

	str = calloc(1, sizeof(kstring_t));
	if (str == NULL) {
		perror("load_bed: calloc");
		goto err1;
	}

	ks = ks_init(fd);

	while (ks_getuntil(ks, KS_SEP_SPACE, str, &dret) >= 0) {
		int tid, start, end;

		lineno++;

		tid = bam_name2id(bam_hdr, str->s);
		if (tid < 0) {
			fprintf(stderr, "%s: line %d: couldn't find RNAME `%s' in sam header\n",
					bed_fn, lineno, str->s);
			goto err2;
		}

		if (dret == '\n') {
			fprintf(stderr, "%s: line %d: expected 3 columns, only got 1\n",
					bed_fn, lineno);
			goto err2;
		}

		if (ks_getuntil(ks, KS_SEP_SPACE, str, &dret) < 0)
			break;

		if (dret == '\n') {
			fprintf(stderr, "%s: line %d: expected 3 columns, only got 2\n",
					bed_fn, lineno);
			goto err2;
		}

		start = atoi(str->s);
		if (ks_getuntil(ks, KS_SEP_SPACE, str, &dret) < 0)
			break;
		end = atoi(str->s);

		// got 3 columns
		tmp = malloc(sizeof(bed_t));
		if (tmp == NULL) {
			perror("load_bed: malloc");
			goto err2;
		}

		tmp->tid = tid;
		tmp->start = start;
		tmp->end = end;
		tmp->next = NULL;

		if (tail == NULL) {
			head = tail = tmp;
		} else {
			tail->next = tmp;
			tail = tmp;
		}

		if (dret != '\n') {
			// consume the rest of the line
			while ((dret = ks_getc(ks)) > 0 && dret != '\n')
				;
		}
		if (dret < 0)
			break;
	}

	goto done;
err2:
	for (tmp=head; tmp!=NULL; ) {
		bed_t *x = tmp;
		tmp = tmp->next;
		free(x);
	}
	head = NULL;
done:
	ks_destroy(ks);
	free(str->s);
	free(str);
err1:
	close(fd);
err0:
	return head;
}


/*
 * Get next alignment.
 */
static int
next_aln(void *data, bam1_t *b)
{
	bam_aux_t *bat = data;
	int ret;

	while (1) {
		if (bat->opt->bed_fn) {
			if (bat->bed_head == NULL) {
				ret = -1;
				break;
			}
			if (bat->bam_iter == NULL) {
				bed_t *x = bat->bed_head;
				bat->bed_head = x->next;
				bat->bam_iter = sam_itr_queryi(bat->bam_idx, x->tid, x->start, x->end);
				free(x);
				if (bat->bam_iter == NULL) {
					ret = -2;
					break;
				}
			}
			ret = sam_itr_next(bat->bam_fp, bat->bam_iter, b);
			if (ret < 0) {
				// iterator exhausted
				hts_itr_destroy(bat->bam_iter);
				bat->bam_iter = NULL;
				continue;
			}
		} else {
			ret = sam_read1(bat->bam_fp, bat->bam_hdr, b);
			if (ret < 0)
				break;
		}

		if (b->core.flag & (BAM_FUNMAP | BAM_FQCFAIL | BAM_FDUP) ||
				b->core.qual < bat->opt->min_mapq)
			// skip these
			continue;

		break;
	}

	return ret;
}

int
mark_5mC(opt_t *opt)
{
	bam_aux_t bat;
	bam_plp_t plpiter;
	const bam_pileup1_t *plp;
	int tid, pos, n;
	int i;
	int ret;
	int depth;
	int nt[256] = {0,};

	bat.opt = opt;

	bat.bam_fp = sam_open(opt->bam_fn, "r");
	if (bat.bam_fp == NULL) {
		fprintf(stderr, "bam_open: %s: %s\n", opt->bam_fn, strerror(errno));
		ret = -1;
		goto err0;
	}

	bat.bam_hdr = sam_hdr_read(bat.bam_fp);
	if (bat.bam_hdr == NULL) {
		fprintf(stderr, "%s: couldn't read header\n", opt->bam_fn);
		ret = -2;
		goto err1;
	}

	bat.bed_head = NULL;
	if (opt->bed_fn) {
		bat.bam_idx = sam_index_load(bat.bam_fp, opt->bam_fn);
		if (bat.bam_idx == NULL) {
			fprintf(stderr, "Couldn't load index for %s\n", opt->bam_fn);
			ret = -3;
			goto err2;
		}
		bat.bed_head = load_bed(opt->bed_fn, bat.bam_hdr);
		if (bat.bed_head == NULL) {
			ret = -4;
			goto err3;
		}
	}

	bat.bam_iter = NULL;
	plpiter = bam_plp_init(next_aln, &bat);

	//printf("track type=\"bedGraph\" color=200,100,0 altColor=0,100,200 description=\"proportion of methylated Cytosines\"\n");

	while ((plp = bam_plp_auto(plpiter, &tid, &pos, &n)) != 0) {

		if (bat.bam_iter && (pos < bat.bam_iter->beg || pos > bat.bam_iter->end))
			continue;

		nt['c'] = nt['C'] = nt['G'] = nt['g'] = 0;
		for (i=0; i<n; i++) {
			const bam_pileup1_t *p = plp + i;
			uint8_t *xf;
			int c;

			if (p->is_del || p->is_refskip)
				continue;
			
			if (bam_get_qual(p->b)[p->qpos] < opt->min_baseq)
				continue;

			xf = bam_aux_get(p->b, "XF");
			if (xf == NULL)
				continue;

			if (bam_is_rev(p->b))
				c = cmap[xf[p->b->core.l_qseq - p->qpos]];
			else
				c = xf[p->qpos];

			nt[c]++;
		}

		if ((depth = nt['c'] + nt['C'])) {
			printf("%s\t%d\t%d\t+\t%d\t%d\t%d\n",
					bat.bam_hdr->target_name[tid],
					pos,
					pos+1,
					(100*nt['c'] + (depth/2)) / depth,
					nt['c'],
					depth
					);
		}
		if ((depth = nt['g'] + nt['G'])) {
			printf("%s\t%d\t%d\t-\t%d\t%d\t%d\n",
					bat.bam_hdr->target_name[tid],
					pos,
					pos+1,
					(100*nt['g'] + (depth/2)) / depth,
					nt['g'],
					depth
					);
		}
	}

	ret = 0;

	bam_plp_destroy(plpiter);

	if (bat.bed_head) {
		bed_t *x, *xx;
		for (x=bat.bed_head; x!=NULL; ) {
			xx = x;
			x = x->next;
			free(xx);
		}
	}
err3:
	hts_idx_destroy(bat.bam_idx);
err2:
	bam_hdr_destroy(bat.bam_hdr);
err1:
	sam_close(bat.bam_fp);
err0:
	return ret;
}

void
usage(char *argv0)
{
	fprintf(stderr, "usage: %s [-b regions.bed] in.bam\n", argv0);
	exit(1);
}

int
main(int argc, char **argv)
{
	opt_t opt;
	int c;

	memset(&opt, 0, sizeof(opt_t));
	opt.min_mapq = 25;
	opt.min_baseq = 10;

	while ((c = getopt(argc, argv, "b:")) != -1) {
		switch (c) {
			case 'b':
				opt.bed_fn = optarg;
				break;
			default:
				usage(argv[0]);
		}
	}

	if (argc-optind != 1) {
		usage(argv[0]);
	}

	opt.bam_fn = argv[optind];

	if (mark_5mC(&opt) < 0) {
		return 1;
	}

	return 0;
}

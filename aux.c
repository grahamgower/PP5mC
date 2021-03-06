/*
 * Copyright (c) 2018 Graham Gower <graham.gower@gmail.com>
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
#include <htslib/sam.h>
	
int
aux2rrqqhp(bam1_t *b, char **r1, char **r2, char **q1, char **q2, char **hp)
{
	uint8_t *aux;

	aux = bam_aux_get(b, "r1");
	if (aux == NULL) {
		fprintf(stderr, "%s: missing auxiliary field r1:Z\n",
				bam_get_qname(b));
		return -1;
	}
	*r1 = bam_aux2Z(aux);
	if (*r1 == NULL) {
		fprintf(stderr, "%s: invalid auxiliary field r1, not r1:Z\n",
				bam_get_qname(b));
		return -1;
	}

	aux = bam_aux_get(b, "r2");
	if (aux == NULL) {
		fprintf(stderr, "%s: missing auxiliary field r2:Z\n",
				bam_get_qname(b));
		return -1;
	}
	*r2 = bam_aux2Z(aux);
	if (*r2 == NULL) {
		fprintf(stderr, "%s: invalid auxiliary field r2, not r2:Z\n",
				bam_get_qname(b));
		return -1;
	}

	aux = bam_aux_get(b, "q1");
	if (aux == NULL) {
		fprintf(stderr, "%s: missing auxiliary field q1:Z\n",
				bam_get_qname(b));
		return -1;
	}
	*q1 = bam_aux2Z(aux);
	if (*q1 == NULL) {
		fprintf(stderr, "%s: invalid auxiliary field q1, not q1:Z\n",
				bam_get_qname(b));
		return -1;
	}

	aux = bam_aux_get(b, "q2");
	if (aux == NULL) {
		fprintf(stderr, "%s: missing auxiliary field q2:Z\n",
				bam_get_qname(b));
		return -1;
	}
	*q2 = bam_aux2Z(aux);
	if (*q2 == NULL) {
		fprintf(stderr, "%s: invalid auxiliary field q2, not q2:Z\n",
				bam_get_qname(b));
		return -1;
	}

	// hp is optional
	*hp = NULL;
	aux = bam_aux_get(b, "hp");
	if (aux != NULL) {
		*hp = bam_aux2Z(aux);
		if (*hp == NULL) {
			fprintf(stderr, "%s: invalid auxiliary field hp, not hp:Z\n",
					bam_get_qname(b));
			return -1;
		}
	}

	return 0;
}


/*-
 * Copyright 1997, 1998-2003 John-Mark Gurney.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 * 1. Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHOR AND CONTRIBUTORS ``AS IS'' AND
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHOR OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
 * SUCH DAMAGE.
 *
 *	$Id: fib.h 300 2011-02-03 18:46:33Z pawel@BIOEXPLORATORIUM.PL $
 *
 */

#ifndef __FIB_H__
#define __FIB_H__

#ifdef __cplusplus
extern "C" {
#endif

struct fibheap {
	int	(*fh_cmp_fnct)(void *, void *);
	int	fh_n;
	int	fh_Dl;
	struct	fibheap_el **fh_cons;
	struct	fibheap_el *fh_min;
	struct	fibheap_el *fh_root;
	void	*fh_neginf;
	int	fh_keys;
	int size;
};


struct fibheap_el {
	int	fhe_degree;
	int	fhe_mark;
	struct	fibheap_el *fhe_p;
	struct	fibheap_el *fhe_child;
	struct	fibheap_el *fhe_left;
	struct	fibheap_el *fhe_right;
	unsigned int fhe_key;
	void	*fhe_data;
};
typedef int (*voidcmp)(void *, void *);

/* functions for key heaps */
void fh_makekeyheap(struct fibheap *, int size);
void fh_insertkey(struct fibheap *, unsigned int, void *, struct fibheap_el *);
unsigned int fh_minkey(struct fibheap *);
unsigned int fh_replacekey(struct fibheap *, struct fibheap_el *, unsigned int);
void *fh_replacekeydata(struct fibheap *, struct fibheap_el *, unsigned int, void *);
void fh_initheap(struct fibheap *, int size);

/* functions for void * heaps */
void fh_makeheap(struct fibheap *, int size);
voidcmp fh_setcmp(struct fibheap *, voidcmp);
void *fh_setneginf(struct fibheap *, void *);
void fh_insert(struct fibheap *, void *, struct fibheap_el *);

/* shared functions */
void *fh_extractmin(struct fibheap *);
void *fh_min(struct fibheap *);
void *fh_replacedata(struct fibheap *, struct fibheap_el *, void *);
void *fh_delete(struct fibheap *, struct fibheap_el *);
void fh_deleteheap(struct fibheap *);
struct fibheap *fh_union(struct fibheap *, struct fibheap *);

#ifdef __cplusplus
};
#endif
#endif /* _FIB_H_ */

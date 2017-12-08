/*
 * Copyright 2003-2010, 2017 Pawel Daniluk
 * 
 * This file is part of PyDesc.
 * 
 * PyDesc is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * PyDesc is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with PyDesc.  If not, see <http://www.gnu.org/licenses/>.
 */
 


#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include"bitops.h"
#include"random.h"

int bits_32(unsigned int v) 
{
	int c;

	v = v - ((v >> 1) & 0x55555555);                    // reuse input as temporary
	v = (v & 0x33333333) + ((v >> 2) & 0x33333333);     // temp
	v =((v + (v >> 4)) & 0x0F0F0F0F); // count
	c = (int) (((v * 0x1010101) >> 24)); // count

	return c;
}

int bits_64(unsigned long long int v) 
{
	int c;

	v = v - ((v >> 1) & 0x5555555555555555ULL);                    // reuse input as temporary
	v = (v & 0x3333333333333333ULL) + ((v >> 2) & 0x3333333333333333ULL);     // temp
	v = ((v + (v >> 4)) & 0x0F0F0F0F0F0F0F0FULL); // count
	c = (int) (((v * 0x0101010101010101ULL) >> 56)); // count

	return c;
}

int bits(t_maskcell m)
{
#ifdef __x86_64__
	return bits_64(m);
#else
	return bits_32(m);
#endif
}

t_bitmask *mask_alloc(int len)
{
	t_bitmask *res=calloc(1, sizeof(t_bitmask));
    res->len = len;
    res->mask = calloc(MASK_LENGTH(len), sizeof(t_maskcell));

	return res;
}

void mask_free(t_bitmask *mask)
{
    free(mask->mask);
	free(mask);
}

void mask_zeroall(t_bitmask *bitmask)
{
    int len=bitmask->len;
    t_maskcell *mask=bitmask->mask;

	memset(mask, 0, MASK_LENGTH(len)*sizeof(t_maskcell));
}

void mask_setall(t_bitmask *bitmask)
{
    int len=bitmask->len;
    t_maskcell *mask=bitmask->mask;

	memset(mask, 255, (MASK_LENGTH(len)-1)*sizeof(t_maskcell));
	for(int i=MASK_CELL_SIZE*(MASK_LENGTH(len)-1); i<len; i++) {
		SET_BIT(bitmask, i);
	}
}

int mask_size(t_bitmask *bitmask)
{
	int res=0;
    int len=bitmask->len;
    t_maskcell *mask=bitmask->mask;

	for(int i=0; i<=len/MASK_CELL_SIZE; i++) res+=bits(mask[i]);

	return res;
}

int mask_empty(t_bitmask *bitmask)
{
    int len=bitmask->len;
    t_maskcell *mask=bitmask->mask;

	for(int i=0; i<=len/MASK_CELL_SIZE; i++) if(mask[i]) return 0;

	return 1;
}

void mask_not(t_bitmask *bitmask)
{
    int len=bitmask->len;
    t_maskcell *mask=bitmask->mask;

	for(int i=0; i<=len/MASK_CELL_SIZE; i++) mask[i]=~mask[i];
}

void mask_or(t_bitmask *bitmask, t_bitmask *bitmask1)
{
    int len=bitmask->len;
    t_maskcell *mask=bitmask->mask;
    t_maskcell *mask1=bitmask1->mask;

	for(int i=0; i<=len/MASK_CELL_SIZE; i++) mask[i]|=mask1[i];
}

void mask_and(t_bitmask *bitmask, t_bitmask *bitmask1)
{
    int len=bitmask->len;
    t_maskcell *mask=bitmask->mask;
    t_maskcell *mask1=bitmask1->mask;

	for(int i=0; i<=len/MASK_CELL_SIZE; i++) mask[i]&=mask1[i];
}

int mask_and_empty(t_bitmask *bitmask, t_bitmask *bitmask1)
{
    int len=bitmask->len;
    t_maskcell *mask=bitmask->mask;
    t_maskcell *mask1=bitmask1->mask;

	for(int i=0; i<=len/MASK_CELL_SIZE; i++) if(mask[i]&mask1[i]) return 0;

	return 1;
}

void mask_sub(t_bitmask *bitmask, t_bitmask *bitmask1)
{
    int len=bitmask->len;
    t_maskcell *mask=bitmask->mask;
    t_maskcell *mask1=bitmask1->mask;

	for(int i=0; i<=len/MASK_CELL_SIZE; i++) mask[i]&=~mask1[i];
}

int mask_sub_empty(t_bitmask *bitmask, t_bitmask *bitmask1)
{
    int len=bitmask->len;
    t_maskcell *mask=bitmask->mask;
    t_maskcell *mask1=bitmask1->mask;

	for(int i=0; i<=len/MASK_CELL_SIZE; i++) if(mask[i]&~mask1[i]) return 0;

	return 1;
}

void mask_cpy(t_bitmask *bitmask, t_bitmask *bitmask1)
{
    int len=bitmask->len;
    t_maskcell *mask=bitmask->mask;
    t_maskcell *mask1=bitmask1->mask;

	memcpy(mask, mask1, MASK_LENGTH(len)*sizeof(t_maskcell));
}

int mask_eq(t_bitmask *bitmask1, t_bitmask *bitmask2)
{
    int len=bitmask1->len;
    t_maskcell *mask1=bitmask1->mask;
    t_maskcell *mask2=bitmask2->mask;

	for(int i=0; i<=len/MASK_CELL_SIZE; i++) {
		if(mask1[i] != mask2[i]) return 0;
	}

	return 1;
}

int mask_comp(t_bitmask *bitmask1, t_bitmask *bitmask2)
{
    int len=bitmask1->len;
    t_maskcell *mask1=bitmask1->mask;
    t_maskcell *mask2=bitmask2->mask;

	int res=MASK_SUBSET | MASK_SUPERSET | MASK_EQUAL;

	for(int i=0; i<=len/MASK_CELL_SIZE; i++) {
		if(mask1[i] & ~mask2[i]) res=res&~(MASK_EQUAL | MASK_SUBSET);
		if(~mask1[i] & mask2[i]) res=res&~(MASK_EQUAL | MASK_SUPERSET);
	}


	return res;
}

int mask_to_list(int *l, int n, t_bitmask *bitmask)
{
    int len=bitmask->len;

	int p=0;

	for(int i=0; i<len; i++) {
		if(BIT_TEST(bitmask, i)) {
			if(p>=n) return -1;
			l[p++]=i;
		}
	}

	return p;
}

int list_to_mask(t_bitmask *bitmask, int *l, int n)
{
	mask_zeroall(bitmask);

	int p=0;

	for(int i=0; i<n; i++) {
		if(l[i]>=0) {
			SET_BIT(bitmask, l[i]);
			p++;
		}
	}

	return p;
}


t_bitmask *random_mask(int len)
{
	t_bitmask *res=mask_alloc(len);

	for(int i=0; i<len; i++) if(int_ran(0,1)) SET_BIT(res, i);

	return res;
}

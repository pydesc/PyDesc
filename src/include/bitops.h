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
 


#ifndef __BITOPS_H__
#define __BITOPS_H__

typedef unsigned long long int t_maskcell;
typedef struct {
    t_maskcell *mask;
    int len;
} t_bitmask;

#define BITFIELD(x) ((t_maskcell)1<<(x))

#define MASK_CELL_SIZE (sizeof(t_maskcell)*8)

#define BIT_TEST(x,n) ({ typeof(n) aa=n;  ((x->mask)[(aa)/MASK_CELL_SIZE]&BITFIELD((aa)%MASK_CELL_SIZE)); })
#define SET_BIT(x,n) ({ typeof(n) aa=n; ((x->mask)[(aa)/MASK_CELL_SIZE]|=BITFIELD((aa)%MASK_CELL_SIZE)); })
#define UNSET_BIT(x,n) ({ typeof(n) aa=n; ((x->mask)[(aa)/MASK_CELL_SIZE]&=~BITFIELD((aa)%MASK_CELL_SIZE)); })

#define MASK_LENGTH(len) ((len)/MASK_CELL_SIZE+1)

#define MASK_SUBSET 1
#define MASK_SUPERSET 2
#define MASK_EQUAL 4
	

//void fill_bit_count(void);

int bits(t_maskcell i);

t_bitmask *mask_alloc(int len);
void mask_free(t_bitmask *mask);
void mask_zeroall(t_bitmask *mask);
void mask_setall(t_bitmask *mask);
int mask_size(t_bitmask *mask);
int mask_comp(t_bitmask *mask1, t_bitmask *mask2);
void mask_or(t_bitmask *mask, t_bitmask *mask1);
void mask_and(t_bitmask *mask, t_bitmask *mask1);
void mask_sub(t_bitmask *mask, t_bitmask *mask1);
void mask_cpy(t_bitmask *mask, t_bitmask *mask1);
int mask_eq(t_bitmask *mask1, t_bitmask *mask2);
void mask_not(t_bitmask *mask);
int mask_empty(t_bitmask *mask);
int mask_and_empty(t_bitmask *mask, t_bitmask *mask1);
int mask_sub_empty(t_bitmask *mask, t_bitmask *mask1);

int mask_to_list(int *l, int n, t_bitmask *mask);
int list_to_mask(t_bitmask *mask, int *l, int n);

t_bitmask *random_mask(int len);

#endif


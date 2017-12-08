/*
 * Copyright 2017 Pawel Daniluk
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
 
#ifndef __CSTRUCTURES_H__
#define __CSTRUCTURES_H__

#include"maps.h"

#ifdef INCLUDE_SIZEOF

#define CYDESC_STRUCT(__n, __fields) \
typedef struct __n __fields __n; \
int sizeof_ ## __n(void) {return sizeof(__n);};

#else

#define CYDESC_STRUCT(__n, __fields) typedef struct __n __fields __n;

#endif

CYDESC_STRUCT(CSeg, {
  // defined in PyDesc:
    int start;
    int end;

  // undefined in PyDesc:
    int c_start;
    int c_end;
})


CYDESC_STRUCT(CPoint, {
  // defined in PyDesc:
	float x;
	float y;
	float z;

  // undefined in PyDesc:
	float crap;
})

struct CMer;

CYDESC_STRUCT(CMer, {
  // defined in PyDesc:
    int ind;
    int type;

    int next_ind;

    int n_points;
    CPoint *points;

    char **point_names;
    char *type_name;

  // undefined in PyDesc:
    int c_ind;
    int next_c_ind;
    struct CMer *next;
})


CYDESC_STRUCT(CStructure, {
  // defined in PyDesc:
    char *name;
    int n_monomers; // number of mers
    CMer *monomers;
    int n_segs; // number of segments
    CSeg *segs;
    int (*adjusted_number_p)(int, int);

  // undefined in PyDesc:
    t_map *ind_to_c_map;
})

CYDESC_STRUCT(CContact, {
  // defined in PyDesc:
    int mer1;
    int mer2;
    int val;
    int type;

  // undefined in PyDesc:
    int mer1_c_ind;
    int mer2_c_ind;
    int mer1_used;
    int mer2_used;
    int pos;
})

CYDESC_STRUCT(CContactMap, {
  // defined in PyDesc:
    CStructure *structure;
    int n_contacts;
    CContact *contacts;

  // undefined in PyDesc:
    int n_used;
    t_map *c_ind_to_used_map;
    int **contact_array;
})

CYDESC_STRUCT(CElement, {
  // defined in PyDesc:
    int center;
    int start;
    int end;
    int type;
    int status;
    int optional;

  // undefined in PyDesc:
    int center_c_ind;
    int start_c_ind;
    int end_c_ind;
    int center_used;
})

CYDESC_STRUCT(CDescriptor, {
  // defined in PyDesc:
    CStructure *structure;
    CContactMap *contact_map;
    int n_elements;
    CElement *elements;
    int central_element;

  // undefined in PyDesc:

    int central_element_c_ind;
    int central_element_used;
})

CStructure* createCStructure(int* seg, int* loc, float* coo);
void printCStructure(CStructure* structure);

#endif

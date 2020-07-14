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
 
#include<stdlib.h>
#include<stdio.h>
#include<string.h>

#include"simple_macros.h"

#include"cstructures.h"
#include"maps.h"

void free_(void *ptr)
/*
 * According to ctypes documentation on some systems each DLL may have its own heap.
 * In such a case memory has to be freed by DLL which has allocated it.
 * In this way we expose free for this particular purpose.
 */
{
    free(ptr);
}


int CPointTest_array(int n, float *xl, float *yl, float *zl, CPoint *point_arr)
{
    for(int i=0; i<n; i++) {
        if(xl[i]!=point_arr[i].x) return 1;
        if(yl[i]!=point_arr[i].y) return 1;
        if(zl[i]!=point_arr[i].z) return 1;
    }
    return 0;
}


CPoint *CPointTest_array_cpy(int n, CPoint *point_arr)
{
    CPoint *res = calloc(n, sizeof(CPoint));

    memcpy(res, point_arr, n*sizeof(CPoint));

    return res;
}

CMer *CMerTest_array_cpy(int n, CMer *mer_arr)
{
    CMer *res = calloc(n, sizeof(CMer));

    memcpy(res, mer_arr, n*sizeof(CMer));

    return res;
}

CSeg *CSegTest_array_cpy(int n, CSeg *seg_arr)
{
    CSeg *res = calloc(n, sizeof(CSeg));

    memcpy(res, seg_arr, n*sizeof(CSeg));

    return res;
}

int CStructureTest_indices(CStructure *str)
{
    for(int i=0; i<str->n_monomers; i++) {
        CMer *mer=&str->monomers[i];
        if(mer->c_ind != i) return 1;
        if(mer->next_ind>0) {
            if(mer->next->ind != mer->next_ind) return 1;
            if(mer->next != &(str->monomers[mer->next_c_ind])) return 1;
        } else {
            if(mer->next) return 1;
        }

        if(mer->c_ind != map_map(str->ind_to_c_map, mer->ind)) return 1;
        if(mer->ind != map_unmap(str->ind_to_c_map, mer->c_ind)) return 1;
	}

    return 0;
}

int CStructureTest_adjusted_number(CStructure *str, int start, int end)
{
    return str->adjusted_number_p(start, end);
}

int CContactMapTest_indices(CContactMap *cmap, int (*callback)(int, int, int))
{
    int cnt=0;
    for(int i=0; i<cmap->n_used; i++) {
        for(int j=0; j<cmap->n_used; j++) {
            if(cmap->contact_array[i][j]>=0) {
                int c_ind1=map_unmap(cmap->c_ind_to_used_map, i);
                int c_ind2=map_unmap(cmap->c_ind_to_used_map, j);

                CContact *con=&cmap->contacts[cmap->contact_array[i][j]];

                int ind1=map_unmap(cmap->structure->ind_to_c_map, c_ind1);
                int ind2=map_unmap(cmap->structure->ind_to_c_map, c_ind2);
                int cb_res=callback(ind1, ind2, con->val);

                if(cb_res) return 1;
                cnt++;
            }
        }
    }

    if(cnt!=cmap->n_contacts) {
        P_S("cnt!=cmap->n_contacts: "); P_INT(cnt); P_INT(cmap->n_contacts); P_NL;
        return 1;
    }

    return 0;
}

int CDescriptorTest_element_map(CDescriptor *cdesc)
{
    int used_unmap(int used_ind) {
        return map_unmap(cdesc->contact_map->c_ind_to_used_map, used_ind);
    }

    for(int i=0; i<cdesc->n_elements; i++) {
        CElement *el=&(cdesc->elements[i]);
        if(el->center_c_ind != used_unmap(el->center_used)) return 1;
        if(el->center_used != i) return 1;
    }

    return 0;
}


#define INCLUDE_SIZEOF
CYDESC_STRUCT(CTest, {
  // defined in PyDesc:
    int f1;
    int f2;

  // undefined in PyDesc:
    int hf1;
    int hf2;
})
#undef INCLUDE_SIZEOF

int init_CTest_calls=0;
int del_CTest_calls=0;

void init_CTest(CTest *test)
{
    test->hf1=2*test->f1;
    test->hf2=2*test->f2;

    init_CTest_calls++;
}

void del_CTest(CTest *test)
{
    del_CTest_calls++;
}

void init_CTest1(CTest *test)
{
    init_CTest(test);
}

void del_CTest1(CTest *test)
{
    del_CTest(test);
}

int get_init_CTest_calls(void)
{
    return init_CTest_calls;
}

int get_del_CTest_calls(void)
{
    return del_CTest_calls;
}

void print_CPoint(CPoint *point) {
    printf("%f %f %f\n", point->x, point->y, point->z);
}

CPoint *new_CPoint(float x, float y, float z) {
    CPoint *res=malloc(sizeof(CPoint));

    res->x=x;
    res->y=y;
    res->z=z;

    printf("point created: ");
    print_CPoint(res);

    return res;
}


void print_CMer(CMer *mer)
{
    P_INT(mer->ind); P_NL;
    P_INT(mer->type); P_NL;
    P_S(mer->type_name); P_NL;

    P_INT(mer->next_ind); P_NL;

    P_INT(mer->n_points); P_NL;
    for(int i=0;i<mer->n_points;i++) {
        P_S(mer->point_names[i]);
        print_CPoint(&(mer->points[i]));
    }
}

void print_CStructure(CStructure* str)
{
    P_INT(str->n_monomers); P_NL;

    for(int i=0; i<str->n_monomers; i++) {
        print_CMer(&str->monomers[i]);
    }

    P_INT(str->n_segs); P_NL;
    for(int i=0; i<str->n_segs; i++) {
        P_INT(str->segs[i].start);
        P_INT(str->segs[i].end);
        P_NL;
    }
}




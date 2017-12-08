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
#include "stdio.h"

#include"simple_macros.h"
#include"maps.h"
#include"arrays.h"
#include"bitops.h"

#define INCLUDE_SIZEOF
#include "cstructures.h"
#undef INCLUDE_SIZEOF

void init_CStructure(CStructure *str)
{
    int max_ind=0;

    for(int i=0; i<str->n_monomers; i++) {
        CMer *mer=&str->monomers[i];
        mer->c_ind=i;

        max_ind=MAX(max_ind, mer->ind);
    }

    str->ind_to_c_map=map_alloc(max_ind+1, str->n_monomers);

    for(int i=0; i<str->n_monomers; i++) {
        CMer *mer=&str->monomers[i];
        map_set(str->ind_to_c_map, mer->ind, i);
        if(mer->next_ind>0) {
            CMer *m1=mer;
            if(mer->next_ind>mer->ind) {
                while(m1->ind != mer->next_ind) m1++;
            } else {
                while(m1->ind != mer->next_ind) m1--;
            }
            mer->next=m1;
            mer->next_c_ind=m1->c_ind;
        }
	}

    for(int i=0; i<str->n_segs; i++) {
        str->segs[i].c_start = map_map(str->ind_to_c_map, str->segs[i].start);
        str->segs[i].c_end = map_map(str->ind_to_c_map, str->segs[i].end);
    }
}

void del_CStructure(CStructure *str)
{
    map_free(str->ind_to_c_map);
}

void init_CContactMap(CContactMap *contact_map)
{
    int n_monomers=contact_map->structure->n_monomers;
    t_bitmask *used=mask_alloc(n_monomers);

    for(int i=0; i<contact_map->n_contacts; i++) {
        contact_map->contacts[i].mer1_c_ind=map_map(contact_map->structure->ind_to_c_map, contact_map->contacts[i].mer1);
        contact_map->contacts[i].mer2_c_ind=map_map(contact_map->structure->ind_to_c_map, contact_map->contacts[i].mer2);

        contact_map->contacts[i].pos=i;

        SET_BIT(used, contact_map->contacts[i].mer1_c_ind);
        SET_BIT(used, contact_map->contacts[i].mer2_c_ind);
    }

    int n_used = mask_size(used);
    contact_map->n_used = n_used;
    contact_map->c_ind_to_used_map = map_alloc(contact_map->structure->n_monomers, n_used);

    int num=0;

    for(int i=0; i<n_monomers; i++) if(BIT_TEST(used, i)) {
        map_set(contact_map->c_ind_to_used_map, i, num++);
    }

    mask_free(used);

    contact_map->contact_array = (int **)array_alloc(n_used, n_used, sizeof(int));

    array_set((void **)contact_map->contact_array, -1);

    for(int i=0; i<contact_map->n_contacts; i++) {
        CContact *contact = &(contact_map->contacts[i]);

        contact->mer1_used = map_map(contact_map->c_ind_to_used_map, contact->mer1_c_ind);
        contact->mer2_used = map_map(contact_map->c_ind_to_used_map, contact->mer2_c_ind);
        contact_map->contact_array[contact->mer1_used][contact->mer2_used] = contact->pos;
    }
}

void del_CContactMap(CContactMap *contact_map)
{
    map_free(contact_map->c_ind_to_used_map);
    array_free((void **)contact_map->contact_array);
}

void init_CDescriptor(CDescriptor *descriptor)
{
    CStructure *structure=descriptor->structure;

    for(int i=0; i<descriptor->n_elements; i++) {
        CElement *el=&(descriptor->elements[i]);

        el->center_c_ind=map_map(structure->ind_to_c_map, el->center);
        el->start_c_ind=map_map(structure->ind_to_c_map, el->start);
        el->end_c_ind=map_map(structure->ind_to_c_map, el->end);
        el->center_used=map_map(descriptor->contact_map->c_ind_to_used_map, el->center_c_ind);
    }

    descriptor->central_element_c_ind=map_map(structure->ind_to_c_map, descriptor->central_element);
    descriptor->central_element_used=map_map(descriptor->contact_map->c_ind_to_used_map, descriptor->central_element_c_ind);
}

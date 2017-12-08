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
 
#include<math.h>

#undef COMPDESC_DEBUG   // Terribly ugly hack to ensure that t_compdesc_debug is included from libdesc.

#include"simple_macros.h"
#include"compdesc.h"
#include"arrays.h"

#define COMPDESC_DEBUG
#include<desc.h>

t_desc *desc1;
t_desc *desc2;

int set_desc(char *name, t_desc **res)
{
    *res=find_desc(name);

    /* P_POINT(*res); P_NL; */

    if(*res == NULL) {
        return 1;
    }

    /* P_INT((*res)->n_el); P_NL; */
    return 0;
}

int set_desc1(char *name)
{
    return set_desc(name, &desc1);
}

int set_desc2(char *name)
{
    return set_desc(name, &desc2);
}

int my_load_domain(char *name)
{
    t_domain *dom=load_domain(name);

    if(dom == NULL) {
        return 1;
    }

    return 0;
}

void unpack_pars(t_compdesc_parameters *pars)
{
    th_n_el = 0;
    th_n_seg_adj = 0;
    th_seg_adj_len = 18.0;
    th_zero_el_RMSD = pars->th_zero_el_RMSD;
    th_el_RMSD = pars->th_el_RMSD;
    th_pair_RMSD = pars->th_pair_RMSD;
    force_order = pars->force_order;
    n_AA_straighten = 1;
    th_el_align = pars->th_el_align;
    th_seg_align = 0;
    th_min_al_cov = pars->th_mer_align;
    th_align_size_equiv = pars->th_align_size_equiv;
    th_overall_RMSD = pars->th_overall_RMSD;
    th_full_cov = 0.9;
}

void average_coords(char *strname, float ***coords_p, int ***seg_adj_p)
{
    t_domain *dom = find_domain(strname);

    float **coords = array_alloc(dom->n_AA, 7, sizeof(float));

    for(int i=1; i<=dom->n_AA; i++) {
        coords[i-1][0]=dom->CA[i].x;
        coords[i-1][1]=dom->CA[i].y;
        coords[i-1][2]=dom->CA[i].z;
        coords[i-1][3]=dom->STR[i].x;
        coords[i-1][4]=dom->STR[i].y;
        coords[i-1][5]=dom->STR[i].z;
        if(i+1<=dom->n_AA) coords[i-1][6]=dist(&(dom->STR[i]), &(dom->STR[i+1]));
        else coords[i-1][6]=NAN;
    }

    *coords_p=coords;

    int **seg_adj = array_alloc(dom->n_AA, dom->n_AA, sizeof(int));

    for(int i=1; i<=dom->n_AA; i++) {
        for(int j=1; j<=dom->n_AA; j++) {
            if(j<i) seg_adj[i-1][j-1]=0;
            else if(j==i) seg_adj[i-1][j-1]=1;
            else seg_adj[i-1][j-1]=adj_seg_seg(dom, i, j);
        }
    }

    *seg_adj_p = seg_adj;
}

int compdesc_ref(int max_res, t_compdesc_parameters *pars, t_compdesc_result **results, t_compdesc_debug **debug)
{

    load_config(0);
    unpack_pars(pars);

    t_group_entry *alignment;

	int n_res=desc_comp_multi_debug(desc1, desc2, &alignment, max_res, debug);


    *results = calloc(n_res, sizeof(t_compdesc_result));

    t_compdesc_result *res = *results;

    for(int i=0; i<n_res; i++) {
        res[i].n_monomers = alignment[i].n_AA;
        res[i].RMSD = alignment[i].RMSD;
        res[i].TR=alignment[i].TR;

        res[i].desc1_monomers = calloc(res[i].n_monomers, sizeof(int));
        res[i].desc2_monomers = calloc(res[i].n_monomers, sizeof(int));

        int pos=0;

        void f(t_group_entry *al, int aa1, int aa2, int contact) {
            res[i].desc1_monomers[pos] = aa1;
            res[i].desc2_monomers[pos] = aa2;
            pos++;
        }

        alignment_iterate_AA(&(alignment[i]), f);
    }

    return n_res;
}

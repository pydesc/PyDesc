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
 
#ifndef __COMPDESC_DEBUG_H__
#define __COMPDESC_DEBUG_H__

typedef struct t_compdesc_debug {
    int max_al_len;

    int n_el1;
    int n_el2;
    int *elts1;
    int *elts2;
    float th_n_el;
    int stage1;

    int n_seg_adj1;
    int n_seg_adj2;
    float th_n_seg_adj;
    int stage2;

    float zero_el_rmsd;
    float th_zero_el_rmsd;
    int stage3;

    float th_el_rmsd;
    float th_pair_rmsd;
    int n_comp;
    int *comp1;
    int *comp2;
    float *el_rmsd;
    float *pair_rmsd;
    int n_edges;
    int *edge_list1;
    int *edge_list2;
    int stage4;

    char **edge_graph;
    int n_sets;
    int **indep_sets;
    int max_allowed;
    int *n_allowed;
    int ***allowed_comb;
    int *n_align;
    int ***allowed_align;
    int stage5;

    int n_full_align;
    int **full_align;
    int **full_align_cov;
    float *full_align_rmsd;
    int stage6;

    float **full_align_crop_val_run1;
    int **full_align_crop_max_run1;
    int **full_align_run1;
    float *full_align_rmsd_run1;

    int n_res_align;
    int **res_align;
    int stage7;
} t_compdesc_debug;

#ifdef __DESC_H__
int desc_comp_multi_debug(t_desc *desc1, t_desc *desc2, t_group_entry **alignment, int max_n_al, t_compdesc_debug **debug);
#endif
#endif

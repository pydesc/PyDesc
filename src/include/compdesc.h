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
 
#ifndef __COMPDESC_H__
#define __COMPDESC_H__

#include "cstructures.h"
#include "overfit.h"

typedef struct t_compdesc_parameters {
    float th_zero_el_RMSD;
    float th_el_RMSD;
    float th_pair_RMSD;
    int force_order;
    float th_mer_align;
    float th_el_align;
    float th_con_align;
    float th_align_size_equiv;
    float th_overall_RMSD;
} t_compdesc_parameters;

typedef struct t_compdesc_result {
    int n_monomers;
    int* desc1_monomers;
    int* desc2_monomers;

    int n_elements;
    int* desc1_elements;
    int* desc2_elements;

    int n_contacts;
    int* desc1_contacts;
    int* desc2_contacts;

    float RMSD;
    t_transrot TR;
} t_compdesc_result;


#ifdef COMPDESC_DEBUG

typedef struct t_compdesc_debug {

    int n_el1;
    int n_el2;
    int *elts1;
    int *elts2;
    int stage0;

    void *prelim_cov;
    int stage1;

    float zero_el_rmsd;
    int n_all_element_comp;
    void *all_element_comp;
    int n_all_contact_comp;
    void *all_contact_comp;
    void *sim_components;
    void *sim_graph;
    int stage2;

    int n_sets;
    int **indep_sets;
    int *n_allowed;
    int ***allowed_comb_per_set;
    int n_full_align;
    int **allowed_comb;
    void *full_align;
    int stage3;

    float **full_align_crop_val_run1;
    int **full_align_crop_max_run1;
    void *full_align_run1;

    int n_res_align;
    void *res_align;
    int stage4;
} t_compdesc_debug;

int compdesc(CDescriptor *desc1, CDescriptor *desc2, int max_res, t_compdesc_parameters *pars, t_compdesc_result **results, t_compdesc_debug **debug);
#else
int compdesc(CDescriptor *desc1, CDescriptor *desc2, int max_res, t_compdesc_parameters *pars, t_compdesc_result **results);
#endif

int has_debug(void);

#endif

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
 
#ifndef __COMPDESC_PRIV_H__
#define __COMPDESC_PRIV_H__

#include "overfit.h"
#include "cstructures.h"
#include "maps.h"
#include "compdesc.h"


#define IN_DESC 1
#define CENTRAL 2
#define ZERO 4

#define NON_OPT 8
#define NON_OPT_EL 16

#define ALIGNED 64
#define ALIGNED_EL 128

extern float th_zero_el_RMSD;
extern float th_el_RMSD;
extern float th_pair_RMSD;
extern int force_order;
extern float th_mer_align;
extern float th_el_align;
extern float th_con_align;
extern float th_align_size_equiv;
extern float th_overall_RMSD;

typedef struct {
    int n_mer;
    int n_mer_all;
    int n_el;
    int n_el_all;
    int n_con;
    int n_con_all;

    int n_AA;
    int n_contacts;

    int *desc_vec;
    int *con_vec;
} t_desc_coverage;

typedef struct {
    t_desc_coverage desc1;
    t_desc_coverage desc2;
} t_al_coverage;

typedef struct {
    t_map *mer_map;
    t_map *el_map;
    t_map *contact_map;

    int *mer_cov_count;
    int *el_cov_count;

    CDescriptor *desc1;
    CDescriptor *desc2;

    int n_AA;

    t_al_coverage *coverage;

    int dud;

    float RMSD;
    t_transrot TR;
} t_alignment;



typedef struct {
    CElement *el1;
    CElement *el2;
    float RMSD;
    t_overfit_sums overfit_sums;
    int consim_cnt; // number of similar contacts covering element similarity
    int core_cnt; // as above but recursively exluding contacts
    int pos;
} t_elsim;

typedef struct {
    CContact *con1;
    CContact *con2;
    t_elsim *elsim1;
    t_elsim *elsim2;
    float RMSD;
    int pos;
} t_consim;

typedef struct {
    t_elsim *elsim_list;
    t_consim *consim_list;

    t_consim ***consim_map;

    int n_elsim;
    int n_consim;
} t_sim_components;

typedef struct {
    t_elsim *elsim;
    t_consim *consim;
} t_sim;

typedef struct {
    t_sim *sim_list;
    char **sim_graph;

    int n_sim;
} t_sim_graph;

void al_coverage_free(t_al_coverage *al_cov);
void copy_al_coverage(t_al_coverage *dest, t_al_coverage *src);
void destroy_al_coverage(t_al_coverage *al_cov);
int al_coverage_valid(t_al_coverage *al_cov);

int coverage_test_prelim(CDescriptor *desc1, CDescriptor *desc2, t_al_coverage **coverage_p);
int coverage_test_component(CDescriptor *desc1, CDescriptor *desc2, t_sim_components *sim_components, t_al_coverage **coverage_p);
int coverage_test_alignment(t_alignment* align_new, t_al_coverage **coverage_p);

float RMSD_element(int token, CElement *el1, CElement *el2, CStructure *struct1, CStructure *struct2);
float RMSD_pair(int token, t_overfit_sums *sums1, t_overfit_sums *sums2);
float RMSD_frag(t_alignment *align_new, int zero_el1_c_ind, int zero_el2_c_ind, CStructure *struct1, CStructure *struct2);
float RMSD_frag_crop(t_alignment *align, float threshold, int max_del, int *act_del, t_al_coverage *al_cov, float *deb_val, int *deb_max);

void destroy_alignment(t_alignment *align);

void alignment_size(t_alignment *align);
int is_subalignment(t_alignment *al1, t_alignment *al2);

void copy_alignment(t_alignment *dest, t_alignment *align);


#endif // __COMPDESC_PRIV_H__

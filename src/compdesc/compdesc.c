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

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


#include "cstructures.h"

#include "compdesc_priv.h"
#include "compdesc.h"

#include "overfit.h"
#include "simple_macros.h"

#include "arrays.h"


#include "01_prelim.h"
#include "02_components.h"
#include "03_build.h"
#include "04_prune.h"

//#define DEBUG_FILL_ALLOWED

float th_zero_el_RMSD;
float th_el_RMSD;
float th_pair_RMSD;
int force_order;
float th_mer_align;
float th_el_align;
float th_con_align;
float th_align_size_equiv;
float th_overall_RMSD;




static void unpack_pars(t_compdesc_parameters *pars)
{
    th_zero_el_RMSD = pars->th_zero_el_RMSD;
    th_el_RMSD = pars->th_el_RMSD;
    th_pair_RMSD = pars->th_pair_RMSD;
    force_order = pars->force_order;
    th_mer_align = pars->th_mer_align;
    th_el_align = pars->th_el_align;
    th_con_align = pars->th_con_align;
    th_align_size_equiv = pars->th_align_size_equiv;
    th_overall_RMSD = pars->th_overall_RMSD;
}

#ifdef COMPDESC_DEBUG

static int unmap(CDescriptor *desc, int c_ind) {
    if(c_ind>=0) {
        return map_unmap(desc->structure->ind_to_c_map, c_ind);
    } else {
        return c_ind;
    }
}

static void save_descs_to_deb(CDescriptor *desc1, CDescriptor *desc2, t_compdesc_debug *deb)
{
    deb->n_el1=desc1->n_elements;
    deb->n_el2=desc2->n_elements;

    void deb_save_elts(int **deb_elts, CDescriptor *desc) {
        *deb_elts=calloc(desc->n_elements, sizeof(int));
        int pos=1;
        for(int i=0; i<desc->n_elements; i++) if(i!=desc->central_element_used) (*deb_elts)[pos++]=unmap(desc, desc->elements[i].center_c_ind);
        (*deb_elts)[0]=unmap(desc, desc->central_element_c_ind);
    }

    deb_save_elts(&(deb->elts1), desc1);
    deb_save_elts(&(deb->elts2), desc2);

    deb->stage0=1;
}
#endif

int compdesc(CDescriptor *desc1, CDescriptor *desc2, int max_res, t_compdesc_parameters *pars, t_compdesc_result **results 
#ifdef COMPDESC_DEBUG
        , t_compdesc_debug **debug
#endif
        )
{
    int n_full_align=0;
    int *al_order=0;
    int token=overfit_get_token();
    t_alignment *full_align=0;

    t_sim_components *sim_components=0;
    t_sim_graph *sim_graph=0;

    int n_res=0;

    unpack_pars(pars);



#ifdef COMPDESC_DEBUG
    t_compdesc_debug *deb;
    t_compdesc_debug fake_deb;
    if(debug) {
        *debug=calloc(1, sizeof(t_compdesc_debug));
        deb=*debug;
    } else {
        deb=&fake_deb;
    }

    save_descs_to_deb(desc1, desc2, deb);

#endif

    // 1. Preliminary checks
#ifdef COMPDESC_DEBUG
    if(!check_sizes(desc1, desc2, deb)) goto fail;
#else
    if(!check_sizes(desc1, desc2)) goto fail;
#endif

    // 2. 0-element RMSD, similar contacts


#ifdef COMPDESC_DEBUG
    sim_components=compare_elements(token, desc1, desc2, deb);
#else
    sim_components=compare_elements(token, desc1, desc2);
#endif

    if(sim_components==0) goto fail;


    /* P_INT(sim_components->n_elsim) P_INT(sim_components->n_consim) P_NL; */

#ifdef COMPDESC_DEBUG
    sim_graph=build_graph(sim_components, deb);
#else
    sim_graph=build_graph(sim_components);
#endif

    /* P_INT(sim_graph->n_sim) P_NL; */

    clean_ion_contacts(token, desc1, desc2, sim_components, sim_graph);

    // 3. Possible alignments and their sizes
#ifdef COMPDESC_DEBUG
    n_full_align=build_alignments(desc1, desc2, &full_align, sim_components, sim_graph, deb);
#else
    n_full_align=build_alignments(desc1, desc2, &full_align, sim_components, sim_graph);
#endif


    if(n_full_align==0) goto fail;

    /* P_INT(n_full_align) P_NL; */

    al_order=calloc(n_full_align, sizeof(int));

    for(int i=0; i<n_full_align; i++) {
        al_order[i]=i;
    }

	// 4. Best alignment with qualifying RMSD
#ifdef COMPDESC_DEBUG
    int cnt=prune_alignments(desc1, desc2, full_align, al_order, n_full_align, max_res, deb);
#endif

#ifdef COMPDESC_DEBUG
    n_res=save_results(desc1, desc2, results, full_align, al_order, n_full_align, max_res, deb);
#else
    n_res=save_results(desc1, desc2, results, full_align, al_order, n_full_align, max_res);
#endif

fail:

#ifndef COMPDESC_DEBUG
    if(sim_components) sim_components_free(sim_components);
    if(sim_graph) sim_graph_free(sim_graph);
#endif

    for(int i=0; i<n_full_align; i++) destroy_alignment(&full_align[i]);
    if(full_align) free(full_align);

    if(al_order) free(al_order);
    overfit_release_token(token);

    return n_res;
}

void free_results(t_compdesc_result* res, int n)
{
    for(int i=0; i<n; i++) {
        free(res[i].desc1_monomers);
        free(res[i].desc1_elements);
        free(res[i].desc1_contacts);
        free(res[i].desc2_monomers);
        free(res[i].desc2_elements);
        free(res[i].desc2_contacts);
    }

    free(res);
}

int has_debug(void)
{
    int res;

#ifdef COMPDESC_DEBUG
    res = 1;
#else
    res = 0;
#endif

    return res;
}


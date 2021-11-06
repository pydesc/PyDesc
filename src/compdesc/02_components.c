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

#include"simple_macros.h"
#include"compdesc_priv.h"
#include"cstructures.h"
#include"arrays.h"
#include"02_components.h"
#include"string.h"

static int element_overlap(CElement *el1, CElement *el2)
{
    if(el1->start_c_ind<=el2->start_c_ind && el2->start_c_ind<=el1->end_c_ind) return 1;
    if(el2->start_c_ind<=el1->start_c_ind && el1->start_c_ind<=el2->end_c_ind) return 1;

    return 0;
}

static int element_offset(CElement *el1, CElement *el2)
{
    return el1->center_c_ind-el2->center_c_ind;
}

void sim_components_free(t_sim_components *sim_components)
{
    free(sim_components->elsim_list);
    free(sim_components->consim_list);

    array_free((void **)sim_components->consim_map);

    free(sim_components);
}

void sim_graph_free(t_sim_graph *sim_graph)
{
    free(sim_graph->sim_list);
    array_free((void **)sim_graph->sim_graph);
    free(sim_graph);
}

t_sim_components *compare_elements(int token, CDescriptor *desc1, CDescriptor *desc2
#ifdef COMPDESC_DEBUG
        , t_compdesc_debug *deb
#endif
        )
    /*
     * TODO: Add support for descriptors not having a central element.
     */

{

    CElement *zero_el1=&desc1->elements[desc1->central_element_used];
    CElement *zero_el2=&desc2->elements[desc2->central_element_used];

    // 3. 0-element RMSD
    float zero_el_RMSD=RMSD_element(token, zero_el1, zero_el2, desc1->structure, desc2->structure);
#ifdef COMPDESC_DEBUG
    deb->zero_el_rmsd=zero_el_RMSD;
#endif

    if(zero_el_RMSD>th_zero_el_RMSD) return 0;

#ifdef COMPDESC_DEBUG
    deb->stage2=1;
#endif

    int **elel_to_elsim=(int **)array_alloc(desc1->n_elements, desc2->n_elements, sizeof(int));

    for(int i=0; i<desc1->n_elements; i++)
        for(int j=0; j<desc2->n_elements; j++)
            elel_to_elsim[i][j]=-1;

    t_elsim *elsim_list=calloc(desc1->n_elements*desc2->n_elements, sizeof(t_elsim));
    t_consim *consim_list=calloc(desc1->contact_map->n_contacts*desc2->contact_map->n_contacts, sizeof(t_consim));

    int n_elsim=0;
    int n_consim=0;

    elsim_list[n_elsim].el1=zero_el1;
    elsim_list[n_elsim].el2=zero_el2;
    elsim_list[n_elsim].RMSD=zero_el_RMSD;
    elel_to_elsim[desc1->central_element_used][desc2->central_element_used]=n_elsim;
    overfit_sumsave_str(token, &(elsim_list[n_elsim].overfit_sums));
    n_elsim++;


    // 4. RMSD of pairs
#ifdef COMPDESC_DEBUG
    deb->n_all_element_comp=1+(desc1->n_elements-1)*(desc2->n_elements-1);
    t_elsim *deb_all_element_comp=calloc(deb->n_all_element_comp, sizeof(t_elsim));
    deb->all_element_comp=deb_all_element_comp;
    int deb_pos=0;
    deb_all_element_comp[deb_pos++]=elsim_list[0];
#endif

    for(int i=0; i<desc1->n_elements; i++) if(i!=desc1->central_element_used) {
        CElement *el1=&desc1->elements[i];
        for(int j=0; j<desc2->n_elements; j++) if(j!=desc2->central_element_used) {
            CElement *el2=&desc2->elements[j];
            float dist;

            if(desc1->structure->monomers[el1->center_c_ind].type!=desc2->structure->monomers[el2->center_c_ind].type) continue;

            if(element_overlap(zero_el1, el1) || element_overlap(zero_el2, el2))
                if(element_offset(zero_el1, el1) != element_offset(zero_el2, el2)) continue;

            if(force_order)
                if(element_offset(zero_el1, el1) && element_offset(zero_el2, el2)) continue;

            dist=RMSD_element(token, el1, el2, desc1->structure, desc2->structure);


            elsim_list[n_elsim].el1=el1;
            elsim_list[n_elsim].el2=el2;
            elsim_list[n_elsim].RMSD=dist;
            overfit_sumsave_str(token, &(elsim_list[n_elsim].overfit_sums));

#ifdef COMPDESC_DEBUG
            deb_all_element_comp[deb_pos++]=elsim_list[n_elsim];
#endif

            if(dist>th_el_RMSD) {
                continue;
            }

            elel_to_elsim[i][j]=n_elsim;
            n_elsim++;
        }
    }


#ifdef COMPDESC_DEBUG
    deb->n_all_element_comp=deb_pos;


    deb->stage2=2;
    deb->n_all_contact_comp=desc1->contact_map->n_contacts*desc2->contact_map->n_contacts/2;
    t_consim *deb_all_contact_comp=calloc(deb->n_all_contact_comp, sizeof(t_consim));
    deb->all_contact_comp=deb_all_contact_comp;
    deb_pos=0;
#endif

    for(int i=0; i<desc1->contact_map->n_contacts; i++) {
        CContact *con1=&desc1->contact_map->contacts[i];
        /*
         * If each contact ocurs twice in contact_map (A-B and B-A),
         * to avoid duplications, pick only one of the pair in the first descriptor.
         * In such case, include this line:
        if(con1->mer1_used>con1->mer2_used) continue;
         */
        for(int j=0; j<desc2->contact_map->n_contacts; j++) {
            CContact *con2=&desc2->contact_map->contacts[j];
            float dist;

            CElement *el11, *el12, *el21, *el22;

            el11=&desc1->elements[con1->mer1_used];
            el12=&desc1->elements[con1->mer2_used];
            el21=&desc2->elements[con2->mer1_used];
            el22=&desc2->elements[con2->mer2_used];

            if(elel_to_elsim[el11->center_used][el21->center_used]<0 || elel_to_elsim[el12->center_used][el22->center_used]<0) continue;
            t_elsim *elsim1=&(elsim_list[elel_to_elsim[el11->center_used][el21->center_used]]);
            t_elsim *elsim2=&(elsim_list[elel_to_elsim[el12->center_used][el22->center_used]]);
            if(force_order)
                if(element_offset(el11, el21) && element_offset(el12, el22)) continue;

            dist=RMSD_pair(token, &(elsim1->overfit_sums), &(elsim2->overfit_sums));

            consim_list[n_consim].con1=con1;
            consim_list[n_consim].con2=con2;
            consim_list[n_consim].elsim1=elsim1;
            consim_list[n_consim].elsim2=elsim2;
            consim_list[n_consim].RMSD=dist;
#ifdef COMPDESC_DEBUG
            deb_all_contact_comp[deb_pos++]=consim_list[n_consim];
#endif

            if(dist>th_pair_RMSD) {
                continue;
            }

            elsim1->consim_cnt++;
            elsim2->consim_cnt++;

            n_consim++;
        }
    }

#ifdef COMPDESC_DEBUG
    deb->n_all_contact_comp=deb_pos;
    deb->stage2=3;
#endif

    /*
     * Remove elements not covered by any contacts.
     * Central element is exempt from this filter.
     *
     *
     * CAUTION: First loop breaks pointers to t_elsim structures in consim_list.
     */

    int pos=1;
    for(int i=0; i<n_elsim; i++) {
        if(i==0) continue;
        if(elsim_list[i].consim_cnt>0) {
            elsim_list[pos]=elsim_list[i];
            elsim_list[pos].pos=pos;
            t_elsim *elsim=&elsim_list[i];
            elel_to_elsim[elsim->el1->center_used][elsim->el2->center_used]=pos;
            pos++;
        }
    }
    n_elsim=pos;

    elsim_list=realloc(elsim_list, n_elsim * sizeof(t_elsim));
    consim_list=realloc(consim_list, n_consim * sizeof(t_consim));
    t_consim ***consim_map=(t_consim ***)array_alloc(n_elsim, n_elsim, sizeof(t_consim *));

    for(int i=0; i<n_consim; i++) {
        CContact *con1=consim_list[i].con1;
        CContact *con2=consim_list[i].con2;

        int elsim1_ind=elel_to_elsim[con1->mer1_used][con2->mer1_used];
        int elsim2_ind=elel_to_elsim[con1->mer2_used][con2->mer2_used];

        consim_list[i].elsim1=&(elsim_list[elsim1_ind]);
        consim_list[i].elsim2=&(elsim_list[elsim2_ind]);

        consim_list[i].pos=i;

        consim_map[elsim1_ind][elsim2_ind]=&(consim_list[i]);
        consim_map[elsim2_ind][elsim1_ind]=&(consim_list[i]);
    }



    array_free((void **)elel_to_elsim);

    t_sim_components *res=calloc(1, sizeof(t_sim_components));

    res->elsim_list=elsim_list;
    res->consim_list=consim_list;
    res->consim_map=consim_map;
    res->n_elsim=n_elsim;
    res->n_consim=n_consim;

#ifdef COMPDESC_DEBUG
    deb->sim_components=(void *)res;
    deb->stage2=3;
#endif

    if(!coverage_test_component(desc1, desc2, res, 0)) {
#ifndef COMPDESC_DEBUG
        sim_components_free(res);
#endif
        res=0;
    }

#ifdef COMPDESC_DEBUG
    deb->stage2=4;
#endif

    return res;
}

t_sim_graph *build_graph(t_sim_components *sim_components
#ifdef COMPDESC_DEBUG
        , t_compdesc_debug *deb
#endif
        )
{
    t_elsim *elsim_list=sim_components->elsim_list;
    t_consim ***consim_map=sim_components->consim_map;
    int n_elsim=sim_components->n_elsim;

    t_sim *sim_list=calloc(n_elsim, sizeof(t_sim));
    int n_sim=0;

    for(int i=0; i<n_elsim; i++) {
        sim_list[n_sim].elsim=&(elsim_list[i]);
        n_sim++;
    }

    char **sim_graph=(char **) array_alloc(n_sim, n_sim, sizeof(char));

    int elsim_overlap(t_elsim *elsim1, t_elsim *elsim2) {
        return element_overlap(elsim1->el1, elsim2->el1) || element_overlap(elsim1->el2, elsim2->el2);
    }

    int elsim_contr(t_elsim *elsim1, t_elsim *elsim2) {
        if(elsim_overlap(elsim1, elsim2)) return element_offset(elsim1->el1, elsim2->el1) != element_offset(elsim1->el2, elsim2->el2);
        else return 0;
    }

    int elsim_order(t_elsim *elsim1, t_elsim *elsim2) {
        return element_offset(elsim1->el1, elsim2->el1) * element_offset(elsim1->el2, elsim2->el2)<0 ? 0 : 1;
    }

    int elsim_connected(t_elsim *elsim1, t_elsim *elsim2) {
        return consim_map[elsim1->pos][elsim2->pos]!=0;
    }

    int elconsim_test(t_elsim *elsim1, t_consim *consim2, int (*elsim_test)(t_elsim *elsim1, t_elsim *elsim2)) {
        return elsim_test(elsim1, consim2->elsim1) || elsim_test(elsim1, consim2->elsim2);
    }

    int consim_test(t_consim *consim1, t_consim *consim2, int (*elsim_test)(t_elsim *elsim1, t_elsim *elsim2)) {
        return elconsim_test(consim1->elsim1, consim2, elsim_test) || elconsim_test(consim1->elsim2, consim2, elsim_test);
    }

    int sim_test(t_sim *sim1, t_sim *sim2, int (*elsim_test)(t_elsim *elsim1, t_elsim *elsim2)) {
        if(sim1->elsim) {
            if(sim2->elsim) return elsim_test(sim1->elsim, sim2->elsim);
            else return elconsim_test(sim1->elsim, sim2->consim, elsim_test);
        } else {
            if(sim2->elsim) return elconsim_test(sim2->elsim, sim1->consim, elsim_test);
            else return consim_test(sim1->consim, sim2->consim, elsim_test);
        }
    }

    for(int i=0; i<n_sim; i++) {
        sim_graph[i][i]=1;
        for(int j=i+1; j<n_sim; j++) {
            if(sim_test(&sim_list[i], &sim_list[j], elsim_contr)) sim_graph[i][j]=1;        //contradictory

            if(force_order) {
                if(!sim_test(&sim_list[i], &sim_list[j], elsim_order)) sim_graph[i][j]=1;   //contradictory
            }

            if(!sim_graph[i][j]) {
                if(sim_test(&sim_list[i], &sim_list[j], elsim_connected)) sim_graph[i][j]=-1;   //connected by a consim
            }

            sim_graph[j][i]=sim_graph[i][j];
        }
    }

    t_sim_graph *res=calloc(1, sizeof(t_sim_graph));

    res->sim_list=sim_list;
    res->sim_graph=sim_graph;
    res->n_sim=n_sim;

#ifdef COMPDESC_DEBUG
    deb->sim_graph=(void *)res;
    deb->stage2=5;
#endif

    return res;
}

void clean_ion_contacts(int token, CDescriptor *desc1, CDescriptor *desc2, t_sim_components *sim_components, t_sim_graph *graph)
{
    /*
     *
     * CAVEAT: This function works for graphs which have only t_elsim nodes
     *
     */


    for(int k=0; k<graph->n_sim; k++) {
        t_elsim *ion_sim=graph->sim_list[k].elsim;

        if(ion_sim==0) {
            printf("clean_ion_contacts used on t_sim_graph with t_consim nodes\n");
            abort();
        }

        if(ion_sim->el1->start_c_ind!=ion_sim->el1->end_c_ind) continue;
        if(desc1->structure->monomers[ion_sim->el1->center_c_ind].n_points!=1) continue;

        for(int i=0; i<graph->n_sim; i++) if(graph->sim_graph[i][k]==-1) {
            for(int j=i+1; j<graph->n_sim; j++) if(graph->sim_graph[j][k]==-1) {
                if(graph->sim_graph[i][j]==0) {

                    float dist=RMSD_pair(token, &(graph->sim_list[i].elsim->overfit_sums), &(graph->sim_list[j].elsim->overfit_sums));

                    if(dist>2*th_overall_RMSD) {
                        graph->sim_graph[i][j]=1;
                    }

                }
            }
        }
    }


    int dfs_path[graph->n_sim+1];

    dfs_path[0]=0;

    int visited[graph->n_sim];

    memset(visited, 0, sizeof(visited));

    void dfs(int start) {
        visited[start]=1;

        dfs_path[++dfs_path[0]]=start;

        for(int i=0; i<graph->n_sim; i++) {
            if(graph->sim_graph[start][i]==-1) {
                int log=1;
                for(int j=1; j<=dfs_path[0]; j++) {
                    int path_el=dfs_path[j];
                    if(i==path_el || graph->sim_graph[i][path_el]==1) {
                        log=0;
                        break;
                    }
                }

                if(log) {
                    dfs(i);
                }
            }
        }

        dfs_path[0]--;
    }

    dfs(0);


    for(int i=0; i<graph->n_sim; i++) {
        if(!visited[i]) {
            for(int j=1; j<graph->n_sim; j++) {
                graph->sim_graph[i][j]=0;
                graph->sim_graph[j][i]=0;
            }
        }
    }

}

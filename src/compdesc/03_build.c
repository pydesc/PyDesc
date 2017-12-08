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

#include<stdio.h>
#include<string.h>

#include"simple_macros.h"

#include"arrays.h"

#include"compdesc_priv.h"
#include"cstructures.h"
#include"02_components.h"
#include"03_build.h"


static int element_length(CElement *el)
{
    return el->end-el->start+1;
}

static void fill_alignment(t_alignment *align)
{
    CDescriptor *desc1=align->desc1;
    CDescriptor *desc2=align->desc2;

    align->mer_map = map_alloc(desc1->structure->n_monomers, desc2->structure->n_monomers);
    align->mer_cov_count = calloc(desc1->structure->n_monomers, sizeof(int));

    void fill_mer_map(int el_used1, int el_used2) {
        CElement *el1 = &(desc1->elements[el_used1]);
        CElement *el2 = &(desc2->elements[el_used2]);

        if(element_length(el1) != element_length(el2)) {
            printf("el1: %d len: %d\n", el1->center, element_length(el1));
            printf("el2: %d len: %d\n", el2->center, element_length(el2));
            printf("Aborting\n");
            abort();
        }

        for(int pos1=el1->start_c_ind, pos2=el2->start_c_ind; pos1<=el1->end_c_ind; pos1++, pos2++) {
            int cm1=map_unmap(align->mer_map, pos2);
            int cm2=map_map(align->mer_map, pos1);

            if(cm1!=-1 && cm1!=pos1) {
                printf("fill_alignment: map_map_error  "); P_INT(cm1); P_INT(pos1); P_NL;
                abort();
            }
            if(cm2!=-1 && cm2!=pos2) {
                printf("fill_alignment: map_map_error  "); P_INT(cm2); P_INT(pos2); P_NL;
                abort();
            }
            map_set(align->mer_map, pos1, pos2);

            align->mer_cov_count[pos1]++;
        }
    }

    map_iter(align->el_map, fill_mer_map);
}

static void convert_comb_align(t_alignment *align, int *allowed_comb, t_sim_components *sim_components, t_sim *sim_list, CDescriptor *desc1, CDescriptor *desc2)
{
    memset(align, 0, sizeof(t_alignment));

    align->el_map = map_alloc(desc1->n_elements, desc2->n_elements);
    align->contact_map = map_alloc(desc1->contact_map->n_contacts, desc2->contact_map->n_contacts);

    align->el_cov_count = calloc(desc1->n_elements, sizeof(int));

    t_elsim *elsim_list[allowed_comb[0]*2];
    int n_el=0;

    /* P_T("convert_comb_align") P_NL; */
    /* P_T("elts1: "); */

    for(int i=0; i<allowed_comb[0]; i++) {
        void add_elsim(t_elsim *elsim) {
            elsim_list[n_el++]=elsim;
            map_set(align->el_map, elsim->el1->center_used, elsim->el2->center_used);
            /* printf("%d ", elsim->el1->center_used); */
        }

        t_sim *sim=&(sim_list[allowed_comb[i+1]]);

        if(sim->elsim) {
            add_elsim(sim->elsim);
        } else if(sim->consim) {
            add_elsim(sim->consim->elsim1);
            add_elsim(sim->consim->elsim2);
        } else {
            printf("both elsim and consim false in convert_comb_align\n");
            abort();
        }
    }

    /* P_NL; */

    for(int i=0; i<n_el; i++) {
        for(int j=0; j<n_el; j++) {
            t_consim *consim=sim_components->consim_map[elsim_list[i]->pos][elsim_list[j]->pos];

            /* P_INT(elsim_list[i]->pos) P_INT(elsim_list[j]->pos) */

            if(!consim) {
                /* P_INT(-1) P_NL; */
                continue;
            }

            align->el_cov_count[consim->elsim1->el1->center_used]++;
            align->el_cov_count[consim->elsim2->el1->center_used]++;

            /* P_INT(consim->pos) P_NL; */
            map_set(align->contact_map, consim->con1->pos, consim->con2->pos);

        }
    }


    align->desc1=desc1;
    align->desc2=desc2;

    fill_alignment(align);
}

static void edge_graph_dfs(char **edge_graph, int n_edges, int start, int val, char *visited, char *reachable, void (*f)(int e))
{
    visited[start]=1;

    f(start);

    for(int i=0; i<n_edges; i++) if(!visited[i] && reachable[i] && edge_graph[start][i]==val) edge_graph_dfs(edge_graph, n_edges, i, val, visited, reachable, f);
}

static int fill_indep_sets(char **edge_graph, int n_edges, int **indep_sets, int max_set)
{
    /*
     *
     *  CAVEAT: It is assumed that zero-element is at position 0.
     *
     */

    int n_sets=0;
    char visited[n_edges];
    char reachable[n_edges];

    memset(reachable, 1, sizeof(reachable));
    memset(visited, 0, sizeof(visited));

    reachable[0]=2;

    int int_comp(const void *a, const void *b) { return *((int *)a)-*((int *)b); }

    void mark_set(int e) {
        indep_sets[n_sets][++indep_sets[n_sets][0]]=e;
    }

    void mark_reachable(int e) {
        reachable[e]=2;
    }

    edge_graph_dfs(edge_graph, n_edges, 0, -1, visited, reachable, mark_reachable);

    for(int i=0; i<n_edges; i++) reachable[i]--;
    memset(visited, 0, sizeof(visited));

    for(int i=0; i<n_edges; i++) if(!visited[i] && reachable[i]) {
        indep_sets[n_sets][0]=0;
        edge_graph_dfs(edge_graph, n_edges, i, 1, visited, reachable, mark_set);
        n_sets++;
    }

    for(int i=0; i<n_sets; i++) {
        qsort(indep_sets[i]+1, indep_sets[i][0], sizeof(int), int_comp);
    }

    return n_sets;
}

static int fill_allowed(int ***allowed_comb_p, int *subset, char **edge_graph, int max_len, int n_edges)
{
    int n_comb=0;

    int n_set=subset[0];
    int elim[n_set+1];

    int res[n_set+1];

    int max_comb=n_edges;

    /* max_comb=1; // TODO: Temporary for testing purposes. Remove when done. */

    int **allowed_comb=(int **) array_alloc(max_comb, max_len+1, sizeof(int));

    for(int i=1; i<=n_set; elim[i++]=1000);

    res[0]=0;

    void do_traverse(int start, int lev)
    {
        int next;
        int log=0;

#ifdef DEBUG_FILL_ALLOWED
        printf("do_traverse(%d): start: %d subset[start]: %d  lev %d\n", lev, start, subset[start], lev);
#endif

        for(int i=1; i<=n_set; i++) if(elim[i]>lev && i>=start) {
            res[++res[0]]=subset[i];

#ifdef DEBUG_FILL_ALLOWED
            printf("do_traverse(%d): add_res: %d (%d)\n", lev, i, subset[i]);
#endif

            next=-1;
            log=0;

            for(int j=n_set; j>0; j--) if(elim[j]>lev) {
                if(edge_graph[subset[i]][subset[j]]==1) {
                    elim[j]=lev;
#ifdef DEBUG_FILL_ALLOWED
                    printf("do_traverse(%d): elim: %d (%d)\n", lev, j, subset[j]);
#endif
                } else {
                    if(j<i) {
                        log=1;
#ifdef DEBUG_FILL_ALLOWED
                        printf("do_traverse(%d): log=1 j: %d(%d)\n", lev, j, subset[j]);
#endif
                    } else next=j;
                }
            }

            if(next==-1) {
                if(n_comb>=max_comb) {
                    max_comb*=2;
                    allowed_comb=(int **)array_realloc((void **)allowed_comb, max_comb, max_len+1, sizeof(int));
                }
                if(!log) {
#ifdef DEBUG_FILL_ALLOWED
                    printf("do_traverse(%d): saving: ", lev);
#endif
                    for(int j=0; j<=res[0]; j++) {
                        allowed_comb[n_comb][j]=res[j];
#ifdef DEBUG_FILL_ALLOWED
                        printf("%d ", res[j]);
#endif
                    }
#ifdef DEBUG_FILL_ALLOWED
                    printf("\n");
#endif
                    n_comb++;
                }
            } else  do_traverse(next, lev+1);

            for(int j=n_set; j>0; j--) if(edge_graph[subset[i]][subset[j]]==1) {
                elim[j]=(elim[j]>=lev)?1000:elim[j];
#ifdef DEBUG_FILL_ALLOWED
                if(elim[j]>=lev) printf("do_traverse(%d): unelim: %d (%d)\n", lev, j, subset[j]);
#endif
            }

            res[0]--;
        }
    }

    do_traverse(1,0);

    *allowed_comb_p=allowed_comb;

    return n_comb;
}

static int enumerate_cliques(char **graph, int n_sim, int ***res_p, int min_el
#ifdef COMPDESC_DEBUG
        , t_compdesc_debug *deb
#endif
)
{

    /* short indep_sets[min_el][n_sim+1]; */
    int **indep_sets=(int **)array_alloc(min_el, (n_sim+1), sizeof(int));

    /* P_POINT(indep_sets) P_NL; */

    int n_sets=fill_indep_sets(graph, n_sim, indep_sets, min_el);

    /* P_INT(n_sets) P_NL; */
#ifdef COMPDESC_DEBUG
    deb->n_sets=n_sets;
    /* P_INT(deb->n_sets) P_NL; */
    deb->indep_sets=indep_sets;
#endif

    /* P_INT(min_el) P_INT(n_sets) P_NL; */

    int *n_align=calloc(n_sets, sizeof(int));
    int ***allowed_comb=calloc(n_sets, sizeof(int **));

    for(int i=0; i<n_sets; i++) {
        n_align[i]=fill_allowed(&allowed_comb[i], indep_sets[i], graph, min_el, n_sim);
    }

#ifdef COMPDESC_DEBUG
    deb->n_allowed=n_align;
    deb->allowed_comb_per_set=allowed_comb;
    deb->stage3=1;
#endif

#ifndef COMPDESC_DEBUG
    array_free((void **)indep_sets);
#endif


    // 6. Possible alignments and their sizes
    int n_full_align=1;
    int max_align_size=0;
    int selected[n_sets];
    for(int i=0; i<n_sets; i++) {
        selected[i]=n_align[i]?0:-1;
        n_full_align*=MAX(n_align[i] , 1);

        int max_size=0;

        for(int j=0; j<n_align[i]; j++) {
            max_size=MAX(max_size, allowed_comb[i][j][0]);
        }

        max_align_size+=max_size;

        /* P_INT(i) P_INT(n_align[i]) P_INT(max_size) P_NL; */
    }

    /* P_INT(n_sim) P_NL; */

    int **res=(int **)array_alloc(n_full_align, max_align_size+1, sizeof(int));

    n_full_align=0;

    int done=0;

    while(!done) {
        for(int i=0; i<n_sets; i++) if(selected[i]>=0) {
            for(int j=1; j<=allowed_comb[i][selected[i]][0]; j++) {
                res[n_full_align][++res[n_full_align][0]]=allowed_comb[i][selected[i]][j];
            }
        }

        n_full_align++;

        int pos=0;
        for(;;) {

            if(n_align[pos]>0) {
                selected[pos]++;
                if(selected[pos]>=n_align[pos]) selected[pos++]=0;
                else break;
            } else pos++;

            if(pos>=n_sets) {
                done=1;
                break;
            }
        }
    }

#ifndef COMPDESC_DEBUG
    for(int i=0; i<n_sets; i++) {
        array_free((void **)allowed_comb[i]);
    }
    free(allowed_comb);
    free(n_align);
#endif

    for(int i=0; i<n_full_align; i++) {
        char visited[n_sim];
        char reachable[n_sim];

        int tmp_comb[max_align_size+1];

        memset(visited, 0, sizeof(visited));
        memset(reachable, 0, sizeof(reachable));
        memset(tmp_comb, 0, sizeof(tmp_comb));

        for(int j=1; j<=res[i][0]; j++) reachable[res[i][j]]=1;

        void mark_reachable(int e) {
            tmp_comb[++tmp_comb[0]]=e;
        }

        edge_graph_dfs(graph, n_sim, 0, -1, visited, reachable, mark_reachable);

        /* P_INT_ARR_L(res[i]) P_NL; */
        /* P_INT_ARR_L(tmp_comb) P_NL; */

        memcpy(res[i], tmp_comb, sizeof(tmp_comb));

    }

#ifdef COMPDESC_DEBUG
    deb->n_full_align=n_full_align;
    deb->allowed_comb=res;
#endif


    *res_p=res;
    return n_full_align;
}

int build_alignments(CDescriptor *desc1, CDescriptor *desc2, t_alignment **full_align_p, t_sim_components *sim_components, t_sim_graph *sim_graph
#ifdef COMPDESC_DEBUG
        , t_compdesc_debug *deb
#endif
        )
{
    CElement *get_zero_el(CDescriptor *desc) {
        return &desc->elements[desc->central_element_used];
    }

    int n_el1=desc1->n_elements;
    int n_el2=desc2->n_elements;

    CElement *zero_el1=get_zero_el(desc1);
    CElement *zero_el2=get_zero_el(desc2);
    int min_el=MIN(n_el1, n_el2)+1;

	// 5. Alignments
    char **graph=sim_graph->sim_graph;
    t_sim *sim_list=sim_graph->sim_list;
    int n_sim=sim_graph->n_sim;

    int **allowed_comb=0;

#ifdef COMPDESC_DEBUG
    int n_allowed=enumerate_cliques(graph, n_sim, &allowed_comb, min_el, deb);
#else
    int n_allowed=enumerate_cliques(graph, n_sim, &allowed_comb, min_el);
#endif

    t_alignment *full_align=calloc(n_allowed, sizeof(t_alignment));

    for(int i=0; i<n_allowed; i++) {
        convert_comb_align(&full_align[i], allowed_comb[i], sim_components, sim_list, desc1, desc2);

        alignment_size(&full_align[i]);
        RMSD_frag(&full_align[i], zero_el1->center_c_ind, zero_el2->center_c_ind, desc1->structure, desc2->structure);

        /* P_INT(full_align[i].n_AA) P_NL; */
    }

    *full_align_p=full_align;

#ifndef COMPDESC_DEBUG
    array_free((void **)allowed_comb);
#endif

#ifdef COMPDESC_DEBUG
    deb->full_align=calloc(n_allowed, sizeof(t_alignment));

    for(int i=0; i<n_allowed; i++) {
        copy_alignment(&(((t_alignment *)deb->full_align)[i]), &full_align[i]);
    }

    deb->stage3=2;
#endif

    return n_allowed;
}

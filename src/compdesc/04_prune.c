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
#include"04_prune.h"

int prune_alignments(CDescriptor *desc1, CDescriptor *desc2, t_alignment *full_align, int *al_order, int n_full_align, int max_n_al
#ifdef COMPDESC_DEBUG
        , t_compdesc_debug *deb
#endif
        )
{
    int max_size=0;
    int al_size_comp(const void *a, const void *b) {

        int n_a=*(int *)a;
        int n_b=*(int *)b;


        int size_dif=-(full_align[n_a].n_AA-full_align[n_b].n_AA);

        int dud_a=full_align[n_a].dud;
        int dud_b=full_align[n_b].dud;

        if(dud_a && dud_b) return n_a-n_b;

        if(dud_a) return 1;
        if(dud_b) return -1;

        double RMSD_dif=full_align[n_a].RMSD-full_align[n_b].RMSD;

        if(max_size)
            if(full_align[n_a].n_AA >= max_size - th_align_size_equiv && full_align[n_b].n_AA >= max_size - th_align_size_equiv)
                size_dif=0;

        if(full_align[n_a].RMSD>th_overall_RMSD || full_align[n_b].RMSD>th_overall_RMSD) RMSD_dif=-RMSD_dif;

        if(size_dif) return size_dif;
        else if(RMSD_dif<0.0) return -1;
        else if(RMSD_dif>0.0) return 1;
        else return -(full_align[n_a].n_AA+-full_align[n_b].n_AA);
    }


#ifdef COMPDESC_DEBUG
    int run1=1;
    deb->full_align_crop_val_run1=(float **)array_alloc(n_full_align, 6, sizeof(float));
    deb->full_align_crop_max_run1=(int **)array_alloc(n_full_align, 6, sizeof(int));
#endif

	// 7. Best alignment with qualifying RMSD

    int max_cnt=0;

    for(int i=0; i<n_full_align; i++) {
        if(!coverage_test_alignment(&(full_align[al_order[i]]), 0)) {
            full_align[al_order[i]].dud=1;
        } else {
            if(full_align[al_order[i]].n_AA>max_size) {
                max_size=full_align[al_order[i]].n_AA;
                max_cnt=1;
            } else if(full_align[al_order[i]].n_AA==max_size) {
                max_cnt++;
            }
        }
    }

    qsort(al_order, n_full_align, sizeof(int), al_size_comp);


    int max=max_n_al==0?n_full_align:MIN(max_n_al, n_full_align);

    int i=0;
    while(i<max) {
        t_al_coverage *al_cov;

        int changed=0;

        if(coverage_test_alignment(&(full_align[al_order[i]]), &al_cov)) {
            int tmp=0;

            RMSD_frag_crop(&full_align[al_order[i]], th_overall_RMSD, 5, &tmp, al_cov, 0, 0);
            alignment_size(&full_align[al_order[i]]);

            int new_pos=i;

            if(tmp) {
                if(max_size==full_align[al_order[i]].n_AA+tmp) max_cnt--;

                if(max_cnt==0) {
                    max_size=0;
                    for(int i=0; i<n_full_align; i++) {
                        if(full_align[al_order[i]].dud) break;
                        if(full_align[al_order[i]].n_AA<max_size-th_align_size_equiv) break;
                        if(full_align[al_order[i]].n_AA>max_size) {
                            max_size=full_align[al_order[i]].n_AA;
                            max_cnt=1;
                        } else if(full_align[al_order[i]].n_AA==max_size) {
                            max_cnt++;
                        }
                    }
                    qsort(al_order, n_full_align, sizeof(int), al_size_comp);
                    i=0;
                }

                for(; new_pos<n_full_align-1; new_pos++) {
                    if(al_size_comp(&al_order[new_pos], &al_order[new_pos+1])>0) {
                        SWAP(al_order[new_pos], al_order[new_pos+1]);
                    } else break;
                }
                changed=1;
            }

            if(!changed) {
                int log=0;
                for(int j=0; j<i; j++) {
                    if(is_subalignment(&(full_align[al_order[i]]), &(full_align[al_order[j]]))) {
                        full_align[al_order[i]].dud=1;
                        log=1;
                        break;
                    }
                    if(is_subalignment(&(full_align[al_order[j]]), &(full_align[al_order[i]]))) {
                        full_align[al_order[j]].dud=1;
                        log=1;
                        SWAP(i, j);
                        break;
                    }
                }
                if(log) {
                    new_pos=i;
                    for(; new_pos<n_full_align-1; new_pos++) {
                        if(al_size_comp(&al_order[new_pos], &al_order[new_pos+1])>0) {
                            SWAP(al_order[new_pos], al_order[new_pos+1]);
                            changed=1;
                        } else break;
                    }
                }
            }
        } else {
            full_align[al_order[i]].dud=1;
            break;
        }

        if(!changed) i++;

        al_coverage_free(al_cov);
    }

#ifdef COMPDESC_DEBUG
    deb->stage4=1;
#endif

    int cnt=0;

    for(int i=0; i<n_full_align; i++) {
        if(!coverage_test_alignment(&(full_align[i]), 0) || full_align[i].RMSD>th_overall_RMSD || full_align[i].dud) continue;
        cnt++;
    }


    return cnt;
}

int save_results(CDescriptor *desc1, CDescriptor *desc2, t_compdesc_result **results, t_alignment *full_align, int *al_order, int n_full_align, int max_n_al
#ifdef COMPDESC_DEBUG
        , t_compdesc_debug *deb
#endif
        )
{
#ifdef COMPDESC_DEBUG
    deb->res_align=calloc(n_full_align, sizeof(t_alignment));
#endif

    *results = calloc(n_full_align, sizeof(t_compdesc_result));

    int n_res=0;

    for(int i=0; i<n_full_align; i++) {
        t_alignment *al=&(full_align[al_order[i]]);

        t_compdesc_result *res=&((*results)[n_res]);

        if(!coverage_test_alignment(al, 0) || al->RMSD>th_overall_RMSD || al->dud) continue;
        res->RMSD=al->RMSD;
        res->TR=al->TR;

        void save_map(t_map *map, int *n, int **vec1, int **vec2) {
            *vec1=calloc(map->len1, sizeof(int));
            *vec2=calloc(map->len1, sizeof(int));

            int pos=0;
            void f(int pos1, int pos2) { 
                (*vec1)[pos] = pos1;
                (*vec2)[pos] = pos2;
                pos++;
            }

            map_iter(map, f);

            *vec1=realloc(*vec1, pos*sizeof(int));
            *vec2=realloc(*vec2, pos*sizeof(int));

            *n=pos;
        }

        save_map(al->mer_map, &res->n_monomers, &res->desc1_monomers, &res->desc2_monomers);
        save_map(al->el_map, &res->n_elements, &res->desc1_elements, &res->desc2_elements);
        save_map(al->contact_map, &res->n_contacts, &res->desc1_contacts, &res->desc2_contacts);

        for(int j=0; j<res->n_monomers; j++) {
            res->desc1_monomers[j] = map_unmap(desc1->structure->ind_to_c_map, res->desc1_monomers[j]);
            res->desc2_monomers[j] = map_unmap(desc2->structure->ind_to_c_map, res->desc2_monomers[j]);
        }


#ifdef COMPDESC_DEBUG
        copy_alignment(&(((t_alignment *)deb->res_align)[i]), &full_align[al_order[i]]);
#endif

        n_res++;
        if(max_n_al && n_res>=max_n_al) break;
    }


#ifdef COMPDESC_DEBUG
    deb->n_res_align=n_res;
    deb->stage4=2;
#endif

    *results = realloc(*results, n_res*sizeof(t_compdesc_result));
    return n_res;
}

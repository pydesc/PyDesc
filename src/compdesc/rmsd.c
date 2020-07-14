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
#include<string.h>

#include "compdesc_priv.h"
#include "simple_macros.h"
#include "cstructures.h"

static void RMSD_add_monomers(int token,  int c_ind1, int c_ind2, CStructure *struct1, CStructure *struct2)
{
    CMer *mer1=&struct1->monomers[c_ind1];
    CMer *mer2=&struct2->monomers[c_ind2];

    if(mer1->type != mer2->type) {
        printf("RMSD_add_monomers: types differ "); P_INT(mer1->type); P_INT(mer2->type); P_NL;
        P_INT(mer1->ind); P_INT(mer2->ind); P_NL;
        P_S(mer1->type_name); P_S(mer2->type_name); P_NL;
        abort();
    }

    int n_points=mer1->n_points;

    for(int i=0; i<n_points; i++) {
        float_al *p1=(float_al *)&(mer1->points[i]);
        float_al *p2=(float_al *)&(mer2->points[i]);
        /* P_FLOAT_ARR(p1, 4); P_FLOAT_ARR(p2, 4); P_NL; */
        overfit_add(token, p1, p2);
    }
}

static void RMSD_add_segment(int token,  short start1, short end1, short start2, short end2, CStructure *struct1, CStructure *struct2)
{
    if(end1-start1 != end2-start2) {
        printf("RMSD_add_segment: lengths differ "); P_INT(start1); P_INT(end1); P_INT(start2); P_INT(end2); P_NL;
        abort();
    }

    for(int pos1=start1, pos2=start2; pos1<=end1 && pos2<=end2; pos1++, pos2++) {
        RMSD_add_monomers(token, pos1, pos2, struct1, struct2);
    }
}

void RMSD_setstart(int token, int c_ind1, int c_ind2, CStructure *struct1, CStructure *struct2)
{
    float_al *p1=(float_al *)&(struct1->monomers[c_ind1].points[0]);
    float_al *p2=(float_al *)&(struct2->monomers[c_ind2].points[0]);
    overfit_setstart(token, p1, p2);
}

float RMSD_element(int token, CElement *el1, CElement *el2, CStructure *struct1, CStructure *struct2)
{
    if(el1->type!=el2->type) return 1000;


    overfit_reset(token);

    RMSD_setstart(token, el1->center_c_ind, el2->center_c_ind, struct1, struct2);

    RMSD_add_segment(token, el1->start_c_ind, el1->end_c_ind, el2->start_c_ind, el2->end_c_ind, struct1, struct2);

    return fast_overfit(token, 0);
}

float RMSD_pair(int token, t_overfit_sums *sum1, t_overfit_sums *sum2)
{
    overfit_reset(token);
    overfit_sumadd_str(token, sum1);
    overfit_sumadd_str(token, sum2);
    return fast_overfit(token, 0);

}

float RMSD_frag(t_alignment *align, int zero_el1_c_ind, int zero_el2_c_ind, CStructure *struct1, CStructure *struct2)
{
    int token=overfit_get_token();

    overfit_reset(token);

    RMSD_setstart(token, zero_el1_c_ind, zero_el2_c_ind, struct1, struct2);

    void add_to_overfit(int c_ind1, int c_ind2) {
        RMSD_add_monomers(token, c_ind1, c_ind2, struct1, struct2);
    }

    map_iter(align->mer_map, add_to_overfit);

    align->RMSD=fast_overfit(token, &(align->TR));
    overfit_release_token(token);

    return align->RMSD;
}

static int remove_from_alignment(t_alignment *align, short worst_el1_used) //, t_al_coverage *al_cov)
    /*
     * CAUTION: This function does not account for elements covered but not explicitely aligned are included in coverage.
     */
{
    CDescriptor *desc1=align->desc1;

    t_al_coverage *al_cov=align->coverage;

    int elements_to_del[desc1->n_elements];
    int n_el_del=0;

    int contacts_to_del[desc1->contact_map->n_contacts*2];  // To avoid checking if the contact has been already deleted we allow dupes. Each contact can be deleted only twice.
    int n_con_del=0;

    elements_to_del[n_el_del++]=worst_el1_used;

    int el_pos=0;
    int con_pos=0;

    int done=0;

    /* P_INT(worst_el1_used) P_NL; */

    /* P_T("remove_from_alignment\n") */

retry:
    while(!done) {
        done=1;
        for(int i=el_pos; i<n_el_del; i++) {
            /* P_INT(i) P_INT(elements_to_del[i]) P_NL; */
            for(int j=0; j<desc1->n_elements; j++) if(j!=elements_to_del[i]) {
                int con=desc1->contact_map->contact_array[elements_to_del[i]][j];

                /* P_INT(elements_to_del[i]) P_INT(j) P_INT(con) */
                /* if(con>=0) P_INT(map_map(align->contact_map, con)); */

                if(con>=0 && map_map(align->contact_map, con)>=0) {
                    contacts_to_del[n_con_del++]=con;
                    /* printf(" ***"); */
                }
                /* P_NL; */

                con=desc1->contact_map->contact_array[j][elements_to_del[i]];

                /* P_INT(j) P_INT(elements_to_del[i]) P_INT(con) */
                /* if(con>=0) P_INT(map_map(align->contact_map, con)); */

                if(con>=0 && map_map(align->contact_map, con)>=0) {
                    contacts_to_del[n_con_del++]=con;
                    /* printf(" ***"); */
                }
                /* P_NL; */
            }
        }

        el_pos=n_el_del;

        for(int i=con_pos; i<n_con_del; i++) {
            CContact *con=&desc1->contact_map->contacts[contacts_to_del[i]];

            void process_el(int el_used) {

                /* P_INT(contacts_to_del[i]) P_INT(el_used) P_INT(align->el_cov_count[el_used]); */

                if(align->el_cov_count[el_used]>0) {
                    align->el_cov_count[el_used]--;
                    if(!align->el_cov_count[el_used]) {
                        elements_to_del[n_el_del++]=el_used;
                        /* printf(" ***"); */
                        done=0;
                    }
                }

                /* P_NL; */
            }

            process_el(con->mer1_used);
            process_el(con->mer2_used);
        }

        con_pos=n_con_del;
    }

    for(int i=0; i<n_con_del; i++) {
        void update_con(int con_pos, t_desc_coverage *desc_cov) {
            if(!(desc_cov->con_vec[con_pos] & NON_OPT)) {
                desc_cov->n_con_all--;
            }

            desc_cov->n_con--;
        }

        update_con(contacts_to_del[i], &al_cov->desc1);
        update_con(map_map(align->contact_map, contacts_to_del[i]), &al_cov->desc2);

        /* CContact *con=&desc1->contact_map->contacts[contacts_to_del[i]]; */
        /* int el1=con->mer1_used; */
        /* int el2=con->mer2_used; */
        /* P_T("removing ") P_INT(el1) P_INT(el2) P_NL; */

        map_unset(align->contact_map, contacts_to_del[i], 1);
    }

    /* P_INT(n_el_del) P_NL; */

    for(int i=0; i<n_el_del; i++) {
        CElement *el=&desc1->elements[elements_to_del[i]];

        /* P_INT(el->center_used); P_NL; */

        if(!(al_cov->desc1.desc_vec[el->center_c_ind] & NON_OPT_EL)) al_cov->desc1.n_el_all--;
        if(!(al_cov->desc2.desc_vec[map_map(align->mer_map,el->center_c_ind)] & NON_OPT_EL)) al_cov->desc2.n_el_all--;

        al_cov->desc1.n_el--;
        al_cov->desc2.n_el--;
    }


    for(int i=0; i<n_el_del; i++) {
        CElement *el=&desc1->elements[elements_to_del[i]];

        for(int j=el->start_c_ind; j<=el->end_c_ind; j++) {
            if(align->mer_cov_count[j]>0) {
                align->mer_cov_count[j]--;
                if(!align->mer_cov_count[j]) {
                    if(!(al_cov->desc1.desc_vec[j] & NON_OPT)) al_cov->desc1.n_mer_all--;
                    if(!(al_cov->desc2.desc_vec[map_map(align->mer_map,j)] & NON_OPT)) al_cov->desc2.n_mer_all--;

                    al_cov->desc1.n_mer--;
                    al_cov->desc2.n_mer--;

                    map_unset(align->mer_map, j, 1);
                }
            }
        }
        /* P_T("removing ") P_INT(el->center_used)  P_NL; */
        map_unset(align->el_map, el->center_used, 1);
    }


    int comp[align->el_map->len1*2];

    for(int i=0; i<align->el_map->len1; i++) {
        comp[i]=i;
    }

    int n_comp=align->el_map->len1;

    void update_comp(int pos) {
        if(comp[pos]==pos) return;

        update_comp(comp[pos]);
        comp[pos]=comp[comp[pos]];
    }


    void process_con(int con1_pos, int con2_pos) {
        CContact *con=&desc1->contact_map->contacts[con1_pos];

        int el1_used=con->mer1_used;
        int el2_used=con->mer2_used;

        /* P_INT(el1_used) P_INT(el2_used) P_NL; */

        update_comp(el1_used);
        update_comp(el2_used);

        if(comp[el1_used]!=comp[el2_used]) {
            comp[comp[el1_used]]=n_comp;
            comp[comp[el2_used]]=n_comp;
            comp[n_comp]=n_comp;
            n_comp++;
        }
    }

    map_iter(align->contact_map, process_con);

    for(int i=0; i<desc1->n_elements; i++) {
        update_comp(i);
    }

    /* P_INT_ARR(comp, desc1->n_elements) P_NL; */

    void remove_not_connected_element(int el1_used, int el2_used) {
        if(comp[el1_used]!=comp[desc1->central_element_used]) {
            elements_to_del[n_el_del++]=el1_used;
        }
    }


    n_el_del=0;
    n_con_del=0;
    el_pos=0;
    con_pos=0;

    map_iter(align->el_map, remove_not_connected_element);

    /* P_INT_ARR(elements_to_del, n_el_del) P_NL; */

    if(n_el_del) {
        /* printf("remove_elements retry\n"); */
        goto retry;
    }




    int old_n_AA=align->n_AA;

    alignment_size(align);

    /* P_INT(old_n_AA) P_INT(align->n_AA) P_NL; */

    /*
     * Old version which marks all elements covered by alignment as aligned.
     */
    al_coverage_free(align->coverage);
    align->coverage=0;
    coverage_test_alignment(align, 0);

    return old_n_AA-align->n_AA;
}


static float RMSD_frag_crop_int(t_alignment *align, float threshold, int max_del, int *act_del, int *set_aside, int n_set_aside, int start_el_used, float *deb_val, int *deb_max)
{
    float RMSD;

    CDescriptor *desc1=align->desc1;
    CDescriptor *desc2=align->desc2;

    int zero_el1_c_ind = desc1->central_element_c_ind;
    int zero_el2_c_ind = desc2->central_element_c_ind;

    int zero_el1_used = desc1->central_element_used;

    if(max_del<0) {
        printf("RMSD_frag_crop_int: "); P_INT(max_del); printf("<0\n");
        abort();
    }
    RMSD = RMSD_frag(align, zero_el1_c_ind, zero_el2_c_ind, desc1->structure, desc2->structure);
    *act_del=0;

    int max[max_del+1];
    float val[max_del+1];
    int saved_n_set_aside_local[max_del+1];

    for(int i=0; i<=max_del; i++) {
        max[i]=-1;
        val[i]=1000;
    }


    /* printf("RMSD_frag_crop: "); P_INT(max_del) P_FLOAT(RMSD) P_NL; */

    int n_set_aside_local=0;

    if(RMSD>threshold && max_del>0) {
        void try_remove(int el1_used, int el2_used, int extend_set_aside) {
            if(el2_used<0) return;
            if(el1_used==zero_el1_used) {
                return;
            }

            t_alignment tmp_align;
            /* t_al_coverage tmp_al_cov; */

            memset(&tmp_align, 0, sizeof(tmp_align));
            /* memset(&tmp_al_cov, 0, sizeof(tmp_al_cov)); */

            copy_alignment(&tmp_align, align);
            /* copy_al_coverage(&tmp_al_cov, al_cov); */

            /* int del_count=remove_from_alignment(&tmp_align, pos2, desc1, desc2); */
            int del_count=remove_from_alignment(&tmp_align, el1_used);//, &tmp_al_cov);



            int act;

            if(del_count==0) {
                if(extend_set_aside) set_aside[n_set_aside+n_set_aside_local++]=el1_used;
            } else if(del_count<=max_del && al_coverage_valid(tmp_align.coverage)) {
                float val1=RMSD_frag_crop_int(&tmp_align, threshold, max_del-del_count, &act, set_aside, n_set_aside+n_set_aside_local, el1_used+1, 0,0);
                if(val1<val[act+del_count]) {
                    max[act+del_count]=el1_used; val[act+del_count]=val1;
                    saved_n_set_aside_local[act+del_count]=n_set_aside_local;
                }
            }

            destroy_alignment(&tmp_align);
        }

        for(int i=0; i<n_set_aside; i++) {
            try_remove(set_aside[i], map_map(align->el_map, set_aside[i]), 0);
        }

        for(int i=start_el_used; i<align->el_map->len1; i++) {
            try_remove(i, map_map(align->el_map, i), 1);
        }

        int min_c=0;
        float min_val=val[0];

        for(int j=1; j<=max_del; j++) if(val[j]<min_val) {
            min_val=val[j];
            min_c=j;
            if(min_val<=threshold) break;
        }

        int act;

        if(max[min_c]>=0) {
            int del_count=remove_from_alignment(align, max[min_c]);//, al_cov);

            RMSD=RMSD_frag_crop_int(align, threshold, max_del-del_count, &act, set_aside, n_set_aside+saved_n_set_aside_local[min_c], max[min_c]+1, 0,0);
            *act_del=act+del_count;
        }
    }
    if(deb_val) {
        for(int j=0; j<=max_del; j++) {
            deb_val[j]=val[j];
            deb_max[j]=max[j];
        }
    }

    /* P_INT(max_del) P_INT(*act_del) P_FLOAT(RMSD) P_NL; */
    return RMSD;
}

float RMSD_frag_crop(t_alignment *align, float threshold, int max_del, int *act_del, t_al_coverage *al_cov, float *deb_val, int *deb_max)
{
    int *set_aside=calloc(align->el_map->len1, sizeof(int));

    float res=RMSD_frag_crop_int(align, threshold, max_del, act_del, set_aside, 0, 0, deb_val, deb_max);
    free(set_aside);
    return res;
}

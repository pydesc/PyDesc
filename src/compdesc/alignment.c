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

#include "simple_macros.h"

#include "compdesc_priv.h"

//#define COVERAGE_CACHE_TEST

void destroy_alignment(t_alignment *align)
{
    if(align->el_map) map_free(align->el_map);
    if(align->mer_map) map_free(align->mer_map);
    if(align->contact_map) map_free(align->contact_map);
    if(align->mer_cov_count) free(align->mer_cov_count);
    if(align->el_cov_count) free(align->el_cov_count);
    if(align->coverage) al_coverage_free(align->coverage);

    memset(align, 0, sizeof(t_alignment));
}

int is_subalignment(t_alignment *al1, t_alignment *al2)
{

    if(!map_is_included(al1->mer_map, al2->mer_map)) return 0;
    if(!map_is_included(al1->el_map, al2->el_map)) return 0;
    if(!map_is_included(al1->contact_map, al2->contact_map)) return 0;

    return 1;
}

void copy_alignment(t_alignment *dest, t_alignment *align)
{
    destroy_alignment(dest);


    if(align->el_map) {
        dest->el_map=map_alloc(align->el_map->len1, align->el_map->len2);
        void update_el(int key, int val) {map_set(dest->el_map, key, val);}
        map_iter(align->el_map, update_el);
    }

    if(align->mer_map) {
        dest->mer_map=map_alloc(align->mer_map->len1, align->mer_map->len2);
        void update_mer(int key, int val) {map_set(dest->mer_map, key, val);}
        map_iter(align->mer_map, update_mer);
    }

    if(align->contact_map) {
        dest->contact_map=map_alloc(align->contact_map->len1, align->contact_map->len2);
        void update_contact(int key, int val) {map_set(dest->contact_map, key, val);}
        map_iter(align->contact_map, update_contact);
    }

    if(align->mer_cov_count) {
        dest->mer_cov_count=calloc(align->mer_map->len1, sizeof(int));
        memcpy(dest->mer_cov_count, align->mer_cov_count, align->mer_map->len1*sizeof(int));
    }

    if(align->el_cov_count) {
        dest->el_cov_count=calloc(align->el_map->len1, sizeof(int));
        memcpy(dest->el_cov_count, align->el_cov_count, align->el_map->len1*sizeof(int));
    }

    if(align->coverage) {
        dest->coverage=calloc(1, sizeof(t_al_coverage));
        copy_al_coverage(dest->coverage, align->coverage);
    }

    dest->desc1=align->desc1;
    dest->desc2=align->desc2;
    dest->n_AA=align->n_AA;
    dest->RMSD=align->RMSD;
    dest->TR=align->TR;
}

void alignment_size(t_alignment *align)
{
    int n_AA=0;
    void count(int c_ind1, int c_ind2) {
        n_AA++;
    }

    map_iter(align->mer_map, count);

    align->n_AA=n_AA;
}

static int desc_coverage_valid(t_desc_coverage *desc_cov)
{
    int check_one(int n, int n_all, float th) {
        if((float)n/n_all>=th || n_all==0) return 0;
        return 1;
    }

    if(check_one(desc_cov->n_mer, desc_cov->n_mer_all, th_mer_align)) return 0;
    if(check_one(desc_cov->n_el, desc_cov->n_el_all, th_el_align)) return 0;
    if(check_one(desc_cov->n_con, desc_cov->n_con_all, th_con_align)) return 0;

    return 1;
}

int al_coverage_valid(t_al_coverage *al_cov)
{
    if(!desc_coverage_valid(&al_cov->desc1)) return 0;
    if(!desc_coverage_valid(&al_cov->desc2)) return 0;

    return 1;
}

static void create_desc_vector(int *desc_vec, CDescriptor *desc)
{
    for(int i=0; i<desc->n_elements; i++) {
        CElement *el=&desc->elements[i];

        int nopt=0;
        int type;
        if(i==desc->central_element_used) {
            type=ZERO;
        } else {
            type=CENTRAL;
            nopt=1-el->optional;
        }

        for(int j=el->start_c_ind; j<=el->end_c_ind; j++) {
            if(j==el->center_c_ind) {
                desc_vec[j] |= type | nopt*(NON_OPT_EL | NON_OPT) | IN_DESC;
            } else {
                desc_vec[j] |= nopt*(NON_OPT) | IN_DESC;
            }
        }
    }
}

static void create_con_vector(int *con_vec, CDescriptor *desc)
{

    for(int i=0; i<desc->contact_map->n_contacts; i++) {
        CContact *con=&desc->contact_map->contacts[i];

        int nopt=con->val==2 ? 1 : 0;

        if(nopt) con_vec[i] |= NON_OPT;
    }
}

static void calculate_coverage(CDescriptor *desc, t_desc_coverage *desc_cov)
{
    int n_mer_u=0, n_mer_all=0;
    int n_el_u=0, n_el_all=0;
    int n_con_u=0, n_con_all=0;

    int n_AA=desc_cov->n_AA;
    int n_contacts=desc_cov->n_contacts;

    int *desc_vec=desc_cov->desc_vec;
    int *con_vec=desc_cov->con_vec;

    if(n_AA<3) {
        P_INT(n_AA); printf("<3\n");
        abort();
    }


    for (int i=0;i<n_AA;i++) {
        if(desc_vec[i] & ALIGNED) {
            n_mer_u++;
            if((desc_vec[i] & ALIGNED_EL) && (desc_vec[i] & (CENTRAL | ZERO))) n_el_u++;

        }

        if(desc_vec[i] & (NON_OPT | ALIGNED)) {
            n_mer_all++;
            if((desc_vec[i] & (NON_OPT_EL | ALIGNED_EL)) && (desc_vec[i] & (CENTRAL | ZERO))) n_el_all++;
        }

    }

    for(int i=0; i<n_contacts; i++) {
        if(con_vec[i] & ALIGNED) n_con_u++;
        if(con_vec[i] & (NON_OPT | ALIGNED)) n_con_all++;
    }


    desc_cov->n_mer=n_mer_u;
    desc_cov->n_mer_all=n_mer_all;
    desc_cov->n_el=n_el_u;
    desc_cov->n_el_all=n_el_all;
    desc_cov->n_con=n_con_u;
    desc_cov->n_con_all=n_con_all;

}

static void fill_desc_coverage(CDescriptor *desc, t_desc_coverage *desc_cov)
{
    int n_AA=desc->structure->n_monomers;
    int n_contacts=desc->contact_map->n_contacts;

    desc_cov->n_AA=n_AA;
    desc_cov->n_contacts=n_contacts;

    desc_cov->desc_vec=calloc(n_AA, sizeof(int));
    desc_cov->con_vec=calloc(n_contacts, sizeof(int));

    create_desc_vector(desc_cov->desc_vec, desc);
    create_con_vector(desc_cov->con_vec, desc);
}

static void copy_desc_coverage(t_desc_coverage *dest, t_desc_coverage *src)
{
    *dest=*src;

    dest->desc_vec=calloc(dest->n_AA, sizeof(int));
    dest->con_vec=calloc(dest->n_contacts, sizeof(int));

    memcpy(dest->desc_vec, src->desc_vec, dest->n_AA*sizeof(int));
    memcpy(dest->con_vec, src->con_vec, dest->n_contacts*sizeof(int));
}

void copy_al_coverage(t_al_coverage *dest, t_al_coverage *src)
{
    copy_desc_coverage(&dest->desc1, &src->desc1);
    copy_desc_coverage(&dest->desc2, &src->desc2);
}

#ifdef COVERAGE_CACHE_TEST
static int is_equal_desc_coverage(t_desc_coverage *dcov1, t_desc_coverage *dcov2)
{
    int res=1;

    if(dcov1->n_AA != dcov2->n_AA) res=0;
    if(dcov1->n_contacts != dcov2->n_contacts) res=0;
    if(dcov1->n_mer != dcov2->n_mer) res=0;
    if(dcov1->n_mer_all != dcov2->n_mer_all) res=0;
    if(dcov1->n_el != dcov2->n_el) res=0;
    if(dcov1->n_el_all != dcov2->n_el_all) res=0;
    if(dcov1->n_con != dcov2->n_con) res=0;
    if(dcov1->n_con_all != dcov2->n_con_all) res=0;

    if(!res) {
        P_INT(dcov1->n_AA) P_INT(dcov2->n_AA) P_NL;
        P_INT(dcov1->n_contacts) P_INT(dcov2->n_contacts) P_NL;
        P_INT(dcov1->n_mer) P_INT(dcov2->n_mer) P_NL;
        P_INT(dcov1->n_mer_all) P_INT(dcov2->n_mer_all) P_NL;
        P_INT(dcov1->n_el) P_INT(dcov2->n_el) P_NL;
        P_INT(dcov1->n_el_all) P_INT(dcov2->n_el_all) P_NL;
        P_INT(dcov1->n_con) P_INT(dcov2->n_con) P_NL;
        P_INT(dcov1->n_con_all) P_INT(dcov2->n_con_all) P_NL;
    }

    return res;
}

static int is_equal_al_coverage(t_al_coverage *cov1, t_al_coverage *cov2)
{
    int r1=is_equal_desc_coverage(&cov1->desc1, &cov2->desc1);
    int r2=is_equal_desc_coverage(&cov1->desc2, &cov2->desc2);

    return r1&&r2;

}
#endif

static void destroy_desc_coverage(t_desc_coverage *desc_cov)
{
    if(desc_cov->desc_vec) free(desc_cov->desc_vec);
    if(desc_cov->con_vec) free(desc_cov->con_vec);
}

void destroy_al_coverage(t_al_coverage *al_cov)
{
    destroy_desc_coverage(&al_cov->desc1);
    destroy_desc_coverage(&al_cov->desc2);
}

void al_coverage_free(t_al_coverage *al_cov)
{
    destroy_al_coverage(al_cov);
    free(al_cov);
}

static int coverage_test_int(CDescriptor *desc1, CDescriptor *desc2, void mark_aligned(t_al_coverage *), t_al_coverage **coverage_p)
{
    t_al_coverage *cov=calloc(1, sizeof(t_al_coverage));

    fill_desc_coverage(desc1, &cov->desc1);
    fill_desc_coverage(desc2, &cov->desc2);

    mark_aligned(cov);

    calculate_coverage(desc1, &(cov->desc1));
    calculate_coverage(desc2, &(cov->desc2));

    int res=al_coverage_valid(cov);

    if(coverage_p) *coverage_p=cov;
    else {
        al_coverage_free(cov);
    }

    return res;
}

int coverage_test_prelim(CDescriptor *desc1, CDescriptor *desc2, t_al_coverage **coverage_p)
{
    void mark_descriptors(t_al_coverage *cov) {
        void mark_desc(t_desc_coverage *dcov) {
            for(int i=0; i<dcov->n_AA; i++) {
                if(dcov->desc_vec[i] & IN_DESC && ! (dcov->desc_vec[i] & NON_OPT)) {
                    dcov->desc_vec[i] |= ALIGNED;
                }

                if(dcov->desc_vec[i] & IN_DESC && ! (dcov->desc_vec[i] & NON_OPT_EL)) {
                    if(dcov->desc_vec[i] & (ZERO | CENTRAL)) dcov->desc_vec[i] |= ALIGNED_EL;
                }
            }

            for(int i=0; i<dcov->n_contacts; i++) {
                if(! (dcov->con_vec[i] & NON_OPT)) {
                    dcov->con_vec[i] |= ALIGNED;
                }
            }
        }

        mark_desc(&(cov->desc1));
        mark_desc(&(cov->desc2));
    }

    t_al_coverage *cov;

    coverage_test_int(desc1, desc2, mark_descriptors, &cov);

    void recalculate_prelim_cov(int *n1, int *n_all1, int *n2, int *n_all2) {
        int tmp_n1 = MIN(*n_all1, *n_all2);
        int tmp_n2 = MIN(*n_all1, *n_all2);
        int tmp_n_all1 = *n_all1 - *n1;
        int tmp_n_all2 = *n_all2 - *n2;

        *n1 = tmp_n1;
        *n2 = tmp_n2;
        *n_all1 = tmp_n_all1;
        *n_all2 = tmp_n_all2;
    }

    recalculate_prelim_cov(&(cov->desc1.n_mer), &(cov->desc1.n_mer_all), &(cov->desc2.n_mer), &(cov->desc2.n_mer_all));
    recalculate_prelim_cov(&(cov->desc1.n_el), &(cov->desc1.n_el_all), &(cov->desc2.n_el), &(cov->desc2.n_el_all));
    recalculate_prelim_cov(&(cov->desc1.n_con), &(cov->desc1.n_con_all), &(cov->desc2.n_con), &(cov->desc2.n_con_all));

    /*
     * Old version which marks all elements covered by alignment as aligned.
     */

    cov->desc1.n_el=cov->desc1.n_el_all;
    cov->desc2.n_el=cov->desc2.n_el_all;

    int res=al_coverage_valid(cov);

    if(coverage_p) *coverage_p=cov;
    else {
        al_coverage_free(cov);
    }

    return res;
}

int coverage_test_component(CDescriptor *desc1, CDescriptor *desc2, t_sim_components *sim_components, t_al_coverage **coverage_p)
{
    void mark_components(t_al_coverage *cov) {
        void mark_elsim(t_elsim *elsim) {
            void mark_el(t_desc_coverage *dcov, CElement *el) {
                for(int j=el->start_c_ind; j<=el->end_c_ind; j++) {
                    dcov->desc_vec[j] |= ALIGNED;
                    if(j==el->center_c_ind) {
                        dcov->desc_vec[j] |= ALIGNED_EL;
                    }
                }
            }

            mark_el(&cov->desc1, elsim->el1);
            mark_el(&cov->desc2, elsim->el2);
        }

        void mark_consim(t_consim *consim) {
            cov->desc1.con_vec[consim->con1->pos] |= ALIGNED;
            cov->desc2.con_vec[consim->con2->pos] |= ALIGNED;
        }

        for(int i=0; i<sim_components->n_elsim; i++) mark_elsim(&sim_components->elsim_list[i]);
        for(int i=0; i<sim_components->n_consim; i++) mark_consim(&sim_components->consim_list[i]);
    }

    t_al_coverage *cov;

    coverage_test_int(desc1, desc2, mark_components, &cov);

    /*
     * Old version which marks all elements covered by alignment as aligned.
     */

    cov->desc1.n_el=cov->desc1.n_el_all;
    cov->desc2.n_el=cov->desc2.n_el_all;

    int res=al_coverage_valid(cov);

    if(coverage_p) *coverage_p=cov;
    else {
        al_coverage_free(cov);
    }

    return res;
}

int coverage_test_alignment(t_alignment* alignment, t_al_coverage **coverage_p)
{
    CDescriptor *desc1=alignment->desc1;
    CDescriptor *desc2=alignment->desc2;

    void mark_alignment(t_al_coverage *cov) {
        void mark_aligned(int pos1, int pos2) {
            cov->desc1.desc_vec[pos1] |= ALIGNED;
            cov->desc2.desc_vec[pos2] |= ALIGNED;
        }

        map_iter(alignment->mer_map, mark_aligned);

        /*
         * Old version which marks all elements covered by alignment as aligned.
         */

        int seg_start=-1;
        int seg_end=-1;
        int el=0;

        /*
        * CAUTION: This code does not work if there exist elements el1, e2
        * such that:
        * el1.central_c_ind<el2.central_c_ind
        * el1.end_c_ind>el2.end_c_ind
        */

        void mark_aligned_el(int used, CDescriptor *desc, int *d_vec) {
            if(desc->elements[used].start_c_ind>seg_end) {
                seg_start=desc->elements[used].start_c_ind;
                seg_end=desc->elements[used].end_c_ind;
            } else {
                seg_end=desc->elements[used].end_c_ind;
            }

            while(el<desc->n_elements) {
                CElement *element=&desc->elements[el];

                if(element->end_c_ind>seg_end) break;

                if(element->start_c_ind>=seg_start) {
                    d_vec[element->center_c_ind] |= ALIGNED_EL;
                }
                el++;
            }
        }

        void mark_aligned_el1(int used1, int used2) {
            mark_aligned_el(used1, desc1, cov->desc1.desc_vec);
        }

        void mark_aligned_el2(int used1, int used2) {
            mark_aligned_el(used2, desc2, cov->desc2.desc_vec);
        }

        map_iter(alignment->el_map, mark_aligned_el1);

        seg_start=-1; seg_end=-1; el=0;

        map_iter_rev(alignment->el_map, mark_aligned_el2);

        void mark_aligned_con(int pos1, int pos2) {
            cov->desc1.con_vec[pos1] |= ALIGNED;
            cov->desc2.con_vec[pos2] |= ALIGNED;
        }

        map_iter(alignment->contact_map, mark_aligned_con);
    }

    int res=0;

#ifdef COVERAGE_CACHE_TEST
    t_al_coverage *old_cov=0;

    if(alignment->coverage) {
        old_cov=alignment->coverage;
        alignment->coverage=0;
    }

    res=coverage_test_int(desc1, desc2, mark_alignment, &alignment->coverage);

    if(old_cov) {
        if(is_equal_al_coverage(old_cov, alignment->coverage)) {
        } else {
            printf("********* coverages differ ************************\n");
            abort();
        }


        al_coverage_free(old_cov);
    }
#else
    if(!alignment->coverage) {
        res=coverage_test_int(desc1, desc2, mark_alignment, &alignment->coverage);
    } else {
        res=al_coverage_valid(alignment->coverage);
    }

#endif


    if(coverage_p) {
        *coverage_p=calloc(1, sizeof(t_al_coverage));

        copy_al_coverage(*coverage_p, alignment->coverage);
    }

    return res;
}

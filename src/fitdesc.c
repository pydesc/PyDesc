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
 
#include "fitdesc.h"
#include "simple_macros.h"
#include<stdlib.h>
#include<stdio.h>
#include<assert.h>
#include<string.h>
#include<math.h>

//#define MD_TEST

static void getMapOrderingSegmentsAccordingToLength(int* seg_map, const CStructure* desc) // checked
{
    int n_segs = desc->n_segs;
    int* lengths = (int*) malloc(2*n_segs*sizeof(int));

    for(int i=0;i<n_segs;i++) {
        lengths[2*i] = desc->segs[i].end - desc->segs[i].start + 1;
        lengths[2*i+1] = i;
    }

    int compareIntsDecreasing(const void* a, const void* b) {
        return( *((int*)b) - *((int*)a) );
    }

    qsort(lengths, n_segs, 2*sizeof(int), (comparison_fn_t)compareIntsDecreasing);

    for(int i=0;i<n_segs;i++) {
        seg_map[i] = lengths[2*i+1];
    }

    free(lengths);
}

static void initializeTree(t_segfit** tree, int* n_branches, int token) // checked
{
    tree[0] = (t_segfit*) malloc( 1 * sizeof(t_segfit) );
    n_branches[0] = 1; 
    // I'm creating an element "tree[0][0]" which doesn't store any segment's fit
    // but it helps to initialize the first (i.e. "t_segfit* tree[1]") "branch" with the 
    // fit of segment desc->segs[seg_map[0]]

    tree[0][0].source = NULL;

    // TODO: czy ponizsze jest w ogole do czegokolwiek potrzebne
    //	tree[0][0].start_in_dom_ind = -1;
    //	tree[0][0].start_in_dom_mer_array = -1;
    //	tree[0][0].length = 0;
    //	tree[0][0].seg_num = 0;

    overfit_reset(token);
    overfit_sumsave_str(token, &(tree[0][0].sums));
}

static int isConflict(int begin1, int length1, int begin2, int length2) // checked
{
    if( begin1>begin2 && (begin1-begin2>=length2))
        return 0; // meanning: there's no conflict, this fit is fine
    if( begin2>begin1 && (begin2-begin1>=length1))
        return 0; // same as above

    return 1; // this fit is wrong
}


/*
   "noClashesWithAlreadyAligned" returns "1" if the segment with first amino-acid "aa" (taken from the domain) has no clashes with the already assigned segments
   "0" otherwise
   */
static int noClashesWithAlreadyAligned(t_segfit** tree, int height, int branch, int seg_len, int aa) // checked
{
    t_segfit tmp = tree[height][branch];

    while( tmp.source!=0 ) { // until it reaches NULL (check initializeTree)
        if( isConflict(aa,seg_len,tmp.start_in_dom_mer_array,tmp.length) == 1 )
            return 0;
        else
            tmp = *tmp.source;
    }

    return 1;
}

static int incompatibleAAtypesOccur(const CStructure* dom, int aa, const CStructure* desc, int seg_num, int seg_len)
{
    int start = desc->segs[seg_num].c_start;
    //printf("start = %d\n",start);

    for(int i=0;i<seg_len;i++)
        if( dom->monomers[aa+i].type != desc->monomers[start+i].type ) {
            //printf("dom->mers[aa+i].ind = %d, desc->mers[start+i].ind = %d\n",dom->monomers[aa+i].ind,desc->monomers[start+i].ind);
            return 0;
        }

    return 1;
}

static void overfitAddMers(int token, CMer* desc_mer, CMer* dom_mer) // checked
{
    assert(desc_mer->type==dom_mer->type); // double check

    for(int i=0;i<desc_mer->n_points;i++)
        overfit_add(token, (float*) &(dom_mer->points[i]), (float*) &(desc_mer->points[i]));
}

#ifdef MD_TEST
static void print_CMer(CMer *mer)
{
    printf("ind: %5d \ttype: %5d \tnext_ind: %5d \tc_ind: %5d \tnext_c_ind: %5d\n",mer->ind,mer->type,mer->next_ind,mer->c_ind,mer->next_c_ind);

}
static void print_CStructure(const CStructure* str)
{
    P_INT(str->n_monomers); P_NL;

    P_INT(str->n_segs); P_NL;
    for(int i=0; i<str->n_segs; i++) {
        P_INT(str->segs[i].start);
        P_INT(str->segs[i].end);
        P_NL;
    }

    for(int i=0; i<str->n_monomers; i++)
        print_CMer(&str->monomers[i]);
}
#endif // MD_TEST

static void addToSums(int seg_num, int seg_len, int aa, const CStructure* desc, const CStructure* dom, int token) // checked
{
    int start = desc->segs[seg_num].c_start;
    for(int j=0;j<seg_len;j++)
        overfitAddMers(token,&(desc->monomers[start+j]),&(dom->monomers[aa+j])); // this is needed by "fast_overfit" in order to find RMSD
}

static void reallocTree(t_segfit** tree, int* n_branches, int new_height) // checked
{
    if( n_branches[new_height] == 0) {
        n_branches[new_height] = 1; 	// so it's easier to realloc twice the current size (otherwise, the current size would be 0)
        tree[new_height] = NULL;		// realloc(NULL,5*sizeof(t_segfit)) (for example) is the same as malloc(5*sizeof(t_segfit))
    }

    t_segfit* tmp = (t_segfit*) realloc( tree[new_height], 2*n_branches[new_height]*sizeof(t_segfit));
    // now "tree[new_height]" no longer points at NULL

    if(!(tmp))	{
        free(tree[new_height]);
        perror("realloc");
        exit(EXIT_FAILURE);
    } else {
        n_branches[new_height] *= 2;
        tree[new_height] = tmp;
    }
}

static int* getLengthMap(int* length_map, const CStructure* dom)
{
    int tmp = 0;
    for(int aa=dom->n_monomers-1;aa>=0;aa--) { // changed from "aa>0" to "aa>=0"
        int prev_number = dom->monomers[aa].next_ind; // artificial, but works
        int curr_number = dom->monomers[aa].ind;
        if(prev_number-curr_number==1)
            tmp += 1;
        else
            tmp = 1;

        length_map[aa] = tmp;
    }

    return length_map;
}

static void makeBranch(t_segfit** tree, int* n_branches, int i, int seg_num, int* length_map, const CStructure* desc, const CStructure* dom, float r, int token) // TODO
{

    int height = i;
    int new_branch_number = 0;	// to know how many branches were made so far (see below), i.e. tree[height][new_branch_number] should be the next created "branch"

    /*
     * "forOneSegfit" takes "tree[height][branch]" and creates an yet unknown number of branches
     * on the level "new_height" of "tree"
     * the size of the level "new_height" is being reallocated during the procedure
     *
     * I defined this function here, bacause it would require a lot of arguments otherwise, and it's being used only here
     */
    int forOneSegfit(int branch) {
        t_transrot TR;
        int new_height = height+1;
        int seg_len = desc->segs[seg_num].end - desc->segs[seg_num].start + 1; 

        /*printf("seg_len = %d\n",seg_len);	
          printf("length_map = "); for(int i=0;i<dom->n_monomers;i++) printf("%d ",length_map[i]); printf("\n");
          printf("incompatibleAAtypesOccur(dom,0,desc,seg_num,seg_len) = %d",incompatibleAAtypesOccur(dom,0,desc,seg_num,seg_len) );*/

        int d1 = tree[height][branch].sums.N + seg_len;
        int d2 = desc->n_monomers;

        float r_adj = r * sqrtf((float)d2 / (float)d1);

        for( int aa = 0; aa < dom->n_monomers-seg_len+1; aa++) {
            // "aa" is the amino-acid number in domain "dom", to which segment "Si" from "desc" will be alligned to
            // so that (aa,aa+1,aa+2) will fit to an examplary segment (length 3)
            if( length_map[aa] < seg_len || incompatibleAAtypesOccur(dom,aa,desc,seg_num,seg_len)==0 || noClashesWithAlreadyAligned(tree, height, branch, seg_len, aa)==0 )
                continue;

            /* printf("elo\n"); */
            overfit_reset(token);
			addToSums(seg_num,seg_len,aa,desc,dom,token); // see the definition of "addToSums"

			float RMSDForOneSeg = fast_overfit(token,&TR); // see "overfit.c"

/* printf("RMSDForOneSeg = %f.3\n",RMSDForOneSeg); */
			if( RMSDForOneSeg > r )
				continue;

			overfit_sumadd_str(token, &(tree[height][branch].sums));

			float RMSDForCurrentSegs = fast_overfit(token, &TR);
			if(RMSDForCurrentSegs > r_adj)
                continue;

            if(seg_len == 1) {
                CMer *mer1 = &(desc->monomers[desc->segs[seg_num].c_start]);
                CMer *mer2 = &(dom->monomers[aa]);

                if(mer1->n_points == 1) {
                    if(overfit_trans_dist((float*) &(mer1->points[0]), (float*) &(mer2->points[0]), TR) > r) {
                        continue;
                    }
                }
            }

            if(new_branch_number>=n_branches[new_height])
                reallocTree(tree,n_branches,new_height);

            overfit_sumsave_str(token, &(tree[new_height][new_branch_number].sums));

            tree[new_height][new_branch_number].start_in_dom_ind = dom->monomers[aa].ind;
            tree[new_height][new_branch_number].start_in_dom_mer_array = aa;
            tree[new_height][new_branch_number].length = seg_len;
            tree[new_height][new_branch_number].seg_num = seg_num;

            tree[new_height][new_branch_number].source = &(tree[height][branch]);

            tree[new_height][new_branch_number].RMSD = RMSDForCurrentSegs;
            tree[new_height][new_branch_number].TR = TR;

            int start_num = desc->segs[seg_num].c_start;
            tree[new_height][new_branch_number].start_in_desc_ind = desc->monomers[start_num].ind;
            tree[new_height][new_branch_number].start_in_desc_mer_array = start_num;

            new_branch_number++;
		}

		return new_branch_number;
	}


	for(int branch=0; branch < n_branches[height]; branch++)
		new_branch_number = forOneSegfit(branch);

	n_branches[height+1] = new_branch_number; 	// height+1 is new_height; newly created level of "tree" has "new_branch_number" number of branches
}

static void getWholeFit(t_fitdesc_result* res, t_segfit* leaf, CStructure* desc)
{
	t_segfit tmp = *leaf;

	res->RMSD = tmp.RMSD;
    res->TR=tmp.TR;

	while(tmp.source) {
		int desc_start_from_zero = tmp.start_in_desc_mer_array;
		for(int i=0;i<tmp.length;i++) {
			res->motif_monomers[desc_start_from_zero+i] = tmp.start_in_desc_ind+i;
			res->structure_monomers[desc_start_from_zero+i] = tmp.start_in_dom_ind+i;
		}

		tmp = *(tmp.source);
	}
}

static void allocResult(t_fitdesc_result* result, int n_mers)
{
	result->n_monomers = n_mers;
	result->motif_monomers = (int*) calloc(n_mers, sizeof(int));
	result->structure_monomers = (int*) calloc(n_mers, sizeof(int));
}

static void storeFits(t_segfit* leaves, int final_n_fits, CStructure* desc,t_fitdesc_result** results)
{
	//printf("final_n_fits = %d\n",final_n_fits);
	t_fitdesc_result* tmp = (t_fitdesc_result*) malloc(final_n_fits * sizeof(t_fitdesc_result));

	for( int leafNum = 0; leafNum < final_n_fits; leafNum++ ) {
		allocResult(&tmp[leafNum],desc->n_monomers);
		getWholeFit(&tmp[leafNum],&leaves[leafNum],desc);
	}

	*results = &tmp[0];
}

static void destroyTree(t_segfit** tree, int tree_height, int* n_branches)//checked
{
	for(int i=0;i<tree_height;i++) {
		if(n_branches[i]!=0)
			free(tree[i]);
	}

	free(tree);
}

#ifdef MD_TEST
static void print_t_fitdesc_result(t_fitdesc_result* result) {
	printf("RMSD = %.3f\n",result->RMSD);
	for(int i=0;i<result->n_monomers;i++) 
		printf("%3d --- %3d\n",result->motif_monomers[i],result->structure_monomers[i]);

}
#endif

/////////////////////////////////////////////////////////////////////////////////////////////
///////////////////  NON - STATIC   FUNCTIONS :  ////////////////////////////////////////////

int fitdesc(CStructure* desc, CStructure* dom, float r, int I, t_fitdesc_result** results)
{
#ifdef MD_TEST
	printf("\nPRINTING DESC:\n"); print_CStructure(desc);
	printf("\nPRINTING DOM:\n"); print_CStructure(dom);
	printf("\n\n");
#endif // MD_TEST
	int token = overfit_get_token();

	int tree_height = desc->n_segs+1; // one extra level of tree is used for initialization of the whole tree (it doesn't store any of the segment's fit)
	t_segfit** tree = (t_segfit**) malloc(tree_height * sizeof(t_segfit*));
	int n_branches[tree_height];
	for(int i=1;i<tree_height;i++) // i==0: see "initializeTree"
		n_branches[i] = 0;
	// "n_branches[i]" will store the number of branches on the i-th level of tree

	initializeTree(tree,n_branches,token);

	int seg_map[desc->n_segs];
	getMapOrderingSegmentsAccordingToLength(seg_map,desc); // for ordering the segments according to their length
/* for(int i=0;i<desc->n_segs;i++) printf("seg_map[%d] = %d\n",i,seg_map[i]); */

	int length_map[dom->n_monomers];
	getLengthMap(length_map,dom); // for amino-acid number "aa", "length_map[aa]" is the maximum length of a segment, which can be aligned starting from "aa" in "dom"


	for(int i=0;i<desc->n_segs;i++) {
		int seg_num = seg_map[i];	//so that fitting starts from the longest segment
		makeBranch(tree,n_branches,i,seg_num,length_map,desc,dom,r,token);	// here "i" is the same as tree's "height"
	}
	// after checking whether any acceptable fits were found (and if so - how many),  we prepare the output:

	// final number of fits found:
	int final_n_fits = n_branches[tree_height-1];
	t_segfit* leaves = tree[tree_height-1];

	if(final_n_fits != 0) {
		///////////////////////BEGIN SORTING(leaves)//////////////////////////////////////
		int compareRMSDs(const t_segfit *a, const t_segfit *b) { // for increasing order
			if((a->RMSD - b->RMSD) > 0)
				return 1;
			if((a->RMSD - b->RMSD) == 0)
				return 0;
			else
				return -1;
		}

		qsort(leaves, final_n_fits, sizeof(t_segfit), (comparison_fn_t)compareRMSDs);
		//////////////////////END SORTING(leaves)/////////////////////////////////////////////

		if(I!=0 && final_n_fits >= abs(I) )
			final_n_fits = I;

		storeFits(leaves,final_n_fits,desc,results); // RESULTS are being produced HERE
	}


	destroyTree(tree,tree_height,n_branches);

	overfit_release_token(token);

#ifdef MD_TEST
	printf("PRINT-OUT:\n");
	printf("dom->n_monomers = %d\n",dom->n_monomers);
	printf("\nPRINTING DESC:\n"); print_CStructure(desc);
	printf("\nPRINTING DOM:\n"); print_CStructure(dom);
	printf("\n\n");
	printf("final_n_fits = %d\n",final_n_fits);
	if( final_n_fits>0 ) {
		printf("PROPER:\n");
		print_t_fitdesc_result(&results[0][0]);
	}
	
	if( final_n_fits>1 ) {
		printf("ERRORNOUS:\n");
		print_t_fitdesc_result(&results[0][1]);
	}
	
	printf("########################################################################################\n");
	printf("########################################################################################\n");
	printf("########################################################################################\n");
#endif // MD_TEST

  return final_n_fits;
}

void free_results(t_fitdesc_result* res, int n)
{
    for(int i=0; i<n; i++) {
        free(res[i].motif_monomers);
        free(res[i].structure_monomers);
    }

    free(res);
}


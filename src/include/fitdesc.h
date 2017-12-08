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
 
#ifndef __FITDESC_H__
#define __FITDESC_H__

#include "cstructures.h"
#include "overfit.h"



/*
 *  "struct t_segfit" stores <<redundant>> information (explained below) about fitting of a single segment;
 *	it's used only by fitdesc and after it's being used to produce output, the memory is being relieved
 *	"t_segfit" is used to form a data structure explained in "fitdescDoc.pdf";
 *	every "t_segfit" has a "t_segfit* source"; meaning: when fitting the i-th segment,
 *  and having previous (i-1) segments already fitted, we need to know how were the previous segments
 *	fitted - knowledge of their position and so-called sums ("t_overfit_sums sums") (below) is required 
 *	in order to find if the fit is correct
 *	"t_overfit_sums sums" is used by Kabsch algorith implemented in "overfit.c" (where I would like to 
 *  refer interested reader :>) and it's being stored in order to find RMSD
 *	"t_transrot TR" is stored only for for output (see "t_group_output" in "desc.h")
 *
 */
typedef struct t_segfit{
  int start_in_dom_ind; 		// starting aa number in domain
  int start_in_dom_mer_array; // as above, but here it's the relative number
  int length;   	// fitted segment's length
  int seg_num; 		// segment's number in seed desc
  int start_in_desc_ind; 	// starting aa number of segment in descriptor
  int start_in_desc_mer_array; // as above, but here it's the relative number
  

  struct t_segfit* source; // pointing at the t_segfit's source
  
  t_overfit_sums sums;
  
  float RMSD;
  t_transrot TR; 		
  
} t_segfit;


typedef struct t_fitdesc_result {
  int n_monomers;
  int* motif_monomers;
  int* structure_monomers;

  float RMSD;
  t_transrot TR;
} t_fitdesc_result;

int fitdesc(CStructure* desc, CStructure* dom, float r, int I, t_fitdesc_result** results);
void free_results(t_fitdesc_result* res, int n);

#endif

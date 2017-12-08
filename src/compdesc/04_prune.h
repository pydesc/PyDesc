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

#ifndef __04_PRUNE_H__
#define __04_PRUNE_H__

#include "compdesc_priv.h"

#ifdef COMPDESC_DEBUG
int prune_alignments(CDescriptor *desc1, CDescriptor *desc2, t_alignment *full_align, int *al_order, int n_full_align, int max_n_al, t_compdesc_debug *deb);
int save_results(CDescriptor *desc1, CDescriptor *desc2, t_compdesc_result **results, t_alignment *full_align, int *al_order, int n_full_align, int max_n_al, t_compdesc_debug *deb);
#else
int prune_alignments(CDescriptor *desc1, CDescriptor *desc2, t_alignment *full_align, int *al_order, int n_full_align, int max_n_al);
int save_results(CDescriptor *desc1, CDescriptor *desc2, t_compdesc_result **results, t_alignment *full_align, int *al_order, int n_full_align, int max_n_al);
#endif

#endif

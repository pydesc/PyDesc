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

#ifndef __03_BUILD_H__
#define __03_BUILD_H__

#include "compdesc_priv.h"
#include "02_components.h"

#ifdef COMPDESC_DEBUG
int build_alignments(CDescriptor *desc1, CDescriptor *desc2, t_alignment **full_align_p, t_sim_components *sim_components, t_sim_graph *sim_graph, t_compdesc_debug *deb);
#else
int build_alignments(CDescriptor *desc1, CDescriptor *desc2, t_alignment **full_align_p, t_sim_components *sim_components, t_sim_graph *sim_graph);
#endif

#endif

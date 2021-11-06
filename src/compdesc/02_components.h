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

#ifndef __02_COMPONENTS_H__
#define __02_COMPONENTS_H__

#include "compdesc_priv.h"



#ifdef COMPDESC_DEBUG
t_sim_components *compare_elements(int token, CDescriptor *desc1, CDescriptor *desc2, t_compdesc_debug *deb);
t_sim_graph *build_graph(t_sim_components *sim_components, t_compdesc_debug *deb);
#else
t_sim_components *compare_elements(int token, CDescriptor *desc1, CDescriptor *desc2);
t_sim_graph *build_graph(t_sim_components *sim_components);
#endif

void clean_ion_contacts(int token, CDescriptor *desc1, CDescriptor *desc2, t_sim_components *sim_components, t_sim_graph *graph);

void sim_components_free(t_sim_components *sim_components);
void sim_graph_free(t_sim_graph *sim_graph);

#endif

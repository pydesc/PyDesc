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

#include"compdesc_priv.h"
#include"cstructures.h"
#include"01_prelim.h"

int check_sizes(CDescriptor *desc1, CDescriptor *desc2
#ifdef COMPDESC_DEBUG
        , t_compdesc_debug *deb
#endif
        )
{
#ifdef COMPDESC_DEBUG
    if(!coverage_test_prelim(desc1, desc2, (t_al_coverage **)&deb->prelim_cov)) return 0;
    deb->stage1=1;
#else
    if(!coverage_test_prelim(desc1, desc2, 0)) return 0;
#endif

    return 1;
}


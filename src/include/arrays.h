/*
 * Copyright 2010, 2017 Pawel Daniluk
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
 


#ifndef __ARRAYS_H__
#define __ARRAYS_H__

#include<stdlib.h>

int memory_available(size_t size, float frac);
size_t array_size(long int n, long int m, long int size);
void **array_alloc(long int n, long int m, long int size);
void **array_realloc(void **arr, long int n, long int m, long int size);
void array_set(void **arr, int val);
void array_cpy(void **dst, void ** src);
void array_free(void **array);

#endif

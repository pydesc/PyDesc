/*
 * Copyright 2014, 2017 Pawel Daniluk
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
 


#ifndef __MAPS_H__
#define __MAPS_H__

#include<stdlib.h>

typedef struct {
    int len1;
    int len2;
    int *map12;
    int *map21;
} t_map;

t_map *map_alloc(int len1, int len2);
void map_free(t_map *map);

void map_set(t_map *map, int val1, int val2);
void map_unset(t_map *map, int key, int dir);

int map_is_included(t_map *map1, t_map *map2);
void map_cpy(t_map *desc, t_map *src);

int map_get(t_map *map, int key, int dir);

void map_iter_dir(t_map *map, void (*f)(int, int), int dir);

inline int __attribute__((always_inline)) map_map(t_map *map, int key)
{
    return map_get(map, key, 1);
}

inline int __attribute__((always_inline)) map_unmap(t_map *map, int key)
{
    return map_get(map, key, 2);
}

inline void __attribute__((always_inline)) map_iter(t_map *map, void (*f)(int, int))
{
    map_iter_dir(map, f, 1);
}

inline void __attribute__((always_inline)) map_iter_rev(t_map *map, void (*f)(int, int))
{
    map_iter_dir(map, f, 2);
}


#endif

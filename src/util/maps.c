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
 
#include<stdlib.h>
#include<stdio.h>
#include<string.h>
#include<unistd.h>



#include"simple_macros.h"
#include"maps.h"

t_map *map_alloc(int len1, int len2)
{
    t_map *res=malloc(sizeof(t_map));

    res->len1=len1;
    res->len2=len2;

    res->map12=calloc(len1, sizeof(int));
    res->map21=calloc(len2, sizeof(int));

    for(int i=0; i<len1; res->map12[i++]=-1);
    for(int i=0; i<len2; res->map21[i++]=-1);

    return res;
}


void map_free(t_map *map)
{
    free(map->map12);
    free(map->map21);

    free(map);
}

void map_set(t_map *map, int val1, int val2)
{
    int oldmap2=map->map12[val1];
    int oldmap1=map->map21[val2];

    if(oldmap2>=0 && oldmap2 != val2) {
        map->map12[val1]=-1;
        map->map21[oldmap2]=-1;
    }

    if(oldmap1>=0 && oldmap1 != val1) {
        map->map21[val2]=-1;
        map->map12[oldmap1]=-1;
    }

    map->map12[val1]=val2;
    map->map21[val2]=val1;
}

void map_unset(t_map *map, int key, int dir)
{
    if(dir==1) {
        int oldmap=map->map12[key];
        if(oldmap>=0) {
            map->map12[key]=-1;
            map->map21[oldmap]=-1;
        }
    } else {
        int oldmap=map->map21[key];
        if(oldmap>=0) {
            map->map21[key]=-1;
            map->map12[oldmap]=-1;
        }
    }
}

int map_get(t_map *map, int key, int dir)
{
    if(key<0) {
        P_INT(key) P_NL
        abort();
    }
    if(dir==1) {
        if(key>=map->len1) {
            P_INT(key) P_INT(map->len1) P_NL;
            abort();
        }
        return map->map12[key];

    } else {
        if(key>=map->len2) {
            P_INT(key) P_INT(map->len2) P_NL;
            abort();
        }
        return map->map21[key];
    }
}

void map_iter_dir(t_map *map, void (*f)(int, int), int dir)
{
    if(dir==1) {
        for(int i=0; i<map->len1; i++) {
            if(map->map12[i] != -1) f(i, map->map12[i]);
        }
    } else {
        for(int i=0; i<map->len2; i++) {
            if(map->map21[i] != -1) f(map->map21[i], i);
        }
    }
}

void map_cpy(t_map *desc, t_map *src)
{
    memcpy(desc->map12, src->map12, src->len1 * sizeof(int));
    memcpy(desc->map21, src->map21, src->len2 * sizeof(int));
}

int map_is_included(t_map *map1, t_map *map2)
{
    for(int i=0; i<map1->len1; i++) {
        int mi=map_map(map1, i);

        if(mi!=-1 && mi!=map_map(map2, i)) return 0;
    }

    return 1;
}

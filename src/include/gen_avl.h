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
 
#ifndef __GEN_AVL_H__
#define __GEN_AVL_H__

//#define NULL 0

#define avl_empty (struct avl_node *) NULL

#ifdef __cplusplus
extern "C" {
#endif

struct avl_node {
	char * key;
	char * data;
	
	short key_size;
	short data_size;
	
	short avl_height;
	struct avl_node * avl_left;
	struct avl_node * avl_right;
};

void avl_find (char * key, short key_size, struct avl_node ** ret_node, struct avl_node ** ptree);
void avl_insert (struct avl_node * new_node, struct avl_node ** ptree);
void avl_remove (struct avl_node * node_to_delete, struct avl_node ** ptree);
void avl_destroy_tree (struct avl_node ** ptree, void destroy_node(struct avl_node *p));

void avl_walk_tree (struct avl_node **ptree, void f(struct avl_node *p));

#ifdef __cplusplus
};
#endif

#endif

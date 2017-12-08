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
 
#include<unistd.h>
#include "gen_avl.h"

/*
 * A struct avl_node has the following fields:
 *   avl_left     left son of a tree node
 *   avl_right    right son of a tree node
 *   avl_height   1+max(heightof(left),heightof(right))
 * The empty tree is represented as NULL.
 */

/* Since the trees are balanced, their height will never be large. */
#define avl_maxheight	100
#define heightof(tree)	((tree) == avl_empty ? 0 : (tree)->avl_height)

/*
 * Consistency and balancing rules:
 * 1. tree->avl_height == 1+max(heightof(tree->avl_left),heightof(tree->avl_right))
 * 2. abs( heightof(tree->avl_left) - heightof(tree->avl_right) ) <= 1
 * 3. foreach node in tree->avl_left: node->avl_key <= tree->avl_key,
 *    foreach node in tree->avl_right: node->avl_key >= tree->avl_key.
 */

/*
 * Rebalance a tree.
 * After inserting or deleting a node of a tree we have a sequence of subtrees
 * nodes[0]..nodes[k-1] such that
 * nodes[0] is the root and nodes[i+1] = nodes[i]->{vm_avl_left|vm_avl_right}.
 */
static void avl_rebalance (struct avl_node *** nodeplaces_ptr, int count)
{
	for ( ; count > 0 ; count--) {
		struct avl_node ** nodeplace = *--nodeplaces_ptr;
		struct avl_node * node = *nodeplace;
		struct avl_node * nodeleft = node->avl_left;
		struct avl_node * noderight = node->avl_right;
		int heightleft = heightof(nodeleft);
		int heightright = heightof(noderight);
		if (heightright + 1 < heightleft) {
			/*                                                      */
			/*                            *                         */
			/*                          /   \                       */
			/*                       n+2      n                     */
			/*                                                      */
			struct avl_node * nodeleftleft = nodeleft->avl_left;
			struct avl_node * nodeleftright = nodeleft->avl_right;
			int heightleftright = heightof(nodeleftright);
			if (heightof(nodeleftleft) >= heightleftright) {
				/*                                                        */
				/*                *                    n+2|n+3            */
				/*              /   \                  /    \             */
				/*           n+2      n      -->      /   n+1|n+2         */
				/*           / \                      |    /    \         */
				/*         n+1 n|n+1                 n+1  n|n+1  n        */
				/*                                                        */
				node->avl_left = nodeleftright; 
				nodeleft->avl_right = node;
				nodeleft->avl_height = 1 + (node->avl_height = 1 + heightleftright);
				*nodeplace = nodeleft;
			} else {
				/*                                                        */
				/*                *                     n+2               */
				/*              /   \                 /     \             */
				/*           n+2      n      -->    n+1     n+1           */
				/*           / \                    / \     / \           */
				/*          n  n+1                 n   L   R   n          */
				/*             / \                                        */
				/*            L   R                                       */
				/*                                                        */
				nodeleft->avl_right = nodeleftright->avl_left;
				node->avl_left = nodeleftright->avl_right;
				nodeleftright->avl_left = nodeleft;
				nodeleftright->avl_right = node;
				nodeleft->avl_height = node->avl_height = heightleftright;
				nodeleftright->avl_height = heightleft;
				*nodeplace = nodeleftright;
			}
		}
		else if (heightleft + 1 < heightright) {
			/* similar to the above, just interchange 'left' <--> 'right' */
			struct avl_node * noderightright = noderight->avl_right;
			struct avl_node * noderightleft = noderight->avl_left;
			int heightrightleft = heightof(noderightleft);
			if (heightof(noderightright) >= heightrightleft) {
				node->avl_right = noderightleft; 
				noderight->avl_left = node;
				noderight->avl_height = 1 + (node->avl_height = 1 + heightrightleft);
				*nodeplace = noderight;
			} else {
				noderight->avl_left = noderightleft->avl_right;
				node->avl_right = noderightleft->avl_left;
				noderightleft->avl_right = noderight;
				noderightleft->avl_left = node;
				noderight->avl_height = node->avl_height = heightrightleft;
				noderightleft->avl_height = heightright;
				*nodeplace = noderightleft;
			}
		}
		else {
			int height = (heightleft<heightright ? heightright : heightleft) + 1;
			if (height == node->avl_height)
				break;
			node->avl_height = height;
		}
	}
}

static int inline key_less_than(char * key1, short key1_size, char * key2, short key2_size)
{
	int i;
	for (i = 0; (i < key1_size) && (i < key2_size); i++) 
		if (key1[i] > key2[i]) return (0);
		else if(key1[i] < key2[i]) return(1);
    
	if (key1_size > key2_size) return (0);

	return (1);
}

static int inline key_equals(char * key1, short key1_size, char * key2, short key2_size)
{
	int i;

	if (key1_size != key2_size) return (0);

	for (i = 0; i < key1_size; i++) 
		if (key1[i] != key2[i]) return (0);
    
	return (1);
}

/* Find a node in the tree */ 
void avl_find (char * key, short key_size, struct avl_node ** ret_node, struct avl_node ** ptree)
{
	struct avl_node ** nodeplace = ptree;

	for (;;) {
		struct avl_node * node = *nodeplace;
		if (node == avl_empty) {
			* ret_node = avl_empty;
			return;
		}
		if (key_equals (key, key_size, node->key, node->key_size))
			break;
		if (key_less_than (key, key_size, node->key, node->key_size))
			nodeplace = &node->avl_left;
		else
			nodeplace = &node->avl_right;
	}
	* ret_node = * nodeplace;
}

/* Insert a node into a tree. */
void avl_insert (struct avl_node * new_node, struct avl_node ** ptree)
{
	char * key = new_node->key;
	short key_size = new_node->key_size;
	struct avl_node ** nodeplace = ptree;
	struct avl_node ** stack[avl_maxheight];
	int stack_count = 0;
	struct avl_node *** stack_ptr = &stack[0]; /* = &stack[stackcount] */
	for (;;) {
		struct avl_node * node = *nodeplace;
		if (node == avl_empty)
			break;
		*stack_ptr++ = nodeplace; stack_count++;
		if (key_less_than (key, key_size, node->key, node->key_size))
			nodeplace = &node->avl_left;
		else
			nodeplace = &node->avl_right;
	}
	new_node->avl_left = avl_empty;
	new_node->avl_right = avl_empty;
	new_node->avl_height = 1;
	*nodeplace = new_node;
	avl_rebalance(stack_ptr,stack_count);
}

/* Removes a node out of a tree. */
void avl_remove (struct avl_node * node_to_delete, struct avl_node ** ptree)
{
	char * key = node_to_delete->key;
	short key_size = node_to_delete->key_size;
	struct avl_node ** nodeplace = ptree;
	struct avl_node ** stack[avl_maxheight];
	int stack_count = 0;
	struct avl_node *** stack_ptr = &stack[0]; /* = &stack[stackcount] */
	struct avl_node ** nodeplace_to_delete;
	for (;;) {
		struct avl_node * node = *nodeplace;
		*stack_ptr++ = nodeplace; stack_count++;
		if (key_equals (key, key_size, node->key, node->key_size))
			break;
		if (key_less_than (key, key_size, node->key, node->key_size))
			nodeplace = &node->avl_left;
		else
			nodeplace = &node->avl_right;
	}
	nodeplace_to_delete = nodeplace;
	/* Have to remove node_to_delete = *nodeplace_to_delete. */
	if (node_to_delete->avl_left == avl_empty) {
		*nodeplace_to_delete = node_to_delete->avl_right;
		stack_ptr--; stack_count--;
	} else {
		struct avl_node *** stack_ptr_to_delete = stack_ptr;
		struct avl_node ** nodeplace = &node_to_delete->avl_left;
		struct avl_node * node;
		for (;;) {
			node = *nodeplace;
			if (node->avl_right == avl_empty)
				break;
			*stack_ptr++ = nodeplace; stack_count++;
			nodeplace = &node->avl_right;
		}
		*nodeplace = node->avl_left;
		/* node replaces node_to_delete */
		node->avl_left = node_to_delete->avl_left;
		node->avl_right = node_to_delete->avl_right;
		node->avl_height = node_to_delete->avl_height;
		*nodeplace_to_delete = node; /* replace node_to_delete */
		*stack_ptr_to_delete = &node->avl_left; /* replace &node_to_delete->vm_avl_left */
	}
	avl_rebalance(stack_ptr, stack_count);
}

void avl_destroy_tree (struct avl_node ** ptree, void destroy_node(struct avl_node *p))
{

	while(*ptree != avl_empty) {
		struct avl_node *node=*ptree;

		avl_remove(node, ptree);

		destroy_node(node);
	}
}

void avl_walk_tree (struct avl_node **ptree, void f(struct avl_node *p))
{
	if(*ptree == avl_empty) return;


	f(*ptree);

	avl_walk_tree(&(*ptree)->avl_left, f);
	avl_walk_tree(&(*ptree)->avl_right, f);
}

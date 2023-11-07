/*    
    Copyright 2013-2023 Onera.

    This file is part of Cassiopee.

    Cassiopee is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Cassiopee is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "proto.h"

tree::tree(E_Int n, E_Int s)
:
	enabled(n),
	level(n),
	children((1+s)*n),
	indir(n),
	parent(n),
	last(0),
	size(n),
  stride(s),
  nleaves(n),
  l2g(n)
{
	for (E_Int i = 0; i < n; i++) {
		enabled[i] = 1;
		level[i] = 0;
		indir[i] = -1;
		parent[i] = -1;
	}

	for (E_Int i = 0; i < (1+stride)*n; i++)
		children[i] = -1;
  
  for (E_Int i = 0; i < n; i++)
    l2g[i] = i;
}

tree::tree(E_Int s)
:
	enabled(),
	level(),
	children(),
	indir(),
	parent(),
	last(0),
	size(0),
  stride(s),
  nleaves(0),
  l2g()
{}

void tree::setSizeAndStride(E_Int n, E_Int s)
{
  enabled.resize(n);
  level.resize(n);
  children.resize((1+s)*n);
  indir.resize(n);
  parent.resize(n);
  last = 0;
  size = n;
  stride = s;
  nleaves = n;
  l2g.resize(n);

  for (E_Int i = 0; i < n; i++) {
		enabled[i] = 1;
		level[i] = 0;
		indir[i] = -1;
		parent[i] = -1;
	}

	for (E_Int i = 0; i < (1+stride)*n; i++)
		children[i] = -1;
  
  for (E_Int i = 0; i < n; i++)
    l2g[i] = i;
}

void tree_insert_children(tree *T, E_Int id, E_Int start, E_Int n)
{
	E_Int i, *pt;

	T->enabled[id] = 0;
	
	pt = &T->enabled[start];
	for (i = 0; i < n; i++)
		pt[i] = 1;
  
  T->nleaves += n - 1;
	
	pt = &T->level[start];
	E_Int lvl = T->level[id]+1;
	for (i = 0; i < n; i++)
		pt[i] = lvl;
	
	pt = &T->children[T->last];
	pt[0] = n;
	pt++;
	for (i = 0; i < n; i++)
		pt[i] = start + i;
	T->indir[id] = T->last;

	pt = &T->parent[start];
	for (i = 0; i < n; i++)
		pt[i] = id;

	T->last += n+1;
}

E_Float tree_memsize(tree *T)
{
  size_t memsize = 0;
  memsize += T->enabled.size();
  memsize += T->level.size();
  memsize += T->parent.size();
  memsize += T->indir.size();
  memsize += T->children.size();

  return memsize*sizeof(E_Int)/1000000.;
}

void tree_resize(tree *T, E_Int increment, E_Int stride)
{
	E_Int new_size = T->size + increment;

	T->enabled.resize(new_size);
	T->level.resize(new_size);
	T->indir.resize(new_size);
	T->parent.resize(new_size);

	E_Int new_size_c = new_size * (stride + 1);
	
	T->children.resize(new_size_c); 

	// init
	for (E_Int i = T->size; i < new_size; i++) {
		T->enabled[i] = 1;
		T->level[i] = 0;
		T->indir[i] = -1;
		T->parent[i] = -1;
	}

	for (E_Int i = T->last; i < new_size_c; i++) {
		T->children[i] = -1;
	}

	T->size = new_size;
}

E_Int tree_get_nchildren(tree *T, E_Int id)
{
	if (T->indir[id] == -1) return 0;
	return T->children[T->indir[id]];
}

E_Int *tree_get_children(tree *T, E_Int id)
{
	if (T->indir[id] == -1) return NULL;
	return &T->children[T->indir[id]+1];
}

void tree_print(tree *T)
{
  printf("size: %d\n", T->size);
  printf("last: %d\n", T->size);
  assert((E_Int)T->enabled.size() == T->size);
  assert((E_Int)T->level.size() == T->size);
  assert((E_Int)T->parent.size() == T->size);
  assert((E_Int)T->indir.size() == T->size);
  printf("enabled:\n");
  for (E_Int i = 0; i < T->size; i++) {printf("%d ", T->enabled[i]);} puts("");
  printf("level:\n");
  for (E_Int i = 0; i < T->size; i++) {printf("%d ", T->level[i]);} puts("");
  printf("parent:\n");
  for (E_Int i = 0; i < T->size; i++) {printf("%d ", T->parent[i]);} puts("");
	for (E_Int i = 0; i < T->size; i++) {
		E_Int where = T->indir[i];
		if (where == -1) continue;
		printf("%d: ", i);
		E_Int nchild = T->children[where];
		E_Int *children = &T->children[where+1];
		for (E_Int j = 0; j < nchild; j++)
			printf("%d ", children[j]);
		puts("");
	}
}

void tree_get_face_leaves(tree *ft, E_Int face, std::vector<E_Int> &leaves)
{
  E_Int nchildren = tree_get_nchildren(ft, face);

  if (nchildren == 0) {
    leaves.push_back(face);
    return;
  }

  E_Int *children = tree_get_children(ft, face);
  for (E_Int i = 0; i < nchildren; i++)
    tree_get_face_leaves(ft, children[i], leaves);
}

#include "proto.h"

tree::tree(E_Int n, E_Int stride)
:
	enabled(n),
	level(n),
	children((1+stride)*n),
	indir(n),
	parent(n),
	last(0),
	size(n)
{
	for (E_Int i = 0; i < n; i++) {
		enabled[i] = 1;
		level[i] = 0;
		indir[i] = -1;
		parent[i] = -1;
	}

	for (E_Int i = 0; i < (1+stride)*n; i++)
		children[i] = -1;
}

void tree_insert_children(tree *T, E_Int id, E_Int start, E_Int n)
{
	E_Int i, *pt;

	T->enabled[id] = 0;
	
	pt = &T->enabled[start];
	for (i = 0; i < n; i++)
		pt[i] = 1;
	
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
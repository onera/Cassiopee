#include "Proto.h"

void tree_init(Element **tree, E_Int nelem)
{
  for (E_Int i = 0; i < nelem; i++) {
    tree[i] = (Element *)XMALLOC(sizeof(Element));
    Element *Elem = tree[i];
    Elem->children = NULL;
    Elem->nchildren = 0;
    Elem->parent = i;
    Elem->position = 0;
    Elem->type = -1;
    Elem->level = 0;
  }
}

void tree_print(Element **tree, E_Int nelem)
{
  puts("");
  for (E_Int i = 0; i < nelem; i++) {
    E_Int nchildren = tree[i]->nchildren;
    E_Int *children = tree[i]->children;
    printf("%d -> ", i);
    for (E_Int j = 0; j < nchildren; j++) {
      printf("%d ", children[j]);
    }
    puts("");
  }
  puts("");
}
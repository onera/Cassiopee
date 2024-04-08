#include "proto.h"

event::event(point *P)
: p(P), s(NULL), left(NULL), right(NULL)
{}

void event_show(event *E)
{
  printf("<%lu,", E->p->id);
  if (E->s) printf("%lu>\n", E->s->id);
  else printf("NULL>\n");
}

void event_print(event *E)
{
  if (E == NULL) return;

  event_print(E->left);
  event_show(E);
  event_print(E->right);
}

/* INSERTION */

event *event_insert(event *&root, point *p)
{
    if (root == NULL) {
        root = new event(p);
        return root;
    }

    int cmp = cmp_points(root->p, p);

    if (cmp == 0) return root;

    if (cmp == -1) return event_insert(root->right, p);

    return event_insert(root->left, p);
}

/* SEARCH */

event *event_lookup(event *root, point *p)
{
    if (root == NULL)
        return root;

    int cmp = cmp_points(root->p, p);

    if (cmp == 0) return root;
    
    if (cmp == -1) return event_lookup(root->right, p);

    return event_lookup(root->left, p);
}

event *event_locate(event *root, E_Float x, E_Float y)
{
  if (root == NULL)
    return root;

  int cmp = cmp_xyz(root->p->x, root->p->y, x, y);

  if (cmp == 0)
    return root;
  else if (cmp == -1)
    return event_locate(root->right, x, y);
  else
    return event_locate(root->left, x, y);
}

event *event_min(event *root)
{
    event *current = root;

    while (current->left != NULL)
        current = current->left;

    return current;
}

/* DELETION */

event *event_delete(event *root, point *p)
{
    if (root == NULL) return root;

    if (root->p == p) {
      // Node item is reached.

      // If one of the children is empty:
      if (root->left == NULL) {
        event *tmp = root->right;
        delete root;
        return tmp;
      } else if (root->right == NULL) {
        event *tmp = root->left;
        delete root;
        return tmp;
      }

      // If both of the children exist:
      else {
        event *succ_parent = root;

        // Find successor: leftmost node in right subtree.
        event *succ = root->right;
        while (succ->left != NULL) {
          succ_parent = succ;
          succ = succ->left;
        }

        if (succ_parent != root)
          succ_parent->left = succ->right;
        else
          succ_parent->right = succ->right;

        // Copy successor data to root.
        root->p = succ->p;
        root->s = succ->s;

        // Delete successor and return root.
        delete succ;
        return root;
      }
    }

    int cmp = cmp_points(root->p, p);

    assert(cmp);

    if (cmp == -1) {
      root->right = event_delete(root->right, p);
      return root;
    } else if (cmp == 1) {
      root->left = event_delete(root->left, p);
      return root;
    }
}

#include "proto.h"

status::status(segment *S)
: s(S), p(NULL), c(NULL), left(NULL), right(NULL)
{}

void status_show(status *S)
{
  printf("<S%lu,",S->s->id);
  if (S->p) printf("P%lu,", S->p->id);
  else printf("NULL,");
  if (S->c) printf("S%lu>\n", S->c->id);
  else printf("NULL>\n");
}

void status_print(status *root)
{
  if (root == NULL) return;

  status_print(root->left);
  status_show(root);
  status_print(root->right);
}

/* SEGMENT ORDER */

static
int _status_cmp_segments_sweep(segment *s1, segment *s2, E_Float xs,
  E_Float ys)
{
    double px1 = s1->p->x;
    double py1 = s1->p->y;
    double qx1 = s1->q->x;
    double qy1 = s1->q->y;
    double px2 = s2->p->x;
    double py2 = s2->p->y;
    double qx2 = s2->q->x;
    double qy2 = s2->q->y;
    double dx1 = qx1 - px1;
    double dy1 = qy1 - py1;
    double dx2 = qx2 - px2;
    double dy2 = qy2 - py2;

    // Test if the segments have the same underlying line.
    double T1 = dy1*dx2 - dy2*dx1;
    int sign1 = get_sign(T1);
    if (sign1 == 0) {
        // Parallel lines.
        // Is p2 collinear to s1?
        int sign2 = get_orient(s1->p, s1->q, s2->p);
        if (sign2 == 0) {
            // Yes.
            // q2 should be collinear too.
            assert(get_orient(s1->p, s1->q, s2->q) == 0);

            return 0;
        }
    }

    // Lines are different.
    // Most one of the lines is vertical.
    // Deal first with the cases where one of the lines is vertical.

    // s1 is vertical, s2 is not vertical.
    // Return -1 is l x s2 is above psweep. 
    if (f_eq(dx1, 0.0)) {
        int i = get_sign(py2*dx2 - px2*dy2 + dy2*xs - ys*dx2);
        return (i <= 0) ? 1 : -1;
    }

    // s2 is vertical, s1 is not vertical.
    // Return -1 if l x s1 is below or equal to psweep.
    if (f_eq(dx2, 0.0)) {
        int i = get_sign(py1*dx1 - px1*dy1 + dy1*xs - ys*dx1);
        return (i <= 0) ? -1 : 1;
    }

    // Neither s1 nor s2 is vertical. Compare l x s1 and l x s2
    double T2 = dx2*(py1*dx1 + dy1*(xs - px1)) - dx1*(py2*dx2 + dy2*(xs - px2));
    int sign2 = get_sign(T2);
    if (sign2 != 0) return sign2;

    // s1 and s2 intersect the sweep line in the same point.
    // We compare slopes:
    // s1 had larger slope than s2 iff T1 > 0

    double T3 = (py1 - ys)*dx1 - dy1*(xs - px1);

    // The common intersection is above psweep iff T3 > 0. In this case
    // we return -sign(T1) and sign(T1) otherwise.

    int sign3 = get_sign(T3);
    return (sign3 <= 0) ? sign1 : -sign1;
}

static
int _status_cmp_segments(segment *s0, segment *s1, E_Float xs, E_Float ys)
{
    int cmp = _status_cmp_segments_sweep(s0, s1, xs, ys);

    if (cmp == 0) {
        cmp = cmp_segments_lexico(s0, s1);
    }

    assert(cmp != 0);

    return cmp;
}

/* INSERTION */

status *status_insert(status *&root, segment *s, E_Float xs, E_Float ys)
{
    if (root == NULL) {
        root = new status(s);
        return root;
    }

    if (root->s == s) {
        return root;
    }

    int cmp = _status_cmp_segments(root->s, s, xs, ys);

    if (cmp == -1)
      return status_insert(root->right, s, xs, ys);

    return status_insert(root->left, s, xs, ys);
}

/* SEARCH */

status *status_lookup(status *root, segment *s, E_Float xs, E_Float ys)
{
    if (root == NULL || root->s == s)
        return root;

    int cmp = _status_cmp_segments(root->s, s, xs, ys);

    if (cmp == -1)
      return status_lookup(root->right, s, xs, ys);

    return status_lookup(root->left, s, xs, ys);
}

status *status_locate(status *root, segment *s, E_Float xs, E_Float ys)
{
    if (root == NULL) return root;

    int cmp = _status_cmp_segments_sweep(root->s, s, xs, ys);

    if (cmp == 0)
        return root;
    else if (cmp == -1)
        return status_locate(root->right, s, xs, ys);
    else
        return status_locate(root->left, s, xs, ys);
}

static
void _status_succ(status *root, segment *s, E_Float xs, E_Float ys,
  status *&succ)
{
    if (root == NULL) return;

    if (root->s == s) {
        // Successor is leftmost child in right subtree
        if (root->right != NULL) {
            status *tmp = root->right;
            while (tmp->left)
                tmp = tmp->left;
            succ = tmp;
        }

        return;
    }

    int cmp = _status_cmp_segments(root->s, s, xs, ys);

    assert(cmp != 0);

    if (cmp == -1) {
        _status_succ(root->right, s, xs, ys, succ);
    } else {
        succ = root;
        _status_succ(root->left, s, xs, ys, succ);
    }
}

status *status_succ(status *root, segment *s, E_Float xs, E_Float ys)
{
    status *succ = NULL;
    _status_succ(root, s, xs, ys, succ);
    return succ;
}

static
void _status_pred(status *root, segment *s, E_Float xs, E_Float ys,
  status *&pred)
{
    if (root == NULL) return;

    if (root->s == s) {
        // Predecessor is rightmost child in left subtree
        if (root->left != NULL) {
            status *tmp = root->left;
            while (tmp->right)
                tmp = tmp->right;
            pred = tmp;
        }

        return;
    }

    int cmp = _status_cmp_segments(root->s, s, xs, ys);

    assert(cmp != 0);

    if (cmp == -1) {
        pred = root;
        _status_pred(root->right, s, xs, ys, pred);
    } else {
        _status_pred(root->left, s, xs, ys, pred);
    }
}

status *status_pred(status *root, segment *s, E_Float xs, E_Float ys)
{
    status *pred = NULL;
    _status_pred(root, s, xs, ys, pred);
    return pred;
}

/* DELETION */

status *status_delete(status *root, segment *s, E_Float xs, E_Float ys)
{
    if (root == NULL) return root;

    if (root->s == s) {
        // Node item is reached.

        // If one of the children is empty:
        if (root->left == NULL) {
            status *tmp = root->right;
            delete root;
            return tmp;
        } else if (root->right == NULL) {
            status *tmp = root->left;
            delete root;
            return tmp;
        }

        // If both of the children exist:
        else {
            status *succ_parent = root;

            // Find successor: leftmost node in right subtree.
            status *succ = root->right;
            while (succ->left != NULL) {
                succ_parent = succ;
                succ = succ->left;
            }

            if (succ_parent != root)
                succ_parent->left = succ->right;
            else
                succ_parent->right = succ->right;

            // Copy successor data to root.
            root->s = succ->s;
            root->p = succ->p;
            root->c = succ->c;

            // Delete successor and return root.
            delete succ;
            return root;
        }
    }

    int cmp = _status_cmp_segments(root->s, s, xs, ys);

    assert(cmp);

    if (cmp == -1) {
        root->right = status_delete(root->right, s, xs, ys);
        return root;
    } else {
        root->left = status_delete(root->left, s, xs, ys);
        return root;
    }
}

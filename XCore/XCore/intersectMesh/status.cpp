#include "proto.h"
#include <cstddef>
#include <cassert>
#include <cstdio>

snode::snode(segment *S)
: s(S), inf(NULL), left(NULL), right(NULL)
{}

status::status()
: root(NULL), xs(0), ys(0)
{}

snode *status::insert(segment *s)
{
    return _insert(root, s);
}

snode *status::_insert(snode *&root, segment *s)
{
    if (root == NULL) {
        root = new snode(s);
        return root;
    }

    int cmp = _segment_cmp(root->s, s);

    if (cmp == 0)
        return root;
    else if (cmp < 0)
        return _insert(root->right, s);
    else
        return _insert(root->left, s);
}

snode *status::locate(segment *s)
{
    return _locate(root, s);
}

snode *status::_locate(snode *root, segment *s)
{
    if (root == NULL)
        return NULL;

    int cmp = _segment_cmp_sweep(root->s, s);

    if (cmp == 0)
        return root;
    else if (cmp < 0)
        return _locate(root->right, s);
    else
        return _locate(root->left, s);
}

snode *status::lookup(segment *s)
{
    return _lookup(root, s);
}

snode *status::_lookup(snode *root, segment *s)
{
    if (root == NULL) {
        root = new snode(s);
        return root;
    }

    int cmp = _segment_cmp(root->s, s);

    if (cmp == 0)
        return root;
    else if (cmp < 0)
        return _insert(root->right, s);
    else
        return _insert(root->left, s);
}

void status::erase(segment *s)
{
    root = _erase(root, s);
}

snode *status::_erase(snode *root, segment *s)
{
    if (!root)
        return root;

    int cmp = _segment_cmp(root->s, s);

    if (cmp < 0) {
        root->right = _erase(root->right, s);
        return root;
    } else if (cmp > 0) {
        root->left = _erase(root->left, s);
        return root;
    }

    assert(root->s == s);

    if (root->left == NULL) {
        snode *tmp = root->right;
        delete root;
        return tmp;
    } else if (root->right == NULL) {
        snode *tmp = root->left;
        delete root;
        return tmp;
    } else {
        snode *succ_parent = root;

        snode *succ = root->right;
        while (succ->left) {
            succ_parent = succ;
            succ = succ->left;
        }

        if (succ_parent != root)
            succ_parent->left = succ->right;
        else
            succ_parent->right = succ->right;

        root->s = succ->s;
        root->inf = succ->inf;

        delete succ;
        return root;
    }
}

snode *status::pred(snode *sit)
{
    return pred(sit->s);
}

snode *status::pred(segment *s)
{
    snode *pre = NULL;
    _pred(root, s, pre);
    return pre;
}

void status::_pred(snode *root, segment *s, snode *&pre)
{
    if (!root)
        return;

    int cmp = _segment_cmp(root->s, s);

    if (cmp == 0) {
        assert(root->s == s);

        if (root->left != NULL) {
            snode *tmp = root->left;
            while (tmp->right)
                tmp = tmp->right;
            pre = tmp;
        }

        return;
    }

    if (cmp > 0) {
        _pred(root->left, s, pre);
    } else {
        pre = root;
        _pred(root->right, s, pre);
    }
}

snode *status::succ(snode *sit)
{
    return succ(sit->s);
}

snode *status::succ(segment *s)
{
    snode *suc = NULL;
    _succ(root, s, suc);
    return suc;
}

void status::_succ(snode *root, segment *s, snode *&suc)
{
    if (!root)
        return;

    int cmp = _segment_cmp(root->s, s);

    if (cmp == 0) {
        assert(root->s == s);

        if (root->right != NULL) {
            snode *tmp = root->right;
            while (tmp->left)
                tmp = tmp->left;
            suc = tmp;
        }

        return;
    }

    if (cmp > 0) {
        suc = root;
        _succ(root->left, s, suc);
    } else {
        _succ(root->right, s, suc);
    }
}


int status::_segment_cmp(segment *s0, segment *s1)
{
    if (s0 == s1)
        return 0;

    int cmp = _segment_cmp_sweep(s0, s1);

    if (cmp == 0) {
        //cmp = sign(s0->id - s1->id);
        cmp = segment_cmp_lexico(s0, s1);
    }

    return cmp;
}

int status::_segment_cmp_sweep(segment *s0, segment *s1)
{
    double px0 = s0->p->x;
    double py0 = s0->p->y;
    double qx0 = s0->q->x;
    double qy0 = s0->q->y;
    double px1 = s1->p->x;
    double py1 = s1->p->y;
    double qx1 = s1->q->x;
    double qy1 = s1->q->y;
    double dx0 = qx0 - px0;
    double dy0 = qy0 - py0;
    double dx1 = qx1 - px1;
    double dy1 = qy1 - py1; 

    double T1 = dy0*dx1 - dy1*dx0;
    int sign1 = sign(T1);
    if (sign1 == 0) {
        int sign2 = vertex_orient(s0->p, s0->q, s1->p);
        if (sign2 == 0) {
            assert(vertex_orient(s0->p, s0->q, s1->q) == 0);
            return 0;
        }
    }

    if (sign(dx0) == 0) {
        double T2 = py1*dx1 - px1*dy1 + dy1*xs - ys*dx1;
        int sign2 = sign(T2);
        return (sign2 <= 0) ? 1 : -1;
    }

    if (sign(dx1) == 0) {
        double T2 = py0*dx0 - px0*dy0 + dy0*xs - ys*dx0;
        int sign2 = sign(T2);
        return (sign2 <= 0) ? -1 : 1;
    }

    double T2 = dx1*(py0*dx0 + dy0*(xs - px0)) -
                dx0*(py1*dx1 + dy1*(xs - px1));
    int sign2 = sign(T2);
    if (sign2 != 0)
        return sign2;

    double T3 = (py0*dx0 - px0*dy0) + (dy0*xs - ys*dx0);
    int sign3 = sign(T3);
    return (sign3 <= 0) ? sign1 : -sign1;
}

void status::update_sweep_position(double x, double y)
{
    xs = x;
    ys = y;
}

void status::print()
{
    snode_print(root);
}

void snode_print(snode *root)
{
    if (!root)
        return;

    snode_print(root->left);
    root->print();
    snode_print(root->right);
}

void snode::print()
{
    printf("<S%d>\n", s->id);
}

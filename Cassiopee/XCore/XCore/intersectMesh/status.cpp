/*    
    Copyright 2013-2024 Onera.

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
#include <cstddef>
#include <cassert>

#include "status.h"
#include "primitives.h"

Status::Status()
: root(NULL), rx(0.0), ry(0.0)
{}

static
Int Segment_cmp(Segment *s0, Segment *s1, Float rx, Float ry)
{
    if (s0 == s1) return 0;

    Int cmp = compare(*s0, *s1, rx, ry);

    if (cmp == 0) {
        cmp = cmp_mySeg(*s0, *s1);
    }

    return cmp;
}

Snode *Status::insert(Segment *key, void *inf)
{
    if (key == NULL) return NULL;
    return insert_(root, key, inf);
}

Snode *Status::insert_(Snode *&root, Segment *key, void *inf)
{
    if (root == NULL) {
        root = new Snode(key, inf);
        return root;
    }

    Int cmp = Segment_cmp(root->key, key, rx, ry);

    if (cmp == 0) {
        root->inf = inf;
        return root;
    } else if (cmp < 0) {
        return insert_(root->right, key, inf);
    } else {
        return insert_(root->left, key, inf);
    }
}

Snode *Status::locate(Segment *seg)
{
    if (seg == NULL) return NULL;
    return locate_(root, seg);
}

Snode *Status::locate_(Snode *root, Segment *seg)
{
    if (root == NULL) return NULL;

    Int cmp = compare(*root->key, *seg, rx, ry);

    if (cmp == 0) return root;
    else if (cmp < 0) return locate_(root->right, seg);
    else return locate_(root->left, seg);
}

Snode *Status::lookup(Segment *seg)
{
    if (seg == NULL) return NULL;
    return lookup_(root, seg);
}

Snode *Status::lookup_(Snode *root, Segment *seg)
{
    if (root == NULL) return NULL;

    Int cmp = Segment_cmp(root->key, seg, rx, ry);

    if (cmp == 0) return root;
    else if (cmp < 0) return lookup_(root->right, seg);
    else return lookup_(root->left, seg);
}

Snode *Status::pred(Segment *seg)
{
    Snode *pre = NULL;
    pred_(root, seg, pre);
    return pre;
}

Snode *Status::succ(Segment *seg)
{
    Snode *suc = NULL;
    succ_(root, seg, suc);
    return suc;
}

Snode *Status::pred(Snode *sit)
{
    return pred(sit->key);
}

Snode *Status::succ(Snode *sit)
{
    return succ(sit->key);
}

void Status::pred_(Snode *root, Segment *seg, Snode *&pre)
{
    if (root == NULL) return;

    Int cmp = Segment_cmp(root->key, seg, rx, ry);
    
    if (cmp == 0) {
        assert(root->key == seg);

        if (root->left != NULL) {
            Snode *tmp = root->left;
            while (tmp->right) tmp = tmp->right;
            pre = tmp;
        }

        return;
    }

    if (cmp > 0) {
        pred_(root->left, seg, pre);
    } else {
        pre = root;
        pred_(root->right, seg, pre);
    }
}

void Status::succ_(Snode *root, Segment *seg, Snode *&suc)
{
    if (root == NULL) return;

    Int cmp = Segment_cmp(root->key, seg, rx, ry);
    
    if (cmp == 0) {
        assert(root->key == seg);

        if (root->right != NULL) {
            Snode *tmp = root->right;
            while (tmp->left) tmp = tmp->left;
            suc = tmp;
        }

        return;
    }

    if (cmp > 0) {
        suc = root;
        succ_(root->left, seg, suc);
    } else {
        succ_(root->right, seg, suc);
    }
}

void Status::erase(Segment *seg)
{
    root = erase_(root, seg);
}

Snode *Status::erase_(Snode *root, Segment *seg)
{
    if (root == NULL) return NULL;

    Int cmp = Segment_cmp(root->key, seg, rx, ry);

    if (cmp < 0) {
        root->right = erase_(root->right, seg);
        return root;
    } else if (cmp > 0) {
        root->left = erase_(root->left, seg);
        return root;
    }

    assert(root->key == seg);

    if (root->left == NULL) {
        Snode *tmp = root->right;
        delete root;
        return tmp;
    } else if (root->right == NULL) {
        Snode *tmp = root->left;
        delete root;
        return tmp;
    } else {
        Snode *succ_parent = root;

        Snode *succ = root->right;
        while (succ->left) {
            succ_parent = succ;
            succ = succ->left;
        }

        if (succ_parent != root) succ_parent->left = succ->right;
        else succ_parent->right = succ->right;

        root->key = succ->key;
        root->inf = succ->inf;

        delete succ;
        return root;
    }
}

void Status::print()
{
    return Snode::print_tree(root);
}

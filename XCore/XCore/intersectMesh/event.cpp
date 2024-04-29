#include "proto.h"
#include <cstddef>
#include <cstdio>
#include <cassert>

event::event(vertex *V)
: key(V), inf(NULL), left(NULL), right(NULL)
{}

event::event(double x, double y, int oid)
{
    key = new vertex(x, y, oid);
    inf = NULL;
    left = NULL;
    right = NULL;
}

void event::inorder(std::vector<vertex *> &V)
{
    if (left) left->inorder(V);
    V.push_back(key);
    if (right) right->inorder(V);
}

queue::queue()
: root(NULL), nelem(0)
{}

void queue::inorder(std::vector<vertex *> &V)
{
    root->inorder(V);
}

event *queue::insert(double x, double y, int oid)
{
    return _insert(root, x, y, oid);
}

event *queue::_insert(event *&root, double x, double y, int oid)
{
    if (root == NULL) {
        root = new event(x, y, oid);
        return root;
    }

    vertex *key = root->key;

    int cmp = xy_cmp(key->x, key->y, x, y);
    
    if (cmp == 0) {
        if (key->oid[0] == -1) key->oid[0] = oid;
        else key->oid[1] = oid;
        return root;
    }

    else if (cmp < 0)
        return _insert(root->right, x, y, oid);

    else
        return _insert(root->left, x, y, oid);
}

event *queue::insert(vertex *v)
{
    return _insert(root, v);
}

event *queue::_insert(event *&root, vertex *v)
{
    if (root == NULL) {
        root = new event(v);
        return root;
    }

    int cmp = vertex_cmp_xy(root->key, v);

    if (cmp == 0)
        return root;
    else if (cmp < 0)
        return _insert(root->right, v);
    else
        return _insert(root->left, v);
}

vertex *queue::locate_v(double x, double y)
{
    return locate(x, y)->key;
}

event *queue::locate(double x, double y)
{
    return _locate(root, x, y);
}

event *queue::_locate(event *root, double x, double y)
{
    if (root == NULL)
        return NULL;

    int cmp = xy_cmp(root->key->x, root->key->y, x, y);

    if (cmp == 0)
        return root;
    else if (cmp < 0)
        return _locate(root->right, x, y);
    else
        return _locate(root->left, x, y);
}

event *queue::lookup(vertex *v)
{
    return _lookup(root, v);
}

event *queue::_lookup(event *root, vertex *v)
{
    if (root == NULL)
        return NULL;

    int cmp = vertex_cmp(root->key, v);

    if (cmp == 0)
        return root;
    else if (cmp < 0)
        return _lookup(root->right, v);
    else
        return _lookup(root->left, v);
}

event *queue::min()
{
    if (root == NULL)
        return NULL;

    event *curr = root;

    while (curr->left != NULL)
        curr = curr->left;

    return curr;
}

void queue::erase(event *E)
{
    erase(E->key);
}

void queue::erase(vertex *v)
{
    root = _erase(root, v);
}

event *queue::_erase(event *root, vertex *v)
{
    if (!root)
        return root;

    int cmp = vertex_cmp(root->key, v);

    if (cmp < 0) {
        root->right = _erase(root->right, v);
        return root;
    } else if (cmp > 0) {
        root->left = _erase(root->left, v);
        return root;
    }

    assert(root->key == v);

    if (root->left == NULL) {
        event *tmp = root->right;
        delete root;
        return tmp;
    } else if (root->right == NULL) {
        event *tmp = root->left;
        delete root;
        return tmp;
    } else {
        event *succ_parent = root;

        event *succ = root->right;
        while (succ->left) {
            succ_parent = succ;
            succ = succ->left;
        }

        if (succ_parent != root)
            succ_parent->left = succ->right;
        else
            succ_parent->right = succ->left;

        root->key = succ->key;
        root->inf = succ->inf;

        delete succ;
        return root;
    }
}

int queue::empty()
{
    return root == NULL;
}

void event::print()
{
    /*
    printf("<P%d,", key->id);
    if (inf) printf("S%d>\n", inf->id);
    else printf("NULL>\n");
    */
}

void queue::print()
{
    event_print(root);
}

void event_print(event *root)
{
    if (root == NULL)
        return;

    event_print(root->left);
    root->print();
    event_print(root->right);
}

#include "snode.h"

Snode::Snode(Segment *Key, void *Inf)
: key(Key), inf(Inf), left(NULL), right(NULL)
{}

void Snode::print_tree(Snode *root)
{
    if (root == NULL) return;

    print_tree(root->left);
    root->print();
    print_tree(root->right);
}

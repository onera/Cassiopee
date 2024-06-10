#include "event.h"
#include "vertex.h"

Event::Event(E_Float x, E_Float y)
{
    key = new Vertex(x, y);
    inf = NULL;
    left = right = NULL;
}

Event::Event(E_Float x, E_Float y, E_Int oid, E_Int color)
{
    key = new Vertex(x, y, oid, color);
    inf = NULL;
    left = right = NULL;
}

void Event::inorder(std::vector<Vertex *> &V) const
{
    if (left) left->inorder(V);
    V.push_back(key);
    if (right) right->inorder(V);
}

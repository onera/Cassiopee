#include "vertex.h"

Vertex::Vertex(E_Float X, E_Float Y, E_Int Oid, E_Int color)
: x(X), y(Y), rep(NULL), id(-1), left(NULL)
{
    oid[color] = Oid;
    oid[(color+1)%2] = -1;
}

Vertex::Vertex(E_Float X, E_Float Y)
: x(X), y(Y), rep(NULL), id(-1), left(NULL)
{
    oid[0] = -1;
    oid[1] = -1;
}

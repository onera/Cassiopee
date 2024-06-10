#include <cstddef>

#include "face.h"

Face::Face()
: rep(NULL)
{
    oid[0] = oid[1] = -1;
}

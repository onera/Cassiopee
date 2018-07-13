#include <stdlib.h>
#include <float.h>
#include "gpc.h"

void share_eps_(double *b)
{
   b[0] = GPC_EPSILON;
   return;
}

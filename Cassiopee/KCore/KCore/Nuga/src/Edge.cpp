/*



--------- NUGA v1.0



*/
//Authors : Sam Landier (sam.landier@onera.fr)

#include "Nuga/include/Edge.h"
#include "Nuga/include/defs.h"
#include "Nuga/include/maths.hxx"
#include <math.h>

//=============================================================================
void
K_MESH::Edge::getBoundary(const Edge& E1, const Edge& E2, boundary_type& b)
{
  b = IDX_NONE;
  const E_Int& n10 = E1._nodes[0];
  const E_Int& n11 = E1._nodes[1];
  const E_Int& n20 = E2._nodes[0];
  const E_Int& n21 = E2._nodes[1];

  if (n11 == n20) b = n11;
  else if (n10 == n21) b = n10;
  else if (n10 == n20) b = n10;
  else if (n11 == n21) b = n11;
}

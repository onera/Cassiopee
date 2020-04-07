/*



NUGA



*/
//Author : Sâm Landier (sam.landier@onera.fr)

#ifndef NUGA_T6_HXX
#define NUGA_T6_HXX

#include "MeshElement/Polygon.h"
#include "Nuga/include/subdivision.hxx"

namespace NUGA
{

  struct T6 : public K_MESH::Polygon
  {
    using parent_type = K_MESH::Polygon;
    
    template <eSUBDIV_TYPE STYPE>
    static void split(ngon_unit& PGs, const E_Int* nodes, E_Int firstChild)
    {
      // set them in _ng.PGs
      E_Int* p = PGs.get_facets_ptr(firstChild);

      p[0] = nodes[0];  p[1] = nodes[3];  p[2] = nodes[5]; 
      p[4] = nodes[3];  p[5] = nodes[1];  p[6] = nodes[4]; 
      p[8] = nodes[5];  p[9] = nodes[4];  p[10] = nodes[2]; 
      p[12] = nodes[3]; p[13] = nodes[4]; p[14] = nodes[5]; 
    }
  };
}

#endif

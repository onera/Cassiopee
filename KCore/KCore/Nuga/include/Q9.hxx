/*
 
 
 
              NUGA 
 
 
 
 */

#ifndef NUGA_Q9_HXX
#define NUGA_Q9_HXX

#include "Nuga/include/subdivision.hxx"

namespace NUGA
{
  struct Q9
  {
    template <eSUBDIV_TYPE STYPE>
    static void split(ngon_unit& PGs, const E_Int* nodes, E_Int firstChild)
    {
      // set them in PGs
      E_Int* q41 = PGs.get_facets_ptr(firstChild);
      E_Int* q42 = PGs.get_facets_ptr(firstChild+1);
      E_Int* q43 = PGs.get_facets_ptr(firstChild+2);
      E_Int* q44 = PGs.get_facets_ptr(firstChild+3);

      q41[0] = nodes[0]; q41[1] = nodes[4]; q41[2] = nodes[8]; q41[3] = nodes[7];
      q42[0] = nodes[4]; q42[1] = nodes[1]; q42[2] = nodes[5]; q42[3] = nodes[8];
      q43[0] = nodes[8]; q43[1] = nodes[5]; q43[2] = nodes[2]; q43[3] = nodes[6];
      q44[0] = nodes[7]; q44[1] = nodes[8]; q44[2] = nodes[6]; q44[3] = nodes[3];
    }
  };
}

#endif

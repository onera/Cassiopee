/*
 
 
 
              NUGA 
 
 
 
 */
//Authors : Sâm Landier (sam.landier@onera.fr), Alexis Rouil (alexis.rouil@onera.fr)

#ifndef __K_MESH_BASIC_H
#define __K_MESH_BASIC_H

#include "Nuga/include/Quadrangle.h"
#include "Nuga/include/Hexahedron.h"
#include "Nuga/include/Tetrahedron.h"
#include "Nuga/include/Pyramid.h"
#include "Nuga/include/Prism.h"

#include "Nuga/include/DynArray.h"
#include "Nuga/include/ArrayAccessor.h"
#include "Nuga/include/ngon_t.hxx"
#include "Nuga/include/defs.h"

namespace K_MESH
{

class Basic {
public: 
    static constexpr E_Int NB_BOUNDS=6;//fixme

    //
    template< typename ngo_t>
    static void reorder_pgs(ngo_t& ng, const K_FLD::IntArray& F2E, E_Int i);
    //
    template <typename ngunit_t>
    static inline void iso_barycenter(const K_FLD::FloatArray& crd, const ngunit_t & PGs, const E_Int* first_pg, E_Int nb_pgs, E_Int index_start, E_Float* G);
    //
    E_Float quality(const K_FLD::FloatArray& crd, E_Float* Vol){return 1;}
};

///
template< typename ngo_t>
void Basic::reorder_pgs(ngo_t& ng, const K_FLD::IntArray& F2E, E_Int i)
{
  if (Hexahedron::is_of_type(ng.PGs, ng.PHs.get_facets_ptr(i), ng.PHs.stride(i)))
    Hexahedron::reorder_pgs(ng, F2E, i);
  else if (Tetrahedron::is_of_type(ng.PGs, ng.PHs.get_facets_ptr(i), ng.PHs.stride(i)))
    Tetrahedron::reorder_pgs(ng, F2E, i);
  else if (Pyramid::is_of_type(ng.PGs, ng.PHs.get_facets_ptr(i), ng.PHs.stride(i)))
    Pyramid::reorder_pgs(ng, F2E, i);
  else if (Prism::is_of_type(ng.PGs, ng.PHs.get_facets_ptr(i), ng.PHs.stride(i)))
    Prism::reorder_pgs(ng, F2E, i);
}

///
template <typename ngunit_t>
void Basic::iso_barycenter(const K_FLD::FloatArray& crd, const ngunit_t & PGs, const E_Int* first_pg, E_Int nb_pgs, E_Int index_start, E_Float* G)
{    
  if (Hexahedron::is_of_type(PGs, first_pg, nb_pgs)){
    Hexahedron::iso_barycenter(crd, PGs, first_pg, nb_pgs, index_start, G);
  }
  else if (Tetrahedron::is_of_type(PGs, first_pg, nb_pgs)){
    Tetrahedron::iso_barycenter(crd, PGs, first_pg, nb_pgs, index_start, G);
  }
  else if (Pyramid::is_of_type(PGs, first_pg, nb_pgs)){
    Pyramid::iso_barycenter(crd, PGs, first_pg, nb_pgs, index_start, G);
  }
  else if (Prism::is_of_type(PGs, first_pg, nb_pgs)){
    Prism::iso_barycenter(crd, PGs, first_pg, nb_pgs, index_start, G);
  }
}

}
#endif /* BASIC_H */


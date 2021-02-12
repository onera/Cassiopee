/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

#ifndef _GEOMALGO_H_
#define _GEOMALGO_H_

#include "Nuga/include/DynArray.h"
#include "Nuga/include/Edge.h"
#include "Nuga/include/macros.h"

namespace NUGA {

template <typename ElementType>
class GeomAlgo {
public:
  ///
  template <typename Coordinate_t, typename Connectivity_t, int DIM>
  static void neighboring (const K_FLD::ArrayAccessor<Coordinate_t>& coords,
                           const K_FLD::ArrayAccessor<Connectivity_t>& conn, Connectivity_t& neighbors);
  ///
  inline static void reversi_chimera_skin (const Vector_t<K_FLD::FloatArray*>& crds, const Vector_t<K_FLD::IntArray*>& cnts, bool outward=true);
  
  /// Based on quality
  inline static void get_swapE (const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect, const K_FLD::IntArray& neighbors, 
                                const std::set<K_MESH::NO_Edge>& hE, E_Float tol, std::vector<std::pair<E_Int, E_Int> >& swapE);
  /// Based on quality (with input har edges)
  //inline static void get_swapE (const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect, const K_FLD::IntArray& neighbors, 
   //                             E_Float tol, std::vector<std::pair<E_Int, E_Int> >& swapE);
  
  ///
  inline static E_Float angle_measure(const E_Float* ni, const E_Float* nj, const E_Float* E0, const E_Float* E1);

  ///
  inline static void min_quality(const K_FLD::FloatArray& crd, const K_FLD::IntArray& cnt, E_Float& minq, E_Int& imin);
  
private:

  GeomAlgo(void){}

  ~GeomAlgo(void){}
};

} // End namespace NUGA

#include "GeomAlgo.cxx"

#endif

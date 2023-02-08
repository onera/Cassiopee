/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

#ifndef _DELAUNAY_TRIANGULATOR_H_
#define _DELAUNAY_TRIANGULATOR_H_

#include "Nuga/include/T3Mesher.h"

namespace DELAUNAY
{

class Triangulator
{
public:
  Triangulator();
  E_Int run(const K_FLD::FloatArray& coord, const E_Int* pNodes, E_Int nb_nodes, E_Int index_start, K_FLD::IntArray& connectM, K_FLD::IntArray& neighbors, bool do_not_shuffle, bool improve_quality) const ;
  E_Int run(const K_FLD::FloatArray& coord, const E_Int* pNodes, E_Int nb_nodes, E_Int index_start, K_FLD::IntArray& connectM, bool do_not_shuffle, bool improve_quality) const ;
  
#ifdef FLAG_STEP
  mutable E_Float tcreate, tconnect, trun, ttra, tcompact;
#endif
  
private:
  static E_Int __set_connectE2 (const E_Int* pNodes, E_Int nb_nodes, K_FLD::IntArray& connectE2, E_Int index_start);
  
  Triangulator(const Triangulator&);
  
  mutable std::vector<std::pair<E_Int, E_Int> > _swapE;

private:
  mutable DELAUNAY::T3Mesher<E_Float> _mesher;
  mutable DELAUNAY::MeshData _data;

  mutable K_FLD::IntArray connectE2, connectE2b;
  mutable K_FLD::FloatArray Wpos;
  mutable Vector_t<E_Int> oldIds;
  
#ifdef DEBUG_TRIANGULATOR
public:
  static bool dbg_enabled;
#endif
};
}

#endif

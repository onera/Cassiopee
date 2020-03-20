

/*
 
 
 
              NUGA 
 
 
 
 */
//Author : Sâm Landier (sam.landier@onera.fr)

#ifndef NUGA_COLLIDER_HXX
#define NUGA_COLLIDER_HXX

#include "Nuga/include/macros.h"
#include <vector>
#include "Fld/FldArray.h"
#include "Fld/ArrayWriter.h"
#include "Nuga/Delaunay/Triangulator.h"
#include "MeshElement/Triangle.h"
#include "Fld/ngon_unit.h"
#ifdef COLLIDE_DBG
#include "IO/io.h"
#endif


namespace NUGA
{
  namespace COLLIDE
  {
    
    enum eOVLPTYPE { NONE, OVERSET, ABUTTING, INTERSECT};

    using collision_func = bool (*) (const E_Float* P0, const E_Float* P1, const E_Float* P2, const E_Float* Q0, const E_Float* Q1, const E_Float* Q2, E_Float tol);
    
    // Valid for any VOL or SURF elements
    template <typename crd_t, int DIM, typename ELT1, typename ELT2>
    bool simplicial_colliding(const crd_t& crd1, ELT1& e1, const crd_t& crd2, ELT2& e2, collision_func COLLIDING_FUNC, E_Float abstol, E_Int& t1, E_Int& t2)
    {
      return false;
    }
    
    ///
    template <typename acrd_t, typename acnt_t>
    bool __reorient(acrd_t& acrd1, acnt_t*& cnT4)
    {
      return true;
    }

    ///  
    template <typename crd_t, typename cnt_t1, typename cnt_t2, typename loc_t, typename ELT1, typename ELT2, int DIM>
    void compute(const crd_t& crd1,  short px1, short py1, short pz1, const cnt_t1& cnt1, 
                  const crd_t& crd2,  short px2, short py2, short pz2, const cnt_t2& cnt2, 
                  std::vector<E_Int>& ids1, std::vector<E_Int>& ids2, E_Float tolerance = EPSILON)
    {
      assert (false);
    }
 
/// fixme : this function must be reworked : logic inside could be wrong and assume frank X (by using fast_intersectT3 that doesnt behave well otherwise) 
template <typename crd_t, typename cnt_t1, typename cnt_t2, typename loc_t, typename ELT1, typename ELT2, int DIM>
void compute(const K_FLD::ArrayAccessor<crd_t>& acrd1,  const cnt_t1& cnt1, const loc_t& loc1,
             const crd_t& crd2,  short px2, short py2, short pz2, const cnt_t2& cnt2, 
             std::vector<E_Int>& is_x1, std::vector<E_Int>& is_x2, E_Float tolerance = EPSILON)
{
}


///
template <typename ELT1, typename ELT2, typename loc_t>
void compute_overlap(const K_FLD::FloatArray& crd1, const ngon_unit& PGs1,
                     const K_FLD::FloatArray& crd2, const ngon_unit& PGs2, const loc_t& loc2, 
                     double ps_min,//overlap criterion
                     std::vector<E_Int>& is_x1, std::vector<E_Int>& is_x2,
                     E_Float RTOL, bool swap, const E_Float* norm2 = nullptr) //norm2 when right direction is known upon entry
{
}

///
template<typename loc_t>
void compute_overlap(const K_FLD::FloatArray& crd1, const K_FLD::IntArray& edges1,
                     const K_FLD::FloatArray& crd2, const K_FLD::IntArray& edges2, const loc_t& loc2, 
                     std::vector<eOVLPTYPE>& is_x1, std::vector<eOVLPTYPE>& is_x2,
                     E_Float RTOL)
{
}
  
} //COLLIDE
}   // NUGA

#endif /* NUGA_COLLIDER_HXX */


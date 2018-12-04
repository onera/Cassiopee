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


namespace NUGA
{
  namespace COLLIDE
  {
    // Valid for any VOL or SURF elements
    template <typename acrd_t, int DIM, typename ELT1, typename ELT2>
    bool simplicial_colliding(const acrd_t& acrd1, ELT1& e1, const acrd_t& acrd2, ELT2& e2)
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
      //assert (false);
    }
 
    ///
    template <typename crd_t, typename cnt_t1, typename cnt_t2, typename loc_t, typename ELT1, typename ELT2, int DIM>
    void compute(const K_FLD::ArrayAccessor<crd_t>& acrd1,  const cnt_t1& cnt1, const loc_t& loc1,
                  const crd_t& crd2,  short px2, short py2, short pz2, const cnt_t2& cnt2, 
                  std::vector<E_Int>& is_x1, std::vector<E_Int>& is_x2, E_Float tolerance = EPSILON)
    {
    }
  }
}

#endif /* NUGA_COLLIDER_HXX */


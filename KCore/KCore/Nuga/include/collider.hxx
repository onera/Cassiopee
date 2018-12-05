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
      //
      E_Int T1[3], T2[3] = {0};
      E_Float P0[DIM], P1[DIM], P2[DIM], Q0[DIM], Q1[DIM], Q2[DIM];
      
      DELAUNAY::Triangulator dt;
      
      e1.triangulate(dt, acrd1);
      e2.triangulate(dt, acrd2);
  
      for (E_Int i=0; i < ELT1::NB_TRIS; ++i)
      {
        e1.triangle(i, T1);
        acrd1.getEntry(T1[0], P0);
        acrd1.getEntry(T1[1], P1);
        acrd1.getEntry(T1[2], P2);
      
        E_Int nb_tris2 = e2.nb_tris();
        for (E_Int j=0; j < nb_tris2; ++j)
        {
          e2.triangle(j, T2);
          acrd2.getEntry(T2[0], Q0);
          acrd2.getEntry(T2[1], Q1);
          acrd2.getEntry(T2[2], Q2);
      
          if (K_MESH::Triangle::fast_intersectT3<DIM>(P0, P1, P2, Q0, Q1, Q2))
            return true;
        }
      }
      
      return false;
    }
    
    ///
    template <typename acrd_t, typename acnt_t>
    bool __reorient(acrd_t& acrd1, acnt_t*& cnT4)
    {
      std::vector<size_t> elts_to_reorient;
      size_t nb_elts(cnT4->size());
      const E_Int NB_NODES = cnT4->stride();
      E_Float Ni[3], Nj[3], Nk[3], Nl[3];
      E_Int t4[NB_NODES];
      
      using cnt_t = typename acnt_t::ArrayType;

      for (size_t i = 0; i < nb_elts; ++i)
      {
        cnT4->getEntry(i, t4);

        acrd1.getEntry(t4[0], Ni);
        acrd1.getEntry(t4[1], Nj);
        acrd1.getEntry(t4[2], Nk);
        acrd1.getEntry(t4[3], Nl);

        if (K_FUNC::zzdet4(Ni, Nj, Nk, Nl) < -E_EPSILON) //must be reoriented
          elts_to_reorient.push_back(i);
      }

      size_t s = elts_to_reorient.size();

      if (s == 0)
        return false;
      //need to create a copy of the input connectivity
      cnt_t* new_connect = new cnt_t(cnT4->array());
      K_FLD::ArrayWriter<cnt_t> writer(*cnT4, *new_connect);
      //permut the first and second nodes
      for (size_t i = 0; i < s; ++i)
        std::swap(writer.getVal(elts_to_reorient[i], 0), writer.getVal(elts_to_reorient[i], 1));

      //reassign the new array to the accessor
      E_Int shift = cnT4->shift();
      delete cnT4;
      cnT4 = new K_FLD::ArrayAccessor<cnt_t>(*new_connect, shift);

      return true;
    }

    ///  
    template <typename crd_t, typename cnt_t1, typename cnt_t2, typename loc_t, typename ELT1, typename ELT2, int DIM>
    void compute(const crd_t& crd1,  short px1, short py1, short pz1, const cnt_t1& cnt1, 
                  const crd_t& crd2,  short px2, short py2, short pz2, const cnt_t2& cnt2, 
                  std::vector<E_Int>& ids1, std::vector<E_Int>& ids2, E_Float tolerance = EPSILON)
    {
      assert (false);
//      using acnt1_t = 
//      
//      acrd_t acrd1(crd1, px1, py1, pz1);
//      acnt1_t acnt1(cnt1, -1);
//      
//      __reorient(acrd1, acnt1);
//      
//      loc_t loc1(acrd1, acnt1);
//      
//      compute(acrd1,  cnt1, loc1, crd2,  px2, py2, pz2, cnt2, ids1, ids2, tolerance);

    }
 
    ///
    template <typename crd_t, typename cnt_t1, typename cnt_t2, typename loc_t, typename ELT1, typename ELT2, int DIM>
    void compute(const K_FLD::ArrayAccessor<crd_t>& acrd1,  const cnt_t1& cnt1, const loc_t& loc1,
                  const crd_t& crd2,  short px2, short py2, short pz2, const cnt_t2& cnt2, 
                  std::vector<E_Int>& is_x1, std::vector<E_Int>& is_x2, E_Float tolerance = EPSILON)
    {
      using acrd_t = K_FLD::ArrayAccessor<crd_t>;
      using acnt1_t = K_FLD::ArrayAccessor<cnt_t1>;
      using acnt2_t = K_FLD::ArrayAccessor<cnt_t2>;
      
      acnt1_t acnt1(cnt1, -1);
      acrd_t acrd2(crd2, px2, py2, pz2);
      acnt2_t acnt2(cnt2, -1);

      size_t sz2((size_t)acnt2.size());

      is_x2.resize(sz2, 0);
      is_x1.resize(acnt1.size(), 0);

      std::vector<E_Int> cands1;
      ELT1 e1;
      ELT2 e2;
      size_t i1;

#ifndef COLLIDE_DBG
#pragma omp parallel for private (cands1, e1, e2, i1)
#endif
      for (size_t i2 = 0; i2 < sz2; ++i2)
      {
        
        acnt2.getEntry(i2, e2);

        loc1.get_candidates(e2, acrd2, cands1);

        for (i1=0; (i1<cands1.size()) && !is_x2[i2]; ++i1)
        {
          acnt1.getEntry(cands1[i1], e1);

          is_x1[i1] = is_x2[i2] = simplicial_colliding<acrd_t, DIM>(acrd1, e1, acrd2, e2);
        }
      }
    }
  }
}

#endif /* NUGA_COLLIDER_HXX */


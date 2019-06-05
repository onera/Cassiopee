

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
#ifdef COLLIDE_DBG
#include "IO/io.h"
#endif


namespace NUGA
{
  namespace COLLIDE
  {
    
    using collision_func = bool (*) (const E_Float* P0, const E_Float* P1, const E_Float* P2, const E_Float* Q0, const E_Float* Q1, const E_Float* Q2, E_Float tol);
    
    // Valid for any VOL or SURF elements
    template <typename acrd_t, int DIM, typename ELT1, typename ELT2>
    bool simplicial_colliding(const acrd_t& acrd1, ELT1& e1, const acrd_t& acrd2, ELT2& e2, collision_func COLLIDING_FUNC, E_Float abstol)
    {
      //
      E_Int T1[3], T2[3] = {0};
      E_Float P0[DIM], P1[DIM], P2[DIM], Q0[DIM], Q1[DIM], Q2[DIM];
      
      DELAUNAY::Triangulator dt;
      
      e1.triangulate(dt, acrd1);
      e2.triangulate(dt, acrd2);

      for (E_Int i=0; i < e1.nb_tris(); ++i)
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
          
#ifdef COLLIDE_DBG

//  K_FLD::IntArray tmpE;
//  for (size_t i=0; i < 3; ++i)
//  {
//    E_Int E[] = {i,(i+1)%3};
//    tmpE.pushBack(E, E+2);
//  }
//  for (size_t i=0; i < 3; ++i)
//  {
//    E_Int E[] = {3+i,3+(i+1)%3};
//    tmpE.pushBack(E, E+2);
//  }
//
//  K_FLD::FloatArray crd;
//  crd.pushBack(P0,P0+3);
//  crd.pushBack(P1,P1+3);
//  crd.pushBack(P2,P2+3);
//  crd.pushBack(Q0,Q0+3);
//  crd.pushBack(Q1,Q1+3);
//  crd.pushBack(Q2,Q2+3);
//
//  MIO::write("conf.mesh", crd, tmpE, "BAR");
          
#endif
      
          if (COLLIDING_FUNC(P0, P1, P2, Q0, Q1, Q2, abstol))
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
      
      using cnt_t = typename acnt_t::array_type;

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
 
    /// fixme : this function must be reworked : logic inside could be wrong and assume frank X (by using fast_intersectT3 that doesnt behave well otherwise) 
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

        for (i1=0; (i1<cands1.size()) && !is_x2[i2]; ++i1) //fixme : logic seems wrong : if i2 is already colliding, it cannot be used to invalidate any i2
        {
          acnt1.getEntry(cands1[i1], e1);

          is_x1[i1] = is_x2[i2] = simplicial_colliding<acrd_t, DIM>(acrd1, e1, acrd2, e2, K_MESH::Triangle::fast_intersectT3<3>, -1.); //do not pass the tol as we use here a predicate
        }
      }
    }


///
template <typename ELT1, typename ELT2, typename loc_t>
void compute_overlap(const K_FLD::FloatArray& crd1, const ngon_unit& PGs1,
                     const K_FLD::FloatArray& crd2, const ngon_unit& PGs2, const loc_t& loc2, 
                     double ps_min,//overlap criterion
                     std::vector<E_Int>& is_x1, std::vector<E_Int>& is_x2,
                     E_Float RTOL, bool swap, const E_Float* norm2 = nullptr) //norm2 when right direction is known upon entry
{
  is_x1.clear();
  is_x2.clear();
  
  E_Int nb_pgs1 = PGs1.size();
  
  is_x1.resize(nb_pgs1, 0);
  is_x2.resize(PGs2.size(), 0);
  
  using acrd_t = K_FLD::ArrayAccessor<K_FLD::FloatArray>;
  //using acnt_t = K_FLD::ArrayAccessor<K_FLD::IntArray>;
  
  acrd_t acrd1(crd1);
  acrd_t acrd2(crd2);
  
  std::vector<E_Int> cands2;
  E_Float n1[3], n2[3];
  size_t j;
  
  if (norm2 != nullptr) //prescribed orthogonal direction
  {
    n2[0] = norm2[0]; n2[1] = norm2[1]; n2[2] = norm2[2];
    K_FUNC::normalize<3>(n2);
  }
  
  // for nodal tolerance
  K_FLD::FloatArray L;
  K_CONNECT::MeshTool::computeIncidentEdgesSqrLengths(crd1, PGs1, L);
  
#ifndef COLLIDE_DBG
//#pragma omp parallel private (cands2, n1, n2, j)
//#pragma omp for
#endif
  for (E_Int i = 0; i < nb_pgs1; ++i)
  {
    const E_Int* nodes = PGs1.get_facets_ptr(i);
    E_Int nb_nodes = PGs1.stride(i);
    
    ELT1 PG1(nodes, nb_nodes, -1);

    loc2.get_candidates(PG1, acrd1, cands2);
    
#ifdef COLLIDE_DBG
  /*if (!cands2.empty())
  {
    ngon_unit pgs;
    for (size_t i=0; i < cands2.size(); ++i)
      pgs.add(PGs2.stride(cands2[i]), PGs2.get_facets_ptr(cands2[i]));
    pgs.updateFacets();
    NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>::draw_PGs("cands", crd2, pgs);
  }*/
#endif
    
    if (cands2.empty()) continue;

    ELT1::template normal<acrd_t, 3>(acrd1, nodes, nb_nodes, 1, n1);
    
    // Tolerance hypothese : constant per PG : take the min over the polygon nodes
    E_Float Lref = K_CONST::E_MAX_FLOAT;
    for (E_Int n=0; n < nb_nodes;++n)
      Lref = MIN(Lref, L(0,nodes[n]-1));

    E_Float abstol = MAX(E_EPSILON, RTOL*::sqrt(Lref));
    //std::cout << abstol << std::endl;
    bool isx = false;
    for (j = 0; (j < cands2.size()) && !isx; ++j)
    {
      E_Int J = cands2[j];

      const E_Int* nodes2 = PGs2.get_facets_ptr(J);
      E_Int nb_nodes2 = PGs2.stride(J);
      
      if (norm2 == nullptr)
        ELT2::template normal<acrd_t, 3>(acrd2, nodes2, nb_nodes2, 1, n2);
      
      double ps = ::fabs(K_FUNC::dot<3>(n1,n2));
      if (ps < ps_min) continue;
      
      // Polygons pairs are now roughly overlapping/parallel
      
      ELT2 PG2(nodes2, nb_nodes2, -1);
      
      isx = NUGA::COLLIDE::simplicial_colliding<acrd_t, 3>(acrd1, PG1, acrd2, PG2, K_MESH::Triangle::overlap, abstol);
      
      if (swap && !isx) //deeper test
      {
        PG1.shuffle_triangulation();
        PG2.shuffle_triangulation();
        
        isx = NUGA::COLLIDE::simplicial_colliding<acrd_t, 3>(acrd1, PG1, acrd2, PG2, K_MESH::Triangle::overlap, abstol);
      }
      //assert (i > -1 && i < is_x1.size());
      //assert (J > -1 && J < is_x2.size());
      is_x1[i] |= isx;
      is_x2[J] |= isx;
    }
  }
}
  } //COLLIDE
}   // NUGA

#endif /* NUGA_COLLIDER_HXX */


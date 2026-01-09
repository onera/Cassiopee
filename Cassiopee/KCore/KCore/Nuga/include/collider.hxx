/*    
    Copyright 2013-2026 ONERA.

    This file is part of Cassiopee.

    Cassiopee is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Cassiopee is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.
*/
//Authors : Sam Landier (sam.landier@onera.fr)

#ifndef NUGA_COLLIDER_HXX
#define NUGA_COLLIDER_HXX

#include "Nuga/include/macros.h"
#include <vector>
#include "Nuga/include/ArrayWriter.h"
#include "Nuga/include/Triangulator.h"
#include "Nuga/include/Triangle.h"
#include "Nuga/include/ngon_unit.h"
#include "Nuga/include/maths.hxx"
#include "Nuga/include/mesh_t.hxx"
#include "Nuga/include/vertex.h"

//#define COLLIDER_DBG
#ifdef COLLIDER_DBG
//#include "IO/io.h"
#include "Nuga/include/medit.hxx"
#endif


namespace NUGA
{
  namespace COLLIDE
  {
    using collision_func = bool (*) (const E_Float* P0, const E_Float* P1, const E_Float* P2, const E_Float* Q0, const E_Float* Q1, const E_Float* Q2, E_Float tol);
    
    // Valid for any VOL or SURF elements
    template <typename crd_t, short DIM, typename ELT1, typename ELT2>
    bool simplicial_colliding(const crd_t& crd1, ELT1& e1, const crd_t& crd2, ELT2& e2, collision_func COLLIDING_FUNC, E_Float abstol, E_Int& t1, E_Int& t2)
    {
      //
      E_Int T1[3], T2[3] = {0};
      typename crd_t::pt_t P0, P1, P2, Q0, Q1, Q2;
            
      t1=t2 = IDX_NONE;
      
      //DELAUNAY::Triangulator dt;
      //e1.triangulate(dt, crd1);
      //e2.triangulate(dt, crd2);
      e1.cvx_triangulate(crd1);
      e2.cvx_triangulate(crd2);

      for (E_Int i=0; i < e1.nb_tris(); ++i)
      {
        e1.triangle(i, T1);
        P0 = crd1.col(T1[0]);
        P1 = crd1.col(T1[1]);
        P2 = crd1.col(T1[2]);

        K_SEARCH::BBox3D b1;
        b1.compute(crd1, T1, 3, 0);
        if (abstol > 0.) b1.enlarge(abstol);
      
        E_Int nb_tris2 = e2.nb_tris();
        for (E_Int j=0; j < nb_tris2; ++j)
        {
          e2.triangle(j, T2);
          Q0 = crd2.col(T2[0]);
          Q1 = crd2.col(T2[1]);
          Q2 = crd2.col(T2[2]);

          K_SEARCH::BBox3D b2;
          b2.compute(crd2, T2, 3, 0);
          if (abstol > 0.) b2.enlarge(abstol);

          // FILTER 1 : faces 3D boxes
          if (!K_SEARCH::BbTree3D::boxesAreOverlapping(&b1, &b2, EPSILON)) // face-box test worth it
          {
            /*ngon_type ng1, ng2;
            K_FLD::FloatArray c1, c2;
            b1.convert2NG(c1, ng1);
            b2.convert2NG(c2, ng2);
            medith::write("b1", c1, ng1);
            medith::write("b2", c2,  ng2);*/
            continue;
          }

          if (COLLIDING_FUNC(P0, P1, P2, Q0, Q1, Q2, abstol))
          {
#ifdef COLLIDER_DBG

            K_FLD::IntArray tmpT(3, 2);
            tmpT(0, 0) = 0;
            tmpT(1, 0) = 1;
            tmpT(2, 0) = 2;
            tmpT(0, 1) = 3;
            tmpT(1, 1) = 4;
            tmpT(2, 1) = 5;

            K_FLD::FloatArray crd;
            crd.pushBack(P0, P0 + 3);
            crd.pushBack(P1, P1 + 3);
            crd.pushBack(P2, P2 + 3);

            K_SEARCH::BBox3D b(crd);
            b.enlarge(0.05);

            K_FLD::FloatArray crd1;
            crd1.pushBack(Q0, Q0 + 3);
            crd1.pushBack(Q1, Q1 + 3);
            crd1.pushBack(Q2, Q2 + 3);

            K_SEARCH::BBox3D b1(crd1);
            b1.enlarge(0.05);

            crd.pushBack(crd1);

            medith::write("conf", crd, tmpT, "TRI");

            K_SEARCH::BBox3D bx;
            bool are_x_box = K_SEARCH::BBox3D::intersection(b, b1, bx);

            //std::cout << "box are colliding ? ! " << are_x_box << std::endl;
#endif
            t1=i; t2=j;
            return true;
          }
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

        if (NUGA::zzdet4(Ni, Nj, Nk, Nl) < -EPSILON) //must be reoriented
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
    template <typename crd_t, typename cnt_t1, typename cnt_t2, typename loc_t, typename ELT1, typename ELT2, short DIM>
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
    template <typename crd_t, typename cnt_t1, typename cnt_t2, typename loc_t, typename ELT1, typename ELT2, short DIM>
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

#ifndef COLLIDER_DBG
#pragma omp parallel for private (cands1, e1, e2, i1)
#endif
      for (size_t i2 = 0; i2 < sz2; ++i2)
      {
        
        acnt2.getEntry(i2, e2);

        loc1.template get_candidates<ELT2, acrd_t>(e2, acrd2, cands1);

        for (i1=0; (i1<cands1.size()) && !is_x2[i2]; ++i1) //fixme : logic seems wrong : if i2 is already colliding, it cannot be used to invalidate any i2
        {
          acnt1.getEntry(cands1[i1], e1);
          E_Int t1,t2;
          is_x1[i1] = is_x2[i2] = simplicial_colliding<acrd_t, DIM>(acrd1, e1, acrd2, e2, K_MESH::Triangle::fast_intersectT3<3>, -1., t1, t2); //do not pass the tol as we use here a predicate
        }
      }
    }

///
template <typename ELT1, typename ELT2, typename loc_t>
void compute_overlap(const K_FLD::FloatArray& crd1, const ngon_unit& PGs1,
                     const K_FLD::FloatArray& crd2, const ngon_unit& PGs2, const loc_t& loc2, 
                     std::vector<E_Int>& is_x1/*1-based w negval for abutt*/,
                     std::vector<E_Int>& is_x2/*1-based w negval for abutt*/,
                     E_Float RTOL, 
                     E_Float ps_min=0.99/*overlap criterion*/, bool swap = true, const E_Float* norm2 = nullptr) //norm2 when right direction is known upon entry
{
  is_x1.clear();
  is_x2.clear();
  
  E_Int nb_pgs1 = PGs1.size();
  E_Int it1, it2;
  
  is_x1.resize(nb_pgs1, IDX_NONE);
  is_x2.resize(PGs2.size(), IDX_NONE);
  
  using crd_t = K_FLD::FloatArray;
    
  std::vector<E_Int> cands2;
  E_Float n1[3], n2[3];
  size_t j;
  
  if (norm2 != nullptr) //prescribed orthogonal direction
  {
    n2[0] = norm2[0]; n2[1] = norm2[1]; n2[2] = norm2[2];
    NUGA::normalize<3>(n2);
  }
  
  // for nodal tolerance
  K_FLD::FloatArray L;
  NUGA::MeshTool::computeIncidentEdgesSqrLengths(crd1, PGs1, L);
  
#ifndef COLLIDER_DBG
//#pragma omp parallel private (cands2, n1, n2, j)
//#pragma omp for
#endif
  for (E_Int i = 0; i < nb_pgs1; ++i)
  {
    const E_Int* nodes = PGs1.get_facets_ptr(i);
    E_Int nb_nodes = PGs1.stride(i);
    
    ELT1 PG1(nodes, nb_nodes, -1);

    loc2.get_candidates(PG1, crd1, cands2);
    
#ifdef COLLIDER_DBG
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

    ELT1::template normal<crd_t, 3>(crd1, nodes, nb_nodes, 1, n1);
    
    // Tolerance hypothese : constant per PG : take the min over the polygon nodes
    E_Float Lref = NUGA::FLOAT_MAX;
    for (E_Int n=0; n < nb_nodes;++n)
      Lref = MIN(Lref, L(0,nodes[n]-1));

    E_Float abstol = MAX(ZERO_M, RTOL*sqrt(Lref));
    for (j = 0; (j < cands2.size()); ++j)
    {
      E_Int J = cands2[j];

      const E_Int* nodes2 = PGs2.get_facets_ptr(J);
      E_Int nb_nodes2 = PGs2.stride(J);
      
      if (norm2 == nullptr)
        ELT2::template normal<crd_t, 3>(crd2, nodes2, nb_nodes2, 1, n2);
      
      E_Float ps = NUGA::dot<3>(n1,n2);
      if (::fabs(ps) < ps_min) continue;
      
      // Polygons pairs are now roughly overlapping/parallel : important in order to have a relevant result when using simplicial_colliding with Triangle::overlap
      
      ELT2 PG2(nodes2, nb_nodes2, -1);
      
      bool isx = NUGA::COLLIDE::simplicial_colliding<crd_t, 3>(crd1, PG1, crd2, PG2, K_MESH::Triangle::overlap, abstol, it1, it2);
      
      if (swap && !isx) //deeper test
      {
        PG1.shuffle_triangulation();
        PG2.shuffle_triangulation();
        
        isx = NUGA::COLLIDE::simplicial_colliding<crd_t, 3>(crd1, PG1, crd2, PG2, K_MESH::Triangle::overlap, abstol, it1, it2);
      }
      //assert (i > -1 && i < is_x1.size());
      //assert (J > -1 && J < is_x2.size());
      if (isx)
      {
        if (ps < 0.) // ABUTTING
        {
          is_x1[i] = -(J+1); // 1-based because of neg val logic
          is_x2[J] = -(i+1);
        }
        else         // OVERSET
        {
          is_x1[i] = (J + 1); // 1-based because of neg val logic
          is_x2[J] = (i + 1);
        }
      }
    }
  }
}

///
template <typename ELT1, typename ELT2, typename loc_t>
void compute_overlap(const K_FLD::FloatArray& crd1, const K_FLD::IntArray& edges1,
                     const K_FLD::FloatArray& crd2, const K_FLD::IntArray& edges2, const loc_t& loc2, 
                     std::vector<E_Int>& is_x1/*1-based w negval for abutt*/, std::vector<E_Int>& is_x2 /*1-based w negval for abutt*/,
                     E_Float RTOL)
{
  is_x1.clear();
  is_x2.clear();
  
  E_Int nb_elts1 = edges1.cols();
  //E_Int it1, it2;
  
  is_x1.resize(nb_elts1, IDX_NONE);
  is_x2.resize(edges2.cols(), IDX_NONE);
  
  using acrd_t = K_FLD::ArrayAccessor<K_FLD::FloatArray>;
  //using acnt_t = K_FLD::ArrayAccessor<K_FLD::IntArray>;
  
  acrd_t acrd2(crd2);
  
  std::vector<E_Int> cands2;
  //E_Float n1[3], n2[3];
  size_t j;

  // for nodal tolerance
  K_FLD::FloatArray L;
  NUGA::MeshTool::computeIncidentEdgesSqrLengths(crd1, edges1, L);
  
#ifndef COLLIDER_DBG
//#pragma omp parallel private (cands2, n1, n2, j)
//#pragma omp for
#endif
  for (E_Int i = 0; i < nb_elts1; ++i)
  {
    const E_Int* nodes = edges1.col(i);
    E_Int nb_nodes = 2;
    
    ELT1 E1(nodes);

    loc2.get_candidates(E1, crd1, cands2);
    
#ifdef COLLIDER_DBG
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

    // Tolerance hypothese : constant per PG : take the min over the polygon nodes
    E_Float Lref = NUGA::FLOAT_MAX;
    for (E_Int n=0; n < nb_nodes;++n)
      Lref = MIN(Lref, L(0,nodes[n]-1));

    E_Float abstol = MAX(EPSILON, RTOL*sqrt(Lref));
    E_Float u00, u01, u10, u11;
    E_Bool overlap = false;
    
    for (j = 0; (j < cands2.size()); ++j)
    {
      E_Int J = cands2[j];

      const E_Int* nodes2 = edges2.col(J);
      //E_Int nb_nodes2 = 2;
      
      // Polygons pairs are now roughly overlapping/parallel : important in order to have a relevant result when using simplicial_colliding with Triangle::overlap

      ELT2 E2(nodes2);

      K_MESH::Edge::intersect<3> (crd1.col(nodes[0]), crd1.col(nodes[1]), 
                               crd2.col(nodes2[0]), crd2.col(nodes2[1]), 
                               abstol, true/*tol_is_absolute*/,
                               u00, u01, u10, u11, overlap);
      //assert (i > -1 && i < is_x1.size());
      //assert (J > -1 && J < is_x2.size());

      if (overlap)
      {
        E_Float d1[3], d2[3];
        NUGA::diff<3>(crd1.col(nodes[1]),  crd1.col(nodes[0]),  d1);
        NUGA::diff<3>(crd2.col(nodes2[1]), crd2.col(nodes2[0]), d2);
        E_Float ps = NUGA::dot<3>(d1,d2);

        if (ps < 0.) // ABUTTING
        {
          is_x1[i] = -(J + 1); // 1-based because of neg val logic
          is_x2[J] = -(i + 1);
        }
        else        // OVERSET
        {
          is_x1[i] = J + 1; // 1-based because of neg val logic
          is_x2[J] = i + 1;
        }
      }
    }
  }
}

///
template <typename aelt_t, typename bound_mesh_t>
bool get_colliding(const aelt_t& e1, const bound_mesh_t& mask_bit, std::vector<E_Int>& cands, E_Int idx_start, E_Float RTOL, bool only_first_found);

/// Polygon vs polyLine impl.
template <> inline
bool get_colliding<NUGA::aPolygon, edge_mesh_t>
(const NUGA::aPolygon& ae1, const edge_mesh_t& mask_bit, std::vector<E_Int>& cands, E_Int idx_start, E_Float RTOL, bool first_found)
{
 
  bool hasX(false);

  // reduce mask to candidates
  edge_mesh_t lmask(mask_bit, cands, idx_start);

  E_Float normal1[3];
  ae1.normal<3>(normal1);
  const E_Float* plane_pt = ae1.m_crd.col(0); // choosing first node

#ifdef DEBUG_XCELLN
  //medith::write("cutter_front_before_projection", lmask.crd, lmask.cnt);
#endif

  E_Float l21 = ae1.Lref2();
  E_Float l22 = lmask.Lref2();//compute it before proj
  E_Float ATOL(RTOL * sqrt(std::min(l21, l22)));

  // project candidates on e1's plane => 2D problem
  STACK_ARRAY(E_Float, lmask.crd.cols(), signed_dists);
  for (int i = 0; i < lmask.crd.cols(); ++i)
  {
    // orthogonal projection on a plane (projection dir is the normal to the plane)
    // overwrite it
    signed_dists[i] = NUGA::project(plane_pt, normal1, lmask.crd.col(i), normal1, lmask.crd.col(i));
  }

  // discard those that are too far
  std::vector<bool> keep(cands.size(), true);

  for (int i = 0; (i < lmask.cnt.cols()); ++i)
  {
    int e1 = lmask.cnt(0, i);
    int e2 = lmask.cnt(1, i);

    E_Float z1 = signed_dists[e1];
    E_Float z2 = signed_dists[e2];

    if (z1*z2 < 0.) continue; // means crossing 

    E_Float mz = std::min(::fabs(z1), ::fabs(z2));
    keep[i] = (mz < ATOL);
  }

  K_CONNECT::keep<bool> pred(keep);
  K_CONNECT::IdTool::compress(cands, pred);

  if (cands.empty()) return false;

  lmask.compress(keep);

  // close enough so go to 2D for real test

  // compute collision between e1 and each candidate until founding one collision

  // a. triangulate e1, check if any lmask point falls into the triangulation
  DELAUNAY::Triangulator dt;
  ae1.triangulate(dt);

  // b. go 2D both crd (a copy of it)  and lmask.crd

  // Computes the fitting box coordinate system optimizing the view over the contour
  K_FLD::FloatArray P(3, 3), iP(3, 3);

  NUGA::computeAFrame(normal1, P);
  iP = P;
  K_FLD::FloatArray::inverse3(iP);

  K_FLD::FloatArray crd2D(ae1.m_crd);
  NUGA::transform(crd2D, iP);// Now we are in the fitting coordinate system.
  crd2D.resize(2, crd2D.cols()); // e1 is 2D now.
  NUGA::transform(lmask.crd, iP);// Now we are in the fitting coordinate system.
  lmask.crd.resize(2, lmask.crd.cols()); // lmask is 2D now.

#ifdef COLLIDER_DBG
  {
    K_FLD::FloatArray c1(crd2D), c2(lmask.crd);
    c1.resize(3, c1.cols(), 0.);
    c2.resize(3, c2.cols(), 0.);
    medith::write("subj2D", c1, ae1.begin(), ae1.nb_nodes(), 0);
    medith::write("cutter_front2D", c2, lmask.cnt);
  }
#endif

  // c. detect collisions
  
  E_Int ncands(cands.size());
  keep.clear();
  keep.resize(ncands, false);

  E_Int T1[3];
  for (int i = 0; (i < ae1.nb_tris()); ++i)
  {
    ae1.triangle(i, T1);
    const E_Float* P1 = crd2D.col(T1[0]);
    const E_Float* Q1 = crd2D.col(T1[1]);
    const E_Float* R1 = crd2D.col(T1[2]);

    for (int j = 0; (j<ncands); ++j)
    {
      if (keep[j]) continue;

      K_MESH::Edge e = lmask.element(j);
      const E_Float* P2 = lmask.crd.col(e.node(0));
      const E_Float* Q2 = lmask.crd.col(e.node(1));

      bool isx = K_MESH::Triangle::fast_intersectE2_2D(P1, Q1, R1, P2, Q2, ATOL);

      if (!isx) //deeper check
      {
        E_Float u00, u01;
        E_Int tx;
        E_Bool overlap;

        isx = K_MESH::Triangle::intersect<2>(P1, Q1, R1, P2, Q2, ATOL, true, u00, u01, tx, overlap);
      }
      if (isx)
      {
        if (first_found)
        {
          cands[0] = cands[j];
          cands.resize(1);
          return true;
        }

        hasX = keep[j] = true;
      }
    }
  }

  if (hasX)
  {
    K_CONNECT::keep<bool> pred1(keep);
    K_CONNECT::IdTool::compress(cands, pred1);
  }
  return hasX;
}

/// Polyhedron vs surface impl.
template <> inline
bool get_colliding<NUGA::aPolyhedron<UNKNOWN>, pg_smesh_t>
(const NUGA::aPolyhedron<UNKNOWN>& ae1, const pg_smesh_t& mask_bit, std::vector<E_Int>& cands, E_Int idx_start, E_Float RTOL, bool first_found)
{
  bool hasX(false);

  // reduce mask to candidates : noy hpc way for above lineic (yet)
  //pg_smesh_t lmask(mask_bit, cands, idx_start);
  const pg_smesh_t& lmask = mask_bit;

#ifdef COLLIDER_DBG
  medith::write("lmask", lmask.crd, lmask.cnt);
#endif

  //E_Float l21 = ae1.Lref2();
  //E_Float l22 = lmask.Lref2(cands, idx_start);
  //E_Float ATOL(RTOL * sqrt(std::min(l21, l22)));

  // compute collision between e1 and each candidate until founding one collision

  // a. triangulate e1, check if any lmask point falls into the triangulation
  //DELAUNAY::Triangulator dt;
  //ae1.triangulate(dt);
  //ae1.cvx_triangulate(ae1.m_crd);

  E_Float Lref21 = ae1.Lref2();

  // (b. projection for 2D)

  // c. detect collisions

  E_Int ncands(cands.size());
  std::vector<bool> keep(ncands, false);

  for (int j = 0; (j<ncands); ++j)
  {
    if (keep[j]) continue;

    K_MESH::Polygon e2 = lmask.element(cands[j] - idx_start);

    E_Float Lref22 = e2.Lref2(lmask.crd);

    E_Float abstol = RTOL * sqrt(std::min(Lref21, Lref22));

    //do not pass the tol as we use here a predicate
    E_Int t1, t2;
    using crd_t = K_FLD::FloatArray;
    bool isx = simplicial_colliding<crd_t, 3>(ae1.m_crd, ae1, lmask.crd, e2, K_MESH::Triangle::fast_intersectT3<3>, abstol, t1, t2);

    if (!isx) //deeper check
    {
      //fixme : required ?
    }
    if (isx)
    {
#ifdef CLASSIFYER_DBG
      //medith::write("xe2", lmask.crd, e2.begin(), e2.nb_nodes(), 1);
#endif
      if (first_found)
      {
        cands[0] = cands[j];
        cands.resize(1);
        return true;
      }

      hasX = keep[j] = true;
    }
  }

  if (hasX)
  {
    K_CONNECT::keep<bool> pred(keep);
    K_CONNECT::IdTool::compress(cands, pred);
  }
  return hasX;
}

/// Vertex vs surface impl.
template <> inline
bool get_colliding<vertex, pg_smesh_t>
(const vertex& pt, const pg_smesh_t& surface, std::vector<E_Int>& cands, E_Int idx_start, E_Float ARTOL, bool first_found_dummy)
{
  cands.clear();
  if (surface.ncells() == 0) return false;

  auto loc = surface.get_localizer();

  if (ARTOL < 0.)// relative
  {
    E_Float Lref = sqrt(pt.val2);
    ARTOL *= -Lref;
  }

  K_SEARCH::BBox3D bb1;
  bb1.minB[0] = pt.vec[0] - ARTOL;
  bb1.minB[1] = pt.vec[1] - ARTOL;
  bb1.minB[2] = pt.vec[2] - ARTOL;
  bb1.maxB[0] = pt.vec[0] + ARTOL;
  bb1.maxB[1] = pt.vec[1] + ARTOL;
  bb1.maxB[2] = pt.vec[2] + ARTOL;
  
  loc->get_candidates(bb1, cands, idx_start, -1./*dummy*/);

  return !cands.empty();
}
/// above wrapper for list of vertices : returns indir 'pt to faces'
inline void get_colliding
(const std::vector<vertex>& pts, const pg_smesh_t& surface, E_Float ARTOL, eMetricType mtype, std::vector<std::vector<E_Int>>& pt_to_elt)
{
  pt_to_elt.clear();
  if (surface.ncells() == 0) return;

  size_t npts = pts.size();

  pt_to_elt.resize(npts);

  if (ARTOL < 0.)
    surface.get_nodal_metric2(mtype); // to compute it if missing

  DELAUNAY::Triangulator dt;
  std::vector<E_Int> cands;
  E_Int T[3];
  E_Float UV[3];

  //
  for (size_t i = 0; i < npts; ++i)
  {
    //const E_Float& x = pts[i].vec[0];
    //const E_Float& y = pts[i].vec[1];
    //const E_Float& z = pts[i].vec[2];

    get_colliding(pts[i], surface, cands, surface.index_start, ARTOL, false/*dummy*/);
    
    if (cands.empty()) continue; // regular
    
    //find sticking PGs
    for (size_t c = 0; c < cands.size(); ++c)
    {
      E_Int PGi = cands[c] - surface.index_start;
      auto PG = surface.element(PGi);

      PG.triangulate(dt, surface.crd);

      for (E_Int t = 0; t < PG.nb_tris(); ++t)
      {
        PG.triangle(t, T);

        const E_Float * P0 = surface.crd.col(T[0]);
        const E_Float * P1 = surface.crd.col(T[1]);
        const E_Float * P2 = surface.crd.col(T[2]);

        bool interfere, inside;
        E_Float d = K_MESH::Triangle::minDistanceToPoint(P0, P1, P2, pts[i].vec, UV, interfere, inside);
        if (!interfere) continue;

        // check PG is not too far
        E_Float TOLi = ARTOL;
        if (ARTOL < 0.) // relative
        {
          E_Float PGLref2 = PG.Lref2(surface.crd);// (surface.nodal_metric2);
          TOLi *= -sqrt(std::min(pts[i].val2, PGLref2));
        }

        if (d >= TOLi) continue;

        // sticking PG found
        pt_to_elt[i].push_back(PGi);
        break;
      }
    }
  }
}

///
template <typename mesh_t1, typename mesh_t2>
void compute(const mesh_t1& m1, const mesh_t2& m2,
  std::vector<bool>& is_x1, std::vector<bool>& is_x2, E_Float RTOL = ZERO_M,
  bool roughly_by_bbox = false)
{
  is_x1.resize(m1.ncells(), 0);
  is_x2.resize(m2.ncells(), 0);

  auto m2_loc = m2.get_localizer();
  if (m2_loc == nullptr) return;

  if (!roughly_by_bbox)
    m1.get_nodal_metric2();//for ae

  std::vector<E_Int> cands;

#ifndef COLLIDER_DBG
#pragma omp parallel for private (cands)
#endif
  for (E_Int i = 0; i < m1.ncells(); ++i)
  {
    cands.clear();

    if (roughly_by_bbox)
    {
      auto e1 = m1.element(i);

      m2_loc->get_candidates(e1, m1.crd, cands, 1, RTOL); //return as 0-based (fixme for volumic, was 1-based)
      if (cands.empty()) continue;
    }
    else
    {
      auto ae1 = m1.aelement(i);

      m2_loc->get_candidates(ae1, ae1.m_crd, cands, 1, RTOL); //return as 0-based (fixme for volumic, was 1-based)
      if (cands.empty()) continue;

#ifdef DEBUG_METRIC
      if (cands.size() > 50)
      {
        std::cout << "found more than 50 cands" << std::endl;
        medith::write("cands", fldm.crd, fldm.cnt, &cands, 1);
        medith::write("ae1", ae1);
      }
#endif

      bool is_x = get_colliding(ae1, m2, cands, 1, RTOL, true/*returns at first found*/);
      if (!is_x) continue;
    }

    is_x1[i] = is_x2[cands[0] - 1] = true;
  }
}

}   // COLLIDE
}   // NUGA

#endif /* NUGA_COLLIDER_HXX */


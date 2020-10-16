

/*
 
 
 
              NUGA 
 
 
 
 */
//Author : Sâm Landier (sam.landier@onera.fr)

#ifndef NUGA_COLLIDER_HXX
#define NUGA_COLLIDER_HXX

#include "Nuga/include/macros.h"
#include <vector>
#include "Fld/FldArray.h"
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
    template <typename crd_t, int DIM, typename ELT1, typename ELT2>
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
      
        E_Int nb_tris2 = e2.nb_tris();
        for (E_Int j=0; j < nb_tris2; ++j)
        {
          e2.triangle(j, T2);
          Q0 = crd2.col(T2[0]);
          Q1 = crd2.col(T2[1]);
          Q2 = crd2.col(T2[2]);
          
#ifdef COLLIDER_DBG

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
          {
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
                     double ps_min = 0.95/*overlap criterion*/,bool swap = true, const E_Float* norm2 = nullptr) //norm2 when right direction is known upon entry
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

    E_Float abstol = MAX(ZERO_M, RTOL*::sqrt(Lref));
    //std::cout << "compute_overlap for VOL : " << abstol << std::endl;
    for (j = 0; (j < cands2.size()); ++j)
    {
      E_Int J = cands2[j];

      const E_Int* nodes2 = PGs2.get_facets_ptr(J);
      E_Int nb_nodes2 = PGs2.stride(J);
      
      if (norm2 == nullptr)
        ELT2::template normal<crd_t, 3>(crd2, nodes2, nb_nodes2, 1, n2);
      
      double ps = NUGA::dot<3>(n1,n2);
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

    E_Float abstol = MAX(EPSILON, RTOL*::sqrt(Lref));
    E_Float u00, u01, u10, u11;
    E_Bool overlap = false;

    //std::cout << abstol << std::endl;
    
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
        double ps = NUGA::dot<3>(d1,d2);

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
bool get_colliding(const aelt_t& e1, const bound_mesh_t& mask_bit, std::vector<E_Int>& cands, E_Int idx_start, double RTOL, bool only_first_found);

/// Polygon vs polyLine impl.
template <> inline
bool get_colliding<NUGA::aPolygon, edge_mesh_t>
(const NUGA::aPolygon& ae1, const edge_mesh_t& mask_bit, std::vector<E_Int>& cands, E_Int idx_start, double RTOL, bool first_found)
{
 
  bool hasX(false);

  // reduce mask to candidates
  edge_mesh_t lmask(mask_bit, cands, idx_start);

  double normal1[3];
  ae1.normal<3>(normal1);
  const double * plane_pt = ae1.m_crd.col(0); // choosing first node

#ifdef DEBUG_XCELLN
  //medith::write("cutter_front_before_projection", lmask.crd, lmask.cnt);
#endif

  double l21 = ae1.Lref2();
  double l22 = lmask.Lref2();//compute it before proj
  double ATOL(RTOL * ::sqrt(std::min(l21, l22)));

  // project candidates on e1's plane => 2D problem
  STACK_ARRAY(double, lmask.crd.cols(), signed_dists);
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

    double z1 = signed_dists[e1];
    double z2 = signed_dists[e2];

    if (z1*z2 < 0.) continue; // means crossing 

    double mz = std::min(::fabs(z1), ::fabs(z2));
    keep[i] = (mz < ATOL);
  }

  K_CONNECT::keep<bool> pred(keep);
  K_CONNECT::IdTool::compress(cands, pred);

  if (cands.empty())
    return false;

  lmask.compress(keep);

  // close enough so go to 2D for real test

#ifdef DEBUG_XCELLN
  //medith::write("cutter_front_projected", lmask.crd, lmask.cnt);
#endif

  // compute collision between e1 and each candidate until founding one collision

  // a. triangulate e1, check if any lmask point falls into the triangulation
  DELAUNAY::Triangulator dt;
  ae1.triangulate(dt);

  // b. go 2D both crd (a copy of it)  and lmask.crd

  // Computes the fitting box coordinate system optimizing the view over the contour
  K_FLD::FloatArray P(3, 3), iP(3, 3);

  FittingBox::computeAFrame(normal1, P);
  iP = P;
  K_FLD::FloatArray::inverse3(iP);

  K_FLD::FloatArray crd2D(ae1.m_crd);
  FittingBox::transform(crd2D, iP);// Now we are in the fitting coordinate system.
  crd2D.resize(2, crd2D.cols()); // e1 is 2D now.
  FittingBox::transform(lmask.crd, iP);// Now we are in the fitting coordinate system.
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

  int T1[3];
  for (E_Int i = 0; (i < ae1.nb_tris()); ++i)
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
(const NUGA::aPolyhedron<UNKNOWN>& ae1, const pg_smesh_t& mask_bit, std::vector<E_Int>& cands, E_Int idx_start, double RTOL, bool first_found)
{
  bool hasX(false);

  // reduce mask to candidates : noy hpc way for above lineic (yet)
  //pg_smesh_t lmask(mask_bit, cands, idx_start);
  const pg_smesh_t& lmask = mask_bit;

#ifdef COLLIDER_DBG
  medith::write("lmask", lmask.crd, lmask.cnt);
#endif

  //double l21 = ae1.Lref2();
  //double l22 = lmask.Lref2(cands, idx_start);
  //double ATOL(RTOL * ::sqrt(std::min(l21, l22)));

  // compute collision between e1 and each candidate until founding one collision

  // a. triangulate e1, check if any lmask point falls into the triangulation
  //DELAUNAY::Triangulator dt;
  //ae1.triangulate(dt);
  ae1.cvx_triangulate(ae1.m_crd);
  

  // (b. projection for 2D)

  // c. detect collisions

  E_Int ncands(cands.size());
  std::vector<bool> keep(ncands, false);

  int T1[3];
  for (E_Int i = 0; (i < ae1.nb_tris()); ++i)
  {
    ae1.triangle(i, T1);
    //const E_Float* P1 = ae1.m_crd.col(T1[0]);
    //const E_Float* Q1 = ae1.m_crd.col(T1[1]);
    //const E_Float* R1 = ae1.m_crd.col(T1[2]);

    for (int j = 0; (j<ncands); ++j)
    {
      if (keep[j]) continue;

      K_MESH::Polygon e2 = lmask.element(cands[j] - idx_start);

      //do not pass the tol as we use here a predicate
      E_Int t1, t2;
      using crd_t = K_FLD::FloatArray;
      bool isx = simplicial_colliding<crd_t, 3>(ae1.m_crd, ae1, lmask.crd, e2, K_MESH::Triangle::fast_intersectT3<3>, -1., t1, t2);

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
(const vertex& pt, const pg_smesh_t& surface, std::vector<E_Int>& cands, E_Int idx_start, double RTOL, bool first_found_dummy)
{
  cands.clear();

  auto loc = surface.get_localizer();

  double Lref = ::sqrt(pt.val2);

  K_SEARCH::BBox3D bb1;
  bb1.minB[0] = pt.vec[0] - RTOL * Lref;
  bb1.minB[1] = pt.vec[1] - RTOL * Lref;
  bb1.minB[2] = pt.vec[2] - RTOL * Lref;
  bb1.maxB[0] = pt.vec[0] + RTOL * Lref;
  bb1.maxB[1] = pt.vec[1] + RTOL * Lref;
  bb1.maxB[2] = pt.vec[2] + RTOL * Lref;
  
  loc->get_candidates(bb1, cands, idx_start, -1./*dummy*/);

  return !cands.empty();
}
/// above wrapper for list of vertices : returns indir 'pt to face'
inline void get_colliding
(const std::vector<vertex>& pts, const pg_smesh_t& surface, double RTOL, std::vector<E_Int>& pt_to_elt)
{
  pt_to_elt.clear();

  size_t npts = pts.size();

  pt_to_elt.resize(npts, IDX_NONE);

  surface.get_nodal_metric2(); // to compute it if missing

  DELAUNAY::Triangulator dt;
  std::vector<E_Int> cands;
  E_Int T[3];
  E_Float UV[3];
  //
  for (size_t i = 0; i < npts; ++i)
  {
    get_colliding(pts[i], surface, cands, surface.index_start, RTOL, false/*dummy*/);
    
    if (cands.empty()) continue; // regular
    
    //find best PG
    double dmin{ NUGA::FLOAT_MAX };
    E_Int PGtarget{ IDX_NONE };
    //
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
        double d = K_MESH::Triangle::minDistanceToPoint(P0, P1, P2, pts[i].vec, UV, interfere, inside);
        if (!interfere) continue;

        if (d < dmin)
        {
          PGtarget = PGi;
          dmin = d;
          break;
        }
      }
    }

    if (PGtarget != IDX_NONE) //check is not too far
    {
      auto PG = surface.element(PGtarget);
      double PGLref2 = PG.Lref2(surface.nodal_metric2);
      double TOLi = RTOL * ::sqrt(std::min(pts[i].val2, PGLref2));
      if (dmin > TOLi)
        PGtarget = IDX_NONE;
      else
        PGtarget = PGtarget;
    }

    pt_to_elt[i] = PGtarget;
  }

#ifdef COLLIDER_DBG
  K_FLD::FloatArray crd = surface.crd;
  std::vector<bool> keep(surface.ncells(), false);
  for (size_t i = 0; i < pt_to_elt.size(); ++i)
  {
    if (pt_to_elt[i] != IDX_NONE)
    {
      crd.pushBack(pts[i].vec, pts[i].vec + 3);
      keep[pt_to_elt[i]] = true;
    }
  }
  K_CONNECT::keep<bool> pred(keep);
  std::vector<E_Int> oids, nids;

  ngon_unit pgs;
  surface.cnt.extract_by_predicate(pred, pgs, oids, nids);

  medith::write("singular", crd, pgs);
  
  ngon_type ng(pgs, true);
  K_FLD::IntArray cnto;
  ng.export_to_array(cnto);
  MIO::write("singular.tp", crd, cnto, "NGON");
  cnto.clear();
  MIO::write("nodes.tp", crd, cnto, "BAR");
  cnto.resize(2, 1, 0);
  cnto(1, 0) = 1;
  medith::write("nodes", crd, cnto, "BAR");
  //K_FLD::IntArray dum;
  //MIO::write("D:\\slandier\\DATA\\tmp\\immersion\\nodes.tp", crd, dum);

#endif
}
  
}   // COLLIDE
}   // NUGA

#endif /* NUGA_COLLIDER_HXX */


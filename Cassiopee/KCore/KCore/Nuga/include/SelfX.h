/*    
    Copyright 2013-2025 Onera.

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

#ifndef SELFX_H
#define SELFX_H


namespace NUGA
{

typedef K_FLD::FloatArray crd_t;
typedef K_FLD::IntArray cnt_t;

#if !defined(NETBEANSZ) && !defined(VISUAL)
inline bool getBoundary(const E_Int* t0, const E_Int* t1, E_Int& i, E_Int& j)
{
  for (i=0; i < 3; ++i)
  {
    const E_Int& t00=*(t0+i);
    const E_Int& t01=*(t0+(i+1)%3);
      
    for (j=0; j < 3; ++j)
    {
      const E_Int& t10=*(t1+j);
      const E_Int& t11=*(t1+(j+1)%3);
      
      if ((t00==t10) && (t01==t11))
        return true;
      if ((t00==t11) && (t01==t10))
        return true;
    }
  }
  
  i=j=IDX_NONE;
  return false;
}
#endif

///
template <short DIM>
bool __fast_discard
(const K_FLD::FloatArray& pos, const E_Int* t0, const E_Int* t1, E_Float tol)
{
  const E_Int& T0 = *t0;
  const E_Int& T1 = *(t0+1);
  const E_Int& T2 = *(t0+2);
  const E_Int& p0 = *t1;
  const E_Int& p1 = *(t1+1);
  const E_Int& p2 = *(t1+2);
  const E_Float* Q0 = pos.col(T0);
  const E_Float* Q1 = pos.col(T1);
  const E_Float* Q2 = pos.col(T2);
  const E_Float* P0 = pos.col(p0);
  const E_Float* P1 = pos.col(p1);
  const E_Float* P2 = pos.col(p2);
  
  // 1. Check if points are on the same side of the t's plane
  E_Float _U1[DIM],_U2[DIM], _U3[DIM];
  NUGA::diff<DIM>(P1, P0, _U1);
  NUGA::diff<DIM>(P2, P0, _U2);
  NUGA::crossProduct<DIM>(_U1,_U2,_U3);
  NUGA::normalize<DIM>(_U3);
        
  bool is_far[] = {false,false, false};
  E_Float h0(0.), h1(0.), h2(0.);
  NUGA::diff<DIM>(Q0, P0, _U1);
  h0 = NUGA::dot<DIM>(_U3, _U1);
  is_far[0] = (h0 >= tol) || (h0 <= -tol);
    
  NUGA::diff<DIM>(Q1, P0, _U1);
  h1 = NUGA::dot<DIM>(_U3, _U1);
  is_far[1] = (h1 >= tol) || (h1 <= -tol);
    
  NUGA::diff<DIM>(Q2, P0, _U1);
  h2 = NUGA::dot<DIM>(_U3, _U1);
  is_far[2] = (h2 >= tol) || (h2 <= -tol);
    
  E_Int s[3];
  s[0]=SIGN(h0);
  s[1]=SIGN(h1);
  s[2]=SIGN(h2);
    
  if (is_far[0] && is_far[1] && is_far[2] && (s[0] == s[1]) && (s[0]==s[2]))
    return true;
  
  bool overlapping = (!is_far[0] && !is_far[1] && !is_far[2]);
  E_Int shn=IDX_NONE;//shared node's rank

  if ( (T0 == p0) || (T0 == p1) || (T0 == p2) ) shn=0;
  if ( (T1 == p0) || (T1 == p1) || (T1 == p2) ) shn=1;
  if ( (T2 == p0) || (T2 == p1) || (T2 == p2) ) shn=2;
    
  // 2. Check if they are not overlapping but they share a point and the 2 others are far and in the same side => discard
  
  if (!overlapping) //nor coplanar
  {
    if (shn != IDX_NONE)
    {
      if ( (s[(shn+1)%3] == s[(shn+2)%3]) && (is_far[(shn+1)%3] && is_far[(shn+2)%3]) )
        return true;  
    
      // 3. Check if they are not overlapping but they share an edge => discard
      E_Int i,j;    
      if (getBoundary(t0,t1, i,j)) //true if an edge is shared and therefore i,j are valued
        return true;
    }
  }
  else //overlapping or just coplanar : COMMENTED FOR NGON BOOLEAN DOINT ONE PASS ONLY (very few case like that occur at first iter)
  {
    //if sharing an edge and opposite nodes are on each side => just coplanar
    if (shn != IDX_NONE) //they must share at least a node
    {
      E_Int i,j;    
      if (getBoundary(t0,t1, i,j)) //true if an edge is shared and therefore i,j are valued
      {
        E_Float shE[3], U[3], n0[3], n1[3];
        E_Int n0op = *(t0+(i+2)%3);
        E_Int n1op = *(t1+(j+2)%3);
        
        NUGA::diff<3>(pos.col(*(t0+(i+1)%3)), pos.col(*(t0+i)), shE);
        
        NUGA::diff<3>(pos.col(n0op), pos.col(*(t0+i)), U);
        NUGA::crossProduct<3>(shE,U,n0);
        
        NUGA::diff<3>(pos.col(n1op), pos.col(*(t0+i)), U);
        NUGA::crossProduct<3>(shE,U,n1);
        
        if (NUGA::dot<3>(n0, n1) < 0.)
          return true;
      }
    }   
  }
  
  return false;
}

///
template <short DIM>
E_Bool __intersect
(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect, E_Int t, E_Int e0, E_Int e1, E_Float tol, E_Bool& coplanar)
{
  
  E_Float u0[2], E[DIM], eps(/*100.*EPSILON*/tol);
  E_Bool  overlap, intersect;
  E_Int tx[2];  
  K_FLD::IntArray::const_iterator pS = connect.col(t);

  if (*pS == e0 && *(pS+1) == e1)
    return false;
  if (*pS == e1 && *(pS+1) == e0)
    return false;

  if (*(pS+1) == e0 && *(pS+2) == e1)
    return false;
  if (*(pS+1) == e1 && *(pS+2) == e0)
    return false;

  if (*pS == e0 && *(pS+2) == e1)
    return false;
  if (*pS == e1 && *(pS+2) == e0)
     return false;
  
  intersect = K_MESH::Triangle::intersect<3>
    (pos, *pS, *(pS+1), *(pS+2), e0, e1, tol, true, u0[0], u0[1], tx, overlap, coplanar);
  
  if (!intersect)
    return false;

  bool one_single_x_point = (u0[1] == NUGA::FLOAT_MAX) || (::fabs(u0[1] - u0[0]) < eps);
  bool share_a_node = (*pS == e0 || *(pS+1) == e0 || *(pS+2) == e0) || (*pS == e1 || *(pS+1) == e1 || *(pS+2) == e1);
  bool one_inside = ((u0[0] > eps) && (u0[0] < 1.-eps));
  
  if (share_a_node && one_single_x_point && !one_inside) //not a real intersection
    return false;
  
  bool add2E, add2T;
  intersect=false;

  NUGA::diff<DIM>(pos.col(e1), pos.col(e0), E);

  for (E_Int n = 0; n < 2; ++n)
  {
    if (u0[n] == NUGA::FLOAT_MAX)
      continue;
    
    add2E = (u0[n] >= eps) && (u0[n] <= (1. - eps));
    add2T = (tx[n] == 0);
          
    if (add2E) // Create new point Ni
      intersect = 1;
    else if (add2T)
      intersect = 3;
  }
  
  if (overlap)
    intersect=2;

  return intersect;
}

///
template <short DIM>
E_Int intersect_from_Conformizer_wo_trace
(const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect, E_Int t1, E_Int t2, E_Float tol)
{
  E_Int  i, r(0);
  E_Int ret(0);
  
//  if (t1.Si == 19494 || t2.Si == 19494)
//  {
//      std::cout << "caught" << std::endl;
//  }
  
//  //////////////////////////////////////////////////////////////////////////////////
//  const E_Float* P1 = pos.col(connect(0,t1.Si));
//  const E_Float* Q1 = pos.col(connect(1,t1.Si)); 
//  const E_Float* R1 = pos.col(connect(2,t1.Si));
//  const E_Float* P2 = pos.col(connect(0,t2.Si));
//  const E_Float* Q2 = pos.col(connect(1,t2.Si));
//  const E_Float* R2 = pos.col(connect(2,t2.Si));
//  
//  if (!K_MESH::Triangle::fast_intersectT3<3>(P1, Q1, R1, P2, Q2, R2))
//  {
//#ifdef DEBUG_TRI_CONFORMIZER
//    ++K_GENERATOR::ConformizerRoot::fastdiscard_counter;
//#endif
//    return false;
//  }
//      
//  //////////////////////////////////////////////////////////////////////////////////
  
  // Fast return : Triangles far from each others  
  if (__fast_discard<DIM>(pos, connect.col(t2), connect.col(t1), tol))
    return 0;
    
  if (__fast_discard<DIM>(pos, connect.col(t1), connect.col(t2), tol))
    return 0;

  //bool insideT1(false), insideT2(false);
  E_Bool coplanarE;
  E_Int nb_coplanarE(0);
  bool are_coplanar=false;
  
  // T2->T1
  for (i = 0; i < 3; ++i)
  {
    const E_Int& e0=connect(i, t2);
    const E_Int& e1=connect((i+1)%3, t2);
    r = __intersect<3>(pos, connect, t1, e0, e1, tol, coplanarE);
    nb_coplanarE = (r==2) ? nb_coplanarE+1 : nb_coplanarE;
    //insideT1 |= (r==3);
    ret |=r;
  }
  
  if (nb_coplanarE > 1)
    are_coplanar=true; // overlap might require a second pass as each triangle split is done separtely so triangulation can be different
  
  nb_coplanarE=0;
  for (i = 0; i < 3; ++i)
  {
    const E_Int& e0=connect(i, t2);
    const E_Int& e1=connect((i+1)%3, t2);
    r = __intersect<3>(pos, connect, t2, e0, e1, tol, coplanarE);
    nb_coplanarE = (r==2) ? nb_coplanarE+1 : nb_coplanarE;
    ret |=r;
    //insideT2 |= (r==3);
  }
  
  if (nb_coplanarE > 1)
    are_coplanar=true; // overlap might require a second pass as each triangle split is done separtely so triangulation can be different
  
  if (are_coplanar)
    ret=2;

  return ret;
}

void concatenate_PG_triangles (E_Int PGi, const ngon_t<cnt_t>& ng, const cnt_t& cntT3, const std::vector<E_Int>& PG_to_T3s, std::vector<E_Int>& pgi_T3s)
{  
  if (PGi+1 > (E_Int)PG_to_T3s.size()-1)
    return;
  
  E_Int nb_T3=PG_to_T3s[PGi+1]-PG_to_T3s[PGi];
  if (nb_T3 <= 0) return;
    
  for (E_Int k=0; k<nb_T3; ++k)
    pgi_T3s.push_back(PG_to_T3s[PGi]+k);
}

void concatenate_PH_triangles (E_Int PHi, const ngon_t<cnt_t>& ng, const cnt_t& cntT3, const std::vector<E_Int>& PG_to_T3s, std::vector<E_Int>& phi_T3s)
{
  phi_T3s.clear();
  
  size_t nb_pgs = ng.PHs.stride(PHi);
  for (size_t j=0; j < nb_pgs; ++j)
  {
    E_Int PGi = ng.PHs.get_facet(PHi, j) -1;
    concatenate_PG_triangles(PGi, ng, cntT3, PG_to_T3s, phi_T3s); //appended
  }
}

template <typename TriangulatorType>
void selfX(const K_FLD::FloatArray& crd, const ngon_t<cnt_t>& ng, std::vector<E_Int>& xlist)
{
  
  xlist.clear();
  
  std::vector<E_Int> T3_to_PG, PG_to_T3s;
  cnt_t cntT3;
  
  std::cout << "selfX : Triangulating ..." << std::endl;
  // Triangulate once all the PGs
  E_Int err = ngon_t<cnt_t>::triangulate_pgs<TriangulatorType>(ng.PGs, crd, cntT3, T3_to_PG, true, true);
  if (err)
  {
    std::cout << "selX ERROR : triangulation failure." << std::endl;
    return;
  }
  
  // Inverse indirection
  PG_to_T3s.reserve(ng.PGs.size()+1);
  int k=-1;
  for (size_t i=0; i < T3_to_PG.size(); ++i)
  {
    if (k != T3_to_PG[i])
    {
      k=T3_to_PG[i];
      PG_to_T3s.push_back(i);
    }
  }
  PG_to_T3s.push_back(T3_to_PG.size());
  
  // Construct the BbTree on the PHs
  std::cout << "selfX : Creating BBTree ..." << std::endl;
  E_Int nb_elts = ng.PHs.size();
  K_SEARCH::BoundingBox<3>* pool = new K_SEARCH::BoundingBox<3>[nb_elts];
  std::vector<K_SEARCH::BoundingBox<3>*> boxes(nb_elts);
  std::vector<E_Int> nodes;
  
  for (E_Int i = 0; i < nb_elts; ++i)
  {
    ng.nodes_ph(i, nodes, true);
    boxes[i] = &pool[i];
    boxes[i]->compute(crd, nodes);
  }
  K_SEARCH::BbTree<3> tree(boxes);
  
  std::set<K_MESH::NO_Edge> processed_PH_pairs, processed_PG_pairs;
  std::vector<E_Int> pgi_T3s, pgj_T3s, candidates;
  
  //std::cout << "selfX : run ..." << std::endl;
  for (E_Int i = 0; i < nb_elts; ++i)
  {  
    if (i > 10000 && i%10000 == 0)
      std::cout << i << " -th element processed over " << nb_elts << std::endl;
    //    
    candidates.clear();
    tree.getOverlappingBoxes(boxes[i]->minB, boxes[i]->maxB, candidates);
    if (candidates.empty())
      continue;
    
    E_Int nb_PGis = ng.PHs.stride(i);
    const E_Int* pPGi = ng.PHs.get_facets_ptr(i);

    // Get all the triangles for this PHi
    //concatenate_triangles(i, ng, cntT3, PG_to_T3s, phi_T3s);
    
    bool is_x=false;
    
    for (size_t j=0; j< candidates.size(); ++j)
    {
      E_Int& jj = candidates[j];
      
      if (!processed_PH_pairs.insert(K_MESH::NO_Edge(i, jj)).second) // if already in
            continue;
      
      E_Int nb_PGjs = ng.PHs.stride(jj);
      const E_Int* pPGj = ng.PHs.get_facets_ptr(jj);
      
      for (E_Int pgi=0; (pgi<nb_PGis) && !is_x; ++pgi)
      {
        E_Int PGi = *(pPGi+pgi)-1;
//        
        pgi_T3s.clear();
        concatenate_PG_triangles(PGi, ng, cntT3, PG_to_T3s, pgi_T3s);
                
        for (E_Int pgj=0; (pgj<nb_PGjs) && !is_x; ++pgj)
        {
          E_Int PGj = *(pPGj+pgj)-1;
          //if (!processed_PG_pairs.insert(K_MESH::NO_Edge(PGi, PGj)).second) // if already in
            //continue;
          
          pgj_T3s.clear();
          concatenate_PG_triangles(PGj, ng, cntT3, PG_to_T3s, pgj_T3s);
          
          // Now do the check T3/T3
          for (size_t i1=0; (i1 < pgi_T3s.size()) && !is_x; ++i1)
          {
            for (size_t i2=0; (i2 < pgj_T3s.size()) && !is_x; ++i2)
            {
              E_Int& I=pgi_T3s[i1];
              E_Int& J=pgj_T3s[i2];
          
              const E_Float* P1 = crd.col(cntT3(0,I));
              const E_Float* P2 = crd.col(cntT3(1,I));
              const E_Float* P3 = crd.col(cntT3(2,I));
          
              const E_Float* Q1 = crd.col(cntT3(0,J));
              const E_Float* Q2 = crd.col(cntT3(1,J));
              const E_Float* Q3 = crd.col(cntT3(2,J));
          
              is_x = K_MESH::Triangle::fast_intersectT3<3>(P1, P2, P3, Q1, Q2, Q3, EPSILON/*dummy tol*/);
              if (is_x)
                is_x = intersect_from_Conformizer_wo_trace<3>(crd, cntT3, I, J, EPSILON);
#ifdef DEBUG_SELFX
              if (is_x)
              {
                std::cout << "intersecting pair I/J : " << I << "/" << J << std::endl;
                K_FLD::IntArray tmp;
                tmp.pushBack(cntT3.col(I), cntT3.col(I)+3);
                tmp.pushBack(cntT3.col(J), cntT3.col(J)+3);
                medith::write("xpair.mesh", crd, tmp, "TRI");
                
              }
#endif
            }
          }
        }
      }
      
      if (is_x){
        xlist.push_back(i); xlist.push_back(jj);
        break;
      }
    }
    if (is_x) break;
  }

  std::cout << "selfX : DONE." ;
  if (xlist.empty())
    std::cout << " No self-intersections found." << std::endl;
  else
  {
    std::cout << " WARNING : intersections found between " << xlist[0] << " and " << xlist[1] << std::endl;
  }
    
  delete [] pool;
}

typedef E_Int(*xfunc)(const K_FLD::FloatArray& crd, const K_FLD::IntArray &cT3, E_Int i, E_Int j, E_Float TOL);

static bool T3intersect(xfunc func, const K_FLD::FloatArray& crd, const K_FLD::IntArray &cT3, E_Float TOL, const E_Int* first_t31, E_Int nb_t31, const E_Int* first_t32, E_Int nb_t32)
{
  //
  for (E_Int i = 0; i < nb_t31; ++i)
  {
    for (E_Int j = 0; j < nb_t32; ++j)
    {
      if (func(crd, cT3, *(first_t31 + i), *(first_t32 + j), TOL) != 0)
        return true;
    }
  }
  return false;
}

template<typename TriangulatorType>
static E_Int PGintersect( xfunc func,
  const K_FLD::FloatArray& crd, const ngon_unit& PGS,
  const E_Int* first_pg1, E_Int nb_pgs1, const E_Int* first_pg2, E_Int nb_pgs2,
  K_FLD::IntArray& cT3, E_Int* xT31, E_Int* nb_T31, E_Int* xT32, E_Int* nb_T32, 
  bool &is_x)
{
  E_Float TOL = EPSILON;

  is_x = false;

  TriangulatorType dt;

  typedef K_SEARCH::BBox3D Bx3D;
  typedef K_SEARCH::BbTree<3, Bx3D> Tree3D;

  bool self_test = (first_pg1 == first_pg2 && nb_pgs1 == nb_pgs2);
  K_SEARCH::BBox3D box1, box2;

  for (E_Int i = 0; i < nb_pgs1; ++i)
  {
    E_Int PG1 = *(first_pg1 + i) - 1;
    const E_Int* nodes1 = PGS.get_facets_ptr(PG1);
    E_Int nb_nodes1 = PGS.stride(PG1);

    box1.compute(crd, nodes1, nb_nodes1, 1/*index_start*/);

    E_Int j = self_test ? i + 1 : 0;
    for (; j < nb_pgs2; ++j)
    {
      E_Int PG2 = *(first_pg2 + j) - 1;
      const E_Int* nodes2 = PGS.get_facets_ptr(PG2);
      E_Int nb_nodes2 = PGS.stride(PG2);

      box2.compute(crd, nodes2, nb_nodes2, 1/*index_start*/);

      // box test
      if (!Tree3D::boxesAreOverlapping(&box1, &box2, TOL))
        continue;

      // need a full test on triangulations
      
      // is PG1 triangulated ?
      if (xT31[i] == IDX_NONE)
      {
        xT31[i] = cT3.cols();
        E_Int err = K_MESH::Polygon::triangulate(dt, crd, nodes1, nb_nodes1, 1, cT3);
        if (err) return 1;
        nb_T31[i] = cT3.cols() - xT31[i];
      }
      // is PG2 triangulated ?
      if (xT32[j] == IDX_NONE)
      {
        xT32[j] = cT3.cols();
        E_Int err = K_MESH::Polygon::triangulate(dt, crd, nodes2, nb_nodes2, 1, cT3);
        if (err) return 1;
        nb_T32[j] = cT3.cols() - xT32[j];
      }

      if (T3intersect(func, crd, cT3, TOL, &xT31[i], nb_T31[i], &xT32[j], nb_T32[j]))
      {
        is_x = true;
        return 0;
      }
    }
  }

  return 0;

}

}

#endif

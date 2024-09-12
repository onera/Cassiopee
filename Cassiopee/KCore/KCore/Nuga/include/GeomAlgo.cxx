/*    
    Copyright 2013-2024 Onera.

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

#ifndef _GEOMALGO_CXX_
#define _GEOMALGO_CXX_

#include "Nuga/include/GeomAlgo.h"
#include "Nuga/include/EltAlgo.h"
#include "Nuga/include/BbTree.h"
#include <assert.h>
#include "Nuga/include/MeshTool.h"
#include "Nuga/include/ArrayWriter.h"
#include "Nuga/include/macros.h"

#ifdef DEBUG_GEOM_ALGO
#include "IO/io.h"
#endif

#define MIN2(a,b) ((a<b) ? a : b)
#define MIN3(a,b,c) ((a<b) ? ((a<c) ? 0 : 2) : ((b<c) ? 1 : 2))

///
template <short DIM>
static void iso_barycenter(const K_FLD::FloatArray& coords, E_Float* isob)
{
  for (size_t i = 0; i < DIM; ++i) isob[i]=0.;
  
  if (coords.cols()==0)
    return;
  size_t sz = coords.getSize();
  for (size_t i = 0; i < sz; ++i)
   for (size_t j = 0; j < DIM; ++j)
     isob[j] += coords(j,i);
  
  E_Float k = 1./(E_Float)coords.cols();
  
  for (size_t i = 0; i < DIM; ++i)
    isob[i] *= k;
}

///
template <typename ElementType>
template <typename Coordinate_t, typename Connectivity_t, short DIM>
void
NUGA::GeomAlgo<ElementType>::neighboring
(const K_FLD::ArrayAccessor<Coordinate_t>& coords, const K_FLD::ArrayAccessor<Connectivity_t>& conn,
 Connectivity_t& neighbors)
{
  size_t ELTS(conn.size()), sz, count(0), nbounds;
  E_Float point[DIM], isob[DIM];
  K_FLD::FloatArray crds;;
  ElementType E;
  typename ElementType::boundary_type B;
  K_FLD::FloatArray barys;
  barys.reserve(DIM, ELTS);
  Vector_t<E_Int> glob_face_id, xglob, facet2elt, new_IDs;
  
  neighbors.clear();

  xglob.push_back(0);

  for (size_t K = 0; K < ELTS; ++K)
  {
    E.set(conn, K);
    nbounds = E.nbounds();
    for (size_t b=0; b < nbounds; ++b)
    {
      E.getBoundary(b, B);
      sz = B.nnodes();
      
      crds.clear();
      crds.reserve(DIM, sz);

      for (size_t n = 0; n < sz; ++n)
      {
        coords.getEntry(B.node(n), point);
        crds.pushBack(point, point+DIM);
      }
      
      iso_barycenter<DIM>(crds, isob);
      barys.pushBack(isob, isob+DIM);
      
      glob_face_id.push_back(count);
      facet2elt.push_back(K);
    } // bound
    xglob.push_back(glob_face_id.size());
  } // Elements

  // Merge
  K_FLD::ArrayAccessor<K_FLD::FloatArray> cab(barys);
  E_Int nmerges = ::merge(cab, EPSILON, new_IDs);

  // Build Neighbour matrix
  std::vector<std::vector<E_Int> > ngh_tmp(ELTS);

  bool fixed_stride=true;
  E_Int stride = 0;
  for (size_t K = 0; K < ELTS; ++K)
  {
    const E_Int & start = xglob[K];
    const E_Int & end = xglob[K+1];
    
    if (stride == 0)
      stride = end-start;
    else if (end-start != stride)
      fixed_stride=false;
    
    for (size_t i = start; i < end; ++i)
    {
      const E_Int & nid = new_IDs[i];
      if (nid== i)
        continue;
      const E_Int& Kn = facet2elt[nid];
      ngh_tmp[K].push_back(Kn);
      ngh_tmp[Kn].push_back(K);
    }
  }
  
  if (fixed_stride)
  {
    K_FLD::ArrayWriter<Connectivity_t> nghA(neighbors);
    nghA.clear();
    nghA.reserve(stride, ELTS);
    
    E_Int* entry(new E_Int(stride));
    E_Int s;
    for (size_t i = 0; i < ngh_tmp.size(); ++i)
    {
      s = ngh_tmp[i].size();
      for (size_t j=0; j < s; ++j)
        entry[j]=ngh_tmp[i][j];
      for (size_t j=s; j < stride; ++j) //fills the tail with nonde index.
        entry[j]=IDX_NONE;
      nghA.push_back(entry, stride);
    }
    delete entry;
  }
  else
  {
    //fixme : PH and PGs
  }

}


inline bool is_a_ref(const K_FLD::FloatArray& crd, const K_FLD::IntArray& cnt, E_Int ith, const K_SEARCH::BbTree3D& tree, bool& outward)
{
  // Check if ith elt can be used as a ref, and how it is oriented
    
  outward=true;
  E_Float G[3], P[3], N[3];
  
  
  //Build a point almost lying on the centroid but just above (EPS distance)
  K_MESH::Triangle T(cnt.col(ith));
      
  T.isoG(crd, T.nodes(), G);
  T.normal(crd, T.nodes(), N);
      
  NUGA::sum<3>(G,N,P); // P is above G at a distance of 1
            
  Vector_t<E_Int> boxes;
  tree.getIntersectingBoxes(P, G, boxes);
  
  if (boxes.empty()) // can this happen ?
    return false;
  
  size_t nb_cands(boxes.size());
  
#ifdef DEBUG_GEOM_ALGO
  /*{
    K_FLD::FloatArray crd1(crd);
    crd1.pushBack(G,G+3);
    crd1.pushBack(P,P+3);
    K_FLD::IntArray cnt1(2,1,0); cnt1(0,0)=crd1.cols()-2;cnt1(1,0)=crd1.cols()-1;
    medith::write("ray.mesh", crd1, cnt1, "BAR");
    {
      K_FLD::IntArray cnt2;
      for (size_t i=0; i < nb_cands; ++i)
      {
        cnt2.pushBack(cnt.col(boxes[i]), cnt.col(boxes[i])+3);
      }
      
      medith::write("cands.mesh", crd1, cnt2, "TRI");
    }
  }*/
#endif
        
  bool first=true;
  E_Int s(0);
  for (size_t b=0; b < nb_cands; ++b)
  {
    if (boxes[b] == ith) continue;
          
    K_MESH::Triangle Ti(cnt.col(boxes[b]));
          
    const E_Float *Q0, *Q1, *Q2;
    Q0 = crd.col(Ti.node(0));
    Q1 = crd.col(Ti.node(1));
    Q2= crd.col(Ti.node(2));
          
    E_Bool parallel, coincident;
    E_Float lambda, min_d, UV[2];
    K_MESH::Triangle::planeLineMinDistance<3>(Q0, Q1, Q2, G, P, EPSILON, true/*tol_is_absolute*/, lambda, UV, parallel, coincident, min_d);
    
    bool intersect= (((1. - UV[0] -  UV[1]) >= -EPSILON)  && (UV[0] >= -EPSILON) && (UV[1] >= -EPSILON));
    if (!intersect)
      continue;
    
    if (first)
    {
      s=(lambda < 0) ? +1 : -1;
      first=false;
      continue;
    }
    
    if (s != ((lambda < 0) ? +1 : -1))
      return false; // cannot rely on this one for the reference
  }
  
  if (s==0)
    return false;
  
  outward = (s==1);
  return true;
}

///
namespace NUGA
{
template <>
inline void GeomAlgo<K_MESH::Triangle>::reversi_chimera_skin
(const Vector_t<K_FLD::FloatArray*>& crds, const Vector_t<K_FLD::IntArray*>& cnts, bool outward)
{
  if (crds.empty() || cnts.empty())
     return;
   
   if (crds.size() != cnts.size())
     return;
   
   size_t nb_zones = crds.size();
   
   //0. Clean by removing duplicates in order to have a correct answer when calling getManifoldNeighbor below
   for (size_t i=1; i < nb_zones; ++i)
   {
     Vector_t<E_Int> dupIds;
     NUGA::MeshTool::removeDuplicated(*cnts[i], dupIds, false);
   }
   
   //1. Gather all in one container : SUPPOSED TO ENCLOSE A VOLUME (CHIMERA SURFACE)
   K_FLD::FloatArray crd(*crds[0]);
   K_FLD::IntArray cnt(*cnts[0]);
   
   std::vector<E_Int> shiftcnt(nb_zones); // to pass from local index to surface index
   shiftcnt[0]=0;
   
   for (size_t i=1; i < nb_zones; ++i)
   {
     E_Int shiftcrd = crd.cols();
     shiftcnt[i]=cnt.cols();
     
     cnt.pushBack(*cnts[i]);
     cnt.shift(shiftcrd, shiftcnt[i]);
     crd.pushBack(*crds[i]);
   }
   
   // 2. Build the BbTree for finding references
   E_Int nb_elts(cnt.cols());
   K_SEARCH::BoundingBox<3>* pool = new K_SEARCH::BoundingBox<3>[nb_elts];
   Vector_t<K_SEARCH::BoundingBox<3>*> boxes(nb_elts);
   for (E_Int i = 0; i < nb_elts; ++i)
   {
     boxes[i] = &pool[i];
     boxes[i]->compute(crd, cnt.col(i), 3/*NB NODES*/);
   }
   K_SEARCH::BbTree3D tree(boxes);
   
   //3. Reorient loop
   Vector_t<E_Int> Krefs, colors, orient;
   K_FLD::IntArray neighbors;
  
#ifndef NETBEANSZ
#pragma omp parallel for private(neighbors, orient, colors)
#endif
   for (size_t i=0; i< nb_zones; ++i)
   {
     //A. Get connex bits, so
     //a. call getManifoldNeighbor
     neighbors.clear();
     NUGA::EltAlgo<K_MESH::Triangle>::getManifoldNeighbours(*cnts[i], neighbors, false);
     //b. call coloring
     NUGA::EltAlgo<K_MESH::Triangle>::coloring (neighbors, colors);
     size_t nb_colors = 1 + *std::max_element(colors.begin(), colors.end());
     
     //B. Get a ref orientation for each bit
     orient.clear();
     orient.resize(cnts[i]->cols(), 1);
     
     Krefs.clear();
     Krefs.resize(nb_colors, IDX_NONE);
     size_t nbref=0;
     bool otwrd;
     //
     for (E_Int K = 0; (K < cnts[i]->cols()) && (nbref != nb_colors); ++K)
     {
       if (Krefs[colors[K]] != IDX_NONE)
         continue;
       
       otwrd = true;
       if(is_a_ref(crd, cnt, K+shiftcnt[i], tree, otwrd))
       {
         ++nbref;
         Krefs[colors[K]]=K;
         if (outward != otwrd)
           orient[K]=-1;
       }
     }
     
     //C. Compute correct orientations
     for (size_t c = 0; c < nb_colors; ++c)
       EltAlgo<K_MESH::Triangle>::reversi_connex(*cnts[i], neighbors, Krefs[c], orient);
     
     //D. Apply reorientation 
     //permut the third and second nodes
      for (size_t k = 0; k < orient.size(); ++k)
        if (orient[k] == -1)
          std::swap((*cnts[i])(1,k), (*cnts[i])(2,k));
   }
   
   delete[] pool;
}

///
template <>
inline void GeomAlgo<K_MESH::Triangle>::get_swapE
 (const K_FLD::FloatArray& pos, const K_FLD::IntArray& connect, const K_FLD::IntArray& neighbors, 
  const std::set<K_MESH::NO_Edge>& hE, E_Float tol, std::vector<std::pair<E_Int, E_Int> >& swapE)
 {
  //
  E_Float dj,d2[3], lambda, q11, q12, q21, q22;
  const E_Int *pS;
  E_Int ROWS(pos.rows()), Sn, bn;
  size_t j;
  K_FLD::IntArray::const_iterator pSn;
  E_Int Ni, Nj, Nk, Nl;
  
  std::vector<bool> frozen(connect.cols(), false); //one swapping edge per element (otherwise swapE entries get inconsistent)
      
  for (E_Int i=0; i < connect.cols(); ++i)
  {
    pS = connect.col(i);
    
    if (frozen[i]) continue;
    
    Ni = *(pS);
    Nj = *(pS+1);
    Nl = *(pS+2);
      
    d2[0]=d2[1]=d2[2]=NUGA::FLOAT_MAX;
      
    q11 = (ROWS == 3) ? K_MESH::Triangle::qualityG<3>(pos.col(Ni), pos.col(Nj), pos.col(Nl)) : K_MESH::Triangle::qualityG<2>(pos.col(Ni), pos.col(Nj), pos.col(Nl));
    if (q11 >= tol)
      continue;
      
    for (j=0; j < 3; ++j)
    {
      Sn = neighbors(j, i);
      if (Sn == IDX_NONE)
        continue;
      
      if (frozen[Sn]) continue;
      
      Ni = *(pS+j);
      Nj = *(pS+(j+1)%3);
      Nl = *(pS+(j+2)%3);
      
      if (hE.find(K_MESH::NO_Edge(Nj,Nl)) != hE.end()) //prevent to swap hard edges
        continue;
      
      dj = (ROWS == 3) ? K_MESH::Edge::edgePointMinDistance2<3>(pos.col(Nj), pos.col(Nl), pos.col(Ni), lambda) : K_MESH::Edge::edgePointMinDistance2<2>(pos.col(Nj), pos.col(Nl), pos.col(Ni), lambda);

      if ((lambda < 0.) || (lambda > 1.))
        continue;
      
//      E_Float L2 = NUGA::sqrDistance(pos.col(*(pS+(j+1)%3)), pos.col(*(pS+(j+2)%3)), 3);
//      if ((L2*lambda*lambda < _tolerance*_tolerance) || (L2 * (1.-lambda) * (1. - lambda) < _tolerance*_tolerance))
//        std::cout << "HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH" << std::endl;
    
      d2[j]=dj;
    }
    
    j = MIN3(d2[0], d2[1], d2[2]);
    
    if (d2[j] == NUGA::FLOAT_MAX)
      continue;
    
    Sn = neighbors(j, i);  
    pSn = connect.col(Sn);
    bn = K_MESH::Triangle::getOppLocalNodeId(i, j, connect, neighbors);

    Ni = *(pS+j);
    Nj = *(pS+(j+1)%3);
    Nl = *(pS+(j+2)%3);
    Nk = *(pSn + bn);
    
    q12 = (ROWS == 3) ? K_MESH::Triangle::qualityG<3>(pos.col(Nj), pos.col(Nk), pos.col(Nl)) : K_MESH::Triangle::qualityG<2>(pos.col(Nj), pos.col(Nk), pos.col(Nl));
    q21 = (ROWS == 3) ? K_MESH::Triangle::qualityG<3>(pos.col(Ni), pos.col(Nj), pos.col(Nk)) : K_MESH::Triangle::qualityG<2>(pos.col(Ni), pos.col(Nj), pos.col(Nk));
    q22 = (ROWS == 3) ? K_MESH::Triangle::qualityG<3>(pos.col(Ni), pos.col(Nk), pos.col(Nl)) : K_MESH::Triangle::qualityG<2>(pos.col(Ni), pos.col(Nk), pos.col(Nl));

    double qmin0 = MIN2(q11,q12); 
    double rmin = (qmin0 == 0.) ? NUGA::FLOAT_MAX : MIN2(q21, q22) / qmin0;

    if (rmin < 1.) //worst cannot get even worst
      continue;

    double qmax0 = MAX(q11,q12); 
    double rmax = (qmax0 == 0.) ? NUGA::FLOAT_MAX : MAX(q21, q22) / qmax0;

    // We want the quality improvement for the total
    // But if the worst improves much than the best deteriorates, it s OK
    if ((q21+q22 < q11+q12) && rmin*rmax < 1.)
      continue; 

    frozen[i]=frozen[Sn]=true;
    
    swapE.push_back(std::make_pair(i,j)); 
  } 
 }

///
template <>
inline void GeomAlgo<K_MESH::Triangle>::min_quality
(const K_FLD::FloatArray& crd, const K_FLD::IntArray& cnt, E_Float& minq, E_Int& imin)
{
  minq = NUGA::FLOAT_MAX;
  imin = IDX_NONE;
  //
  for (E_Int i = 0; i < cnt.cols(); ++i)
  {
    K_FLD::IntArray::const_iterator pS = cnt.col(i);

    E_Int Ni = *pS;
    E_Int Nj = *(pS + 1);
    E_Int Nk = *(pS + 2);

    E_Float q = K_MESH::Triangle::qualityG<3>(crd.col(Ni), crd.col(Nj), crd.col(Nk));

    if (q < minq)
    {
      minq = q;
      imin = i;
    }
  }
}

}
#endif

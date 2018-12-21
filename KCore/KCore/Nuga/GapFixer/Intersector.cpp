/*    
    Copyright 2013-2019 Onera.

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
//Author : Sâm Landier (sam.landier@onera.fr)

#include "Intersector.h"
#include "Connect/MeshTool.h"
#ifdef WIN32
#ifdef E_DEBUG
#include "meshIO/meshIO.h"
#endif
#endif


///
void
Intersector::getXPairs
(const K_FLD::FloatArray& pos, const std::vector<K_FLD::IntArray*>& Cs,
 std::vector<XPair>& pairs)
{
  E_Float ENLARGE_FACTOR = 0.5;
  size_t nb_comps = Cs.size();
  
  // Global component bounding boxes.
  std::vector<K_SEARCH::BBox3D> BBOX(nb_comps);
  std::vector<E_Int> unodes;
  for (size_t s = 0; s < nb_comps; ++s)
  {
    Cs[s]->uniqueVals(unodes);
    BBOX[s].compute(pos, unodes);
  }
  
  // Get Intersection Matrix (coarse with bounding boxes)
  K_FLD::IntArray is_x(nb_comps,nb_comps, 0);
  std::vector<bool> used(nb_comps, false);
  K_SEARCH::BBox3D res;
  for (size_t s1 = 0; s1 < nb_comps; ++s1)
  {
    for (size_t s2 = s1+1; s2 < nb_comps; ++s2) 
    {
      if (K_SEARCH::BBox3D::intersection(BBOX[s1], BBOX[s2], res))
      {
        is_x(s1,s2) = 1;
        used[s1]=used[s2]=true;
      }
    }
  }

  // Box tolerances and min length of each triangle.
  std::vector<std::vector<E_Float> > MinT3Lengths(nb_comps);
  K_FLD::FloatArray L2;
  std::vector<E_Float> minLs;
  for (size_t c = 0; c < nb_comps; ++c)
  {
    if (!used[c])
      continue;
    
    L2.clear();
    K_CONNECT::MeshTool::computeEdgesSqrLengths<3>(pos, *Cs[c], L2);
    MinT3Lengths[c].resize(L2.cols());

    for (E_Int i = 0; i < L2.cols(); ++i)
      MinT3Lengths[c][i] = ::sqrt(L2(0, i));
  }

  // Bounding Boxes Trees.
  std::vector<TreeType* > trees(nb_comps);
  std::vector< std::vector<BBoxType*> > boxes(nb_comps);
  
  BBoxType** pool = new BBoxType*[nb_comps];
  
  //
  for (size_t s = 0; s < nb_comps; ++s)
  {
    pool[s]=0; //init required before deleting
    
    if (!used[s])
      continue;
    
    __create_boxes(pos, *Cs[s], boxes[s], pool[s], BBOX[s], ENLARGE_FACTOR);
        
    // Build the box tree.
    trees[s] = new TreeType(boxes[s]);
  }

  K_FLD::IntArray cc;
  // Intesection test.
  for (size_t s1 = 0; s1 < nb_comps; ++s1)
  {
    for (size_t s2 = s1+1; s2 < nb_comps; ++s2)
    {
      if (is_x(s1,s2))
        __getXPairs(pos, Cs, s1, s2, trees, boxes, MinT3Lengths, ENLARGE_FACTOR, false, pairs);
    }
  }

  // Final cleaning.
  for (size_t s = 0; s < nb_comps; ++s)
  {
    delete [] pool[s];
    // Destroy the box tree.
    delete trees[s];
  }
}

///
void Intersector::__create_boxes
(const K_FLD::FloatArray& crd, const K_FLD::IntArray& cnt, 
 std::vector<K_SEARCH::BBox3D*>& boxes, K_SEARCH::BBox3D*& pool, K_SEARCH::BBox3D& GBbox, E_Float ENLARGE_FACTOR)
{
  size_t nb_elts(cnt.cols());
  
  boxes.reserve(nb_elts);
  pool = new K_SEARCH::BBox3D[nb_elts];
  
  for (E_Int i = 0; i < 3; ++i){GBbox.minB[i] = K_CONST::E_MAX_FLOAT; GBbox.maxB[i] = -K_CONST::E_MAX_FLOAT;}

  for (size_t i = 0; i < nb_elts; ++i)
  {
    K_SEARCH::BBox3D* box = &pool[i];
    box->compute(crd, cnt.col(i), 3);
    boxes.push_back(box);
    
    for (E_Int j = 0; j < 3; ++j)
    {
      GBbox.minB[j] = (GBbox.minB[j] > box->minB[j]) ? box->minB[j] : GBbox.minB[j];
      GBbox.maxB[j] = (GBbox.maxB[j] < box->maxB[j]) ? box->maxB[j] : GBbox.maxB[j];
    }
  }
  
  if (ENLARGE_FACTOR <= 0.)
    return;
  
  // now extend the boxes
  E_Float mL[3], mL0, dx;
  for (size_t i = 0; i < nb_elts; ++i)
  {
    mL0 = -K_CONST::E_MAX_FLOAT;
    for (size_t k = 0; k < 3; ++k)
    {
      mL[k] = boxes[i]->maxB[k] - boxes[i]->minB[k];
      mL0 = std::max(mL0, mL[k]);
    }
    
    //enlarge it
    dx = ENLARGE_FACTOR * mL0;
    for (size_t k = 0; k < 3; ++k)
    {
      boxes[i]->minB[k] -= dx;
      boxes[i]->maxB[k] += dx;
    }
  }
}

inline bool 
Intersector::__overlap
(const K_FLD::FloatArray& pos, K_FLD::IntArray::const_iterator pS1,
 K_FLD::IntArray::const_iterator pS2,
 E_Float tol, E_Bool tol_is_absolute, E_Float Lc)
{
   bool ovlap = __pointIsInside(pos, pS1, pS2, tol, tol_is_absolute, Lc);
   if (!ovlap)
     ovlap = __edgeOverlaps(pos, pS1, pS2, tol, tol_is_absolute, Lc);
   return ovlap;
}

///
bool
Intersector::__pointIsInside
(const K_FLD::FloatArray& pos, K_FLD::IntArray::const_iterator pS1,
 K_FLD::IntArray::const_iterator pS2, E_Float tol, E_Bool tol_is_absolute,
 E_Float Lc)
{
  E_Int Ni;
  bool ret = false, inside;
  E_Float d, ztol(tol), U[2];

  if (!tol_is_absolute)
    ztol *= Lc;
  
  // T2 -> T1
  for (E_Int n = 0; (n < 3) && !ret; ++n)
  {
    Ni = *(pS2+n);
    d = K_MESH::Triangle::minDistanceToPoint(pos, pS1, pos.col(Ni), U, inside);
    ret = inside && (d < ztol);
  }
  
  // T1 -> T2
  for (E_Int n = 0; (n < 3) && !ret; ++n)
  {
    Ni = *(pS1+n);
   
    d = K_MESH::Triangle::minDistanceToPoint(pos, pS2, pos.col(Ni), U, inside);    
    ret = inside && (d < ztol);
  }
  return (E_Bool)ret;
}

///
bool
Intersector::__edgeOverlaps
(const K_FLD::FloatArray& pos, K_FLD::IntArray::const_iterator pS1,
 K_FLD::IntArray::const_iterator pS2, E_Float tol, E_Bool tol_is_absolute,
 E_Float Lc)
{
   E_Int Ni, Nj;
  bool ret = false;
  E_Float lambda, UV[2], min_d;
  E_Bool parallel, coincident;
  E_Float eps(E_EPSILON), abstol(tol), edge_tol_rel(0.1);

  if (!tol_is_absolute)
    abstol *= Lc;
  eps /= Lc;

  E_Float Qi[3], normal1[3], normal2[3];
  K_MESH::Triangle::normal(pos.col(*pS1), pos.col(*(pS1+1)), pos.col(*(pS1+2)), normal1);
  K_FUNC::normalize<3>(normal1);
  K_MESH::Triangle::normal(pos.col(*pS2), pos.col(*(pS2+1)), pos.col(*(pS2+2)), normal2);
  K_FUNC::normalize<3>(normal2);
  for (E_Int k = 0; k < 3; ++k)
  {
    normal1[k] *= abstol;
    normal2[k] *= abstol;
  }

  // T2 -> T1
  
  for (E_Int n = 0; (n < 3) && !ret; ++n)
  {
    Ni = *(pS2+n);
    Nj = *(pS2+(n+1)%3);
    
    for (E_Int m = 0; (m < 3) && !ret; ++m)
    {
      for (E_Int k = 0; k < 3; ++k)
        Qi[k] = pos(k, *(pS1+m)) + normal1[k];

      K_MESH::Triangle::planeLineMinDistance<3>(pos.col(*(pS1+m)), pos.col(*(pS1+(m+1)%3)), Qi, pos.col(Ni), pos.col(Nj), tol, tol_is_absolute, lambda, UV, parallel, coincident, min_d);

      ret = (UV[0] > -eps) && (UV[0] < 1. + eps) && (::fabs(UV[1]) < 1. + eps)  && (lambda > -edge_tol_rel) && (lambda < 1. + edge_tol_rel);
    }
  }

  // T1 -> T2
  
  for (E_Int n = 0; (n < 3) && !ret; ++n)
  {
    Ni = *(pS1+n);
    Nj = *(pS1+(n+1)%3);

    for (E_Int m = 0; (m < 3) && !ret; ++m)
    {
      for (E_Int k = 0; k < 3; ++k)
        Qi[k] = pos(k, *(pS2+m)) + normal2[k];

      K_MESH::Triangle::planeLineMinDistance<3>(pos.col(*(pS2+m)), pos.col(*(pS2+(m+1)%3)), Qi, pos.col(Ni), pos.col(Nj), tol, tol_is_absolute, lambda, UV, parallel, coincident, min_d);

      ret = (UV[0] > -eps) && (UV[0] < 1. + eps) && (::fabs(UV[1]) < 1. + eps)  && (lambda > -edge_tol_rel) && (lambda < 1. + edge_tol_rel);
    }
  }
    
    
 
  return (E_Bool)ret;

}

void
Intersector::__getXPairs
(const K_FLD::FloatArray& pos, const std::vector<K_FLD::IntArray*>& Cs,
 E_Int s1, E_Int s2,
 const std::vector<TreeType* >& trees, 
 const std::vector< std::vector<BBoxType*> >& boxes,
 const std::vector<std::vector<E_Float> >& Lengths,
 E_Float tol, bool tol_is_absolute,
 std::vector<XPair>& pairs)
{
  E_Float                       minL1, minL2, Lc;
  std::vector<E_Int>            bxs;
  K_FLD::IntArray::const_iterator pS1, pS2;
  size_t                        j;
  
  // query intensively a small tree is better
  if (boxes[s1].size() < boxes[s2].size())
    std::swap(s1,s2);
  
  E_Int cols1 = Cs[s1]->cols(), c1;
  //E_Int cols2 = Cs[s2]->cols();
  E_Int c2;

  for (c1 = 0; c1 < cols1; ++c1)
  {
    pS1 = Cs[s1]->col(c1);

    bxs.clear();
    trees[s2]->getOverlappingBoxes(boxes[s1][c1]->minB, boxes[s1][c1]->maxB, bxs);

    minL1 = Lengths[s1][c1];

    for (j = 0; j < bxs.size(); ++j)
    {
      c2 = bxs[j];
      
      minL2 = Lengths[s2][c2];

      Lc = std::max(minL2, minL1);

      pS2 = Cs[s2]->col(c2);
         
      if (__overlap(pos, pS1, pS2, tol, tol_is_absolute, Lc))
        pairs.push_back(XPair(s1, c1, s2, c2));
    }
  }
}

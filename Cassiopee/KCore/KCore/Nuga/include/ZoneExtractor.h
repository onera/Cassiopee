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

#ifndef __ZONEEXTRACTOR_H__
#define	__ZONEEXTRACTOR_H__

#include <vector>
//#include <algorithm>
#define Vector_t std::vector

#include "Nuga/include/ArrayAccessor.h"
#include "BbTree.h"
#include "Nuga/include/ngon_t.hxx"

#define TEMPLATE_COORD_CONNECT template <typename Coordinate_t, typename Connectivity_t>

#define PT_INSIDE_BOX(px,py,pz, minB, maxB) ((minB[0]-EPSILON < px) && (px < maxB[0]+EPSILON) && (minB[1]-EPSILON < py) && (py < maxB[1]+EPSILON) && (minB[2]-EPSILON < pz) && (pz < maxB[2]+EPSILON))

namespace K_CONNECT
{

// Extract 
TEMPLATE_COORD_CONNECT
class ZoneExtractor
{
public:
  
  typedef K_FLD::ArrayAccessor<Coordinate_t> ACoordinate_t;
  typedef K_FLD::ArrayAccessor<Connectivity_t> AConnectivity_t;
    
  ///
  ZoneExtractor(){}
  ///
  ~ZoneExtractor(){}
  
  /// Retrieves any cell in mesh (crd,cnt) inside box(mB,MB).
  template <short DIM, typename ELT_t>
  void getInBox(const E_Float* mB, const E_Float* MB, E_Float tol, const ACoordinate_t& crd, const AConnectivity_t& cnt, Vector_t<E_Int>& xCells_ids);
  template <short DIM>
  void getInBox
  (const E_Float* mB, const E_Float* MB, E_Float tol, const ACoordinate_t& crd, const ngon_unit& ng, Vector_t<E_Int>& xCells_ids);
  template <short DIM>
  void getInBox_omp(const E_Float* mB, const E_Float* MB, E_Float tol, const ACoordinate_t& crd, const AConnectivity_t& cnt, Vector_t<E_Int>& xCells_ids);
  
  /// Retrieves any point cloud crd inside box(mB,MB).
  template <short DIM>
  void getInBox(const E_Float* mB, const E_Float* MB, E_Float tol, const ACoordinate_t& crd, Vector_t<E_Int>& xCells_ids);
  template <short DIM>
  void getInBox_omp(const E_Float* mB, const E_Float* MB, E_Float tol, const ACoordinate_t& crd, Vector_t<E_Int>& xCells_ids);

  /// Retrieves any box in boxes inside box(mB,MB).
  template <short DIM>
  void getInBox(const E_Float* mB, const E_Float* MB, E_Float tol, const Vector_t<K_SEARCH::BoundingBox<DIM>*>& boxes, Vector_t<E_Int>& xnode_ids);
  /// same as above but with a flag instead of a list
  template <short DIM>
  void getInBox(const E_Float* mB, const E_Float* MB, E_Float tol, const Vector_t<K_SEARCH::BoundingBox<DIM>*>& boxes, Vector_t<bool>& is_in);
  //omp version
  template <short DIM>
  void getInBox_omp(const E_Float* mB, const E_Float* MB, E_Float tol, const Vector_t<K_SEARCH::BoundingBox<DIM>*>& boxes, Vector_t<bool>& is_in);
};

typedef ZoneExtractor<K_FLD::FldArrayF, K_FLD::FldArrayI> FldZoneExtractor;
typedef ZoneExtractor<K_FLD::FloatArray, K_FLD::IntArray> DynZoneExtractor;


///
TEMPLATE_COORD_CONNECT
template <short DIM, typename ELT_t>
void ZoneExtractor<Coordinate_t, Connectivity_t>::getInBox
(const E_Float* mB, const E_Float* MB, E_Float tol, const ACoordinate_t& crd, const AConnectivity_t& cnt, Vector_t<E_Int>& xCells_ids)
{
  xCells_ids.clear();
  
  K_SEARCH::BoundingBox<DIM> BB(mB, MB);
    
  for (E_Int i = 0; i < cnt.size(); ++i)
  {
    ELT_t e;
    cnt.getEntry(i, e);
    K_SEARCH::BoundingBox<DIM> bb;
    e.bbox(crd, bb);
        
    if (K_SEARCH::BbTree<DIM>::boxesAreOverlapping(&BB, &bb, tol))
      xCells_ids.push_back(i);
  }
}

///
TEMPLATE_COORD_CONNECT
template <short DIM>
void ZoneExtractor<Coordinate_t, Connectivity_t>::getInBox
(const E_Float* mB, const E_Float* MB, E_Float tol, const ACoordinate_t& crd, const ngon_unit& PGs, Vector_t<E_Int>& xCells_ids)
{
  xCells_ids.clear();
  
  K_SEARCH::BoundingBox<DIM> BB(mB, MB);
  
  E_Int s;
  K_FLD::IntArray e;
  
  for (E_Int i = 0; i < PGs.size(); ++i)
  {
    s = PGs.stride(i);
    e.resize(1, s);
    const E_Int* ptr = PGs.get_facets_ptr(i);
    for (size_t k=0; k < s; ++k)e[k]=*(ptr+k)-1;
    K_SEARCH::BoundingBox<DIM> bb(crd, e.begin(), s);
    
    if (K_SEARCH::BbTree<DIM>::boxesAreOverlapping(&BB, &bb, tol))
      xCells_ids.push_back(i);
  }
}

///
TEMPLATE_COORD_CONNECT
template <short DIM>
void ZoneExtractor<Coordinate_t, Connectivity_t>::getInBox
(const E_Float* mB, const E_Float* MB, E_Float tol, const ACoordinate_t& crd, Vector_t<E_Int>& xnode_ids)
{
  xnode_ids.clear();
  
  size_t sz = crd.size();
  
  const E_Float* px = crd.array().begin(crd.posX(0));
  const E_Float* py = crd.array().begin(crd.posX(1));
  const E_Float* pz = crd.array().begin(crd.posX(2));
  
  for (size_t i = 0; i < sz; ++i)
  {
    if (PT_INSIDE_BOX(px[i],py[i],pz[i],mB,MB))
      xnode_ids.push_back(i);
  }
}

///
TEMPLATE_COORD_CONNECT
template <short DIM>
void ZoneExtractor<Coordinate_t, Connectivity_t>::getInBox
(const E_Float* mB, const E_Float* MB, E_Float tol, const Vector_t<K_SEARCH::BoundingBox<DIM>*>& boxes, Vector_t<E_Int>& xCells_ids)
{
  xCells_ids.clear();

  K_SEARCH::BoundingBox<DIM> BB(mB, MB);

  for (E_Int i = 0; i < boxes.size(); ++i)
  {
    if (K_SEARCH::BbTree<DIM>::boxesAreOverlapping(&BB, boxes[i], tol))
      xCells_ids.push_back(i);
  }
}

///
TEMPLATE_COORD_CONNECT
template <short DIM>
void ZoneExtractor<Coordinate_t, Connectivity_t>::getInBox
(const E_Float* mB, const E_Float* MB, E_Float tol, const Vector_t<K_SEARCH::BoundingBox<DIM>*>& boxes, Vector_t<bool>& is_in)
{
  is_in.clear();
  is_in.resize(boxes.size());

  K_SEARCH::BoundingBox<DIM> BB(mB, MB);

  for (size_t i = 0; i < boxes.size(); ++i)
    is_in[i]=K_SEARCH::BbTree<DIM>::boxesAreOverlapping(&BB, boxes[i], tol);
}

///////////////// OMP ///////////////////////////

///
TEMPLATE_COORD_CONNECT
template <short DIM>
void ZoneExtractor<Coordinate_t, Connectivity_t>::getInBox_omp
(const E_Float* mB, const E_Float* MB, E_Float tol, const ACoordinate_t& crd, const AConnectivity_t& cnt, Vector_t<E_Int>& is_in)
{
  is_in.clear();
  
  E_Int sz = cnt.size();
  
  K_SEARCH::BoundingBox<DIM> BB(mB, MB);
  
  K_FLD::IntArray e;
  E_Int s;
  
  is_in.resize(sz);
  
#pragma omp parallel for private(e,s) 
  for (E_Int i = 0; i < sz; ++i)
  {
    s = cnt.stride(i);
    e.reserve(1, s);
    cnt.getEntry(i, e.begin());
    
    K_SEARCH::BoundingBox<DIM> bb(crd, e.begin(), s);
    
    is_in[i]=K_SEARCH::BbTree<DIM>::boxesAreOverlapping(&BB, &bb, tol);
  }
}

TEMPLATE_COORD_CONNECT
template <short DIM>
void ZoneExtractor<Coordinate_t, Connectivity_t>::getInBox_omp
(const E_Float* mB, const E_Float* MB, E_Float tol, const ACoordinate_t& crd, Vector_t<E_Int>& is_in)
{
  is_in.clear();
  
  size_t sz = crd.size();
  
  is_in.resize(sz);
  
  const E_Float* px = crd.array().begin(crd.posX(0));
  const E_Float* py = crd.array().begin(crd.posX(1));
  const E_Float* pz = crd.array().begin(crd.posX(2));

#pragma omp parallel for
  for (E_Int i = 0; i < sz; ++i)
    is_in[i]=PT_INSIDE_BOX(px[i],py[i],pz[i],mB,MB);
}

///
TEMPLATE_COORD_CONNECT
template <short DIM>
void ZoneExtractor<Coordinate_t, Connectivity_t>::getInBox_omp
(const E_Float* mB, const E_Float* MB, E_Float tol, const Vector_t<K_SEARCH::BoundingBox<DIM>*>& boxes, Vector_t<bool>& is_in)
{
  is_in.clear();
  is_in.resize(boxes.size());

  K_SEARCH::BoundingBox<DIM> BB(mB, MB);

#pragma omp parallel for
  for (E_Int i = 0; i < boxes.size(); ++i)
    is_in[i]=K_SEARCH::BbTree<DIM>::boxesAreOverlapping(&BB, boxes[i], tol);
}

}

#endif

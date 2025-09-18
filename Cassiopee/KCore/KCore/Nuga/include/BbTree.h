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

#ifndef __KCORE_SEARCH_BBTREE_H__
#define __KCORE_SEARCH_BBTREE_H__

#include "Nuga/include/DynArray.h"
#include "Nuga/include/ArrayAccessor.h"
#include "Nuga/include/defs.h"
#include "Nuga/include/Triangle.h"
#include "Nuga/include/macros.h"
#include "Nuga/include/ngon_unit.h"

namespace K_SEARCH
{

///
template <short DIM>
class BoundingBox{

public:
  /// Constructors
  BoundingBox(const E_Float* mB, const E_Float* MB)
  {for (E_Int i = 0; i < DIM; ++i){minB[i] = *(mB++); maxB[i] = *(MB++);}}

  ///
  BoundingBox(const K_FLD::FloatArray& pos, const E_Int* nodes, size_t n)
  {
    const E_Float* Pi;
    for (E_Int i = 0; i < DIM; ++i)
    {minB[i] = NUGA::FLOAT_MAX; maxB[i] = -NUGA::FLOAT_MAX;}

    for (size_t i = 0; i < n; ++i)
    {
      Pi = pos.col(nodes[i]);
      for (E_Int j = 0; j < DIM; ++j)
      {
        minB[j] = (minB[j] > Pi[j]) ? Pi[j] : minB[j];
        maxB[j] = (maxB[j] < Pi[j]) ? Pi[j] : maxB[j];
      }
    }
  }

  ///
  BoundingBox(const K_FLD::FloatArray& pos, E_Int N0=-1)
  {
    const E_Float* Pi;
    for (E_Int i = 0; i < DIM; ++i)
    { minB[i] = NUGA::FLOAT_MAX; maxB[i] = -NUGA::FLOAT_MAX; }

    if (N0 < 0 || N0 >= pos.cols()) N0 = pos.cols(); //if max id is not (or badly) specified

    for (E_Int i = 0; i < N0; ++i)
    {
      Pi = pos.col(i);
      for (E_Int j = 0; j < DIM; ++j)
      {
        minB[j] = (minB[j] > Pi[j]) ? Pi[j] : minB[j];
        maxB[j] = (maxB[j] < Pi[j]) ? Pi[j] : maxB[j];
      }
    }
  }
  
  ///
  template<typename CoordinateArray_t>
  BoundingBox(const K_FLD::ArrayAccessor<CoordinateArray_t>& pos, const E_Int* nodes, size_t n)
  {
    this->compute(pos, nodes, n);
  }

  inline const E_Float& getMinB(E_Int axis) const {return minB[axis];}
  inline const E_Float& getMaxB(E_Int axis) const {return maxB[axis];}

  // Warning : strict comparison
  bool operator==(const BoundingBox& rhs)
  {
    for (E_Int i = 0; i < DIM; ++i)
    {
      if (minB[i] != rhs.minB[i]) return false;
      if (maxB[i] != rhs.maxB[i]) return false;
    }
    return true;
  }
  
  // default impl for 3D. Specialization below after the class scope
  BoundingBox& operator=(const BoundingBox& rhs)
  {
  	minB[0] = rhs.minB[0]; minB[1] = rhs.minB[1]; minB[2] = rhs.minB[2];
  	maxB[0] = rhs.maxB[0]; maxB[1] = rhs.maxB[1]; maxB[2] = rhs.maxB[2];

    return *this;
  }

  bool is_included(const BoundingBox& rhs) const
  {
    for (E_Int i = 0; i < DIM; ++i)
    {
      if (minB[i] + ZERO_M < rhs.minB[i]) return false;
      if (maxB[i] > rhs.maxB[i] + ZERO_M) return false;
    }
    
    return true;
  }

  // default impl for 3D. Specialization below after the class scope
  void enlarge(double RTOL)
  {
    // compute Lref (min & non-null box size)
    double dx = (maxB[0] - minB[0]);
    double dy = (maxB[1] - minB[1]);
    double dz = (maxB[2] - minB[2]);

    double Lref = std::max(dx, std::max(dy, dz)); // ensure to get a non-null val

    if (dx > ZERO_M) Lref = std::min(Lref, dx);
    if (dy > ZERO_M) Lref = std::min(Lref, dy);
    if (dz > ZERO_M) Lref = std::min(Lref, dz);

    Lref *= RTOL;

    minB[0] -= Lref;
    minB[1] -= Lref;
    minB[2] -= Lref;
    maxB[0] += Lref;
    maxB[1] += Lref;
    maxB[2] += Lref;
  }

  // default impl for 3D. Specialization below after the class scope
  void merge(const BoundingBox& b)
  {
    minB[0] = std::min(minB[0], b.minB[0]);
    minB[1] = std::min(minB[1], b.minB[1]);
    minB[2] = std::min(minB[2], b.minB[2]);
    maxB[0] = std::max(maxB[0], b.maxB[0]);
    maxB[1] = std::max(maxB[1], b.maxB[1]);
    maxB[2] = std::max(maxB[2], b.maxB[2]);
  }

  static bool intersection(const BoundingBox& b1, const BoundingBox& b2, BoundingBox& b)
  {
    std::vector<E_Float> palma(4);
        
    for (E_Int i = 0; i < DIM; ++i)
    {
      if (b1.minB[i] > b2.maxB[i]) return false;
      if (b2.minB[i] > b1.maxB[i]) return false;
    }  

    for (E_Int i = 0; i < DIM; ++i)
    {
      palma[0] = b1.minB[i]; palma[1] = b1.maxB[i]; palma[2] = b2.minB[i]; palma[3] = b2.maxB[i];
      std::sort(palma.begin(), palma.end());
      b.minB[i] = palma[1];
      b.maxB[i] = palma[2];
    }

    return true;
  }
  
  template<typename ng_t>
  void convert2NG(K_FLD::FloatArray& crd, ng_t& ng)
  {
    ng.clear();
    crd.clear();
    crd.resize(3, 8);
    
    //P0
    crd(0,0) = minB[0];
    crd(1,0) = minB[1];
    crd(2,0) = minB[2];
    
    //P1
    crd(0,1) = maxB[0];
    crd(1,1) = minB[1];
    crd(2,1) = minB[2];
    
    //P2
    crd(0,2) = maxB[0];
    crd(1,2) = maxB[1];
    crd(2,2) = minB[2];
    
    //P3
    crd(0,3) = minB[0];
    crd(1,3) = maxB[1];
    crd(2,3) = minB[2];
    
    //P4
    crd(0,4) = minB[0];
    crd(1,4) = minB[1];
    crd(2,4) = maxB[2];
    
    //P5
    crd(0,5) = maxB[0];
    crd(1,5) = minB[1];
    crd(2,5) = maxB[2];
    
    //P6
    crd(0,6) = maxB[0];
    crd(1,6) = maxB[1];
    crd(2,6) = maxB[2];
    
    //P7
    crd(0,7) = minB[0];
    crd(1,7) = maxB[1];
    crd(2,7) = maxB[2];
    
    
    E_Int Q4bot[] =  {1,4,3,2};
    E_Int Q4top[] =  {5,6,7,8};
    E_Int Q4left[]=  {1,5,8,4};
    E_Int Q4right[]= {2,3,7,6};
    E_Int Q4front[]= {1,2,6,5};
    E_Int Q4back[]=  {3,4,8,7};
    
    ng.PGs.add(4, Q4bot);
    ng.PGs.add(4, Q4top);
    ng.PGs.add(4, Q4left);
    ng.PGs.add(4, Q4right);
    ng.PGs.add(4, Q4front);
    ng.PGs.add(4, Q4back);
    
    E_Int PH[] = {1,2,3,4,5,6};
    ng.PHs.add(6, PH);
    
    ng.PGs.updateFacets();
    ng.PHs.updateFacets();
  }

  ~BoundingBox()
  {
  }
  
  template <template <typename ELEM, typename ALLOC = std::allocator<ELEM> > class Vector>
  void compute(const K_FLD::FloatArray& pos, const Vector<E_Int> & indices, E_Int idx_start=0)
  {
    const E_Float* Pi;
    for (E_Int i = 0; i < DIM; ++i)
    {minB[i] = NUGA::FLOAT_MAX; maxB[i] = -NUGA::FLOAT_MAX;}

    for (size_t i = 0; i < indices.size(); ++i)
    {
      Pi = pos.col(indices[i] - idx_start);
      for (E_Int j = 0; j < DIM; ++j)
      {
        minB[j] = (minB[j] > Pi[j]) ? Pi[j] : minB[j];
        maxB[j] = (maxB[j] < Pi[j]) ? Pi[j] : maxB[j];
      }
    }
  }
  
  ///
  template<typename CoordinateArray_t, template <typename ELEM, typename ALLOC = std::allocator<ELEM> > class Vector>
  void compute
  (const K_FLD::ArrayAccessor<CoordinateArray_t>& pos, const Vector<E_Int> & indices, E_Int idx_start=0)
  {
    E_Float Pi[DIM];
    for (E_Int i = 0; i < DIM; ++i) Pi[i] = -1.;
    for (E_Int i = 0; i < DIM; ++i)
    { minB[i] = NUGA::FLOAT_MAX; maxB[i] = -NUGA::FLOAT_MAX; }

    for (size_t i = 0; i < indices.size(); ++i)
    {
      pos.getEntry(indices[i] - idx_start, Pi);
      for (E_Int j = 0; j < DIM; ++j)
      {
        minB[j] = (minB[j] > Pi[j]) ? Pi[j] : minB[j];
        maxB[j] = (maxB[j] < Pi[j]) ? Pi[j] : maxB[j];
      }
    }
  }
  
  ///
  template<typename CoordinateArray_t>
  void compute
  (const K_FLD::ArrayAccessor<CoordinateArray_t>& pos)
  {
    E_Float Pi[DIM];
    for (E_Int i = 0; i < DIM; ++i)
    {minB[i] = NUGA::FLOAT_MAX; maxB[i] = -NUGA::FLOAT_MAX;}

    for (E_Int i = 0; i < pos.size(); ++i)
    {
      pos.getEntry(i, Pi);
      for (E_Int j = 0; j < DIM; ++j)
      {
        minB[j] = (minB[j] > Pi[j]) ? Pi[j] : minB[j];
        maxB[j] = (maxB[j] < Pi[j]) ? Pi[j] : maxB[j];
      }
    }
  }
  
  ///
  void compute
  (const K_FLD::FloatArray& pos)
  {
    const E_Float* Pi;
    for (E_Int i = 0; i < DIM; ++i)
    {minB[i] = NUGA::FLOAT_MAX; maxB[i] = -NUGA::FLOAT_MAX;}

    for (E_Int i = 0; i < pos.cols(); ++i)
    {
      Pi = pos.col(i);
      for (E_Int j = 0; j < DIM; ++j)
      {
        minB[j] = (minB[j] > Pi[j]) ? Pi[j] : minB[j];
        maxB[j] = (maxB[j] < Pi[j]) ? Pi[j] : maxB[j];
      }
    }
  }
  
  ///
  template<typename CoordinateArray_t>
  void compute
  (const CoordinateArray_t& pos, const E_Int* nodes, size_t n, E_Int index_start = 0)
  {
    const E_Float* Pi;
    for (E_Int i = 0; i < DIM; ++i)
    {minB[i] = NUGA::FLOAT_MAX; maxB[i] = -NUGA::FLOAT_MAX;}

    for (size_t i = 0; i < n; ++i)
    {
      Pi = pos.col(nodes[i] - index_start);
      for (E_Int j = 0; j < DIM; ++j)
      {
        minB[j] = (minB[j] > Pi[j]) ? Pi[j] : minB[j];
        maxB[j] = (maxB[j] < Pi[j]) ? Pi[j] : maxB[j];
      }
    }
  }
  
  ///
  template<typename CoordinateArray_t>
  void compute
  (const K_FLD::ArrayAccessor<CoordinateArray_t>& pos, const E_Int* nodes, size_t n, E_Int index_start = 0)
  {
    E_Float Pi[DIM];
    for (E_Int i = 0; i < DIM; ++i)
    {minB[i] = NUGA::FLOAT_MAX; maxB[i] = -NUGA::FLOAT_MAX;}

    for (size_t i = 0; i < n; ++i)
    {
      pos.getEntry(nodes[i] - index_start, Pi);
      for (E_Int j = 0; j < DIM; ++j)
      {
        minB[j] = (minB[j] > Pi[j]) ? Pi[j] : minB[j];
        maxB[j] = (maxB[j] < Pi[j]) ? Pi[j] : maxB[j];
      }
    }
  }

public://fixme //private:
  BoundingBox(){for (E_Int i = 0; i < DIM; ++i)
    {minB[i] = NUGA::FLOAT_MAX; maxB[i] = -NUGA::FLOAT_MAX;}}

  
public://fixme //private:
  E_Float minB[DIM];
  E_Float maxB[DIM];
};

typedef BoundingBox<2> BBox2D;
typedef BoundingBox<3> BBox3D;


// 2D - specializations
template <>
inline BBox2D& BBox2D::operator=(const BBox2D& rhs)
{
  minB[0] = rhs.minB[0]; minB[1] = rhs.minB[1];
  maxB[0] = rhs.maxB[0]; maxB[1] = rhs.maxB[1];

  return *this;
}

template <>
inline void BBox2D::enlarge(double RTOL)
{
    // compute Lref (min & non-null box size)
    double dx = (maxB[0] - minB[0]);
    double dy = (maxB[1] - minB[1]);

    double Lref = std::max(dx, dy);

    if (dx > ZERO_M) Lref = std::min(Lref, dx);
    if (dy > ZERO_M) Lref = std::min(Lref, dy);

    Lref *= RTOL;

    minB[0] -= Lref;
    minB[1] -= Lref;
    
    maxB[0] += Lref;
    maxB[1] += Lref;
}

template <>
inline void BBox2D::merge(const BBox2D& b)
{
    minB[0] = std::min(minB[0], b.minB[0]);
    minB[1] = std::min(minB[1], b.minB[1]);
    
    maxB[0] = std::max(maxB[0], b.maxB[0]);
    maxB[1] = std::max(maxB[1], b.maxB[1]);
}

///
template <short DIM, typename BBoxType = BoundingBox<DIM> >
class BbTree {

  public: /** Typedefs */

    typedef BbTree                  self_type;
    typedef NUGA::size_type         size_type;
    typedef K_FLD::IntArray         tree_array_type;
    typedef BBoxType                box_type;
   

#define BBTREE_ROWS 2 // the first row contains the node id, the second(third) contains the left(right) child column id.

  public: /** Constructors and Destructor */

    /// Builds a tree and inserts the boxes from begin to end.
    BbTree(const std::vector<BBoxType*>& boxes, E_Float tolerance=EPSILON);
    
    template <typename array_t>
    BbTree(const K_FLD::FloatArray& crd, const array_t& ng, E_Float tolerance=EPSILON);
    //fixme : hack overload rather than specialization to distinghuish ngon_unit (viso icc 15 compil issue)
    inline BbTree(const K_FLD::FloatArray& crd, const ngon_unit& pgs, E_Float tolerance=EPSILON);
    //fixme : hack overload rather than specialization to distinghuish ngon_unit (viso icc 15 compil issue)
    inline BbTree(const K_FLD::FloatArray& crd, const K_FLD::IntArray& cnt);


    /// Destructor.
    ~BbTree()
    {
      if (_root_id != 0)
        for (size_t i = _root_id; i < _boxes.size(); ++i) delete _boxes[i];
      
      if (_owes_boxes) delete [] _pool;
    };

  public: /** Query methods */

    /// Returns all the boxes which overlap the input box.
    /** Warning: out is not cleared upon entry.*/
    void getOverlappingBoxes(const E_Float* minB, const E_Float* maxB, std::vector<size_type>& out) const;
    
    /// Returns true if an overlapping box is found (false otherwise).
    bool hasAnOverlappingBox(const E_Float* minB, const E_Float* maxB) const;

    /// Returns all the boxes which intersects the input ray defined by P0 and P1.
    /** Warning: out is not cleared upon entry.*/
    void getIntersectingBoxes(const E_Float* P0, const E_Float* P1, std::vector<size_type>& out, E_Float tolerance=EPSILON, bool strict=false) const;
    
    /// Returns the global bounding box : the tree's bounding box
    const BBoxType* getGlobalBox(){return _boxes[_root_id];}
    
    ///
    static E_Bool boxesAreOverlapping(const BBoxType * bb1, const BBoxType * bb2, const E_Float& tol);
    
    ///
    static E_Bool box1IsIncludedinbox2(const BBoxType * bb1, const BBoxType * bb2, const E_Float& tol);

  private: /** Insertion method */

    /** Inserts the boxes in the tree when constructing the BbTree.
     *  Returns the id of the boxes bounding box.*/
    size_type __insert(const std::vector<size_type>::const_iterator& begin,
                       const std::vector<size_type>::const_iterator& end);

    /// Computes the bounding box of the input boxes.
    void __getBoxesBox(const std::vector<BBoxType*>&  boxes, const std::vector<size_type>::const_iterator& begin,
                       const std::vector<size_type>::const_iterator& end, BBoxType& box) const;

    ///
    inline E_Int __getLongestSideAxis(const BBoxType& box) const;

    ///
    void __getOverlappingBoxes(BBoxType* box, E_Int node, std::vector<size_type>& out) const;
    
    ///
    bool __hasAnOverlappingBox(BBoxType* box, E_Int node) const;

    ///
    void __getIntersectingBoxes(const E_Float* P0, const E_Float* P1, E_Int node, const E_Float& abstol, std::vector<size_type>& out) const;
    ///
    void __getIntersectingBoxesStrict(const E_Float* P0, const E_Float* P1, E_Int node, const E_Float& abstol, std::vector<size_type>& out, const BBoxType& boxSeg) const;

    ///
    inline E_Bool __boxIntersectRay(const BBoxType& box, const E_Float* P0, const E_Float* P1, const E_Float& abstol) const;
    ///
    inline E_Bool __boxIntersectSeg(const BBoxType& box, const E_Float* P0, const E_Float* P1, const E_Float& abstol, const BBoxType& boxSeg) const;

#ifdef E_TIME1
public:
  static double _append_tree;
  static double _get_box_boxes;
  static double _get_longest_axis;
  static double _build_vectors;
  static double _fill_vectors;
#endif

 // private:
  public:
    /// BB tree
    tree_array_type          _tree;

    /// boxes
    std::vector<BBoxType*>   _boxes;

    /// tolerance
    E_Float                  _tolerance;

    /// Empty column for the tree.
    size_type                _newL[BBTREE_ROWS];

    /// 
    size_type                _root_id;
    
    ///
    bool                    _owes_boxes;
    BBoxType*               _pool;
    
}; // End class BbTree

typedef BbTree<2, BBox2D> BbTree2D;
typedef BbTree<3, BBox3D> BbTree3D;

} // end namespace


#include "Nuga/include/BbTree.cxx"




#endif /* BBTREE_H_ */

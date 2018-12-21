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

#ifndef KDTREE_H_
#define KDTREE_H_

#include "Fld/DynArray.h"
#include "Def/DefCplusPlusConst.h"
#include "Fld/ArrayAccessor.h"
#include <functional>
#include <algorithm>
#include <vector>

namespace K_SEARCH
{
// ============================================================================
/// Predicate for sorting a container of points according to a given axis.
// ============================================================================
template <typename CoordArrayType>
struct sortingPredicate : public std::binary_function <E_Int, E_Int, bool> 
{
  /// Required by the median balancing method.

  sortingPredicate (const K_FLD::ArrayAccessor<CoordArrayType>& posAcc, E_Int axis)
    : _posAcc(posAcc), _axis(axis){}

  inline bool operator() (E_Int Ni, E_Int Nj) const  
  { return ( _posAcc.getVal(Ni, _axis) <  _posAcc.getVal(Nj, _axis) ); } 

  inline void setAxis(E_Int a){_axis=a;}

  const K_FLD::ArrayAccessor<CoordArrayType>& _posAcc;
  E_Int _axis;
};

///
template <typename ArrayType = K_FLD::FloatArray>
class KdTree {

  public: /** Typedefs */

    typedef           KdTree                            self_type;
    typedef           K_FLD::ArrayAccessor<ArrayType>   coord_access_type;
    typedef           E_Int                             size_type;
    typedef           E_Float                           value_type;
    typedef           K_FLD::IntArray                   tree_array_type;
    typedef           const value_type*                 const_coord_iterator;
    typedef           tree_array_type::iterator         tree_iterator;
    typedef           tree_array_type::const_iterator   const_tree_iterator;

  public: /** Constructors and Destructor */

    /// Build a tree and inserts all the nodes of the input coordinate array pos.
    KdTree(const coord_access_type& posAcc, E_Float tolerance=E_EPSILON);

    /// Builds a tree and inserts only the valid nodes from begin to end.
    KdTree(const coord_access_type& posAcc, 
           std::vector<size_type> indices/*passed by value to preserve input*/,
            E_Float tolerance=E_EPSILON);

    /// Destructor.
    ~KdTree(){};
    
  public: /*Accessors */
    
    E_Int nb_nodes(){ return _tree_sz;}

  public: /** Insertion methods */

    /** Insert a node in the tree. */
    void insert(size_type N);

  public: /** Query methods */

    /// Returns a close node to a given point. It is not necessarily the closest one.
    size_type getClose(const E_Float* pt) const;
    size_type getClose(const E_Float *pt, E_Float& d2) const;

    /// Returns a close node to node N. It is not necessarily the closest one.
    size_type getClose(size_type N) const;
    size_type getClose(size_type N, E_Float& d2) const;

    /// Returns the global closest node.
    size_type getClosest(const E_Float* pt) const;
    size_type getClosest(const E_Float *pt, E_Float& d2) const;

    /// Returns the global closest node to node of index N.
    size_type getClosest(E_Int N) const;
    size_type getClosest(size_type N, E_Float& d2) const;

    /// Returns all the nodes in the input box by appending the vector 'out'.
    /** Warning: out is not cleared upon entry.*/
    void getInBox(const E_Float* minB, const E_Float* maxB, std::vector<size_type>& out) const;

     /// Returns all the nodes in the input sphere centered on C by appending the vector 'out'.
    /** Warning: out is not cleared upon entry.*/
    void getInSphere(const E_Float* C, E_Float radius, std::vector<size_type>& out) const;
    void getInSphere(const E_Float* C, E_Float radius, std::vector<size_type>& out, std::vector<E_Float>& d2) const;
    /// Returns all the nodes in the input sphere centered on node N by appending the vector 'out'.
    /** Warning: out is not cleared upon entry.*/
    void getInSphere(E_Int N, E_Float radius, std::vector<size_type>& out) const;
    void getInSphere(E_Int N, E_Float radius, std::vector<size_type>& out, std::vector<E_Float>& d2) const;

  private:
    /// Underneath method to insert a node in an existing tree (assuming that _tree has the right size).
    void __insert(size_type N);

    /// Underneath algorithm for creating a tree with a set of nodes (balanced construction).
    template <typename InputIterator>
    E_Int __insert(InputIterator begin,InputIterator end, size_type depth);

    /// Underneath algorithm for the getClosest method.
    void __seek_closest(const E_Float *pt, size_type cur_col, size_type axis, E_Float& d2, size_type& m) const;
    
    /// Underneath algorithm for the getClosest method (closest to n).
    void __seek_closest(size_type n, const E_Float *Xn, size_type cur_col, size_type axis, E_Float& d2, size_type& m) const ;
    
    /// Underneath algorithm for the getClose method.
    void __getClosest_through_path(const E_Float* pt, size_type & m, E_Float& d2) const;
    
    /// Underneath algorithm for the getClose method (close to node).
    void __getClosest_through_path(size_type n, const E_Float *Xn, size_type & m, E_Float& d2) const;
    
    /// Underneath algorithm for the getInBox method.
    void __getInBox(size_type ci, size_type axis, const E_Float* mBox, const E_Float* MBox, std::vector<size_type>& out) const;   

  private:

    ///
    const coord_access_type& _posAcc;

    /// kd tree
    tree_array_type	_tree;

    /// current tree size
    size_type       _tree_sz;

    /// dimension (nb of coordinates)
    size_type       _dim;

    /// tolerance
    E_Float         _tolerance;

    /// to extract the coordinates.
    mutable E_Float _Xn[3], _mB[3], _MB[3];

    sortingPredicate<ArrayType> _pred;
    
}; // End class KdTree

} // end namespace

#include "Search/KdTree.cxx"

#endif /* KDTREE_H_ */

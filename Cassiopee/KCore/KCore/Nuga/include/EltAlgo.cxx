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

#ifndef _ELTALGO_CXX_
#define _ELTALGO_CXX_

#include "EltAlgo.h"
#include "Nuga/include/merge.h"
#include "Nuga/include/Triangle.h"
#include "Nuga/include/Polygon.h"
#include "Nuga/include/ArrayAccessor.h"
#include <assert.h>
#include <stack>

///
template <typename ElementType>
template <typename Connectivity_t>
E_Int
NUGA::EltAlgo<ElementType>::getBoundToElements
(const K_FLD::ArrayAccessor<Connectivity_t>& ELTContainer, BoundToEltType& bound_to_elts)
{
  ElementType                     E;
  BoundaryType                    B;
  E_Int*                          pN;
  E_Int                           ROWS(ELTContainer.stride()), COLS(ELTContainer.size()), maxN(0);
  
  pN = new E_Int[ROWS];

  for (E_Int Si = 0; Si < COLS; ++Si)
  {
    ELTContainer.getEntry(Si, pN);
    E.setNodes(pN);

    for (E_Int n = 0; n < ElementType::NB_NODES; ++n)
    {
      E.getBoundary(n, B);
      bound_to_elts[B].push_back(Si);
      maxN = std::max(E_Int(bound_to_elts[B].size()), maxN);
    }
  }
  
  delete [] pN;
  
  return maxN;
}

///
template <typename ElementType>
void
NUGA::EltAlgo<ElementType>::getNeighbours
(const K_FLD::IntArray& ELTContainer, NeighbourType& neighbors)
{
  BoundToEltType                          bound_to_elt;
  typename BoundToEltType::const_iterator itB;
  E_Int                                   i,j;
  size_type                               Si, Sj, sz, COLS(ELTContainer.cols());
  
  neighbors.clear();

  if (COLS == 0)
    return;
  
  assert (ElementType::NB_NODES == ELTContainer.rows());

  // Map the boundary to the connected elements.
  K_FLD::ArrayAccessor<K_FLD::IntArray> actv(ELTContainer);
  getBoundToElements(actv, bound_to_elt);

  neighbors.resize(COLS);
  //
  for (itB = bound_to_elt.begin(); itB != bound_to_elt.end(); ++itB)
  {
    const std::vector<E_Int>& elts = itB->second;
    sz = elts.size();
    for (i = 0; i < sz; ++i)
    {
      Si = elts[i];     
      for (j = i+1; j < sz; ++j)
      {
        Sj = elts[j];
        neighbors[Si].push_back(Sj);
        neighbors[Sj].push_back(Si);
      }
    }
  }
}


///
template <typename ElementType>
inline bool
NUGA::EltAlgo<ElementType>::getManifoldNeighbours
(const K_FLD::IntArray& ELTContainer, K_FLD::IntArray& neighbors, bool strict)
{
  BoundToEltType bound_to_elt;
  typename BoundToEltType::const_iterator itB;
  BoundaryType                            b;
  K_FLD::IntArray::const_iterator         pS;
  E_Int ROWS(ELTContainer.rows()), COLS(ELTContainer.cols()), sz;

  assert (ElementType::NB_NODES == ROWS);
  
  neighbors.clear();
  neighbors.resize(ROWS, COLS, IDX_NONE);
  
  // Map the boundary to the connected elements.
  K_FLD::ArrayAccessor<K_FLD::IntArray> actv(ELTContainer);
  E_Int maxNgh = getBoundToElements(actv, bound_to_elt);

  if (strict && maxNgh != 2) // Non manifold.
    return true;

  for (E_Int i = 0; i < COLS; ++i)
  {
    pS = ELTContainer.col(i);
    for (E_Int n = 0; n < ROWS; ++n)
    {
      ElementType::getBoundary(pS, (n+2)%ROWS, b);
      itB = bound_to_elt.find(b);
      sz = itB->second.size();
      if (sz != 2)
        continue;
      for (E_Int e = 0;  e < sz; ++e)
      {
        if (itB->second[e] != i)
          neighbors((n+2)%ROWS, i) = itB->second[e];
      }
    }
  }
  return false;
}

///
template <typename ElementType>
inline E_Int
NUGA::EltAlgo<ElementType>::getNeighbours (const K_FLD::IntArray& ELTContainer, K_FLD::IntArray& neighbors, bool non_manifold_as_free)
{
  BoundToEltType bound_to_elt;
  typename BoundToEltType::const_iterator itB;
  BoundaryType                            b;
  K_FLD::IntArray::const_iterator         pS;
  E_Int ROWS(ELTContainer.rows()), COLS(ELTContainer.cols()), sz;

  assert (ElementType::NB_NODES == ROWS);
  
  E_Int UNSET = non_manifold_as_free ? IDX_NONE : NON_MANIFOLD_COL;
  
  neighbors.clear();
  neighbors.resize(ROWS, COLS, UNSET);
  
  // Map the boundary to the connected elements.
  K_FLD::ArrayAccessor<K_FLD::IntArray> actv(ELTContainer);
  /*E_Int maxNgh = */getBoundToElements(actv, bound_to_elt);

  for (E_Int i = 0; i < COLS; ++i)
  {
    pS = ELTContainer.col(i);
    for (E_Int n = 0; n < ROWS; ++n)
    {
      ElementType::getBoundary(pS, (n+2)%ROWS, b);
      itB = bound_to_elt.find(b);
      sz = itB->second.size();
      
      if (sz < 2 && !non_manifold_as_free){
        neighbors((n+2)%ROWS, i) = IDX_NONE;
        continue;
      }
      else if (sz > 2) continue;
      
      for (E_Int e = 0;  e < sz; ++e)
      {
        if (itB->second[e] != i)
          neighbors((n+2)%ROWS, i) = itB->second[e];
      }
    }
  }
  return 0;
}

///
template <typename ElementType>
inline void
NUGA::EltAlgo<ElementType>::coloring
(const K_FLD::IntArray& ELTContainer, const NeighbourType& neighbors, const std::set<BoundaryType> & color_bounds, int_vector_type& colors)
{
  size_type                        cols(ELTContainer.cols()), i, sz;
  size_type                        seed(IDX_NONE), S, Sn, color(0), colored(0);
  int_vector_type                  cpool;
  int_vector_type::const_iterator  itC;
  //NeighbourType::const_iterator    itV;
  ElementType                      E, En;
  BoundaryType                     B;
  K_FLD::IntArray::const_iterator  pS, pSn;
  
  colors.clear();
  colors.resize(cols, IDX_NONE);
  
  while (1/*colored < cols*/)
  {
    itC = std::find(colors.begin(), colors.end(), IDX_NONE);
    
    if (itC == colors.end())// Error
      return;

    seed = itC-colors.begin();
    cpool.push_back(seed);

    while (!cpool.empty())
    {
      S = cpool.back();
      cpool.pop_back();

      if (colors[S] != IDX_NONE) // Already colored.
        continue;

      colors[S] = color;
      ++colored;

      assert (S >= 0 && S < (size_type)neighbors.size()); //Error: the neighbors object is not consistent with ELTContainer.

      const std::vector<E_Int>& neigh = neighbors[S];
      sz = neigh.size();

      pS = ELTContainer.col(S);
      E.setNodes(pS);

      for (i = 0; i < sz; ++i)
      {
        Sn = neigh[i];
        
        if (colors[Sn] != IDX_NONE) // Already colored.
          continue;

        pSn = ELTContainer.col(Sn);
        En.setNodes(pSn);

        ElementType::getBoundary(E, En, B);

        if (color_bounds.find(B) != color_bounds.end())
          continue;

        cpool.push_back(Sn);
      }
    }

    ++color;
  }

  //assert ((E_Int)colors.size() == ELTContainer.cols()); //unreachable
}

/// for NGONs
template <typename ElementType>
inline void
NUGA::EltAlgo<ElementType>::coloring (const ngon_unit& neighbors, int_vector_type& colors)
{
  // WARNING : colors are not cleared and suppose that contains coloured elements.

  size_t K, Kn, Kseed(0), color(0), sz, NB_ELTS(neighbors.size());
  int_vector_type                  cpool;
  K_FLD::IntArray neighs;
  
  colors.clear();
  colors.resize(NB_ELTS, IDX_NONE);

  while (1)
  {
    while ((Kseed < NB_ELTS) && (colors[Kseed] != IDX_NONE)) ++Kseed;
    if (NB_ELTS-1 < Kseed)
      return;
    
    cpool.push_back(Kseed);

    while (!cpool.empty())
    {
      K = cpool.back();
      cpool.pop_back();
      
      if (colors[K] != IDX_NONE)
        continue;

      colors[K] = color;

      sz = neighbors.stride(K);
      neighs.reserve(1, sz);
      neighbors.getEntry(K, neighs.begin());
      E_Int* ptr = neighs.begin();
      for (size_t i = 0; i < sz; ++i)
      {
        Kn = *(ptr++);

        if ((Kn != IDX_NONE) && (colors[Kn] == IDX_NONE)) // Not colored.
          cpool.push_back(Kn);
      }
    }

    ++color;
  }
}

// color iff the fontier is of the same color : stop at first inconsistency and return that "bad color"
template <typename ElementType>
template<typename T>
inline bool
NUGA::EltAlgo<ElementType>::coloring_one_connex_homogeneous (const ngon_unit& neighbors, std::vector<T>& colors, size_t Kseed, T UNSET_COL, T color, T FRONT_COL, T& bad_col)
{
  std::vector<T> colors_cpy(colors);//in case of invalid domain (i.e multiple color fronteer)
  bool good_dom = true;
  bool has_boundary_other_than_none = false;
  
  int_vector_type cpool;
  cpool.push_back(Kseed);
  
  K_FLD::IntArray neighs;
  
  while (!cpool.empty() && good_dom)
  {
    E_Int K = cpool.back();
    cpool.pop_back();

    if (colors[K] != UNSET_COL)
      continue;

    colors[K] = color;

    E_Int sz = neighbors.stride(K);
    neighs.reserve(1, sz);
    neighbors.getEntry(K, neighs.begin());
    E_Int* ptr = neighs.begin();
    for (E_Int i = 0; i < sz; ++i)
    {
      E_Int Kn = *(ptr++);
      
      if ( (Kn != IDX_NONE) && (colors[Kn] != UNSET_COL) && (colors[Kn] != color) && (colors[Kn] != FRONT_COL) )
      {
        good_dom = false; //other fronteer than the one expected (FRONT_COL)
        bad_col = colors[Kn];
        break;
      }

      // arbitrary choice to make classifyer::__flag_isolated working.
      // IDX_NONE neighbors are ignored providing there is ate least one elt connected to FRONTCOL.
      // To findout connex island, use another function.
      if ((Kn != IDX_NONE) && (colors[Kn] == FRONT_COL))
        has_boundary_other_than_none = true;

      if ((Kn != IDX_NONE) && (colors[Kn] == UNSET_COL)) // Not colored.
        cpool.push_back(Kn);
    }
  }

  good_dom &= has_boundary_other_than_none;

  if (!good_dom) //revert colors
  {
    colors = colors_cpy;
  }
  return good_dom;
}

// color iff the fontier is of the same color : stop at first inconsistency and return that "bad color"
template <typename ElementType>
template<typename T>
inline bool
NUGA::EltAlgo<ElementType>::coloring_one_connex_homogeneous (const ngon_unit& neighbors, std::vector<T>& colors, size_t Kseed, T UNSET_COL, T color, T FRONT_COL)
{
  E_Int nb_front_col{ 0 };
  
  int_vector_type cpool;
  cpool.push_back(Kseed);
  
  K_FLD::IntArray neighs;
  
  while (!cpool.empty())
  {
    E_Int K = cpool.back();
    cpool.pop_back();

    if (colors[K] != UNSET_COL)
      continue;

    colors[K] = color;

    E_Int sz = neighbors.stride(K);
    neighs.reserve(1, sz);
    neighbors.getEntry(K, neighs.begin());
    E_Int* ptr = neighs.begin();
    for (size_t i = 0; i < sz; ++i)
    {
      E_Int Kn = *(ptr++);

      if (Kn == IDX_NONE) continue;
      
      if ((colors[Kn] != UNSET_COL) && (colors[Kn] != color) && (colors[Kn] != FRONT_COL))
        return false; //good_dom = false; //other fronteer than the one expected (FRONT_COL)

      if (colors[Kn] == UNSET_COL) // Not colored.
        cpool.push_back(Kn);

      if (colors[Kn] == FRONT_COL)
        ++nb_front_col;
    }
  }
  
  //std::cout << "nb of FRONT_COL : " << nb_front_col << std::endl;
  return (nb_front_col > 0); // here either pure none or mixed with front_col : ok iff has front_col
}

// color a connex part surrounded by multiple color frontier
template <typename ElementType>
template<typename T>
inline void
NUGA::EltAlgo<ElementType>::coloring_one_connex_heterogeneous (const ngon_unit& neighbors, std::vector<T>& colors, size_t Kseed, T UNSET_COL, T color)
{
  int_vector_type cpool;
  cpool.push_back(Kseed);
  
  K_FLD::IntArray neighs;
  
  //std::set<E_Int> encountered_cols;

  while (!cpool.empty())
  {
    E_Int K = cpool.back();
    cpool.pop_back();
    
    //encountered_cols.insert(colors[K]);

    if (colors[K] != UNSET_COL)
      continue;

    colors[K] = color;

    E_Int sz = neighbors.stride(K);
    neighs.reserve(1, sz);
    neighbors.getEntry(K, neighs.begin());
    E_Int* ptr = neighs.begin();
    for (E_Int i = 0; i < sz; ++i)
    {
      E_Int Kn = *(ptr++);

      if ((Kn != IDX_NONE) && (colors[Kn] == UNSET_COL)) // Not colored.
        cpool.push_back(Kn);
      //encountered_cols.insert(colors[K]);
    }
  }
  
//  std::cout << "ENCOUNTERED COLORS : " ;
//  for (auto i = encountered_cols.begin(); i != encountered_cols.end(); ++i) std::cout << *i << "/" ;
//  std::cout << std::endl;
}

// color iff the fontier is of the same color : stop at first inconsistency and return that "bad color"
template <typename ElementType>
template<typename T>
inline bool
NUGA::EltAlgo<ElementType>::coloring_one_connex_heterogeneous_if_cols_in_bound(const ngon_unit& neighbors, std::vector<T>& colors, size_t Kseed, T UNSET_COL, T color, E_Int OK_BCOL)
{
  std::vector<T> colors_cpy(colors);//in case of invalid domain (i.e multiple color fronteer)
  
  int_vector_type cpool;
  cpool.push_back(Kseed);

  K_FLD::IntArray neighs;

  std::set<E_Int> boundcols;

  while (!cpool.empty())
  {
    E_Int K = cpool.back();
    cpool.pop_back();

    if (colors[K] != UNSET_COL)
      continue;

    colors[K] = color;

    E_Int sz = neighbors.stride(K);
    neighs.reserve(1, sz);
    neighbors.getEntry(K, neighs.begin());
    E_Int* ptr = neighs.begin();
    for (E_Int i = 0; i < sz; ++i)
    {
      E_Int Kn = *(ptr++);

      if ((Kn != IDX_NONE) && (colors[Kn] == UNSET_COL)) // Not colored.
        cpool.push_back(Kn);
      else if (Kn == IDX_NONE)
        boundcols.insert(IDX_NONE);
      else if (colors[Kn] != UNSET_COL)
        boundcols.insert(colors[Kn]);
    }
  }

  bool good_dom = (boundcols.find(OK_BCOL) != boundcols.end());
  return good_dom;
}

template <typename ElementType>
template<typename T>
inline void
NUGA::EltAlgo<ElementType>::coloring (const ngon_unit& neighbors, std::vector<T>& colors, T UNSET_COL, T FIRST_COL)
{
  // WARNING : colors are not cleared and suppose that contains coloured elements.

  size_t Kseed(0), NB_ELTS(neighbors.size());
  T color(FIRST_COL);
  
  assert (colors.size() >= NB_ELTS);
  
  while (1)
  {
    while ((Kseed < NB_ELTS) && (colors[Kseed] != UNSET_COL)) ++Kseed;
    if (NB_ELTS-1 < Kseed)
      return;
    
    coloring_one_connex_heterogeneous (neighbors, colors, Kseed, UNSET_COL, color++);
  }
}

template <typename ElementType>
template<typename T>
inline void
NUGA::EltAlgo<ElementType>::coloring_homogeneous_zones_only (const ngon_unit& neighbors, std::vector<T>& colors, T UNSET_COL, T FIRST_COL, T FRONT_COL)
{
  // WARNING : colors are not cleared and suppose that contains coloured elements.

  size_t Kseed(0), NB_ELTS(neighbors.size());
  T color(FIRST_COL);
  
  assert (colors.size() >= NB_ELTS);
  
  std::set<E_Int> bad_colors;
  
  while (1)
  {
    while ((Kseed < NB_ELTS) && (colors[Kseed] != UNSET_COL)) ++Kseed;
    if (NB_ELTS-1 < Kseed) break; //nothing to color anymore

    bool good_dom = coloring_one_connex_homogeneous (neighbors, colors, Kseed, UNSET_COL, color++, FRONT_COL);
    if (!good_dom) bad_colors.insert(color-1);//tells if color is ok or not (fully surrounded by the same color)
  }
  
  // reset any bad color
  if (!bad_colors.empty())
  {
    for (size_t i=0; i < colors.size(); ++i)
    {
      if (bad_colors.find(colors[i]) != bad_colors.end())colors[i] = UNSET_COL;
    }
  }
}

template <typename ElementType>
template<typename T>
inline void
NUGA::EltAlgo<ElementType>::coloring_attached_to_col_zones_only(const ngon_unit& neighbors, std::vector<T>& colors, T UNSET_COL, T FIRST_COL, T OK_BCOL)
{
  // WARNING : colors are not cleared and suppose that contains coloured elements.

  size_t Kseed(0), NB_ELTS(neighbors.size());
  T color(FIRST_COL);

  assert(colors.size() >= NB_ELTS);

  std::set<E_Int> bad_colors;

  while (1)
  {
    while ((Kseed < NB_ELTS) && (colors[Kseed] != UNSET_COL)) ++Kseed;
    if (NB_ELTS - 1 < Kseed) break; //nothing to color anymore

    bool good_dom = coloring_one_connex_heterogeneous_if_cols_in_bound(neighbors, colors, Kseed, UNSET_COL, color++, OK_BCOL);
    if (!good_dom) bad_colors.insert(color - 1);//tells if color is ok or not (fully surrounded by the same color)
  }

  // reset any bad color
  if (!bad_colors.empty())
  {
    for (size_t i = 0; i < colors.size(); ++i)
    {
      if (bad_colors.find(colors[i]) != bad_colors.end())colors[i] = UNSET_COL;
    }
  }
}

template <typename ElementType>
inline void
NUGA::EltAlgo<ElementType>::coloring (const K_FLD::IntArray& neighbors, int_vector_type& colors)
{
  // WARNING : colors are not cleared and suppose that contains coloured elements.

  size_t K, Kseed(0), color(0), NB_ELTS(neighbors.cols()), ROWS(neighbors.rows());
  int_vector_type                  cpool;
  K_FLD::IntArray::const_iterator pK;
  
  assert (ROWS == (size_t)ElementType::NB_NODES);
  
  colors.clear();
  colors.resize(NB_ELTS, IDX_NONE);

  while (1)
  {
    while ((Kseed < NB_ELTS) && (colors[Kseed] != IDX_NONE)) ++Kseed;
    if (NB_ELTS-1 < Kseed)
      return;
    
    cpool.push_back(Kseed);

    while (!cpool.empty())
    {
      K = cpool.back();
      cpool.pop_back();
      
      if (colors[K] != IDX_NONE)
        continue;

      colors[K] = color;

      pK = neighbors.col(K);
      for (size_t i = 0; i < ROWS; ++i)
      {
        const E_Int& Kn = *(pK+i);

        if ((Kn != IDX_NONE) && (colors[Kn] == IDX_NONE)) // Not colored.
          cpool.push_back(Kn);
      }
    }

    ++color;
  }
}

///
template <typename ElementType>
void
NUGA::EltAlgo<ElementType>::coloring_pure (const K_FLD::IntArray& neighbors, int_vector_type& colors)
{
  // WARNING : colors are not cleared and suppose that contains coloured elements.

  size_t stride (neighbors.rows()), K, Kseed(0),  color(0), NB_ELTS(neighbors.cols());
  
  colors.clear();
  if (NB_ELTS == 0) return;

  int_vector_type                  cpool;
  K_FLD::IntArray::const_iterator pK;
  
  colors.resize(NB_ELTS, IDX_NONE);

  while (1)
  {
    while ((Kseed < NB_ELTS) && (colors[Kseed] != IDX_NONE)) ++Kseed;
    if (NB_ELTS-1 < Kseed)
      return;

    cpool.push_back(Kseed);

    while (!cpool.empty())
    {
      K = cpool.back();
      cpool.pop_back();
      
      if (colors[K] != IDX_NONE)
        continue;

      colors[K] = color;
      
      pK = neighbors.col(K);
      for (size_t i=0; i < stride; ++i)
      {
        const E_Int& Kn = *(pK+i);
        if ((Kn != IDX_NONE) && (colors[Kn] == IDX_NONE)) // Not colored.
          cpool.push_back(Kn);
      }
    }

    ++color;
  }
}

/// for NGONs
template <typename ElementType>
inline void
NUGA::EltAlgo<ElementType>::extrapolate (const ngon_unit& neighbors, K_FLD::IntArray& properties)
{
  // WARNING : colors are not cleared and suppose that contains coloured elements.

  K_FLD::IntArray neighs;
  
  // Count the number of empty color
  E_Int count(0), NB_ELTS(properties.cols()), ROWS(properties.rows());
  for (E_Int i=0; i < NB_ELTS; ++i) 
    for (E_Int p=0; p < ROWS; ++p)
      if (properties(p,i) == IDX_NONE) ++count;

  while (count)
  {
    for (E_Int i=0; i < NB_ELTS; ++i)
    {
      E_Int nb_neighs = neighbors.stride(i);
      neighs.reserve(1, nb_neighs);
      neighbors.getEntry(i, neighs.begin());
      E_Int* ptr = neighs.begin();
      
      for (E_Int p=0; p < ROWS; ++p)
      {
        if (properties(p,i) != IDX_NONE) continue;

        for (E_Int n = 0; n < nb_neighs; ++n)
        {
          E_Int j = *(ptr++);
          if (j < IDX_NONE && properties(p,j) != IDX_NONE)
          {
            properties(p,i)=properties(p,j);
            --count;
            break;
          }
        }
      }
    } 
  }
}

///
template <>
template <typename Connectivity_t>
void
NUGA::EltAlgo<K_MESH::Triangle>::reversi_connex (const Connectivity_t& cnt, const Connectivity_t& neighbors, E_Int Kseed, int_vector_type& orient)
{

  if (cnt.cols()==0)
    return;
  if (Kseed == IDX_NONE)
    return;

  int_vector_type                  cpool, processed;
  E_Int K, Kngh, i, j, Ni, Nj;
  K_FLD::IntArray::const_iterator pNeigh;
  E_Bool reverse;
  K_FLD::ArrayAccessor<Connectivity_t > conn(cnt);
  
  orient.resize(cnt.cols(), 1);//not cleared on purpose : can be called several time to append orient..
  
  cpool.push_back(Kseed);
  processed.resize(cnt.cols(), false);
  
  while (!cpool.empty())
  {
    K = cpool.back();
    cpool.pop_back();
    processed[K]=true;
    
    K_MESH::Triangle Tk(conn, K);
    
    if (orient[K]==-1)
      Tk.reverse();
    
    pNeigh = neighbors.col(K); // neighbor ptr;
    
    for (size_t k = 0; k < 3; ++k)
    {
      Kngh = *(pNeigh++); 
      
      if ((Kngh == IDX_NONE) || processed[Kngh])
        continue;
      
      K_MESH::Triangle Tkngh(conn, Kngh);
      K_MESH::Triangle::getBoundary(Tk, Tkngh, i, j); //get the shared edge local position
      
      Ni=Tk.node((i+1)%3);//*(pK+(i+1)%3);
      Nj=Tk.node((i+2)%3);//*(pK+(i+2)%3);
      K_MESH::Triangle::getOrientation(Tkngh, Ni, Nj, reverse);
      
      if (reverse)
        orient[Kngh] = -1;
      /*{
        std::swap(*(pKngh+1), *(pKngh+2));   //reorient
        std::swap(*(pNeigh+1), *(pNeigh+2)); //update neighboring consistently.
      }*/
      
      cpool.push_back(Kngh);
      
    }
  }
}

///
template <>
template <typename Connectivity_t>
void
NUGA::EltAlgo<K_MESH::Polygon>::reversi_connex
(const Connectivity_t& PGs, const Connectivity_t& neighbors, E_Int Kseed, int_vector_type& orient)
{
  int_vector_type                  cpool, processed;
  E_Int K, Kngh, i, j, Ni, Nj, s1, s3, nb_neighs;
  K_FLD::IntArray::const_iterator pKngh;
  K_MESH::NO_Edge noE;
  E_Bool reverse;

  typedef K_MESH::Polygon ElementType;
  
  orient.resize(PGs.size(), 1);//not cleared on purpose : can be called several time to append orient..
  
  cpool.push_back(Kseed);
  processed.resize(PGs.size(), false);
  
  while (!cpool.empty())
  {
    K = cpool.back();
    cpool.pop_back();
    processed[K]=true;
    
    s1 = PGs.stride(K);
    ElementType PGk(PGs.get_facets_ptr(K), s1);
    
    nb_neighs = neighbors.stride(K); // neighbor ptr;
    
    for (E_Int n = 0; n < nb_neighs; ++n)
    {
      Kngh = neighbors.get_facet(K,n); 
      
      if ((Kngh == IDX_NONE) || processed[Kngh])
        continue;
      
      pKngh = PGs.get_facets_ptr(Kngh);
      s3 = PGs.stride(Kngh);
      ElementType PGkngh(pKngh, s3);
      ElementType::getBoundary(PGk, PGkngh, i,j); //get the shared edge
      
      Ni=PGk.node(i);
      Nj = PGk.node((i+1)%s1);
      
      ElementType::getOrientation(PGkngh, Ni, Nj, reverse);

      if (orient[K]==-1) reverse = !reverse;
      
      if (reverse) orient[Kngh] = -1;
      
      cpool.push_back(Kngh);    
    }
  }
}

template <>
template <typename Connectivity_t>
E_Int
NUGA::EltAlgo<K_MESH::Polygon>::check_orientation_consistency
(const Connectivity_t& PGs, const Connectivity_t& neighbors, std::set<E_Int>& faultys)
{
  typedef K_MESH::Polygon ElementType;

  int_vector_type                  cpool, processed;
  E_Int Kngh, i, j, Ni, Nj, s1, s3, nb_neighs;
  K_FLD::IntArray::const_iterator pKngh;
  K_MESH::NO_Edge noE;
  E_Bool reverse;
  typedef K_MESH::Polygon ElementType;

  for (E_Int K = 0; K < PGs.size(); ++K)
  {
    s1 = PGs.stride(K);
    ElementType PGk(PGs.get_facets_ptr(K), s1);
    
    nb_neighs = neighbors.stride(K); // neighbor ptr;
    
    for (E_Int n = 0; n < nb_neighs; ++n)
    {
      Kngh = neighbors.get_facet(K,n); 
            
      if (Kngh == IDX_NONE)
        continue;
      
      pKngh = PGs.get_facets_ptr(Kngh);
      s3 = PGs.stride(Kngh);
      ElementType PGkngh(pKngh, s3);
      ElementType::getBoundary(PGk, PGkngh, i,j); //get the shared edge
      
      Ni=PGk.node(i);
      Nj = PGk.node((i+1)%s1);
      
      ElementType::getOrientation(PGkngh, Ni, Nj, reverse);
      
      if (reverse)
      {
        faultys.insert(K);
        faultys.insert(Kngh);
        //return true;
      }
    }
  }
  return !faultys.empty();
}

namespace NUGA // gcc complaint
{
///
template <>
inline E_Int
EltAlgo<K_MESH::Triangle>::fast_swap_edges
(std::vector<std::pair<E_Int, E_Int> >& swapE, K_FLD::IntArray& connect, K_FLD::IntArray& neighbors)
{
  E_Int S, b, Sn, bn, Ni, Nk, S1, S2, b1, b2;
  K_FLD::IntArray::iterator pS, pSn;
      
  for (size_t i=0; i < swapE.size(); ++i)
  {
    std::pair<E_Int, E_Int>& E = swapE[i];

    S = E.first;
    b = E.second;

    pS = connect.col(S);
    Ni = *(pS + b);
    //Nj = *(pS + (b+1) % 3);
    //Nl = *(pS + (b+2) % 3);

    Sn = neighbors(b, S);
    if (Sn == IDX_NONE)
      continue;

    pSn = connect.col(Sn);
    bn = K_MESH::Triangle::getOppLocalNodeId(S, b, connect, neighbors);

    Nk = *(pSn + bn);

    // Neighbors to modify
    S1 = neighbors((b+1) % 3, S);
    S2 = neighbors((bn+1) % 3, Sn);
    b1 = K_MESH::Triangle::getOppLocalNodeId(S, (b+1) % 3, connect, neighbors);
    b2 = K_MESH::Triangle::getOppLocalNodeId(Sn, (bn+1) % 3, connect, neighbors);

    // Update elements S and Sn (connect)
    *(pS  + (b+2)  % 3) = Nk;
    *(pSn + (bn+2) % 3) = Ni;
    
    // Update the neighboring (neighbors)
    neighbors((b+1)  % 3, S)  = Sn;
    neighbors((bn+1) % 3, Sn) = S;
    if ((S1 != IDX_NONE) && (b1 != IDX_NONE))
      neighbors(b1, S1)              = Sn;
    neighbors(bn, Sn)                = S1;
    if ((S2 != IDX_NONE) && (b2 != IDX_NONE))
      neighbors(b2, S2)              = S;
    neighbors(b, S)                  = S2;
  }

  return 0;
}
}

template <typename ElementType>
inline void
NUGA::EltAlgo<ElementType>::smoothd1(const ngon_unit& PHs, const K_FLD::IntArray& noF2E, const std::vector<E_Int>& PHlevel, std::vector<E_Int>& adap_incr)
{
  std::stack<E_Int> stck;

  for (size_t i = 0; i < PHs.size();  i++){
    if (adap_incr[i] != 0){
      stck.push(i);
    }
  }

  while (!stck.empty()){

    E_Int ind_PHi = stck.top(); // index of ith PH
    stck.pop();

    E_Int nb_edges = PHs.stride(ind_PHi); // number of edges of ith PH
    const E_Int* pPHi = PHs.get_facets_ptr(ind_PHi); // find the PG of the ith PH

    for (E_Int j = 0; j < nb_edges; j++){ 

      E_Int ind_PGi = *(pPHi+j)-1; // index of jth PG
      
      E_Int PHG = noF2E(0,ind_PGi);
      E_Int PHD = noF2E(1,ind_PGi);

      if ( (PHD == IDX_NONE) || (PHG == IDX_NONE) )
        continue; // external PG

     
      E_Int incrG = adap_incr[PHG] + PHlevel[PHG];
      E_Int incrD = adap_incr[PHD] + PHlevel[PHD];
      
      if (abs(incrG-incrD) <= 1) // 2:1 rule
        continue;
      
      if (incrG > incrD){
        adap_incr[PHD] += 1;
        stck.push(PHD);
      }
      else{
        adap_incr[PHG] += 1;
        stck.push(PHG);
      }
    }
  }
}

#endif

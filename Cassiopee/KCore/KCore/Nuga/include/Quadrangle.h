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

#ifndef __K_MESH_QUADRANGLE_H__
#define __K_MESH_QUADRANGLE_H__

#include "Nuga/include/defs.h"
#include "Nuga/include/DynArray.h"
#include "Nuga/include/Edge.h"
#include "Nuga/include/IdTool.h"

namespace K_MESH
{

  class Quadrangle
  {
  public:

    using nob_t = K_MESH::NO_Edge;
    
    static constexpr E_Int NB_NODES = 4;
    static constexpr E_Int NB_TRIS = 2;
    
    Quadrangle(){};
    
    template <typename NodeIterator>
    explicit Quadrangle(const NodeIterator pN, E_Int idx_start = 0)
    {_nodes[0] = *pN - idx_start; _nodes[1] = *(pN+1) - idx_start; _nodes[2] = *(pN+2) - idx_start; _nodes[3] = *(pN+3) - idx_start;}

    void setNodes(E_Int n0, E_Int n1, E_Int n2, E_Int n3){_nodes[0]=n0;_nodes[1]=n1;_nodes[2]=n2;_nodes[3]=n3;}
    //void setNodes(E_Int* nodes){for (size_t i = 0; i< 4; ++i)_nodes[i]=*(nodes++);}
       
    E_Int nb_nodes() const {return NB_NODES;}
    E_Int nb_tris() const { return NB_TRIS;}
    E_Int nbounds() const { return 4;}
    
    ///
    inline E_Int* nodes(){return _nodes;} 
    inline const E_Int* nodes() const {return _nodes;}
    
    /// Gets the i-th node.
    inline const E_Int& node(const E_Int& i) const {return _nodes[i];}
    
    template <typename TriangulatorType, typename acrd_t>
    void triangulate (const TriangulatorType& dt, const acrd_t& crd){}
    
    inline void triangle(E_Int i, E_Int* target)
    {
      assert (i >= 0 && i < 2);
      
      if (i==0)
      {
        target[0]=_nodes[0]; target[1]=_nodes[1]; target[2]=_nodes[2];
      }
      else
      {
        target[0]=_nodes[0]; target[1]=_nodes[2]; target[2]=_nodes[3];
      }
    }
    
    static inline void triangulate(const E_Int* nodes, E_Int* target)//WARNING : triangulation is put at target address contiguously
    {
      E_Int t[3];
      t[0] = nodes[0]; t[1] = nodes[1]; t[2] = nodes[3];
      std::copy(&t[0], &t[0]+3, target);
      std::copy(&nodes[0]+1, &nodes[0]+4, target+3);
    }
    
    template<typename box_t, typename CoordAcc>
    void bbox(const CoordAcc& acrd, box_t&bb) const
    {
      for (E_Int i = 0; i < 3; ++i)
      {bb.minB[i] = NUGA::FLOAT_MAX; bb.maxB[i] = -NUGA::FLOAT_MAX;}

      box_t b;
      b.compute(acrd, _nodes, NB_NODES, 0/*idx start*/);
      for (E_Int i = 0; i < NB_NODES; ++i)
      {
        bb.minB[i] = std::min(bb.minB[i], b.minB[i]);
        bb.maxB[i] = std::max(bb.maxB[i], b.maxB[i]);
      }
    }
    
    template <typename CoordAcc>
    void iso_barycenter(const CoordAcc& coord, E_Float* G)
    { 
      //
      for (size_t d=0; d < 3; ++d) G[d]=0.;

      for (E_Int i=0; i < NB_NODES; ++i)
      {
        for (size_t d=0; d < 3; ++d)
        {
          //std::cout << "v : " << coord.getVal(node(i), d) << std::endl;
          G[d] += coord.getVal(_nodes[i], d);
        }
      }
    }

    bool operator==(const K_MESH::Quadrangle& q) const
    {
      return K_CONNECT::IdTool::equal(this->nodes(),q.nodes(),4,true,false);
    }

    static bool need_a_reorient(E_Int PGi, E_Int PHi, bool oriented_if_R, const K_FLD::IntArray & F2E)
    {
      if (F2E(1, PGi) == PHi && oriented_if_R == true) return false;
      else if (F2E(0, PGi) == PHi && oriented_if_R == false) return false;
      else return true;
    }

    static void order_nodes(E_Int *nodes, const E_Int *pN, bool reorient, E_Int i0)
    {
      for (int i = 0; i < 4; i++)
        nodes[i] = pN[i];

      K_CONNECT::IdTool::right_shift<4>(nodes, i0);

      if (reorient) std::swap(nodes[1], nodes[3]);
    }
    
  private:
    E_Int _nodes[4];
  };
}

static uint32_t fcValue(const E_Int *nodes, int len, int pivot)
{
  return (pivot == len-1) ? nodes[0] : nodes[pivot+1];
}

static
int fcIndex(int len, int pivot)
{
  return (pivot == len-1) ? 0 : pivot+1;
}

static
uint32_t rcValue(const E_Int *nodes, int len, int pivot)
{
  return (pivot == 0) ? nodes[len-1] : nodes[pivot-1];
}

static
int rcIndex(int len, int pivot)
{
  return (pivot == 0) ? len-1 : pivot-1;
}

struct Quadrangle_Hash
{
  uint32_t hash(uint32_t val, uint32_t seed) const
  {
    uint32_t HASH = seed;
    HASH += val;
    HASH += HASH << 10;
    HASH ^= HASH >> 6;
    return HASH;
  }

  uint32_t operator()(const K_MESH::Quadrangle& q) const
  {
    const auto& nodes = q.nodes();
    
    // find location of min vertex
    int pivot = 0;
    int len = 4;
    for (int i = 1; i < len; ++i) {
      if (nodes[pivot] > nodes[i]) {
        pivot = i;
      }
    }

    uint32_t seed = 0;

    if (fcValue(nodes, len, pivot) < rcValue(nodes, len, pivot))
    {
      // Forward circulate
      while (len--)
      {
          seed = hash(nodes[pivot], seed);
          pivot = fcIndex(4, pivot);
      }
    }
    else
    {
      // Reverse circulate
      while (len--)
      {
          seed = hash(nodes[pivot], seed);
          pivot = rcIndex(4, pivot);
      }
    }
 
    seed += seed << 3;
    seed ^= seed >> 11;
    seed += seed << 15;
    return seed;
  }
};

#endif

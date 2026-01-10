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

#ifndef _DELAUNAY_MESHDATA_H_
#define _DELAUNAY_MESHDATA_H_

#include "Nuga/include/defs.h"
#include "Nuga/include/DefContainers.h"

namespace DELAUNAY
{
  class MeshData
  {
  public:
    typedef E_Int size_type;
    typedef NUGA::int_vector_type  int_vector_type;
    typedef NUGA::bool_vector_type bool_vector_type;
    
    MeshData():pos(0), connectB(0), unsync_nodes(false), mono_connex(false){hnids.clear();};

    MeshData(K_FLD::FloatArray& p, const K_FLD::IntArray& cB):pos(&p),connectB(&cB), unsync_nodes(false)
    {
      hnids.clear();
    }

    void set(K_FLD::FloatArray& p, const K_FLD::IntArray& cB)
    {
      clear();
      pos = &p;
      connectB = &cB;
    }
    
    void clear()
    {
      hardNodes.clear();
      hardEdges.clear();
      connectM.clear();
      neighbors.clear();
      ancestors.clear();
      colors.clear();
      metrics.clear();
      mask.clear();
      hnids.clear();
      unsync_nodes = false;
      mono_connex = false;
    }
    
    ///
    void sync_hards()
    {
      std::vector<E_Int> hN;
      hN.reserve(hardNodes.size());
      for (auto& Ni : hardNodes)
      {
        if (hnids[Ni] == Ni)hN.push_back(Ni);
      }
      hardNodes=hN;
    
      NUGA::non_oriented_edge_set_type nHE;// = data.hardEdges;
      for (auto& Ei : hardEdges)
      { 
        E_Int Ni = Ei.node(0);
        Ni = (hnids[Ni] == Ni ) ? Ni : hnids[Ni];
        E_Int Nj = Ei.node(1);
        Nj = (hnids[Nj] == Nj ) ? Nj : hnids[Nj];
      
        if (Ni != Nj) nHE.insert(K_MESH::NO_Edge(Ni, Nj));
      }
      hardEdges = nHE;
    }

    K_FLD::FloatArray* pos;
    const K_FLD::IntArray* connectB;
    int_vector_type hardNodes;
    bool unsync_nodes;
    int_vector_type hnids;
    bool mono_connex;
  
    NUGA::non_oriented_edge_set_type hardEdges;
    K_FLD::IntArray connectM;
    K_FLD::IntArray neighbors;
    int_vector_type ancestors;
    int_vector_type colors;
    K_FLD::FloatArray metrics;
    bool_vector_type mask;
  };

  template <typename SurfaceType>
  struct SurfaceMeshData : public MeshData
  {
    SurfaceMeshData(K_FLD::FloatArray& p2D, const K_FLD::FloatArray& p3D, const K_FLD::IntArray& cB, const SurfaceType& s)
      : MeshData(p2D, cB), surface(s), pos3D(p3D), exportUV(false)
    {
    }

    const SurfaceType& surface;
    K_FLD::FloatArray pos3D;
    E_Bool exportUV;
  };
}

#endif


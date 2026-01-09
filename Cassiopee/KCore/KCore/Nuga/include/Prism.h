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

#ifndef NUGA_PRISM_H
#define NUGA_PRISM_H
#include "Nuga/include/defs.h"
#include "Nuga/include/DynArray.h"
#include "Nuga/include/Triangle.h"
#include "Nuga/include/ngon_t.hxx"
#include "Nuga/include/ArrayAccessor.h"

namespace K_MESH
{
class Prism {
public:
    static constexpr E_Int NB_BOUNDS=5;

    Prism(){}
    Prism(const Prism& orig){}

    template <typename ngunit_t>
    Prism(const ngunit_t & PGs, const E_Int* first_pg){}

    ~Prism();

    template< typename ngo_t>
    static void reorder_pgs(ngo_t& ng, const K_FLD::IntArray& F2E, E_Int i);

    template <typename ngunit_t>
    static E_Int get_opposite(const ngunit_t & PGs, const E_Int* first_pg, E_Int k)
    {
      E_Int PGk  = first_pg[k]-1;
      if (PGs.stride(PGk) != 3) return IDX_NONE;//means something only for 

      for (size_t i=0; i < 5; ++i)
      {
        if (i == (size_t)k) continue;
        int PGi  = first_pg[i]-1;
        if (PGs.stride(PGi) == 3) return i;
      }
      return IDX_NONE;
    }

    template <typename ngunit_t>
    static bool is_of_type(const ngunit_t & PGs, const E_Int* first_pg, E_Int nb_pgs)
    {
      if (nb_pgs != 5) return false;

      E_Int s1(0), s2(0);

      for (int i = 0; i<5; i++)
      {
        if (PGs.stride(*(first_pg + i) - 1) == 3) ++s1;
        else if (PGs.stride(*(first_pg + i) - 1) == 4) ++s2;
        else return false;
      }

      return ((s1 == 2) && (s2 == 3));
    }
    ///
    template <typename ngunit_t>
    static inline void iso_barycenter(const K_FLD::FloatArray& crd, const ngunit_t & PGs, const E_Int* first_pg, E_Int nb_pgs, E_Int index_start, E_Float* G);

    E_Float quality(const K_FLD::FloatArray& crd, E_Float* Vol){return 1;}

private:

};

template< typename ngo_t>
void Prism::reorder_pgs(ngo_t& ng, const K_FLD::IntArray& F2E, E_Int i)
{
  //std::cout << "reorder_pgs" << std::endl;  
  std::map<E_Int,E_Int> glmap; // crd1 to 0-26 indexes
  E_Int nb_faces = ng.PHs.stride(i); 
  E_Int* faces = ng.PHs.get_facets_ptr(i);
  E_Int PGi = faces[0] - 1;
  // but convention, first face is bottom, first node is 0 in local numbering (0 to 26)

  E_Int l(0);
  for (int i=0; i< nb_faces; i++){
    if (ng.PGs.stride(faces[i]-1)!=4){
        PGi= faces[i] - 1;
        l=i;
        break;
    }
  }

//  std::cout << "after break" << std::endl;
//  std::cout << "l= " << l << std::endl;
//    for (int l1=0; l1<ng.PHs.stride(i); l1++){
//        E_Int* face = ng.PHs.get_facets_ptr(i);
//        for (int m=0; m<ng.PGs.stride(face[l1]-1); m++){
//            E_Int* pN= ng.PGs.get_facets_ptr(face[l1]-1);
//            std::cout << "pN= " << pN[m]-1 << std::endl;
//        }
//        std::cout << "/////////////" << std::endl;
//    }        
//  std::cout << "l= " << l << std::endl;
//  std::cout << "stride face BOT= " << ng.PGs.stride(PGi) << std::endl;
  

  E_Int* pN = ng.PGs.get_facets_ptr(PGi);

  
  glmap[*pN] = 0; // PHi(0,0) -> 0  
  glmap[*(pN+1)] = 1;
  glmap[*(pN+2)] = 2;

  //std::cout << "F2E = " << F2E(1,PGi) << std::endl;
  if (F2E(1,PGi) != i) // for BOT, PH is the right element. if not, wrong orientation => swap of 1 and 3
  { 
    glmap[*(pN+1)] = 2;
    glmap[*(pN+2)] = 1;
  }
//  for (int i=0; i<ng.PGs.stride(PGi); i++){
//      std::cout << "pN= " << *(pN+i)-1;
//      std::cout << "         it= " << glmap[*(pN+i)] << std::endl;
//  }
  E_Int F1Id(IDX_NONE), F2Id(IDX_NONE), F3Id(IDX_NONE), TOPId(IDX_NONE);

  bool commonNodes[3];

  for (int k = 1; k < nb_faces; ++k)
  {
    int count = 0;
    commonNodes[0] = commonNodes[1] = commonNodes[2] = false;
    int ki = (k+l) % nb_faces;

    E_Int testedPG = faces[ki]-1;

    E_Int* pNode = ng.PGs.get_facets_ptr(testedPG);
    E_Int nb_nodes = ng.PGs.stride(testedPG);
//    std::cout << "testedPG= " << testedPG << "stride= " << ng.PGs.stride(testedPG) << std::endl; 
    for (int j = 0; j < nb_nodes; ++j)
    {
      auto it = glmap.find(pNode[j]);
//      std::cout << "pNode= " << pNode[j];
//      std::cout << "         it= " << it->first << "    local= "<< it->second << std::endl;
      if (it != glmap.end()){
        // found
        count++;
        commonNodes[it->second] = true;
      }
    }
    if (commonNodes[0] && commonNodes[1])
      F1Id = ki;
    else if (commonNodes[1] && commonNodes[2])
      F2Id = ki;
    else if (commonNodes[2] && commonNodes[0])
      F3Id = ki;
    else if (count == 0)
      TOPId = ki;   
    }

  E_Int mol[5];

  mol[0] = faces[l];
  mol[1] = faces[F1Id];
  mol[2] = faces[F2Id];
  mol[3] = faces[F3Id];
  mol[4] = faces[TOPId];

  for (int i = 0; i < nb_faces; ++i)
    faces[i] = mol[i];
   
  //
  assert(ng.PGs.stride(faces[0]-1) == 3);
  assert(ng.PGs.stride(faces[1]-1) == 4);
  assert(ng.PGs.stride(faces[2]-1) == 4);
  assert(ng.PGs.stride(faces[3]-1) == 4);
  assert(ng.PGs.stride(faces[4]-1) == 3);

  assert (l != IDX_NONE && l != F1Id && l != F2Id && l != F3Id && l != TOPId);
  assert (F1Id != IDX_NONE && F1Id != l && F1Id != F2Id && F1Id != F3Id && F1Id != TOPId);
  assert (F2Id != IDX_NONE && F2Id != l && F2Id != F1Id && F2Id != F3Id && F2Id != TOPId);
  assert (F3Id != IDX_NONE && F3Id != l && F3Id != F1Id && F3Id != F2Id && F1Id != TOPId);
  assert (TOPId != IDX_NONE && TOPId != l && TOPId != F1Id && TOPId != F2Id && TOPId != F3Id);
  
//  std::cout <<"after reorder" << std::endl;
//  std::cout <<"/////////////" << std::endl;
//    for (int l1=0; l1<ng.PHs.stride(i); l1++){
//      E_Int* face = ng.PHs.get_facets_ptr(i);
//      for (int m=0; m<ng.PGs.stride(face[l1]-1); m++){
//          E_Int* pN= ng.PGs.get_facets_ptr(face[l1]-1);
//          std::cout << "pN= " << pN[m]-1 << std::endl;
//      }
//      std::cout << "/////////////" << std::endl;
//  }    
}

  template <typename ngunit_t>
  inline void Prism::iso_barycenter(const K_FLD::FloatArray& crd, const ngunit_t & PGs, const E_Int* first_pg, E_Int nb_pgs, E_Int index_start, E_Float* G)
  {
    //WARNING : assuming reodrederd pgs : first is bottom, second is top
#ifdef DEBUG_2019        
    std::cout << "iso_barycenter prism" << std::endl;
#endif
    E_Int new_bary[6];

    for (int i = 0; i < 2; ++i) // 6 points : bot and top nodes
    {
      const E_Int* nodes = PGs.get_facets_ptr(first_pg[i+3*i]-index_start);
#ifdef DEBUG_2019  
      E_Int nb_nodes = PGs.stride(first_pg[i+3*i]-index_start);
      assert (nb_nodes == 3); //fixme : checkme
#endif
      for (int k = 0; k  < 3; ++k)
        new_bary[3*i+k] = nodes[k];   
    }
    
    K_MESH::Polyhedron<STAR_SHAPED>::iso_barycenter(crd, new_bary, 6, 1, G);
  }

}
#endif /* PRISMHEDRON_H */


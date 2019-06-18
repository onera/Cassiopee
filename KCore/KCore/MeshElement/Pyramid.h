/*
 
 
 
              NUGA 
 
 
 
 */

#ifndef __K_MESH_PYRAMID_H__
#define __K_MESH_PYRAMID_H__

#include "Def/DefTypes.h"
#include "Fld/DynArray.h"
#include "MeshElement/Triangle.h"
#include "Fld/ngon_t.hxx"
#include "Fld/ArrayAccessor.h"


namespace K_MESH
{

class Pyramid {

  public:
    static const E_Int NB_NODES = 5;
    static const E_Int NB_TRIS = 6;
    static const E_Int NB_BOUNDS = 5;
      
  public:
    Pyramid():_shift(0){}
    ~Pyramid(){}
    
    Pyramid(const E_Int* nodes, E_Int shift=0):_shift(shift){ for (size_t i = 0; i< NB_NODES; ++i)_nodes[i]=*(nodes++) + shift;}
    
    inline E_Int node(E_Int i){return _nodes[i]+_shift;}
    
    E_Int* nodes() { return _nodes;}
    // const E_Int* nodes() const { return _nodes;}
    
    // E_Int nb_nodes() const {return NB_NODES;}
    E_Int nb_tris() const {return NB_TRIS;}
    
    //void setNodes(E_Int* nodes){for (size_t i = 0; i< 4; ++i)_nodes[i]=*(nodes++);}
    
    template <typename CoordAcc> inline
    void iso_barycenter(const CoordAcc& coord, E_Float* G);
    template <typename ngunit_t> inline
    static void iso_barycenter(const K_FLD::FloatArray& crd, const ngunit_t & PGs, const E_Int* first_pg, E_Int nb_pgs, E_Int index_start, E_Float* G);
    
    template< typename ngo_t>
    static void reorder_pgs(ngo_t& ng, const K_FLD::IntArray& F2E, E_Int i);
    
    
    
    ///
    template <typename TriangulatorType, typename acrd_t>
    void triangulate (const TriangulatorType& dt, const acrd_t& acrd) {} // dummy since it is a basic element
        
    inline void triangle(E_Int i, E_Int* target)
    {
      assert (i >= 0 && i < NB_TRIS);
      
      switch (i)
      {
        case 0 : target[0] = _nodes[0]; target[1] = _nodes[1]; target[2] = _nodes[4]; break;  //014
        case 1 : target[0] = _nodes[1]; target[1] = _nodes[2]; target[2] = _nodes[4]; break;  //124
        case 2 : target[0] = _nodes[2]; target[1] = _nodes[3]; target[2] = _nodes[4]; break;  //234
        case 3 : target[0] = _nodes[3]; target[1] = _nodes[0]; target[2] = _nodes[4]; break;  //304
        
        case 4 : target[0] = _nodes[2]; target[1] = _nodes[1]; target[2] = _nodes[0]; break;  //210
        case 5 : target[0] = _nodes[2]; target[1] = _nodes[0]; target[2] = _nodes[3]; break;  //203
        
        default:break;
      }
    }

    ///
    template<typename box_t, typename CoordAcc>
    void bbox(const CoordAcc& acrd, box_t&bb) const
    {
      for (E_Int i = 0; i < 3; ++i)
        {bb.minB[i] = K_CONST::E_MAX_FLOAT; bb.maxB[i] = -K_CONST::E_MAX_FLOAT;}

      bb.compute(acrd, _nodes, NB_NODES, _shift/*idx start*/);
    }

private:
    E_Int _shift;
    E_Int _nodes[5];
};

///  
template <typename CoordAcc> inline
void Pyramid::iso_barycenter(const CoordAcc& coord, E_Float* G)
{ 
  //
  for (size_t d=0; d < 3; ++d) G[d]=0.;

  for (E_Int i=0; i < NB_NODES; ++i)
  {
    for (size_t d=0; d < 3; ++d)
    {
      //std::cout << "v : " << coord.getVal(node(i), d) << std::endl;
      G[d] += coord.getVal(node(i), d);
    }
  }

  E_Float k = 1./(E_Float)NB_NODES;

  for (size_t i = 0; i < 3; ++i) G[i] *= k;
  //std::cout << "G : " << G[0] << "/" << G[1] << "/" << G[2] << std::endl;

}

///
template <typename ngunit_t>
inline void Pyramid::iso_barycenter(const K_FLD::FloatArray& crd, const ngunit_t & PGs, const E_Int* first_pg, E_Int nb_pgs, E_Int index_start, E_Float* G)
{
  //WARNING : assuming reodrederd pgs : first is bottom, second is top

  E_Int new_bary[5];

  // face BOT
  const E_Int* nodes = PGs.get_facets_ptr(first_pg[0]-index_start);
  E_Int nb_nodes = PGs.stride(first_pg[0]-index_start);

  for (int k = 0; k  < nb_nodes; ++k)
    new_bary[k] = nodes[k];   

  // Summit
  const E_Int* F1 = PGs.get_facets_ptr(first_pg[1]-index_start); // face F1
  E_Int nb_F1 = PGs.stride(first_pg[1]-index_start);
  
  for (int i=0; i< nb_F1; i++){
      E_Int n(0);
      for (int j=0; j < nb_nodes; j++){
          if ( *(F1+i) != *(nodes+j) ){
              n++;
          }
      }
      if ( n == (nb_nodes-1))
        new_bary[4] = *(F1+i);
  }
  
  K_MESH::Polyhedron<STAR_SHAPED>::iso_barycenter(crd, new_bary, 5, 1, G);
}

///
template< typename ngo_t>
void Pyramid::reorder_pgs(ngo_t& ng, const K_FLD::IntArray& F2E, E_Int i)
{
  std::map<E_Int,E_Int> glmap; 
  E_Int nb_faces = ng.PHs.stride(i); 
  E_Int* faces = ng.PHs.get_facets_ptr(i);
  E_Int PGi = faces[0] - 1;

  E_Int l(0);
  for (int i=0; i< nb_faces; i++){
    if (ng.PGs.stride(PGi)!=4){
        PGi= faces[i] - 1;
        l=i;
    }
  }
  
  E_Int* pN = ng.PGs.get_facets_ptr(PGi);

  
  glmap[*pN] = 0; // PHi(0,0) -> 0  
  glmap[*(pN+1)] = 1;
  glmap[*(pN+2)] = 2;
  glmap[*(pN+3)] = 3;


  if (F2E(1,PGi) != i) // for BOT, PH is the right element. if not, wrong orientation => swap of 1 and 3
  { 
    glmap[*(pN+1)] = 3;
    glmap[*(pN+3)] = 1;
  }

  E_Int F1Id(E_IDX_NONE), F2Id(E_IDX_NONE), F3Id(E_IDX_NONE), F4Id(E_IDX_NONE);

  for (int k = 1; k < nb_faces; ++k)
  {
    int count = 0;
    std::vector<bool> commonNodes(4,false);
    E_Int testedPG = faces[(k+l) % nb_faces]-1;
    E_Int* pNode = ng.PGs.get_facets_ptr(testedPG);
    for (int j = 0; j < 3; ++j)
    {
      auto it = glmap.find(pNode[j]);
      if (it != glmap.end()){
        // found
        count++;
        commonNodes[it->second] = true;
      }
    }
    if (commonNodes[0] && commonNodes[1])
      F1Id = k;
    else if (commonNodes[1] && commonNodes[2])
      F2Id = k;
    else if (commonNodes[2] && commonNodes[3])
      F3Id = k;
    else if (commonNodes[3] && commonNodes[0])
      F4Id = k;    
    }

  E_Int mol[5];

  mol[0] = faces[0];
  mol[1] = faces[F1Id];
  mol[2] = faces[F2Id];
  mol[3] = faces[F3Id];
  mol[4] = faces[F4Id];

  for (int i = 0; i < nb_faces; ++i)
    faces[i] = mol[i];
}

}
#endif	/* __K_MESH_PYRAMID_H__ */


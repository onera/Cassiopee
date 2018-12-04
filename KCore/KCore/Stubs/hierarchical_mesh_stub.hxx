/*
 
 
 
              NUGA 
 
 
 
 */
//Authors : Sâm Landier (sam.landier@onera.fr), Alexis Gay (alexis.gay@onera.fr)

#ifndef NUGA_HIERACHICAL_MESH_HXX
#define NUGA_HIERACHICAL_MESH_HXX

#if defined (DEBUG_HIERARCHICAL_MESH) || defined (OUTPUT_ITER_MESH)
#include "IO/io.h"
#include "Nuga/Boolean/NGON_debug.h"
using NGDBG = NGON_debug<K_FLD::FloatArray,K_FLD::IntArray>;
#endif

#include "Nuga/include/tree.hxx"
#include "Nuga/include/q9.hxx"
#include "Nuga/include/h27.hxx"
#include "Connect/IdTool.h"
#include "Nuga/Delaunay/Triangulator.h"

#define NEIGHBOR(PHi, _F2E, PGi) ( (_F2E(0,PGi) == PHi) ? _F2E(1,PGi) : _F2E(0,PGi) )

namespace NUGA
{

template <typename ELT_t>
class tree_trait;

template<>
class tree_trait<K_MESH::Hexahedron>
{
  public:
  using arr_type = K_FLD::IntArray;
  using face_type = K_MESH::Quadrangle;
};

template<>
class tree_trait<K_MESH::Polyhedron<UNKNOWN> >
{
  public:
  using arr_type = ngon_unit;
  using face_type = K_MESH::Polygon;
};


template <typename ELT_t, typename ngo_t = ngon_type, typename crd_t = K_FLD::FloatArray>
class hierarchical_mesh
{
  
  public:
  using elt_type = ELT_t;
  using arr_t = typename tree_trait<elt_type>::arr_type;
  using face_t = typename tree_trait<elt_type>::face_type;
  
  crd_t&                    _crd;
  ngo_t&                    _ng;
  tree<arr_t>               _PGtree, _PHtree;
  K_FLD::IntArray           _F2E;
  bool                      _initialized;

  hierarchical_mesh(crd_t& crd, ngo_t & ng):_crd(crd), _ng(ng), _PGtree(ng.PGs, face_t::NB_NODES), _PHtree(ng.PHs, elt_type::NB_NODES), _initialized(false){}

  E_Int init(){ return 0;}
  
  E_Int adapt(Vector_t<E_Int>& adap_incr, bool do_agglo){ return 0;}
  
  /// face-conformity
  void conformize(){}
  /// Keep only enabled PHs
  void filter_ngon(ngon_type& filtered_ng){}
  
  void refine_PGs(const Vector_t<E_Int> &PHadap){}
  
  void refine_PHs(const Vector_t<E_Int> &PHadap){}
  
  void get_nodes_PHi(E_Int* nodes, E_Int PHi, E_Int centroid_id, E_Int* BOT, E_Int* TOP, E_Int* LEFT, E_Int* RIGHT, E_Int* FRONT, E_Int* BACK){}
  
  void retrieve_ordered_data(E_Int PGi, E_Int i0, bool reorient, E_Int* four_childrenPG, E_Int* LNODES){}
  
  bool need_a_reorient(E_Int PGi, E_Int PHi, bool oriented_if_R){}
  
  E_Int get_i0(E_Int* pFace, E_Int common_node, E_Int* nodes, E_Int nb_edges_face) {return 0;}
  
  void update_F2E(E_Int PHi, E_Int PHchildr0, E_Int* INT, E_Int* BOT, E_Int* TOP, E_Int* LEFT, E_Int* RIGHT, E_Int* FRONT, E_Int* BACK){}
  
  void update_children_F2E(E_Int PGi, E_Int side){}
  
  void get_cell_center(E_Int PHi, E_Float* center){}
  
  void get_enabled_neighbours(E_Int PHi, E_Int* neighbours, E_Int& nb_neighbours){}
  
  void get_higher_level_neighbours(E_Int PHi, E_Int PGi, E_Int* neighbours, E_Int& nb_neighbours){}
  
  bool enabled_neighbours(E_Int PHi){return false;}
  
  void smooth(Vector_t<E_Int>& adap_incr){}
  
  private:
    std::map<K_MESH::NO_Edge,E_Int> _ecenter;
    Vector_t<E_Int>                 _fcenter;
    
    void __compute_edge_center(const Vector_t<E_Int> &PGlist, K_FLD::FloatArray& crd, std::map<K_MESH::NO_Edge,E_Int> & ecenter){}
    void __compute_face_centers(K_FLD::FloatArray& crd, const typename ngo_t::unit_type& pgs, const Vector_t<E_Int> &PGlist, Vector_t<E_Int>& fcenter){}
    inline void __compute_face_center(const crd_t& crd, const E_Int* nodes, E_Int nb_nodes, E_Float* C){}
    void __compute_cell_center(const crd_t& crd, const E_Int* nodes27, E_Float* C){}
  
};

}

#endif

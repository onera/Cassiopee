/*
 
 
 
              NUGA 
 
 
 
 */
//Authors : Sâm Landier (sam.landier@onera.fr), Alexis Gay (alexis.gay@onera.fr)

#ifndef NUGA_GEOM_SENSOR_HXX
#define NUGA_GEOM_SENSOR_HXX

#include "Connect/EltAlgo.h"
#include "MeshElement/Hexahedron.h"
#include "Fld/ngon_t.hxx"
#include <limits.h>
using ngon_type = ngon_t<K_FLD::IntArray>; 

//#define NB_PTS_TO_INCR(nb_pts) (nb_pts/4)
//#define NB_PTS_TO_INCR(nb_pts) ((nb_pts > 0 ) ? nb_pts -1 : 0) // static_sensor : n-1 subdivision where n is the number of source pts in a cell

namespace NUGA
{

/// Octal geometric sensor
template <typename mesh_t, typename crd_t = K_FLD::FloatArray> //ngu for surfacic (PGs) or ngon_t for volumic
class geom_sensor
{
  public:
    
    using data_type = crd_t; //point cloud
    
    geom_sensor(mesh_t& mesh, E_Int max_pts_per_cell = 1, E_Int itermax = -1)
            :_hmesh(mesh), _is_init(false), _refine(true), _has_agglo(false), _agglo(true), _max_pts_per_cell(max_pts_per_cell), _iter_max((itermax <= 0) ? INT_MAX : itermax), _iter(0){}
    
    virtual E_Int init(data_type& data);
    
    bool compute(data_type& data, Vector_t<E_Int>& adap_incr, bool do_agglo);

    void fill_adap_incr(Vector_t<E_Int>& adap_incr, bool do_agglo);
    
    void locate_points(K_SEARCH::BbTree3D& tree, data_type& data);
    
    E_Int get_highest_lvl_cell(E_Float* p, E_Int& PHi);
    
    E_Int get_higher_lvl_cell(E_Float* p, E_Int PHi);
    
    E_Int detect_child(E_Float* p, E_Int PHi, E_Int* children);
    
  protected:
    mesh_t& _hmesh;
    bool _is_init;
    Vector_t<E_Int> _points_to_cell;
    bool _refine; // continue the adaptation
    bool _has_agglo; // one cell might be able to be agglomerated
    bool _agglo; // we have at least 1 agglomeration
    E_Int _max_pts_per_cell;
    E_Int _iter_max, _iter;
};

/// 
template <typename mesh_t, typename crd_t>
E_Int geom_sensor<mesh_t, crd_t>::init(data_type& data)
{ 
  return 0;
}

/// 
template <typename mesh_t, typename crd_t>
bool geom_sensor<mesh_t, crd_t>::compute(data_type& data, Vector_t<E_Int>& adap_incr, bool do_agglo)
{
  
  return true;
}

///
template <typename mesh_t, typename crd_t>
void geom_sensor<mesh_t, crd_t>::fill_adap_incr(std::vector<E_Int>& adap_incr, bool do_agglo)
{
}

///
template <typename mesh_t, typename crd_t>
void geom_sensor<mesh_t, crd_t>::locate_points(K_SEARCH::BbTree3D& tree, data_type& data)
{
}

///
template <typename mesh_t, typename crd_t>
E_Int geom_sensor<mesh_t, crd_t>::get_higher_lvl_cell(E_Float* p, E_Int PHi)
{

  return E_IDX_NONE; // never reached : to silent compilers
}

///
template <typename mesh_t, typename crd_t>
E_Int geom_sensor<mesh_t, crd_t>::get_highest_lvl_cell(E_Float* p, E_Int& PHi)
{ 
  return E_IDX_NONE; // never reached : to silent compilers   
}

///
template <typename mesh_t, typename crd_t>
E_Int geom_sensor<mesh_t, crd_t>::detect_child(E_Float* p, E_Int PHi, E_Int* children)
{
  return 0;
}

}

#endif

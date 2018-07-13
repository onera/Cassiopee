/*
 
 
 
              NUGA 
 
 
 
 */
//Authors : Sâm Landier (sam.landier@onera.fr), Alexis Gay (alexis.gay@onera.fr)

#ifndef NUGA_GEOM_SENSOR_HXX
#define NUGA_GEOM_SENSOR_HXX

#include "Connect/EltAlgo.h"
#include "MeshElement/Hexahedron.h"
#include "Fld/ngon_t.hxx"
using ngon_type = ngon_t<K_FLD::IntArray>; 

//#define NB_PTS_TO_INCR(nb_pts) (nb_pts/4)
#define NB_PTS_TO_INCR(nb_pts) ((nb_pts > 0 ) ? nb_pts -1 : 0)

namespace NUGA
{

/// Octal geometric sensor
template <typename mesh_t, typename crd_t = K_FLD::FloatArray> //ngu for surfacic (PGs) or ngon_t for volumic
class geom_sensor
{      
  public:
    
    using data_type = crd_t; //point cloud
    
    geom_sensor(mesh_t& mesh):_hmesh(mesh), _is_init(false), _continue(true){}
    
    E_Int init(data_type& data);
    
    bool compute(data_type& data, Vector_t<E_Int>& adap_incr);
    
    void locate_points(K_SEARCH::BbTree3D& tree, data_type& data);
    
    E_Int get_highest_lvl_cell(E_Float* p, E_Int& PHi, E_Bool exit);
    
    E_Int get_higher_lvl_cell(E_Float* p, E_Int PHi);
    
    E_Int detect_child(E_Float* p, E_Int PHi, E_Int* children);
    
  private:
    mesh_t& _hmesh;
    bool _is_init;
    Vector_t<E_Int> _points_to_cell;
    bool _continue;
};

/// 
template <typename mesh_t, typename crd_t>
E_Int geom_sensor<mesh_t, crd_t>::init(data_type& data)
{ 
  
  {// brackets to kill the tree asap
    K_SEARCH::BbTree3D tree(_hmesh._crd, _hmesh._ng);
  
    // Count the number of points per cell
    locate_points(tree, data);
  }
  
  return 0;
}

///
template <typename mesh_t, typename crd_t>
bool geom_sensor<mesh_t, crd_t>::compute(data_type& data, Vector_t<E_Int>& adap_incr)
{
  using ELT_t = typename mesh_t::elt_type;

#ifdef FLAG_STEP
  std::cout << _hmesh._ng.PHs.size() << std::endl;
#endif
  
  E_Int nb_elts = _hmesh._ng.PHs.size();
  E_Int nb_pts = _points_to_cell.size();  
  if (_is_init) // first iteration : init is done instead
  {
    _continue=false;

    for (int i = 0; i < nb_pts; ++i)
    {
      E_Int PHi = _points_to_cell[i];

      if (PHi == E_IDX_NONE) continue; //external point
      if (_hmesh._tree._PHtree.children(PHi) != nullptr) // we detect the subdivised PH
      {
        E_Float* p = data.col(i); // in which children of PHi is the point
        E_Int cell = get_higher_lvl_cell(p,PHi);
        _points_to_cell[i] = cell;
      }
    }
  }
      
  _is_init = true; 
  
  //points_to_cell gives the number of points per cell
  {
    Vector_t<E_Int> nb_pts_per_cell(nb_elts,0);
    for (int i = 0; i < nb_pts; ++i)
    {
        if (_points_to_cell[i] != E_IDX_NONE) 
        {
            if (nb_pts_per_cell[_points_to_cell[i]] == 1) _continue = true;
            nb_pts_per_cell[_points_to_cell[i]] += 1;
        }
    }
    
    //adap_incr
    adap_incr.clear();
    adap_incr.resize(nb_elts, 0);
    
    for (int i = 0; i < nb_elts; ++i)
    {
        if (nb_pts_per_cell[i] >= 2)
            adap_incr[i] = 1;
    }
  }
  if (_continue == false ) return false;

  K_CONNECT::EltAlgo<typename mesh_t::elt_type>::smoothd1(_hmesh._ng.PHs, _hmesh._F2E, _hmesh._tree._PHtree.level(), adap_incr); // 2-1 smoothing

  return true;
}

///
template <typename mesh_t, typename crd_t>
void geom_sensor<mesh_t, crd_t>::locate_points(K_SEARCH::BbTree3D& tree, data_type& data)
{
  using ELT_t = typename mesh_t::elt_type;
  
  E_Int nb_src_pts = data.cols();
  
  _points_to_cell.resize(nb_src_pts, E_IDX_NONE);
  
  Vector_t<E_Int> ids;
  
  for (int i = 0; i < nb_src_pts; i++)
  {
    E_Float minB[3];
    E_Float maxB[3];
    E_Float* p = data.col(i);
    
    for (int j = 0; j < 3;j++)
    {
      minB[j] = p[j]-E_EPSILON;
      maxB[j] = p[j]+E_EPSILON;
    }  
    ids.clear();
    tree.getOverlappingBoxes(minB,maxB,ids);
    
    for (size_t j = 0; j < ids.size(); j++)
    {
      E_Int PHi_orient[6];
      ELT_t::get_orient(_hmesh._ng, _hmesh._F2E, ids[j], PHi_orient);
      if (ELT_t::pt_is_inside(_hmesh._ng, _hmesh._crd, ids[j], PHi_orient, p))
      {
        if (_hmesh._tree._PHtree.children(ids[j]) == nullptr)
          _points_to_cell[i] = ids[j];
        else
        {   
          E_Int cell = get_highest_lvl_cell(p, ids[j], false);
          _points_to_cell[i] = cell;
        }

        break;
      }
    }
  }
}

///
template <typename mesh_t, typename crd_t>
E_Int geom_sensor<mesh_t, crd_t>::get_highest_lvl_cell(E_Float* p, E_Int& PHi, E_Bool exit)
{
    if (exit) return PHi;
    
    using ELT_t = typename mesh_t::elt_type;
    
    E_Int* q = _hmesh._tree._PHtree.children(PHi); // the children of PHi
    E_Int nb_children = _hmesh._tree._PHtree.nb_children(PHi);
    E_Bool found = false;
    for (int j = 0; j < nb_children; ++j)
    {
        E_Int PHi_orient[6];
        ELT_t::get_orient(_hmesh._ng, _hmesh._F2E, q[j], PHi_orient);
        if (ELT_t::pt_is_inside(_hmesh._ng, _hmesh._crd, q[j], PHi_orient, p)) // true when the point is located in a child
        {
            found = true;
            if (_hmesh._tree._PHtree.children(q[j]) != nullptr)
                get_highest_lvl_cell(p, q[j], false);
            else
                get_highest_lvl_cell(p, q[j], true);
        }
    }
    if (!found)
    {
        E_Int id = detect_child(p, PHi, q); // should be the most subdivised one, because we move the source point when it happens in order to iterate more
        get_highest_lvl_cell(p, id, true);
    }
}

///
template <typename mesh_t, typename crd_t>
E_Int geom_sensor<mesh_t, crd_t>::get_higher_lvl_cell(E_Float* p, E_Int PHi)
{
    using ELT_t = typename mesh_t::elt_type;
    
    E_Int* q = _hmesh._tree._PHtree.children(PHi); // the children of PHi
    E_Int nb_children = _hmesh._tree._PHtree.nb_children(PHi);
    E_Bool found = false;
    for (int j = 0; j < nb_children; ++j)
    {
        E_Int PHi_orient[6];
        ELT_t::get_orient(_hmesh._ng, _hmesh._F2E, q[j], PHi_orient);
        if (ELT_t::pt_is_inside(_hmesh._ng, _hmesh._crd, q[j], PHi_orient, p)) // true when the point is located in a child
        {
            found = true;
            return q[j];
        }
    }
    if (!found) return detect_child(p,PHi,q); // we have to allocate the closest child, and move the source point (be careful it's still const atm)
    
}

///
template <typename mesh_t, typename crd_t>
E_Int geom_sensor<mesh_t, crd_t>::detect_child(E_Float* p, E_Int PHi, E_Int* children)
{
    using ELT_t = typename mesh_t::elt_type;
    E_Int nb_children = _hmesh._tree._PHtree.nb_children(PHi);
    
    // find the closest child
    E_Int closest_child = 0;
    E_Float center[3];
    _hmesh.get_cell_center(children[0], center);
    E_Float min = K_FUNC::sqrDistance(p, center,3);

    for (int i = 1; i < nb_children; ++i)
    {
        E_Float tmp[3];
        _hmesh.get_cell_center(children[i], tmp);
        E_Float d = K_FUNC::sqrDistance(p,tmp,3);
        if (d < min) 
        {
            min = d;
            closest_child = i;
            for (int i = 0; i < 3; ++i)
                center[i] = tmp[i];
        }
    }
    // move the source point into the closest child
    
    //vector
    E_Float u0 = p[0] - center[0];
    E_Float u1 = p[1] - center[1];
    E_Float u2 = p[2] - center[2];
    E_Float coef = 1;
    
    p[0] = center[0] + coef*u0;
    p[1] = center[1] + coef*u1;
    p[2] = center[2] + coef*u2;
    
    return children[closest_child];
    
}




///////////////////// COARSE SENSOR //////////////////////////////////////////////
//
/////
//template <typename mesh_t, typename crd_t = K_FLD::FloatArray> //ngu for surfacic (PGs) or ngon_t for volumic
//class geom_static_sensor
//{
//  public:
//    
//    using data_type = crd_t; //point cloud
//    
//    geom_static_sensor(const mesh_t& mesh):_hmesh(mesh), _is_done(false){}
//    
//    E_Int init(const data_type& data, Vector_t<E_Int>& adap_incr) { return 0; } //nothing to do
//
//    E_Bool compute(const data_type& data, Vector_t<E_Int>& adap_incr);
//    
//    void locate_points(K_SEARCH::BbTree3D& tree, const data_type& data, Vector_t<E_Int>& pts_per_cell);
//    
//  private :
//      const mesh_t& _hmesh;
//      E_Bool _is_done;
//};
//
/////
//template <typename mesh_t, typename crd_t>
//E_Bool geom_static_sensor<mesh_t, crd_t>::compute(const data_type& data, Vector_t<E_Int>& adap_incr)
//{  
//  
//  if (_is_done) return false; // to compute once only
//  _is_done = true;
//  
//  
//  E_Int nb_elts = _hmesh._ng.PHs.size();
//  
//  adap_incr.clear();
//  adap_incr.resize(nb_elts, 0);// now 0 should be the "no adaptation" value (instead of -1)
//  
//  {// brackets to kill the tree asap
//    K_SEARCH::BbTree3D tree(_hmesh._crd, _hmesh._ng);
//  
//    // Count the number of data per cell
//    locate_points(tree, data, adap_incr);
//  } 
//  
//  // Translate it as an adaptation quantity
//  for (E_Int i=0; i < nb_elts; ++i)
//    adap_incr[i] = NB_PTS_TO_INCR(adap_incr[i]);
//
//  K_CONNECT::EltAlgo<typename mesh_t::elt_type>::smoothd1(_hmesh._ng.PHs, _hmesh._F2E, _hmesh._tree._PHtree.level(), adap_incr); // 2-1 smoothing
//
//  return true;
//  
//}
//
//template <typename mesh_t, typename crd_t>
//void geom_static_sensor<mesh_t, crd_t>::locate_points(K_SEARCH::BbTree3D& tree, const data_type& data, Vector_t<E_Int>& pts_per_cell)
//{
//  E_Int nb_src_pts = data.cols();
//  
//  Vector_t<E_Int> ids;
//  
//  for (int i = 0; i < nb_src_pts; i++)
//  {
//    
//    
//    E_Float minB[3];
//    E_Float maxB[3];
//    const E_Float* p = data.col(i);
//    
//    for (int j = 0; j < 3;j++)
//    {
//        minB[j] = p[j]-E_EPSILON;
//        maxB[j] = p[j]+E_EPSILON;
//    }  
//    ids.clear();
//    tree.getOverlappingBoxes(minB,maxB,ids);
//    
//    for (size_t j = 0; j < ids.size(); j++)
//        pts_per_cell[ids[j]] += 1;
//  }
//
//}
//
//
}

#endif

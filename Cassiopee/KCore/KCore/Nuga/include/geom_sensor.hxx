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

#ifndef NUGA_GEOM_SENSOR_HXX
#define NUGA_GEOM_SENSOR_HXX

#include "Nuga/include/sensor.hxx"
#include "Nuga/include/EltAlgo.h"
#include "Nuga/include/Hexahedron.h"
#include "Nuga/include/ngon_t.hxx"
#include "V1_smoother.hxx"
#include "shell_smoother.hxx"
#include "Nuga/include/BbTree.h"

#include <limits.h>
#ifdef DEBUG_2019
#include <chrono>
#include <ctime>
#endif
using ngon_type = ngon_t<K_FLD::IntArray>; 


namespace NUGA
{
  
  using crd_t = K_FLD::FloatArray;

/// Geometric sensor
template <typename mesh_t, typename sensor_input_t = crd_t>
class geom_sensor : public sensor<mesh_t, sensor_input_t>
{
  public:

    using parent_t = sensor<mesh_t, sensor_input_t>;
    using output_t   = typename mesh_t::output_t;
    
    geom_sensor(mesh_t& mesh, eSmoother smoo_type, E_Int max_pts_per_cell, E_Int itermax)
      :parent_t(mesh,nullptr),
      _max_pts_per_cell(max_pts_per_cell), _iter_max((itermax <= 0) ? INT_MAX : itermax), _iter(0)
    {
      _bbtree = new K_SEARCH::BbTree3D(parent_t::_hmesh._crd, parent_t::_hmesh._ng, EPSILON);

      if (smoo_type == eSmoother::V1_NEIGH)
        parent_t::_smoother = new V1_smoother<mesh_t>();
      else if (smoo_type == eSmoother::SHELL)
        parent_t::_smoother = new shell_smoother<mesh_t>();

    }
    
    virtual E_Int assign_data(const sensor_input_t& data) override;

    bool fill_adap_incr(output_t& adap_incr, bool do_agglo) override;

    bool update() override;

    bool stop() override;

    virtual ~geom_sensor() {delete _bbtree;}

  protected:
    
    void locate_points();

  private:
    
    E_Int get_highest_lvl_cell(const E_Float* p, E_Int PHi) const ;
    
    E_Int get_higher_lvl_cell(const E_Float* p, E_Int PHi) const ;
    
    E_Int detect_child(const E_Float* p, E_Int PHi, const E_Int* children) const ;

  protected:
    
    K_SEARCH::BbTree3D* _bbtree;
    Vector_t<E_Int> _points_to_cell;
    E_Int _max_pts_per_cell;
    E_Int _iter_max, _iter;
    E_Int _cur_nphs;
};

/// 
template <typename mesh_t, typename sensor_input_t>
E_Int geom_sensor<mesh_t, sensor_input_t>::assign_data(const sensor_input_t& data) 
{

  parent_t::assign_data(data);
   
  _iter = 0; //reset
  // fill in _points_to_cell with the initial mesh (giving for each source point the highest level cell containing it)
  _points_to_cell.clear();

  // Count the number of data per cell

  // Check if hmesh has been initialized
  if (not parent_t::_hmesh.is_initialised())
  {
    std::cout << "Error in geom_sensor.hxx / function assign_data: hmesh has not been initialised." << std::endl; 
    return 1;
  }
  
  locate_points();
  
  return 0;
}

///
template <typename mesh_t, typename sensor_input_t>
bool geom_sensor<mesh_t, sensor_input_t>::fill_adap_incr(output_t& adap_incr, bool do_agglo)
{
  bool filled{ false };

  E_Int nb_faces = parent_t::_hmesh._ng.PGs.size();
  E_Int nb_pts = _points_to_cell.size();

  adap_incr.face_adap_incr.clear();
  adap_incr.cell_adap_incr.clear();

  using cell_incr_t = typename output_t::cell_incr_t;
  using face_incr_t = typename output_t::face_incr_t;

  adap_incr.cell_adap_incr.resize(_cur_nphs, cell_incr_t(0));
  adap_incr.face_adap_incr.resize(nb_faces, face_incr_t(0));

  //points_to_cell gives the number of points per cell
  Vector_t<E_Int> nb_pts_per_cell(_cur_nphs,0);
  
  for (E_Int i = 0; i < nb_pts; ++i)
  {
    E_Int PHi = _points_to_cell[i];
    if (PHi != IDX_NONE)
      nb_pts_per_cell[PHi] += 1;
  }

  E_Int maxnbpercell = *std::max_element(ALL(nb_pts_per_cell));

  if (maxnbpercell > _max_pts_per_cell)
  {
    for (E_Int i = 0; i < _cur_nphs; ++i)
    {
      const E_Int* faces = parent_t::_hmesh._ng.PHs.get_facets_ptr(i);
      E_Int nb_faces = parent_t::_hmesh._ng.PHs.stride(i);
      bool admissible_elt = K_MESH::Polyhedron<0>::is_basic(parent_t::_hmesh._ng.PGs, faces, nb_faces);
      if (admissible_elt && (nb_pts_per_cell[i] > _max_pts_per_cell) && (parent_t::_hmesh._PHtree.is_enabled(i))) {// can be and has to be subdivided
        adap_incr.cell_adap_incr[i] = 1;
        filled = true;
      }
    }
  }

  if (do_agglo)
  {
    for (E_Int i = 0; i < _cur_nphs; ++i)
    {
      if (adap_incr.cell_adap_incr[i] == -1) continue; //already processed by brothorhood

      if ((parent_t::_hmesh._PHtree.get_level(i) > 0) && (parent_t::_hmesh._PHtree.is_enabled(i))) // this cell may be agglomerated : check its brothers
      {        
        E_Int father = parent_t::_hmesh._PHtree.parent(i);
        const E_Int* p = parent_t::_hmesh._PHtree.children(father);
        E_Int nb_children = parent_t::_hmesh._PHtree.nb_children(father);
        E_Int sum = 0;
                
        for (E_Int k = 0; k < nb_children; ++k)
        {
          sum += nb_pts_per_cell[p[k]];
          if (!parent_t::_hmesh._PHtree.is_enabled(p[k])) sum = _max_pts_per_cell + 1;
          if (sum > _max_pts_per_cell) break; // number of source points in the father > criteria chosen : cannot agglomerate
        }
        if (sum <= _max_pts_per_cell) // number of source points in the father <= criteria chosen : might be agglomerated : -1 for every child
        {
          for (E_Int k = 0; k < nb_children; ++k) {
            adap_incr.cell_adap_incr[p[k]] = -1;
            filled = true;
          }
        }
      }
    }
  }
  
  return filled;
}

///
template <typename mesh_t, typename sensor_input_t>
void geom_sensor<mesh_t, sensor_input_t>::locate_points() 
{
  // locate every source points in a given mesh
 
  E_Float tol = 10. * ZERO_M;
  
  E_Int nb_src_pts =parent_t::_data.getSize();
  // for each source points, give the highest cell containing it if there is one, otherwise IDX_NONE
  _points_to_cell.resize(nb_src_pts, IDX_NONE);
  
  Vector_t<E_Int> ids;

  bool found=false;
  E_Float minB[3];
  E_Float maxB[3];
  //
  for (E_Int i = 0; i < nb_src_pts; i++)
  {
    const E_Float* p =parent_t::_data.get(i);
    
    for (E_Int j = 0; j < 3;j++)
    {
      minB[j] = p[j]-EPSILON; // small box over the source point
      maxB[j] = p[j]+EPSILON;
    }  
    ids.clear();
    _bbtree->getOverlappingBoxes(minB,maxB,ids); // ids contains the list of PH index partially containing the small box around the source point
    
    if (ids.empty()) continue; // exterior point : IDX_NONE
    
    for (size_t j = 0; j < ids.size(); j++) // which of these boxes has the source point ?
    {
      
      if (K_MESH::Polyhedron<UNKNOWN>::pt_is_inside(ids[j], parent_t::_hmesh._ng.PGs, parent_t::_hmesh._ng.PHs.get_facets_ptr(ids[j]), parent_t::_hmesh._ng.PHs.stride(ids[j]), parent_t::_hmesh._crd, parent_t::_hmesh._F2E, p, tol))
      {
        _points_to_cell[i] = get_highest_lvl_cell(p, ids[j]);
        found = true;
        break;
      }
    }
  }

  if (!found) std::cout << "WARNING : All the source points are outside the mesh. Nothing to do." << std::endl;

}

///
template <typename mesh_t, typename sensor_input_t>
E_Int geom_sensor<mesh_t, sensor_input_t>::get_higher_lvl_cell(const E_Float* p, E_Int PHi) const 
{
  // returns in the PHi's sons where a given source point is
  
  const E_Int* q = parent_t::_hmesh._PHtree.children(PHi); // the children of PHi
  E_Int nb_children = parent_t::_hmesh._PHtree.nb_children(PHi);
  bool found = false;
    
  for (E_Int j = 0; j < nb_children; ++j)
  {
    if (K_MESH::Polyhedron<UNKNOWN>::pt_is_inside(q[j], parent_t::_hmesh._ng.PGs, parent_t::_hmesh._ng.PHs.get_facets_ptr(q[j]), parent_t::_hmesh._ng.PHs.stride(q[j]), parent_t::_hmesh._crd, parent_t::_hmesh._F2E, p, 1.e-14)) // true when the point is located in a child
    {
      found = true;
      return q[j];
    }
  }
  if (!found) return detect_child(p,PHi,q); // if we can't find the child, we have to locate the closest child  (almost never occurs)

  return IDX_NONE; // never reached : to silent compilers

}

///
template <typename mesh_t, typename sensor_input_t>
E_Int geom_sensor<mesh_t, sensor_input_t>::get_highest_lvl_cell(const E_Float* p, E_Int PHi) const 
{
  // will find, in the tree (starting at PHi cell), the highest level cell (ie smallest cell) cointaining a given source point

  // Do not test children if PHi is enabled
  if (parent_t::_hmesh._PHtree.is_enabled(PHi)) return PHi;

  E_Int nb_children = parent_t::_hmesh._PHtree.nb_children(PHi);
  assert(nb_children != 0);
  
  const E_Int* q = parent_t::_hmesh._PHtree.children(PHi); // the children of PHi

  bool found = false;
  for (E_Int j = 0; j < nb_children; ++j)
  {
    if (K_MESH::Polyhedron<UNKNOWN>::pt_is_inside(q[j], parent_t::_hmesh._ng.PGs, parent_t::_hmesh._ng.PHs.get_facets_ptr(q[j]), parent_t::_hmesh._ng.PHs.stride(q[j]), parent_t::_hmesh._crd, parent_t::_hmesh._F2E, p, 1.e-14)) // true when the point is located in a child
    {
      found = true;
      if (parent_t::_hmesh._PHtree.children(q[j]) != nullptr)
        return get_highest_lvl_cell(p, q[j]);
      else
        return q[j];
    }
  }
  if (!found) return detect_child(p, PHi, q);
  
  return IDX_NONE; // never reached : to silent compilers
    
}

///
template <typename mesh_t, typename sensor_input_t>
E_Int geom_sensor<mesh_t, sensor_input_t>::detect_child(const E_Float* p, E_Int PHi, const E_Int* children) const
{
  E_Int nb_children = parent_t::_hmesh._PHtree.nb_children(PHi);
    
  // find the closest child
  E_Int closest_child = 0;
  E_Float center[3];
  parent_t::_hmesh.get_cell_center(children[0], center);
  E_Float min = NUGA::sqrDistance(p, center,3);

  E_Float tmp[3];
  tmp[0] = K_CONST::E_MAX_FLOAT; 
  tmp[1] = K_CONST::E_MAX_FLOAT; 
  tmp[2] = K_CONST::E_MAX_FLOAT; 
  for (E_Int i = 1; i < nb_children; ++i)
  {
    parent_t::_hmesh.get_cell_center(children[i], tmp);
    E_Float d = NUGA::sqrDistance(p, tmp, 3);

    if (d < min)
    {
      min = d;
      closest_child = i;
      for (E_Int i = 0; i < 3; ++i) center[i] = tmp[i];
    }
  }
    
  return children[closest_child];
}

///
template <typename mesh_t, typename sensor_input_t>
bool geom_sensor<mesh_t, sensor_input_t>::update()
{
  //std::cout << "redistrib data geom sensor. : nb pts" << _points_to_cell.size() << std::endl;
  
  // if no subdivision has occured => no reason for a new point dispatch
  if (_cur_nphs == parent_t::_hmesh._ng.PHs.size()) return false;
  
  E_Int nb_pts = _points_to_cell.size();

  for (E_Int i = 0; i < nb_pts; ++i)
  {
    E_Int PHi = _points_to_cell[i];

    if (PHi == IDX_NONE) continue; //external point
    if (parent_t::_hmesh._PHtree.children(PHi) != nullptr) // we detect the subdivised PH : one of his children has this source point
    {
      E_Float* p = parent_t::_data.get(i); // in which children of PHi is the point
      E_Int cell = get_higher_lvl_cell(p,PHi);
      _points_to_cell[i] = cell;
    }
  }
  return true;
}

template <typename mesh_t, typename sensor_input_t>
bool geom_sensor<mesh_t, sensor_input_t>::stop()
{
  //the following must be done here because this has to be refreshed at each iter to be sync for sensor.update
  _cur_nphs = parent_t::_hmesh._ng.PHs.size();

#ifdef FLAG_STEP
  std::cout << "iter : " << _iter << ". nb of PHs : " << parent_t::_hmesh._ng.PHs.size() << std::endl;
#endif

  if (++_iter > _iter_max) return true;
  return false;
}
  

/*************************** COARSE SENSOR ****************************************/
//#define NB_PTS_TO_INCR(nb_pts) (nb_pts/4)
//#define NB_PTS_TO_INCR(nb_pts) ((nb_pts > 0 ) ? nb_pts -1 : 0) // static_sensor : n-1 subdivision where n is the number of source pts in a cell

//
/////
//template <typename mesh_t, typename crd_t = K_FLD::FloatArray> //ngu for surfacic (PGs) or ngon_t for volumic
//class geom_static_sensor
//{
//  public:
//    
//    using data_t = crd_t; //point cloud
//    
//    geom_static_sensor(const mesh_t& mesh):_hmesh(mesh), _is_done(false){}
//    
//    E_Int init(const data_t& data, Vector_t<E_Int>& adap_incr) { return 0; } //nothing to do
//
//    bool compute(const data_t& data, Vector_t<E_Int>& adap_incr);
//    
//    void locate_points(K_SEARCH::BbTree3D& tree, const data_t& data, Vector_t<E_Int>& pts_per_cell);
//    
//  private :
//      const mesh_t& parent_t::_hmesh;
//      bool _is_done;
//};
//
/////
//template <typename mesh_t>
//bool geom_static_sensor<mesh_t, crd_t>::compute(const data_t& data, Vector_t<E_Int>& adap_incr)
//{  
//  
//  if (_is_done) return false; // to compute once only
//  _is_done = true;
//  
//  
//  E_Int nb_elts = parent_t::_hmesh._ng.PHs.size();
//  
//  adap_incr.clear();
//  adap_incr.resize(nb_elts, 0);// now 0 should be the "no adaptation" value (instead of -1)
//  
//  {// brackets to kill the tree asap
//    K_SEARCH::BbTree3D tree(parent_t::_hmesh._crd, parent_t::_hmesh._ng);
//  
//    // Count the number of data per cell
//    locate_points(tree, data, adap_incr);
//  } 
//  
//  // Translate it as an adaptation quantity
//  for (E_Int i=0; i < nb_elts; ++i)
//    adap_incr[i] = NB_PTS_TO_INCR(adap_incr[i]);
//
//  NUGA::EltAlgo<typename mesh_t::elt_type>::smoothd1(parent_t::_hmesh._ng.PHs, parent_t::_hmesh._F2E, parent_t::_hmesh._PHtree.level(), adap_incr); // 2-1 smoothing
//
//  return true;
//  
//}
//
//template <typename mesh_t>
//void geom_static_sensor<mesh_t, crd_t>::locate_points(K_SEARCH::BbTree3D& tree, const data_t& data, Vector_t<E_Int>& pts_per_cell)
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
//        minB[j] = p[j]-EPSILON;
//        maxB[j] = p[j]+EPSILON;
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

} // NUGA
#endif

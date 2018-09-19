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
//#define NB_PTS_TO_INCR(nb_pts) ((nb_pts > 0 ) ? nb_pts -1 : 0) // static_sensor : n-1 subdivision where n is the number of source pts in a cell

namespace NUGA
{

/// Octal geometric sensor
template <typename mesh_t, typename crd_t = K_FLD::FloatArray> //ngu for surfacic (PGs) or ngon_t for volumic
class geom_sensor
{
  public:
    
    using data_type = crd_t; //point cloud
    
    geom_sensor(mesh_t& mesh, E_Int max_pts_per_cell = 1):_hmesh(mesh), _is_init(false), _refine(true), _has_agglo(false), _agglo(true), _max_pts_per_cell(max_pts_per_cell){}
    
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
};

/// 
template <typename mesh_t, typename crd_t>
E_Int geom_sensor<mesh_t, crd_t>::init(data_type& data)
{ 
  // fill in _points_to_cell with the initial mesh (giving for each source point the highest level cell containing it)
  _points_to_cell.clear();
  
  K_SEARCH::BbTree3D tree(_hmesh._crd, _hmesh._ng);
  
  // Count the number of data per cell
  locate_points(tree, data);
  
  return 0;
}

/// 
template <typename mesh_t, typename crd_t>
bool geom_sensor<mesh_t, crd_t>::compute(data_type& data, Vector_t<E_Int>& adap_incr, bool do_agglo)
{

#ifdef FLAG_STEP
  std::cout << _hmesh._ng.PHs.size() << std::endl;
#endif
  
  E_Int nb_elts = _hmesh._ng.PHs.size();
  E_Int nb_pts = _points_to_cell.size();
  _refine=false;
  _agglo = false;
  _has_agglo = false;
  
  if (_is_init) // first iteration : init is done instead
  {
    //update points_to_cell : replace the PH that have been subdivided in adapt by the child(ren) containing the source point(s)
    for (int i = 0; i < nb_pts; ++i)
    {
      E_Int PHi = _points_to_cell[i];

      if (PHi == E_IDX_NONE) continue; //external point
      if (_hmesh._PHtree.children(PHi) != nullptr) // we detect the subdivised PH : one of his children has this source point
      {
        E_Float* p = data.col(i); // in which children of PHi is the point
        E_Int cell = get_higher_lvl_cell(p,PHi);
        _points_to_cell[i] = cell;
      }
    }
  }

  _is_init = true;

  //fill in adap incr thanks to points to cell
  fill_adap_incr(adap_incr, do_agglo);
  
  //apply the 2:1 rule
  _hmesh.smooth(adap_incr);
  
  //detect if at least one agglomeration can be done
  if (_has_agglo)
  {
    for (int i = 0; i < nb_elts; ++i)
    {
      if (adap_incr[i] == -1)
      {
        _agglo = true;
        break;
      }
    }
  }
  //check if both subdivisions and agglomerations are done
  if (_refine == false && _agglo == false)
    return false;
  
  return true;
}

///
template <typename mesh_t, typename crd_t>
void geom_sensor<mesh_t, crd_t>::fill_adap_incr(std::vector<E_Int>& adap_incr, bool do_agglo)
{
  E_Int nb_elts = _hmesh._ng.PHs.size();
  E_Int nb_pts = _points_to_cell.size();
  //points_to_cell gives the number of points per cell
  Vector_t<E_Int> nb_pts_per_cell(nb_elts,0);
  
  for (int i = 0; i < nb_pts; ++i)
  {
    E_Int PHi = _points_to_cell[i];
    if (PHi != E_IDX_NONE)
      nb_pts_per_cell[PHi] += 1;
  }

  //adap_incr
  adap_incr.clear();
  adap_incr.resize(nb_elts, 0);

  for (int i = 0; i < nb_elts; ++i)
  {
    if ( (nb_pts_per_cell[i] >= _max_pts_per_cell + 1) && (_hmesh._PHtree.is_enabled(i))) // has to be subdivided
    {
      adap_incr[i] = 1;
      _refine = true;
    }
  }

  if (do_agglo)
  {
    for (int i = 0; i < nb_elts; ++i)
    {
      if ( (_hmesh._PHtree.get_level(i) > 0) && (_hmesh._PHtree.is_enabled(i)) ) // this cell may be agglomerated : check its brothers
      {
        if (adap_incr[i] == -1) continue;
                
        E_Int father = _hmesh._PHtree.parent(i);
        E_Int* p = _hmesh._PHtree.children(father);
        E_Int nb_children = _hmesh._PHtree.nb_children(father);
        E_Int sum = 0;
                
        for (int k = 0; k < nb_children; ++k)
        {
          sum += nb_pts_per_cell[p[k]];
          if (!_hmesh._PHtree.is_enabled(p[k])) sum = _max_pts_per_cell + 1;
          if (sum > _max_pts_per_cell) break; // number of source points in the father > criteria chosen : cannot agglomerate
        }
        if (sum <= _max_pts_per_cell) // number of source points in the father <= criteria chosen : might be agglomerated : -1 for every child
        {
          for (int k = 0; k < nb_children; ++k)
            adap_incr[p[k]] = -1;
          _has_agglo = true; // true : means that we might have 1 agglomeration
        }
      }
    }
  }
}

///
template <typename mesh_t, typename crd_t>
void geom_sensor<mesh_t, crd_t>::locate_points(K_SEARCH::BbTree3D& tree, data_type& data)
{
  // locate every source points in a given mesh
  using ELT_t = typename mesh_t::elt_type;

  E_Float tol = 10. * ZERO_M;
  
  E_Int nb_src_pts = data.cols();
  // for each source points, give the highest cell containing it if there is one, otherwise E_IDX_NONE
  _points_to_cell.resize(nb_src_pts, E_IDX_NONE);
  
  Vector_t<E_Int> ids;
  
  for (int i = 0; i < nb_src_pts; i++)
  {
    E_Float minB[3];
    E_Float maxB[3];
    E_Float* p = data.col(i);
    
    for (int j = 0; j < 3;j++)
    {
      minB[j] = p[j]-E_EPSILON; // small box over the source point
      maxB[j] = p[j]+E_EPSILON;
    }  
    ids.clear();
    tree.getOverlappingBoxes(minB,maxB,ids); // ids contains the list of PH index partially containing the small box around the source point

    if (ids.empty()) continue; // exterior point : E_IDX_NONE
    
    for (size_t j = 0; j < ids.size(); j++) // which of these boxes has the source point ?
    {
      E_Int PHi_orient[6];
      ELT_t::get_orient(_hmesh._ng, _hmesh._F2E, ids[j], PHi_orient);
      if (ELT_t::pt_is_inside(_hmesh._ng, _hmesh._crd, ids[j], PHi_orient, p, tol))
      {
        if (_hmesh._PHtree.children(ids[j]) == nullptr) // if the cell has no children : source point i is in this cell
          _points_to_cell[i] = ids[j];
        else // if the cell has children : find the highest level (i.e. smallest) cell containing the source point
        {
          E_Int cell = get_highest_lvl_cell(p, ids[j]);
          _points_to_cell[i] = cell;
        }

        break;
      }
    }
  }
}

///
template <typename mesh_t, typename crd_t>
E_Int geom_sensor<mesh_t, crd_t>::get_higher_lvl_cell(E_Float* p, E_Int PHi)
{
  // returns in the PHi's sons where a given source point is
  using ELT_t = typename mesh_t::elt_type;
    
  E_Int* q = _hmesh._PHtree.children(PHi); // the children of PHi
  E_Int nb_children = _hmesh._PHtree.nb_children(PHi);
  bool found = false;
    
  for (int j = 0; j < nb_children; ++j)
  {
    E_Int PHi_orient[6];
    ELT_t::get_orient(_hmesh._ng, _hmesh._F2E, q[j], PHi_orient);
    if (ELT_t::pt_is_inside(_hmesh._ng, _hmesh._crd, q[j], PHi_orient, p, 1.e-14)) // true when the point is located in a child
    {
      found = true;
      return q[j];
    }
  }
  if (!found) return detect_child(p,PHi,q); // if we can't find the child, we have to locate the closest child  (almost never occurs)

  return E_IDX_NONE; // never reached : to silent compilers

}

///
template <typename mesh_t, typename crd_t>
E_Int geom_sensor<mesh_t, crd_t>::get_highest_lvl_cell(E_Float* p, E_Int& PHi)
{
  // will find, in the tree (starting at PHi cell), the highest level cell (ie smallest cell) cointaining a given source point 
  using ELT_t = typename mesh_t::elt_type;

  E_Int nb_children = _hmesh._PHtree.nb_children(PHi);
  
  if (nb_children == 0) return PHi;
  
  E_Int* q = _hmesh._PHtree.children(PHi); // the children of PHi
  
  E_Int PHi_orient[6];
  bool found = false;
  for (int j = 0; j < nb_children; ++j)
  {
    ELT_t::get_orient(_hmesh._ng, _hmesh._F2E, q[j], PHi_orient);
    if (ELT_t::pt_is_inside(_hmesh._ng, _hmesh._crd, q[j], PHi_orient, p, 1.e-14)) // true when the point is located in a child
    {
      found = true;
      if (_hmesh._PHtree.children(q[j]) != nullptr)
        return get_highest_lvl_cell(p, q[j]);
      else
        return q[j];
    }
  }
  if (!found) return detect_child(p, PHi, q);
  
  return E_IDX_NONE; // never reached : to silent compilers
    
}

///
template <typename mesh_t, typename crd_t>
E_Int geom_sensor<mesh_t, crd_t>::detect_child(E_Float* p, E_Int PHi, E_Int* children)
{
  E_Int nb_children = _hmesh._PHtree.nb_children(PHi);
    
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
//  // move the source point into the closest child ( coef ]0,1] : 1 = doesn't move the source point )
//
//  //vector
//  E_Float u0 = p[0] - center[0];
//  E_Float u1 = p[1] - center[1];
//  E_Float u2 = p[2] - center[2];
//  E_Float coef = 1;
//    
//  p[0] = center[0] + coef*u0;
//  p[1] = center[1] + coef*u1;
//  p[2] = center[2] + coef*u2;
    
  return children[closest_child];
}



/*************************** COARSE SENSOR ****************************************/
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
//    bool compute(const data_type& data, Vector_t<E_Int>& adap_incr);
//    
//    void locate_points(K_SEARCH::BbTree3D& tree, const data_type& data, Vector_t<E_Int>& pts_per_cell);
//    
//  private :
//      const mesh_t& _hmesh;
//      bool _is_done;
//};
//
/////
//template <typename mesh_t, typename crd_t>
//bool geom_static_sensor<mesh_t, crd_t>::compute(const data_type& data, Vector_t<E_Int>& adap_incr)
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
//  K_CONNECT::EltAlgo<typename mesh_t::elt_type>::smoothd1(_hmesh._ng.PHs, _hmesh._F2E, _hmesh._PHtree.level(), adap_incr); // 2-1 smoothing
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

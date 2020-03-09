/*
 
 
 
              NUGA 
 
 
 
 */
//Authors : Sâm Landier (sam.landier@onera.fr), Alexis Gay (alexis.gay@onera.fr), Alexis Rouil (alexis.rouil@onera.fr)

#ifndef NUGA_GEOM_SENSOR_HXX
#define NUGA_GEOM_SENSOR_HXX

#include "Connect/EltAlgo.h"
#include "MeshElement/Hexahedron.h"
#include "Fld/ngon_t.hxx"
#include <limits.h>
#ifdef DEBUG_2019
#include <chrono>
#include <ctime>
#endif
using ngon_type = ngon_t<K_FLD::IntArray>; 

//#define NB_PTS_TO_INCR(nb_pts) (nb_pts/4)
//#define NB_PTS_TO_INCR(nb_pts) ((nb_pts > 0 ) ? nb_pts -1 : 0) // static_sensor : n-1 subdivision where n is the number of source pts in a cell
#define NEIGHBOR(PHi, _F2E, PGi) ( (_F2E(0,PGi) == PHi) ? _F2E(1,PGi) : _F2E(0,PGi) )

namespace NUGA
{
  struct dir_incr_type
  {
    std::vector<E_Int>      _adap_incr;
    std::vector<NUGA::eDIR> _ph_dir, _pg_dir;
    
    E_Int& operator[](E_Int i) {return _adap_incr[i]; };
    
    void clear() {_adap_incr.clear();}
    void resize(E_Int sz, E_Int val){_adap_incr.resize(sz, val);}
    std::vector<E_Int>::iterator begin() {return _adap_incr.begin();}
    std::vector<E_Int>::iterator end() {return _adap_incr.end();}
  };
  
  template <eSUBDIV_TYPE STYPE>
  struct adap_incr_type
  {};
  
  template <>
  struct adap_incr_type<ISO>
  {
    using output_type = Vector_t<E_Int>;
  };
  
  template <>
  struct adap_incr_type<DIR>
  {
    using output_type = dir_incr_type;
  };
  
  
  

/// Octal geometric sensor
template <typename mesh_t, typename crd_t = K_FLD::FloatArray> //ngu for surfacic (PGs) or ngon_t for volumic
class geom_sensor
{
  public:
    
    using data_type = crd_t; //point cloud
    static constexpr eSUBDIV_TYPE sub_type = mesh_t::sub_type;
    using output_type   = typename mesh_t::output_type;
    using BbTree3D_type = typename K_SEARCH::BbTree3D;
    
    geom_sensor(mesh_t& mesh, E_Int max_pts_per_cell = 1, E_Int itermax = -1, const void* dummy=nullptr)
      :_hmesh(mesh), _is_init(false), _refine(true), _has_agglo(false), _agglo(true), _max_pts_per_cell(max_pts_per_cell), _iter_max((itermax <= 0) ? INT_MAX : itermax), _iter(0)
    {
      _bbtree = new K_SEARCH::BbTree3D(*_hmesh._crd, *_hmesh._ng, E_EPSILON);
    }
    
    virtual E_Int assign_data(data_type& data);
    
    virtual bool compute(output_type& adap_incr, bool do_agglo);

    void fill_adap_incr(output_type& adap_incr, bool do_agglo);
    
    void fix_adap_incr(dir_incr_type& adap_incr); //DIR : output_type is 
    void fix_adap_incr(std::vector<E_Int>& adap_incr); //ISO
    
    void locate_points();
    
    E_Int get_highest_lvl_cell(const E_Float* p, E_Int& PHi);
    
    E_Int get_higher_lvl_cell(const E_Float* p, E_Int PHi);
    
    E_Int detect_child(const E_Float* p, E_Int PHi, E_Int* children);

    virtual E_Int redistrib_data();

    ~geom_sensor()
    {
      delete _bbtree;
    }

#ifdef DEBUG_2019
    void verif();
    E_Int verif2();
#endif
  
  protected:
    mesh_t& _hmesh;
    data_type* _data;
    BbTree3D_type* _bbtree;
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
E_Int geom_sensor<mesh_t, crd_t>::assign_data(data_type& data) 
{

  _data = &data ;
    
  // fill in _points_to_cell with the initial mesh (giving for each source point the highest level cell containing it)
  _points_to_cell.clear();
  
  // Count the number of data per cell

  // Check if hmesh has been initialized
  if (not _hmesh.is_initialised())
  {
    std::cout << "Error in geom_sensor.hxx / function assign_data: hmesh has not been initialised." << std::endl; 
    return 1;
  }
  
  locate_points();
  
  return 0;
}

#ifdef DEBUG_2019
template <typename mesh_t, typename crd_t>
void geom_sensor<mesh_t, crd_t>::verif()
{
  E_Int err(0);  
  E_Int nb_elts = _hmesh._ng->PHs.size();
  E_Int nb_pts = _points_to_cell.size();
  //points_to_cell gives the number of points per cell
  Vector_t<E_Int> nb_pts_per_cell(nb_elts,0);
  
  for (int i = 0; i < nb_pts; ++i)
  {
    E_Int PHi = _points_to_cell[i];
    if (PHi != E_IDX_NONE)
      nb_pts_per_cell[PHi] += 1;
  }

  E_Int pyra(0);
  for (int i=0; i< nb_elts; i++){
    if (nb_pts_per_cell[i]>1){
      err += 1;
      if (K_MESH::Polyhedron<0>::is_PY5(_hmesh._ng->PGs, _hmesh._ng->PHs.get_facets_ptr(i), _hmesh._ng->PHs.stride(i))){
        pyra += 1;
      }
    }
  }
  
  std::cout << "err= " << err <<std::endl;
  std::cout << "pyra= " << pyra <<std::endl;
}
///
template <typename mesh_t, typename crd_t>
E_Int geom_sensor<mesh_t, crd_t>::verif2()
{
  E_Int err(0);  
  E_Int nb_elts = _hmesh._ng->PHs.size();
  E_Int nb_pts = _points_to_cell.size();
  //points_to_cell gives the number of points per cell
  Vector_t<E_Int> nb_pts_per_cell(nb_elts,0);
  
  for (int i = 0; i < nb_pts; ++i)
  {
    E_Int PHi = _points_to_cell[i];
    if (PHi != E_IDX_NONE)
      nb_pts_per_cell[PHi] += 1;
  }

  for (int i=0; i< nb_elts; i++){
      if (nb_pts_per_cell[i]>1) err += 1; 
  }
  
  std::cout << "err= " << err <<std::endl;
  return err;
}
#endif

/// 
template <typename mesh_t, typename crd_t>
bool geom_sensor<mesh_t, crd_t>::compute(output_type& adap_incr, bool do_agglo)
{

#ifdef FLAG_STEP
  std::cout << "iter : " << _iter << ". nb of PHs : " <<  _hmesh._ng->PHs.size() << std::endl;
#endif
  
  if (++_iter > _iter_max) return false;
  
  E_Int nb_elts = _hmesh._ng->PHs.size();
  E_Int nb_pts = _points_to_cell.size();
  _refine=false;
  _agglo = false;
  _has_agglo = false;
  _is_init = true;

  //fill in adap incr thanks to points to cell
  fill_adap_incr(adap_incr, do_agglo);
  
  //apply the 2:1 rule
  _hmesh.smooth(adap_incr);
  
  // fix inconsistencies (does something only in DIR mode)
  fix_adap_incr(adap_incr);
  
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
void geom_sensor<mesh_t, crd_t>::fill_adap_incr(output_type& adap_incr, bool do_agglo)
{
  E_Int nb_elts = _hmesh._ng->PHs.size();
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
    const E_Int* faces = _hmesh._ng->PHs.get_facets_ptr(i);
    E_Int nb_faces = _hmesh._ng->PHs.stride(i);
    bool admissible_elt = K_MESH::Polyhedron<0>::is_HX8(_hmesh._ng->PGs, faces, nb_faces) || K_MESH::Polyhedron<0>::is_TH4(_hmesh._ng->PGs, faces, nb_faces)
                        || K_MESH::Polyhedron<0>::is_PR6(_hmesh._ng->PGs, faces, nb_faces) || K_MESH::Polyhedron<0>::is_PY5(_hmesh._ng->PGs, faces, nb_faces);

    if ( admissible_elt && (nb_pts_per_cell[i] >= _max_pts_per_cell + 1) && (_hmesh._PHtree.is_enabled(i)) ) // can be and has to be subdivided
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

/// ISO
template <typename mesh_t, typename crd_t>
void geom_sensor<mesh_t, crd_t>::fix_adap_incr(std::vector<E_Int>& adap_incr)
{
  E_Int nb_phs = _hmesh._ng->PHs.size();
  // disable unhandled_elements
  for (E_Int i = 0; i < nb_phs; ++i)
  {
    if (adap_incr[i] == 0) continue;

    E_Int nb_faces = _hmesh._ng->PHs.stride(i); 
    E_Int* faces = _hmesh._ng->PHs.get_facets_ptr(i);
    bool admissible_elt = K_MESH::Polyhedron<0>::is_HX8(_hmesh._ng->PGs, faces, nb_faces) || K_MESH::Polyhedron<0>::is_TH4(_hmesh._ng->PGs, faces, nb_faces)
                           || K_MESH::Polyhedron<0>::is_PR6(_hmesh._ng->PGs, faces, nb_faces) || K_MESH::Polyhedron<0>::is_PY5(_hmesh._ng->PGs, faces, nb_faces);
    
    if (!admissible_elt || !_hmesh._PHtree.is_enabled(i))
      adap_incr[i] = 0;
  }
}
///
template <typename mesh_t, typename crd_t>
void geom_sensor<mesh_t, crd_t>::fix_adap_incr(dir_incr_type& adap_incr)
{
  E_Int nb_phs = _hmesh._ng->PHs.size();
  
  // disable unhandled_elements
  for (E_Int i = 0; i < nb_phs; ++i)
  {
    if (adap_incr[i] == 0) continue;

    E_Int nb_faces = _hmesh._ng->PHs.stride(i); 
    E_Int* faces = _hmesh._ng->PHs.get_facets_ptr(i);
    bool admissible_elt = K_MESH::Polyhedron<0>::is_HX8(_hmesh._ng->PGs, faces, nb_faces) || K_MESH::Polyhedron<0>::is_TH4(_hmesh._ng->PGs, faces, nb_faces)
                           || K_MESH::Polyhedron<0>::is_PR6(_hmesh._ng->PGs, faces, nb_faces) || K_MESH::Polyhedron<0>::is_PY5(_hmesh._ng->PGs, faces, nb_faces);
    
    if (!admissible_elt || !_hmesh._PHtree.is_enabled(i))
      adap_incr[i] = 0;
  }
      
  //solve inconsistencies
  //alexis : todo
  
  
  
  // premiere version : desactiver l'adaptation dans les cellules XYZ (iso) qui sont connectées à une cellule "layer" par un QUAD lateral
  adap_incr._ph_dir.clear();
  adap_incr._ph_dir.resize(_hmesh._ng->PHs.size(), XYZ);
  
  for (E_Int i=0; i < nb_phs; ++i)
  {
    // if type is layer => adap_incr._ph_dir[i] = XY
  }
  
  E_Int nb_pgs = _hmesh._ng->PGs.size();
  adap_incr._pg_dir.clear();
  adap_incr._pg_dir.resize(nb_pgs, NONE);
  
  
  // boucle sur les layer
  
    // appel à get_local
  
    // remplissage approprié de adap_incr._pg_dir pour les 6 faces : X, Y, ou XY
  
  // boucle sur le reste : 
  
    // appel à get_local
  
    // remplissage de adap_incr._pg_dir ou desactivation de adap_incr._ph_dir
  
       // si NONE => remplissage
       // sinon, si valeur différente => desactivation
}
///
template <typename mesh_t, typename crd_t>
void geom_sensor<mesh_t, crd_t>::locate_points()
{
  // locate every source points in a given mesh
  using ELT_t = typename mesh_t::elt_type;

  E_Float tol = 10. * ZERO_M;
  
  E_Int nb_src_pts = _data->cols();
  // for each source points, give the highest cell containing it if there is one, otherwise E_IDX_NONE
  _points_to_cell.resize(nb_src_pts, E_IDX_NONE);
  
  Vector_t<E_Int> ids;

  bool found=false;
  E_Float minB[3];
  E_Float maxB[3];
  //
  for (int i = 0; i < nb_src_pts; i++)
  {
    const E_Float* p = _data->col(i);
    
    for (int j = 0; j < 3;j++)
    {
      minB[j] = p[j]-E_EPSILON; // small box over the source point
      maxB[j] = p[j]+E_EPSILON;
    }  
    ids.clear();
    _bbtree->getOverlappingBoxes(minB,maxB,ids); // ids contains the list of PH index partially containing the small box around the source point
    
    if (ids.empty()) continue; // exterior point : E_IDX_NONE
    
    for (size_t j = 0; j < ids.size(); j++) // which of these boxes has the source point ?
    {
      
      if (K_MESH::Polyhedron<UNKNOWN>::pt_is_inside(ids[j], _hmesh._ng->PGs, _hmesh._ng->PHs.get_facets_ptr(ids[j]), _hmesh._ng->PHs.stride(ids[j]), *_hmesh._crd, _hmesh._F2E, p, tol))
      {
        if (_hmesh._PHtree.children(ids[j]) == nullptr) // if the cell has no children : source point i is in this cell
          _points_to_cell[i] = ids[j];
        else // if the cell has children : find the highest level (i.e. smallest) cell containing the source point
        {
          E_Int cell = get_highest_lvl_cell(p, ids[j]);
          _points_to_cell[i] = cell;
        }

        found = true;
        break;
      }
    }
  }

  if (!found) std::cout << "WARNING : All the source points are outside the mesh. Nothing to do." << std::endl;

}

///
template <typename mesh_t, typename crd_t>
E_Int geom_sensor<mesh_t, crd_t>::get_higher_lvl_cell(const E_Float* p, E_Int PHi)
{
  // returns in the PHi's sons where a given source point is
  using ELT_t = typename mesh_t::elt_type;
  
  E_Int* q = _hmesh._PHtree.children(PHi); // the children of PHi
  E_Int nb_children = _hmesh._PHtree.nb_children(PHi);
  bool found = false;
    
  for (int j = 0; j < nb_children; ++j)
  {
    if (K_MESH::Polyhedron<UNKNOWN>::pt_is_inside(q[j], _hmesh._ng->PGs, _hmesh._ng->PHs.get_facets_ptr(q[j]), _hmesh._ng->PHs.stride(q[j]), *_hmesh._crd, _hmesh._F2E, p, 1.e-14)) // true when the point is located in a child
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
E_Int geom_sensor<mesh_t, crd_t>::get_highest_lvl_cell(const E_Float* p, E_Int& PHi)
{
  // will find, in the tree (starting at PHi cell), the highest level cell (ie smallest cell) cointaining a given source point 
  using ELT_t = typename mesh_t::elt_type;

  // Do not test children if PHi is enabled
  if (_hmesh._PHtree.is_enabled(PHi)) return PHi;

  E_Int nb_children = _hmesh._PHtree.nb_children(PHi);
  assert(nb_children != 0);
  
  E_Int* q = _hmesh._PHtree.children(PHi); // the children of PHi

  bool found = false;
  for (int j = 0; j < nb_children; ++j)
  {
    if (K_MESH::Polyhedron<UNKNOWN>::pt_is_inside(q[j], _hmesh._ng->PGs, _hmesh._ng->PHs.get_facets_ptr(q[j]), _hmesh._ng->PHs.stride(q[j]), *_hmesh._crd, _hmesh._F2E, p, 1.e-14)) // true when the point is located in a child
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
E_Int geom_sensor<mesh_t, crd_t>::detect_child(const E_Float* p, E_Int PHi, E_Int* children)
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

///
template <typename mesh_t, typename crd_t>
E_Int geom_sensor<mesh_t, crd_t>::redistrib_data()
{
  std::cout << "redistrib data geom sensor. " << std::endl;
  E_Int nb_pts = _points_to_cell.size();

  for (int i = 0; i < nb_pts; ++i)
  {
    E_Int PHi = _points_to_cell[i];

    if (PHi == E_IDX_NONE) continue; //external point
    if (_hmesh._PHtree.children(PHi) != nullptr) // we detect the subdivised PH : one of his children has this source point
    {
      E_Float* p = _data->col(i); // in which children of PHi is the point
      E_Int cell = get_higher_lvl_cell(p,PHi);
      _points_to_cell[i] = cell;
    }
  }
  return 0;
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
//  E_Int nb_elts = _hmesh._ng->PHs.size();
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
//  K_CONNECT::EltAlgo<typename mesh_t::elt_type>::smoothd1(_hmesh._ng->PHs, _hmesh._F2E, _hmesh._PHtree.level(), adap_incr); // 2-1 smoothing
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

template <typename mesh_t, typename crd_t = K_FLD::FloatArray>
class geom_sensor2 : public geom_sensor<mesh_t, crd_t>
{
    public:
        using parent_t = geom_sensor<mesh_t, crd_t>;
        Vector_t<E_Int> _Ln;       
        
        
        
        //////////
        geom_sensor2(mesh_t& mesh, E_Int max_pts_per_cell = 1, E_Int itermax = -1, const void* dummy=nullptr): parent_t(mesh, max_pts_per_cell, itermax){}
        void init2();
        bool compute(Vector_t<E_Int>& adap_incr, bool do_agglo);
        void update(E_Int i);
};

template <typename mesh_t, typename crd_t>
void geom_sensor2<mesh_t, crd_t>::init2()
{
    E_Int n_nodes= parent_t::_hmesh._crd->cols();
    _Ln.resize(n_nodes);
}

template <typename mesh_t, typename crd_t>
bool geom_sensor2<mesh_t, crd_t>::compute(Vector_t<E_Int>& adap_incr, bool do_agglo)
{
#ifdef DEBUG_2019    
  auto start = std::chrono::system_clock::now();
#endif

  init2();
    
  if (++parent_t::_iter > parent_t::_iter_max) return false;
  
  E_Int nb_pts = parent_t::_points_to_cell.size();
  parent_t::_refine=false;
  parent_t::_agglo = false;
  parent_t::_has_agglo = false;
  
  if (parent_t::_is_init) // first iteration : init is done instead
  {
    //update points_to_cell : replace the PH that have been subdivided in adapt by the child(ren) containing the source point(s)
    for (int i = 0; i < nb_pts; ++i)
    {
      E_Int PHi = parent_t::_points_to_cell[i];

      if (PHi == E_IDX_NONE) continue; //external point
      if (parent_t::_hmesh._PHtree.children(PHi) != nullptr) // we detect the subdivised PH : one of his children has this source point
      {
        E_Float* p = parent_t::_data->col(i); // in which children of PHi is the point
        E_Int cell = parent_t::get_higher_lvl_cell(p,PHi);
        parent_t::_points_to_cell[i] = cell;
      }
    }
  }

  parent_t::_is_init = true;

  //fill in adap incr thanks to points to cell
  parent_t::fill_adap_incr(adap_incr, do_agglo);
  
  for (int i=0; i<adap_incr.size(); i++){
      if (adap_incr[i]==1)
      update(i);
  }
  
  // NODAL SMOOTHING LOOP
  E_Int Lc;
  E_Int n(0);
  bool carry_on=true;
  while (carry_on)
  {
    n=n+1;  
    carry_on=false;
    for (int i=0; i<parent_t::_hmesh._ng->PHs.size(); i++){
      E_Int Mnodes(0);
      if (adap_incr[i]) continue;
      if (parent_t::_hmesh._PHtree.is_enabled(i)){
      Lc= parent_t::_hmesh._PHtree.get_level(i);
      E_Int* pPHi= parent_t::_hmesh._ng->PHs.get_facets_ptr(i);
      
      // Mnodes definition
      for (int j=0; j<parent_t::_hmesh._ng->PHs.stride(i); j++){
        E_Int PGi = *(pPHi +j) - 1;
        E_Int PHn = NEIGHBOR(i, parent_t::_hmesh._F2E, PGi);
        if (!(PHn == E_IDX_NONE) && !(parent_t::_hmesh._PHtree.is_enabled(PHn)) && !(parent_t::_hmesh._PHtree.get_level(PHn)==0) && !(parent_t::_hmesh._PHtree.is_enabled(parent_t::_hmesh._PHtree.parent(PHn)))){
          E_Int nb_ch= parent_t::_hmesh._PGtree.nb_children(PGi);
          for (int k=0; k< nb_ch; k++){
            E_Int* enf= parent_t::_hmesh._PGtree.children(PGi);
            E_Int stride= parent_t::_hmesh._ng->PGs.stride(*(enf+k)); 
            E_Int* ind = parent_t::_hmesh._ng->PGs.get_facets_ptr(*(enf+k));
            for (int l=0; l<stride; l++){
              Mnodes = std::max(_Ln[*(ind+l)-1],Mnodes);
            }
          }
        }
        else {
          E_Int stride = parent_t::_hmesh._ng->PGs.stride(PGi);
          E_Int *ind = parent_t::_hmesh._ng->PGs.get_facets_ptr(PGi);
          for (int l=0; l< stride; l++){
            Mnodes = std::max(_Ln[*(ind+l)-1],Mnodes);
          }
        }
      }
      // fin Mnodes definition
      if (Lc < Mnodes-1)
      {
        const E_Int* faces = parent_t::_hmesh._ng->PHs.get_facets_ptr(i);
        E_Int nb_faces = parent_t::_hmesh._ng->PHs.stride(i);
        bool admissible_elt = K_MESH::Polyhedron<0>::is_HX8(parent_t::_hmesh._ng->PGs, faces, nb_faces) || K_MESH::Polyhedron<0>::is_TH4(parent_t::_hmesh._ng->PGs, faces, nb_faces)
                             || K_MESH::Polyhedron<0>::is_PR6(parent_t::_hmesh._ng->PGs, faces, nb_faces) || K_MESH::Polyhedron<0>::is_PY5(parent_t::_hmesh._ng->PGs, faces, nb_faces);
    
        if (!admissible_elt)
          continue;
        
        adap_incr[i] = 1;
        update(i);        
        carry_on= true;
      }
    }
  }
}
#ifdef DEBUG_2019
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> t1 = end-start;
  std::cout << "tps compute geom2= " << t1.count() << "s"<< std::endl;
#endif

  for (int i=0; i< adap_incr.size(); i++){
      if (adap_incr[i]==1){
          return true;
      }
  }  
  return false;
}

template <typename mesh_t, typename crd_t>
void geom_sensor2<mesh_t, crd_t>::update(E_Int i)
{         
  E_Int PHi = i;
  E_Int Lc= parent_t::_hmesh._PHtree.get_level(PHi);
  E_Int nb_faces = parent_t::_hmesh._ng->PHs.stride(PHi); 
  for (int j=0; j< nb_faces; j++){    
    E_Int* PGj= parent_t::_hmesh._ng->PHs.get_facets_ptr(PHi);
    E_Int PG = PGj[j] - 1;
    E_Int* pN = parent_t::_hmesh._ng->PGs.get_facets_ptr(PG);
    E_Int nb_nodes = parent_t::_hmesh._ng->PGs.stride(PG);
    for (int l=0; l< nb_nodes; l++){
      _Ln[*(pN+l)-1]=std::max(_Ln[*(pN+l)-1], Lc+1);
    }
  }
}

//////////////////////////////////////////////
template <typename mesh_t, typename crd_t = K_FLD::FloatArray>
class geom_sensor3 : public geom_sensor<mesh_t, crd_t>
{
    public:
        using parent_t = geom_sensor<mesh_t, crd_t>;
        Vector_t<E_Int> _Ln;       
        using data_type = Vector_t<E_Int>;
        
        
        //////////
        geom_sensor3(mesh_t& mesh, E_Int max_pts_per_cell = 1, E_Int itermax = -1, const void* dummy=nullptr): parent_t(mesh, max_pts_per_cell, itermax){}
        E_Int assign_data(data_type level);
        bool compute(Vector_t<E_Int>& adap_incr, bool do_agglo);
        E_Int redistrib_data();
};

template <typename mesh_t, typename crd_t>
E_Int geom_sensor3<mesh_t, crd_t>::assign_data(data_type level)
{
    E_Int n_nodes= level.size(); 
    _Ln.resize(n_nodes);

    for (int i=0; i< n_nodes; i++){
        _Ln[i]= level[i];
    }
    return 0;
    
}      

template <typename mesh_t, typename crd_t>
bool geom_sensor3<mesh_t, crd_t>::compute(Vector_t<E_Int>& adap_incr, bool do_agglo)
{
  // E_Int n_nodes= parent_t::_hmesh._crd->size();
  // _Ln.resize(n_nodes, 0);
  E_Int nb_elt= parent_t::_hmesh._ng->PHs.size();
  bool flag(false);
  adap_incr.clear();
  adap_incr.resize(nb_elt, 0);

  for (int i=0; i< nb_elt; i++){
    if (parent_t::_hmesh._PHtree.is_enabled(i)){
      E_Int* faces= parent_t::_hmesh._ng->PHs.get_facets_ptr(i);
      E_Int n_faces= parent_t::_hmesh._ng->PHs.stride(i);
      for (int k=0; k< n_faces; k++){
        E_Int PGk= *(faces+k)-1;
        E_Int* pN= parent_t::_hmesh._ng->PGs.get_facets_ptr(PGk);
        E_Int n_face_nodes = parent_t::_hmesh._ng->PGs.stride(PGk);
        
        for (int l=0; l< n_face_nodes; l++){
          E_Int nodes_faces= *(pN+l)-1;
          if (_Ln[nodes_faces]>0){
            adap_incr[i]= 1;
            flag=true;
            break;
          }
        }
        if (flag==true) break;
      }
    }
    flag=false;
  }
    
  for (int i=0; i< adap_incr.size(); i++){
    if (adap_incr[i]==1){
      return true;
    }
  }  
  return false;
}


///
template <typename mesh_t, typename crd_t>
E_Int geom_sensor3<mesh_t, crd_t>::redistrib_data()
{
  E_Int n_nodes= parent_t::_hmesh._crd->size();
  _Ln.resize(n_nodes, 0);
  
  for (int i=0; i< n_nodes; i++)
  {
    if (_Ln[i]>0)
      _Ln[i]--;
  }
  return 0;
}
  





}
#endif

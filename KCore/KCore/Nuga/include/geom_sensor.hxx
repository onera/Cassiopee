/*
 
 
 
              NUGA 
 
 
 
 */
//Authors : Sâm Landier (sam.landier@onera.fr), Alexis Gay (alexis.gay@onera.fr), Alexis Rouil (alexis.rouil@onera.fr)

#ifndef NUGA_GEOM_SENSOR_HXX
#define NUGA_GEOM_SENSOR_HXX

#include "Nuga/include/sensor.hxx"
#include "Connect/EltAlgo.h"
#include "MeshElement/Hexahedron.h"
#include "Fld/ngon_t.hxx"
#include <limits.h>
#ifdef DEBUG_2019
#include <chrono>
#include <ctime>
#endif
using ngon_type = ngon_t<K_FLD::IntArray>; 

#define NEIGHBOR(PHi, _F2E, PGi) ( (_F2E(0,PGi) == PHi) ? _F2E(1,PGi) : _F2E(0,PGi) )

namespace NUGA
{
  
  using crd_t = K_FLD::FloatArray;

/// Geometric sensor
template <typename mesh_t>
class geom_sensor : public sensor<mesh_t, crd_t>
{
  public:

    using sensor_data_t = crd_t;
    using parent_t = sensor<mesh_t, crd_t>;
    using sensor_output_t   = typename mesh_t::sensor_output_t;
 
    using BbTree3D_type = typename K_SEARCH::BbTree3D;
    
    geom_sensor(mesh_t& mesh, E_Int max_pts_per_cell = 1, E_Int itermax = -1)
      :parent_t(mesh), _is_init(false), _refine(true), _has_agglo(false), _agglo(true), _max_pts_per_cell(max_pts_per_cell), _iter_max((itermax <= 0) ? INT_MAX : itermax), _iter(0)
    {
      _bbtree = new K_SEARCH::BbTree3D(*parent_t::_hmesh._crd, *parent_t::_hmesh._ng, E_EPSILON);
    }
    
    virtual E_Int assign_data(sensor_data_t& data) override;
    
    bool compute(sensor_output_t& adap_incr, bool do_agglo) override;

    void fill_adap_incr(sensor_output_t& adap_incr, bool do_agglo) ;
    
    void locate_points();
    
    E_Int get_highest_lvl_cell(const E_Float* p, E_Int PHi) const ;
    
    E_Int get_higher_lvl_cell(const E_Float* p, E_Int PHi) const ;
    
    E_Int detect_child(const E_Float* p, E_Int PHi, const E_Int* children) const ;

    virtual void update() override;

    ~geom_sensor()
    {
      delete _bbtree;
    }

  protected:
    
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
template <typename mesh_t>
E_Int geom_sensor<mesh_t>::assign_data(sensor_data_t& data) 
{

  parent_t::assign_data(data);
    
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
template <typename mesh_t>
bool geom_sensor<mesh_t>::compute(sensor_output_t& adap_incr, bool do_agglo)
{

#ifdef FLAG_STEP
  std::cout << "iter : " << _iter << ". nb of PHs : " <<  parent_t::_hmesh._ng->PHs.size() << std::endl;
#endif
  
  if (++_iter > _iter_max) return false;
  
  E_Int nb_elts = parent_t::_hmesh._ng->PHs.size();
  E_Int nb_pts = _points_to_cell.size();
  _refine=false;
  _agglo = false;
  _has_agglo = false;
  _is_init = true;

  //fill in adap incr thanks to points to cell
  fill_adap_incr(adap_incr, do_agglo);
  
  //apply the 2:1 rule
  parent_t::_hmesh.smooth(adap_incr);
  
  // fix inconsistencies
  fix_adap_incr(parent_t::_hmesh, adap_incr);
  
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
template <typename mesh_t>
void geom_sensor<mesh_t>::fill_adap_incr(sensor_output_t& adap_incr, bool do_agglo)
{
  E_Int nb_elts = parent_t::_hmesh._ng->PHs.size();
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
    const E_Int* faces = parent_t::_hmesh._ng->PHs.get_facets_ptr(i);
    E_Int nb_faces = parent_t::_hmesh._ng->PHs.stride(i);
    bool admissible_elt = K_MESH::Polyhedron<0>::is_basic(parent_t::_hmesh._ng->PGs, faces, nb_faces);
    if ( admissible_elt && (nb_pts_per_cell[i] >= _max_pts_per_cell + 1) && (parent_t::_hmesh._PHtree.is_enabled(i)) ) // can be and has to be subdivided
    {
      adap_incr[i] = 1;
      _refine = true;
    }
  }

  if (do_agglo)
  {
    for (int i = 0; i < nb_elts; ++i)
    {
      if ((parent_t::_hmesh._PHtree.get_level(i) > 0) && (parent_t::_hmesh._PHtree.is_enabled(i))) // this cell may be agglomerated : check its brothers
      {
        if (adap_incr[i] == -1) continue;
                
        E_Int father = parent_t::_hmesh._PHtree.parent(i);
        const E_Int* p = parent_t::_hmesh._PHtree.children(father);
        E_Int nb_children = parent_t::_hmesh._PHtree.nb_children(father);
        E_Int sum = 0;
                
        for (int k = 0; k < nb_children; ++k)
        {
          sum += nb_pts_per_cell[p[k]];
          if (!parent_t::_hmesh._PHtree.is_enabled(p[k])) sum = _max_pts_per_cell + 1;
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
template <typename mesh_t>
void geom_sensor<mesh_t>::locate_points() 
{
  // locate every source points in a given mesh
 
  E_Float tol = 10. * ZERO_M;
  
  E_Int nb_src_pts =parent_t::_data->cols();
  // for each source points, give the highest cell containing it if there is one, otherwise E_IDX_NONE
  _points_to_cell.resize(nb_src_pts, E_IDX_NONE);
  
  Vector_t<E_Int> ids;

  bool found=false;
  E_Float minB[3];
  E_Float maxB[3];
  //
  for (int i = 0; i < nb_src_pts; i++)
  {
    const E_Float* p =parent_t::_data->col(i);
    
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
      
      if (K_MESH::Polyhedron<UNKNOWN>::pt_is_inside(ids[j], parent_t::_hmesh._ng->PGs, parent_t::_hmesh._ng->PHs.get_facets_ptr(ids[j]), parent_t::_hmesh._ng->PHs.stride(ids[j]), *parent_t::_hmesh._crd, parent_t::_hmesh._F2E, p, tol))
      {
        if (parent_t::_hmesh._PHtree.children(ids[j]) == nullptr) // if the cell has no children : source point i is in this cell
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
template <typename mesh_t>
E_Int geom_sensor<mesh_t>::get_higher_lvl_cell(const E_Float* p, E_Int PHi) const 
{
  // returns in the PHi's sons where a given source point is
  
  const E_Int* q = parent_t::_hmesh._PHtree.children(PHi); // the children of PHi
  E_Int nb_children = parent_t::_hmesh._PHtree.nb_children(PHi);
  bool found = false;
    
  for (int j = 0; j < nb_children; ++j)
  {
    if (K_MESH::Polyhedron<UNKNOWN>::pt_is_inside(q[j], parent_t::_hmesh._ng->PGs, parent_t::_hmesh._ng->PHs.get_facets_ptr(q[j]), parent_t::_hmesh._ng->PHs.stride(q[j]), *parent_t::_hmesh._crd, parent_t::_hmesh._F2E, p, 1.e-14)) // true when the point is located in a child
    {
      found = true;
      return q[j];
    }
  }
  if (!found) return detect_child(p,PHi,q); // if we can't find the child, we have to locate the closest child  (almost never occurs)

  return E_IDX_NONE; // never reached : to silent compilers

}

///
template <typename mesh_t>
E_Int geom_sensor<mesh_t>::get_highest_lvl_cell(const E_Float* p, E_Int PHi) const 
{
  // will find, in the tree (starting at PHi cell), the highest level cell (ie smallest cell) cointaining a given source point

  // Do not test children if PHi is enabled
  if (parent_t::_hmesh._PHtree.is_enabled(PHi)) return PHi;

  E_Int nb_children = parent_t::_hmesh._PHtree.nb_children(PHi);
  assert(nb_children != 0);
  
  const E_Int* q = parent_t::_hmesh._PHtree.children(PHi); // the children of PHi

  bool found = false;
  for (int j = 0; j < nb_children; ++j)
  {
    if (K_MESH::Polyhedron<UNKNOWN>::pt_is_inside(q[j], parent_t::_hmesh._ng->PGs, parent_t::_hmesh._ng->PHs.get_facets_ptr(q[j]), parent_t::_hmesh._ng->PHs.stride(q[j]), *parent_t::_hmesh._crd, parent_t::_hmesh._F2E, p, 1.e-14)) // true when the point is located in a child
    {
      found = true;
      if (parent_t::_hmesh._PHtree.children(q[j]) != nullptr)
        return get_highest_lvl_cell(p, q[j]);
      else
        return q[j];
    }
  }
  if (!found) return detect_child(p, PHi, q);
  
  return E_IDX_NONE; // never reached : to silent compilers
    
}

///
template <typename mesh_t>
E_Int geom_sensor<mesh_t>::detect_child(const E_Float* p, E_Int PHi, const E_Int* children) const
{
  E_Int nb_children = parent_t::_hmesh._PHtree.nb_children(PHi);
    
  // find the closest child
  E_Int closest_child = 0;
  E_Float center[3];
  parent_t::_hmesh.get_cell_center(children[0], center);
  E_Float min = K_FUNC::sqrDistance(p, center,3);

  for (int i = 1; i < nb_children; ++i)
  {
    E_Float tmp[3];
    parent_t::_hmesh.get_cell_center(children[i], tmp);
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
template <typename mesh_t>
void geom_sensor<mesh_t>::update()
{
  //std::cout << "redistrib data geom sensor. " << std::endl;
  E_Int nb_pts = _points_to_cell.size();

  for (int i = 0; i < nb_pts; ++i)
  {
    E_Int PHi = _points_to_cell[i];

    if (PHi == E_IDX_NONE) continue; //external point
    if (parent_t::_hmesh._PHtree.children(PHi) != nullptr) // we detect the subdivised PH : one of his children has this source point
    {
      E_Float* p = parent_t::_data->col(i); // in which children of PHi is the point
      E_Int cell = get_higher_lvl_cell(p,PHi);
      _points_to_cell[i] = cell;
    }
  }
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
//  E_Int nb_elts = parent_t::_hmesh._ng->PHs.size();
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
//  K_CONNECT::EltAlgo<typename mesh_t::elt_type>::smoothd1(parent_t::_hmesh._ng->PHs, parent_t::_hmesh._F2E, parent_t::_hmesh._PHtree.level(), adap_incr); // 2-1 smoothing
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

template <typename mesh_t>
class geom_sensor2 : public geom_sensor<mesh_t>
{
    public:
        using parent_t = geom_sensor<mesh_t>;
        Vector_t<E_Int> _Ln;       
        
        //////////
        geom_sensor2(mesh_t& mesh, E_Int max_pts_per_cell = 1, E_Int itermax = -1, const void* dummy=nullptr): parent_t(mesh, max_pts_per_cell, itermax){}
        void init2();
        bool compute(Vector_t<E_Int>& adap_incr, bool do_agglo) override;
        void update_nodal_levels(E_Int PHi);
};

template <typename mesh_t>
void geom_sensor2<mesh_t>::init2()
{
    E_Int n_nodes= parent_t::_hmesh._crd->cols();
    _Ln.resize(n_nodes);
}

template <typename mesh_t>
bool geom_sensor2<mesh_t>::compute(Vector_t<E_Int>& adap_incr, bool do_agglo)
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
        E_Float* p = parent_t::parent_t::_data->col(i); // in which children of PHi is the point
        E_Int cell = parent_t::get_higher_lvl_cell(p,PHi);
       parent_t::_points_to_cell[i] = cell;
      }
    }
  }

 parent_t::_is_init = true;

  //fill in adap incr thanks to points to cell
 parent_t::fill_adap_incr(adap_incr, do_agglo);
  
  for (int i=0; i<adap_incr.size(); i++)
  {
    if (adap_incr[i]==1) update_nodal_levels(i);
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
      const E_Int* pPHi= parent_t::_hmesh._ng->PHs.get_facets_ptr(i);
      
      // Mnodes definition
      for (int j=0; j<parent_t::_hmesh._ng->PHs.stride(i); j++){
        E_Int PGi = *(pPHi +j) - 1;
        E_Int PHn = NEIGHBOR(i, parent_t::_hmesh._F2E, PGi);
        if (!(PHn == E_IDX_NONE) && !(parent_t::_hmesh._PHtree.is_enabled(PHn)) && !(parent_t::_hmesh._PHtree.get_level(PHn)==0) && !(parent_t::_hmesh._PHtree.is_enabled(parent_t::_hmesh._PHtree.parent(PHn)))){
          E_Int nb_ch= parent_t::_hmesh._PGtree.nb_children(PGi);
          for (int k=0; k< nb_ch; k++){
            const E_Int* enf= parent_t::_hmesh._PGtree.children(PGi);
            E_Int stride= parent_t::_hmesh._ng->PGs.stride(*(enf+k)); 
            const E_Int* ind = parent_t::_hmesh._ng->PGs.get_facets_ptr(*(enf+k));
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
        bool admissible_elt = K_MESH::Polyhedron<0>::is_basic(parent_t::_hmesh._ng->PGs, faces, nb_faces);
 
        if (!admissible_elt)
          continue;
        
        adap_incr[i] = 1;
        update_nodal_levels(i);        
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

template <typename mesh_t>
void geom_sensor2<mesh_t>::update_nodal_levels(E_Int PHi)
{         
  E_Int Lc= parent_t::_hmesh._PHtree.get_level(PHi);
  E_Int nb_faces = parent_t::_hmesh._ng->PHs.stride(PHi); 
  for (int j=0; j< nb_faces; j++){    
    const E_Int* PGj= parent_t::_hmesh._ng->PHs.get_facets_ptr(PHi);
    E_Int PG = PGj[j] - 1;
    const E_Int* pN = parent_t::_hmesh._ng->PGs.get_facets_ptr(PG);
    E_Int nb_nodes = parent_t::_hmesh._ng->PGs.stride(PG);
    for (int l=0; l< nb_nodes; l++){
      _Ln[*(pN+l)-1]=std::max(_Ln[*(pN+l)-1], Lc+1);
    }
  }
}

} // NUGA
#endif

/*
 
 
 
              NUGA 
 
 
 
 */
//Authors : Sâm Landier (sam.landier@onera.fr), Alexis Gay (alexis.gay@onera.fr), Alexis Rouil (alexis.rouil@onera.fr)

#ifndef NUGA_HIERACHICAL_MESH_HXX
#define NUGA_HIERACHICAL_MESH_HXX

#if defined (DEBUG_HIERARCHICAL_MESH) || defined (OUTPUT_ITER_MESH)
#include "IO/io.h"
#include "Nuga/Boolean/NGON_debug.h"
using NGDBG = NGON_debug<K_FLD::FloatArray,K_FLD::IntArray>;
#endif

#include "Nuga/include/tree.hxx"
#include "Connect/IdTool.h"
#include "Nuga/Delaunay/Triangulator.h"
#include "MeshElement/Basic.h"
#include "MeshElement/Prism.h"
#include "Nuga/include/refiner.hxx"

#define NEIGHBOR(PHi, _F2E, PGi) ( (_F2E(0,PGi) == PHi) ? _F2E(1,PGi) : _F2E(0,PGi) )

namespace NUGA
{

using crd_t = K_FLD::FloatArray;

// general case : varying stride (nb of children in hierarchy)
template <typename ELT_t, eSUBDIV_TYPE STYPE>
struct hmesh_trait
{
  using arr_t = typename NBC<ELT_t, STYPE>::arr_t;
  //using arr_t2 = typename NBC<typename ELT_t::boundary_type, STYPE>::arr_t;

  static const E_Int PHNBC = NBC<ELT_t, STYPE>::nbc;
  static const E_Int PGNBC = NBC<typename ELT_t::boundary_type, STYPE>::nbc;
};

// PRISM : fixed nbc but different type of boundaries (T3 or Q4)
template <>
struct hmesh_trait<K_MESH::Prism, ISO>
{
  using arr_t = K_FLD::IntArray;

  static const E_Int PHNBC = 8;
  static const E_Int PGNBC = 4;
};

// BASIC : variable nbc and different type of boundaries (T3 or Q4)
template <>
struct hmesh_trait<K_MESH::Basic, ISO>
{
  using arr_t = ngon_unit;

  static const E_Int PHNBC = -1;
  static const E_Int PGNBC = -1;
};

///
template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t = ngon_type>
class hierarchical_mesh
{
  public:

    using htrait_t = hmesh_trait<ELT_t, STYPE>;
    using arr_t = typename htrait_t::arr_t; //ngon_unit (general case) or IntArray (fixed nb of children : e.g. all basic but pyra in ISO)
    using sensor_output_t = typename sensor_output_data<STYPE>::type;

    static constexpr  eSUBDIV_TYPE SUBTYPE = STYPE;

    crd_t*                    _crd;             // Coordinates
    ngo_t*                    _ng;              // NGON mesh
    tree<arr_t>               _PGtree, _PHtree; // Polygons/Polyhedra hierarchy
    K_FLD::IntArray           _F2E;             // neighborhood data structure : 2 X (nb_pgs) array : for each PG gives the left and right cells ids that share this PG
    bool                      _initialized;     // flag to avoid initialization more than once.
    refiner<ELT_t, STYPE>     _refiner;         //refinement method must stay static, so this object is here to store _ecenter (not relevant in hmesh)

    ///
    hierarchical_mesh(crd_t& crd, ngo_t & ng):_crd(&crd), _ng(&ng), _PGtree(ng.PGs, htrait_t::PGNBC), _PHtree(ng.PHs, htrait_t::PHNBC), _initialized(false){}

    ///
    E_Int init();
    ///
    E_Int relocate (crd_t& crd, ngo_t & ng) {
      _crd = &crd;
      _ng = &ng;
      _PGtree.set_entities(ng.PGs);
      _PHtree.set_entities(ng.PHs);

      if (ng.PGs.size() != _PGtree.size())
        return 1; //trying to relocate on the wrong mesh
      if (ng.PHs.size() != _PHtree.size())
        return 1; //trying to relocate on the wrong mesh
      
      return 0;
    }
    ///
    E_Int adapt(sensor_output_t& adap_incr, bool do_agglo);
  
    /// face-conformity
    void conformize(Vector_t<E_Int>& pgoids);
    /// Keep only enabled PHs
    void extract_enabled(ngon_type& filtered_ng) const ;
    
    ///
    void get_cell_center(E_Int PHi, E_Float* center) const ;
    ///
    void get_enabled_neighbours(E_Int PHi, E_Int* neighbours, E_Int& nb_neighbours) const ;
    ///
    void get_higher_level_neighbours(E_Int PHi, E_Int PGi, E_Int* neighbours, E_Int& nb_neighbours) const ;
    ///
    bool enabled_neighbours(E_Int PHi) const ; // return true if the PHi-th PH only has enabled neighbours
    ///
    void smooth(sensor_output_t& adap_incr) const ;
    ///
    bool is_initialised() const ;
  
    ///
    void __init(ngo_t& ng, const K_FLD::IntArray& F2E);

    ///
    void __conformize_next_lvl(Vector_t<E_Int>& molec, E_Int PGi, E_Int i) const ;

};

// default implementation for any basci element in ISO mode
template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t>
void hierarchical_mesh<ELT_t, STYPE, ngo_t>::__init(ngo_t& ng, const K_FLD::IntArray& F2E)
{
  // sort the PGs of the initial NGON
  E_Int nb_phs = _ng->PHs.size();
  for (E_Int i = 0; i < nb_phs; ++i)
  { 
    ELT_t::reorder_pgs(*_ng,_F2E,i);
  }
}

template <>
void hierarchical_mesh<K_MESH::Hexahedron, DIR, ngon_type>::__init(ngon_type& ng, const K_FLD::IntArray& F2E)
{
  // alexis : todo
  
  E_Int nb_phs = _ng->PHs.size();
  
  //type init
  ng.PHs._type.clear();
  ng.PHs._type.resize(nb_phs, (E_Int)K_MESH::Polyhedron<0>::eType::HEXA);
  
  // first reorder : by opposite pair
  E_Int generators[2], HX6opposites[6], *pg(generators), *qopp(HX6opposites);
  for (E_Int i=0; i < nb_phs; ++i)
  {
    K_MESH::Polyhedron<0>::is_prismN(ng.PGs, ng.PHs.get_facets_ptr(i), ng.PHs.stride(i), pg, qopp);
    // alexis : utiliser pg et qopp pour réordonner
  }
  
  // promote eventually to layer type + 2nd reorder
  for (E_Int i=0; i < _F2E.cols(); ++i)
  {
    // alexis : si i n'est pas une frontière => continue
    
    E_Int PHi = (_F2E(0,i) != E_IDX_NONE) ? _F2E(0,i) : _F2E(1,i);
    
    E_Int PHcur = PHi;
//    E_Int Basecur = i;
    
    while (true) // climb over the layer
    {
      if (_ng->PHs._type[PHcur] == K_MESH::Polyhedron<0>::eType::LAYER) //already set
        break;
      
      // alexis :mettre la paire de i en premier, i en premier dans la pair
      
      // alexis : verifier l'anisotropie
      
      // alexis : si pas anisotrope => break
      
      // alexis : faire Basecur := 2nd de la pair (top)
      // alexis :marquer ng.PHs.type[PHcur] := LAYER 
      // alexis :faire PHcur = le voisin qui partage le top
    }; 
  }
  
  // third reorder 
  for (E_Int i = 0; i < nb_phs; ++i)
  {
    K_MESH::Hexahedron::reorder_pgs(*_ng,_F2E,i);
  }
  
  //detect now any layer cell
//  double aniso_ratio = 0.2; //fixme : parametre a externaliser
  
  //to do : si aniso => faire  ng.PHs._type =  K_MESH::Polyhedron<0>::eType::LAYER
  
}

///
template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t>
E_Int hierarchical_mesh<ELT_t, STYPE, ngo_t>::init()
{
  if (_initialized) return 0;
  
  E_Int err(0);
  
  // We reorient the PG of our NGON
  _ng->flag_externals(1);
  DELAUNAY::Triangulator dt;
  bool has_been_reversed;
  err = ngon_type::reorient_skins(dt, *_crd, *_ng, has_been_reversed);
  if (err)
    return 1;
  
  //F2E
  ngon_unit neighbors;
  _ng->build_ph_neighborhood(neighbors);
  _ng->build_F2E(neighbors, _F2E);
  
  __init(*_ng, _F2E); // for pure type == reorder_pgs

  _initialized = true;
  
  return err;
}

#ifdef OUTPUT_ITER_MESH
 static E_Int iter = 1;
#endif

template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t>
E_Int hierarchical_mesh<ELT_t, STYPE, ngo_t>::adapt(sensor_output_t& adap_incr, bool do_agglo)
{
  E_Int err(0);

  // identify the sensor
  E_Int max = *std::max_element(ALL(adap_incr));
  E_Int min = (!do_agglo) ? 0 : *std::min_element(ALL(adap_incr));

  if (max == 0 && min == 0) return 0; // no adaptation required
  
  // adaptation has only one sub iteration ? => yes, break at the end
  bool one_generation = ((abs(min) <= 1) && (abs(max) <= 1)) ? true : false;
    
  while (!err)
  {
    if (max > 0) // there is a need to refine
    {
      // refine Faces : create missing children (those PGi with _PGtree.nb_children(PGi) == 0)
      refiner<ELT_t, STYPE>::refine_Faces(adap_incr, *_ng, _PGtree, *_crd, _F2E, _refiner._ecenter);
      // refine Cells with missing children (those PHi with _PHtree.nb_children(PHi) == 0)
      refiner<ELT_t, STYPE>::refine_PHs(adap_incr, *_ng, _PGtree, _PHtree, *_crd, _F2E);
    }
        
    // enable the right PHs & their levels, disable subdivided PHs
    adap_incr.resize(_ng->PHs.size(),0);
    for (E_Int PHi = 0; PHi < _ng->PHs.size(); ++PHi)
    {
      if (!_PHtree.is_enabled(PHi)) continue;
      
      E_Int& adincrPHi = adap_incr[PHi];
      if (adincrPHi == 0) continue;

      if (adincrPHi > 0) // refinement : activate the children, transfer the adapincr & set their level
      {
        E_Int nb_child = _PHtree.nb_children(PHi);
        const E_Int* children = _PHtree.children(PHi);
        for (E_Int j = 0; j < nb_child; ++j)
        {
          _PHtree.enable(*(children+j));
          adap_incr[*(children+j)] = adincrPHi - 1;

          E_Int lvl_p1 = _PHtree.get_level(PHi) + 1;
          _PHtree.set_level(*(children+j), lvl_p1); //fixme : do it somewhere else ?
        }
        adincrPHi = 0;//reset incr
      }
      else // agglomeration : activate the father, transfer that adap incr
      {
        E_Int father = _PHtree.parent(PHi);
        _PHtree.enable(father);
        adap_incr[father] = adincrPHi + 1;
        // reset incr on children
        E_Int nb_child = _PHtree.nb_children(PHi);
        const E_Int* children = _PHtree.children(PHi);
        for (E_Int j = 0; j < nb_child; ++j)
          adap_incr[*(children + j)] = 0;
      }
    }
        
    _ng->PGs.updateFacets();
    _ng->PHs.updateFacets();

#ifdef DEBUG_HIERARCHICAL_MESH  
    //std::cout << "is hybrid at the end of adap ? " << is_hybrid(_ng) << std::endl;
    //    debug_draw(0, _ng, _PHtree, _crd);   
    if (! _ng->attributes_are_consistent()) return false;
#endif
    
#ifdef OUTPUT_ITER_MESH
    ngon_type filtered_ng;
    extract_enabled(filtered_ng);

    std::ostringstream o;
    o << "NGON_it_" << iter << ".plt"; // we create a file at each iteration
    K_FLD::IntArray cnto;
    filtered_ng.export_to_array(cnto);
    MIO::write(o.str().c_str(), *_crd, cnto, "NGON");
    ++iter;
#endif

    if (one_generation) break;
  }

  return err;
}

template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t>
void hierarchical_mesh<ELT_t, STYPE, ngo_t>::conformize(Vector_t<E_Int>& pgoids)
{
  Vector_t<E_Int> old_pgoids;
  _PGtree.get_oids(old_pgoids);
  
  pgoids.clear();//for history (BC and Join preserving)
  
  ngon_unit new_phs;
  Vector_t<E_Int> molec;

  E_Int nb_phs = _ng->PHs.size();
  for (E_Int i = 0; i < nb_phs; ++i)
  {
    if (!_PHtree.is_enabled(i)) continue;
    
    molec.clear();
    E_Int s = _ng->PHs.stride(i);
    const E_Int* pPGi = _ng->PHs.get_facets_ptr(i);

    for (E_Int j = 0; j < s; ++j)
    {
      E_Int PGi = *(pPGi +j) - 1;
      E_Int PHn = NEIGHBOR(i, _F2E, PGi);
      
      if(PHn == E_IDX_NONE)
        molec.push_back(PGi+1);
      else if (_PHtree.is_enabled(PHn))
        molec.push_back(PGi+1);
      else // father or children ?
      {
        E_Int PHf = _PHtree.parent(i);
        if ((PHf != E_IDX_NONE) && _PHtree.is_enabled(PHf))
          molec.push_back(PGi+1);
        else // append the children
        {
          E_Int nbc = _PGtree.nb_children(PGi);
          for (E_Int c=0; c < nbc; ++c){
            E_Int PG_f= *(_PGtree.children(PGi)+c);
            __conformize_next_lvl(molec, PG_f, PHn);
          }
        }        
      }
    }
    
    new_phs.add(molec.size(), &molec[0]);  //alexis : set _type for children ??
  }

  _ng->PHs = new_phs;
  _ng->PHs.updateFacets();
  
  std::vector<E_Int> pgnids, phnids;
  _ng->remove_unreferenced_pgs(pgnids, phnids);
      
  //
  K_CONNECT::IdTool::init_inc(pgoids, _ng->PGs.size());
  for (size_t i=0; i <  pgnids.size(); ++i)
    if (pgnids[i] != E_IDX_NONE)pgoids[pgnids[i]] = old_pgoids[i]; //old_pgoids cannot have E_IDX_NONE : new entities must be self referring
}

template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t>
void hierarchical_mesh<ELT_t, STYPE, ngo_t>::__conformize_next_lvl(Vector_t<E_Int> &molec, E_Int PGi, E_Int i) const 
{
  E_Int sid =  (_PHtree.parent(_F2E(0,PGi)) == i) ? 0 : 1 ;
  E_Int PHj = _F2E(sid,PGi);

  //E_Int PHj = NEIGHBOR(i, _F2E, PGi);

  if(PHj == E_IDX_NONE)
    molec.push_back(PGi+1);
  else if (_PHtree.is_enabled(PHj))
    molec.push_back(PGi+1);
  else // children 
  {
    E_Int nbc = _PGtree.nb_children(PGi);
    for (E_Int c=0; c < nbc; ++c){
      //molec.push_back(*(_PGtree.children(PGi)+c) + 1);
      E_Int PG_f= *(_PGtree.children(PGi)+c);
      __conformize_next_lvl(molec, PG_f, PHj);
    }
  }
}

///
template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t>
void hierarchical_mesh<ELT_t, STYPE, ngo_t>::extract_enabled(ngon_type& filtered_ng) const
{
  filtered_ng.PGs = _ng->PGs;
        
  //E_Int sz = _PGtree.get_parent_size();
  E_Int sz = _ng->PHs.size();
        
  for (int i = 0; i < sz; i++)
  {
    if (_PHtree.is_enabled(i) == true)
    {
      const E_Int* p = _ng->PHs.get_facets_ptr(i);
      E_Int s = _ng->PHs.stride(i);
      filtered_ng.PHs.add(s,p);//alexis : set _type for children ??
    }
  }
    
  filtered_ng.PGs.updateFacets();
  filtered_ng.PHs.updateFacets();
    
  Vector_t<E_Int> pgnids, phnids;
  filtered_ng.remove_unreferenced_pgs(pgnids, phnids);
}


///
template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t>
void hierarchical_mesh<ELT_t, STYPE, ngo_t>::get_higher_level_neighbours
(E_Int PHi, E_Int PGi, E_Int* neighbours, E_Int& nb_neighbours) const
{
  // Warning : assume same orientation for all the descendant of a given PG

  const E_Int* children = _PGtree.children(PGi);
  E_Int nb_children = _PGtree.nb_children(PGi);
  
  for (int i = 0; i < nb_children; ++i)
  {
    E_Int PH = (_F2E(0,PGi) == PHi) ? _F2E(1,children[i]) : _F2E(0,children[i]);
    neighbours[nb_neighbours++] = PH;
  }
}

///
template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t>
void hierarchical_mesh<ELT_t, STYPE, ngo_t>::get_enabled_neighbours
(E_Int PHi, E_Int* neighbours, E_Int& nb_neighbours) const
{
  // fill in up to 24 enabled neighbours and gives the number of enabled neighbours
  const E_Int* p = _ng->PHs.get_facets_ptr(PHi);
  E_Int nb_faces = _ng->PHs.stride(PHi);
    
#ifdef DEBUG_HIERARCHICAL_MESH
  assert (nb_neighbours == 0);
#endif
  for (int i = 0; i < nb_faces; ++i)
  {
    E_Int PGi = p[i] - 1;
      
    E_Int PH = NEIGHBOR(PHi, _F2E, PGi);

    if (PH == E_IDX_NONE) // returns only the enabled neighbours
      continue;
        
    if ( (_PHtree.is_enabled(PH)) || ((_PHtree.get_level(PH) > 0) && (_PHtree.is_enabled(_PHtree.parent(PH)))) )
      neighbours[nb_neighbours++] = PH; // no children : add the PH
    else
      get_higher_level_neighbours(PHi, PGi, neighbours, nb_neighbours); // add the 4 PH of higher level
  }
}

///
template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t>
void hierarchical_mesh<ELT_t, STYPE, ngo_t>::smooth(sensor_output_t& adap_incr) const
{
  std::stack<E_Int> stck;

  for (E_Int i = 0; i < _ng->PHs.size();  i++){
    if (adap_incr[i] != 0){
      stck.push(i);
    }
  }
  
  while (!stck.empty()){

    E_Int ind_PHi = stck.top(); // index of ith PH
    stck.pop();
    E_Int s = _ng->PHs.stride(ind_PHi);

    
    E_Int* neighbours= new E_Int[4*s];//fixme
    E_Int nb_neighbours = 0;
    get_enabled_neighbours(ind_PHi, neighbours, nb_neighbours); // get the effective neighbours (enabled ones) to assert the 2:1 rule

    E_Int incr = adap_incr[ind_PHi] + _PHtree.get_level(ind_PHi);

    for (int i = 0; i < nb_neighbours; ++i)
    {
      E_Int incr_neigh = adap_incr[neighbours[i]] + _PHtree.get_level(neighbours[i]);

      if (abs(incr-incr_neigh) <= 1) // 2:1 rule respected
        continue;

      // not respected : increase by 1 the adap incr of the lowest incr (ind_PHi or neighbours[i])
      E_Int PH_to_mod = (incr > incr_neigh) ? neighbours[i] : ind_PHi;

      if (adap_incr[PH_to_mod] >=0) // increase by 1 to refine : +=1 the adap incr and add the cell in the stack
      {
        adap_incr[PH_to_mod] += 1;
        stck.push(PH_to_mod);
      }
      else // increase by 1 an agglomeration : +=1 the adap incr for the cell and its siblings, and add all of them in the stack
      {
        E_Int father = _PHtree.parent(PH_to_mod);
        const E_Int* p = _PHtree.children(father);
        E_Int nb_children = _PHtree.nb_children(father);
        for (int j = 0; j < nb_children; ++j)
        {
          adap_incr[p[j]] += 1;
          stck.push(p[j]);
        }
      }
    }
    delete[] neighbours;
  }

}

///
template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t>
bool hierarchical_mesh<ELT_t, STYPE, ngo_t>::is_initialised() const
{
  if (_initialized) return true;
  else return false;
}


///
template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t>
void hierarchical_mesh<ELT_t, STYPE, ngo_t>::get_cell_center(E_Int PHi, E_Float* center) const
{    
  ELT_t::iso_barycenter(*_crd, _ng->PGs, _ng->PHs.get_facets_ptr(PHi), _ng->PHs.stride(PHi), 1, center);
}


}

#endif

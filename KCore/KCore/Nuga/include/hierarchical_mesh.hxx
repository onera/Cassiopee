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
#include "Nuga/include/q9.hxx"
#include "Nuga/include/h27.hxx"
#include "Connect/IdTool.h"
#include "Nuga/Delaunay/Triangulator.h"
#include "MeshElement/Basic.h"
#include "MeshElement/Pyramid.h"
#include "MeshElement/Prism.h"
#include "Nuga/include/refiner.hxx"



#define NEIGHBOR(PHi, _F2E, PGi) ( (_F2E(0,PGi) == PHi) ? _F2E(1,PGi) : _F2E(0,PGi) )

namespace NUGA
{
  
// general case : varying number of children, varying type of children
// HYBRID (ISO OR ANISO), 
template <typename ELT_t, eSUBDIV_TYPE STYPE> 
class hmesh_trait
{
  public:
    using arr_type = ngon_unit;
    
    static const E_Int PHNBC = -1;
    static const E_Int PGNBC = -1;
};

// isotropic HEXA subdivision => 8 HEXA children => fixed stride array
template <> 
class hmesh_trait<K_MESH::Hexahedron, ISO>
{
  public:
    using arr_type = K_FLD::IntArray;
    
    static const E_Int PHNBC = 8;
    static const E_Int PGNBC = 4;
};

// isotropic TETRA subdivision => 8 TETRA children => fixed stride array
template <> 
class hmesh_trait<K_MESH::Tetrahedron, ISO>
{
  public:
    using arr_type = K_FLD::IntArray;
    
    static const E_Int PHNBC = 8;
    static const E_Int PGNBC = 4;
};

// isotropic PRISM subdivision => 8 PRISM children => fixed stride array
template <> 
class hmesh_trait<K_MESH::Prism, ISO>
{
  public:
    using arr_type = K_FLD::IntArray;
    
    static const E_Int PHNBC = 8;
    static const E_Int PGNBC = 4;
  
};

///
template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t = ngon_type, typename crd_t = K_FLD::FloatArray>
class hierarchical_mesh
{
  public:

    using elt_type = ELT_t;
    using htrait = hmesh_trait<ELT_t, STYPE>;
    using arr_t = typename htrait::arr_type; //InttArray or ngon_unit
    using self_type = hierarchical_mesh<ELT_t, STYPE, ngo_t, crd_t>;
    using output_type = typename adap_incr_type<STYPE>::output_type;
    static constexpr eSUBDIV_TYPE sub_type = STYPE;
    
    crd_t*                    _crd;             // Coordinates
    ngo_t*                    _ng;              // NGON mesh
    tree<arr_t>               _PGtree, _PHtree; // Polygons/Polyhedra hierarchy
    K_FLD::IntArray           _F2E;             // neighborhood data structure : 2 X (nb_pgs) array : for each PG gives the left and right cells ids that share this PG
    bool                      _initialized;     // flag to avoid initialization more than once.

    ///
    hierarchical_mesh(crd_t& crd, ngo_t & ng):_crd(&crd), _ng(&ng), _PGtree(ng.PGs, htrait::PGNBC), _PHtree(ng.PHs, htrait::PHNBC), _initialized(false)
    {
      init();
    }

    ///
    E_Int init();
    ///
    void relocate (crd_t& crd, ngo_t & ng) {
      _crd = &crd;
      _ng = &ng;
      _PGtree.set_entities(ng.PGs);
      _PHtree.set_entities(ng.PHs);
    }
    ///
    E_Int adapt(output_type& adap_incr, bool do_agglo);
  
    /// face-conformity
    void conformize(Vector_t<E_Int>& pgoids);
    /// Keep only enabled PHs
    void filter_ngon(ngon_type& filtered_ng);
    
    ///
    void get_cell_center(E_Int PHi, E_Float* center);
    ///
    void get_enabled_neighbours(E_Int PHi, E_Int* neighbours, E_Int& nb_neighbours);
    ///
    void get_higher_level_neighbours(E_Int PHi, E_Int PGi, E_Int* neighbours, E_Int& nb_neighbours);
    ///
    bool enabled_neighbours(E_Int PHi); // return true if the PHi-th PH only has enabled neighbours
    ///
    void smooth(output_type& adap_incr);
    ///
    bool is_initialised();
  
#ifdef DEBUG_2019    
    ///
    static bool is_hybrid(const ngon_type& ng);
    ///
    void quality_measure();
    ///
    void check_vol(E_Float Vol_init, bool is_conformize);
    ///
    E_Float vol_init();
    ///
    E_Int control_children(E_Float tol);
    ///
    void control_tree(E_Int k);
    ///
    void vol_enf(E_Int k, E_Float &v);
    ///
    void verif3();
    void verif4();
    void verif5();
    static void debug_ph_actif(E_Int j, const ngon_type& ng, tree<arr_t> PHtree, ngon_unit& ph);
    static void debug_draw(E_Int j, const ngon_type& ng, tree<arr_t> PHtree, const K_FLD::FloatArray& crd);
#endif

  private:
    std::map<K_MESH::NO_Edge,E_Int> _ecenter;
    
    ///
    void __init(ngo_t& ng, const K_FLD::IntArray& F2E);

    ///
    void __conformize_next_lvl(Vector_t<E_Int>& molec, E_Int PGi, E_Int i);

};

// default implementation for any basci element in ISO mode
template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t, typename crd_t>
void hierarchical_mesh<ELT_t, STYPE, ngo_t, crd_t>::__init(ngo_t& ng, const K_FLD::IntArray& F2E)
{
  // sort the PGs of the initial NGON
  for (int i = 0; i < (int)_ng->PHs.size(); ++i)
  {
    ELT_t::reorder_pgs(*_ng,_F2E,i);
  }
}

template <>
void hierarchical_mesh<K_MESH::Hexahedron, DIR, ngon_type, K_FLD::FloatArray>::__init(ngo_t& ng, const K_FLD::IntArray& F2E)
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
template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t, typename crd_t>
E_Int hierarchical_mesh<ELT_t, STYPE, ngo_t, crd_t>::init()
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

template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t, typename crd_t>
E_Int hierarchical_mesh<ELT_t, STYPE, ngo_t, crd_t>::adapt(output_type& adap_incr, bool do_agglo)
{
  E_Int err(0);

  Vector_t<E_Int> PHagglo, PHref, PGref;

  // identify the sensor
  E_Int max = *std::max_element(ALL(adap_incr));
  E_Int min;
  
  if (!do_agglo) min = 0;
  else min = *std::min_element(ALL(adap_incr));
  
  bool loop = true;
  
  if ( (abs(min) <= 1) && (abs(max) <= 1) )
    loop = false;// geom sensor has only one iteration in the while => break at the end if false
  
  while (!err)
  {
    PHagglo.clear();
    if (min < 0) // there is a need to agglomerate
    {
      for (size_t i = 0; i < _ng->PHs.size(); ++i)
        if (adap_incr[i] == -1) PHagglo.push_back(i);
    }

    if (max > 0) // there is a need to refine
    {
      // refine PGs : create missing children (those PGi with _PGtree.nb_children(PGi) == 0)
      refiner<ELT_t, STYPE>::refine_PGs(adap_incr, *_ng, _PGtree, *_crd, _F2E, _ecenter);
      // refine PHs with missing children (those PHi with _PHtree.nb_children(PHi) == 0)
      refiner<ELT_t, STYPE>::refine_PHs(adap_incr, *_ng, _PGtree, _PHtree, *_crd, _F2E);
    }
    
    if (max == 0 && min == 0) break; // adapted
    
#ifdef DEBUG_2019    
//    NGDBG::draw_PH("s1.plt",_crd, _ng, 1);
//    NGDBG::draw_PH("s2.plt",_crd, _ng, 2);
//    NGDBG::draw_PH("s3.plt",_crd, _ng, 3);
//    NGDBG::draw_PH("s4.plt",_crd, _ng, 4);
//    NGDBG::draw_PH("s5.plt",_crd, _ng, 5);
//    NGDBG::draw_PH("s6.plt",_crd, _ng, 6);
//    NGDBG::draw_PH("s7.plt",_crd, _ng, 7);
//    NGDBG::draw_PH("s8.plt",_crd, _ng, 8);
#endif 
    
    // enable the right PHs & their levels, disable subdivided PHs
    adap_incr.resize(_ng->PHs.size(),0);
    for (E_Int i = 0; i < _ng->PHs.size(); ++i)
    {
      E_Int& PHi = i;

      if (!_PHtree.is_enabled(PHi)) continue;

      if (adap_incr[PHi] > 0) // activate the children
      {
        E_Int lvl_p1 = _PHtree.get_level(PHi)+ 1;
        E_Int nb_child = _PHtree.nb_children(PHi);
        const E_Int* children = _PHtree.children(PHi);
        for (E_Int j = 0; j < nb_child; ++j)
        {
          _PHtree.enable(*(children+j));
          //adap_incr[*(children+j)] = adap_incr[PHi] - 1;
          _PHtree.set_level(*(children+j), lvl_p1);
        }
      }
      adap_incr[PHi] = 0;
    }
    // agglomerate
    for (size_t i = 0; i < PHagglo.size(); ++i)
    {
      E_Int father = _PHtree.parent(PHagglo[i]);
      _PHtree.agglomerate(father);
    }
    // reset the agglomerated sons' adap_incr
    for (size_t i = 0; i < PHagglo.size(); ++i)
    {
      E_Int father = _PHtree.parent(PHagglo[i]);
      E_Int* p = _PHtree.children(father);
      E_Int nb_children = _PHtree.nb_children(father);
      for (int j = 0; j < nb_children; ++j)
        adap_incr[p[j]] = 0;
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
    filter_ngon(filtered_ng);

    std::ostringstream o;
    o << "NGON_it_" << iter << ".plt"; // we create a file at each iteration
    K_FLD::IntArray cnto;
    filtered_ng.export_to_array(cnto);
    MIO::write(o.str().c_str(), *_crd, cnto, "NGON");
    ++iter;
#endif

    if (!loop) break;
  }

  return err;
}

template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t, typename crd_t>
void hierarchical_mesh<ELT_t, STYPE, ngo_t, crd_t>::conformize(Vector_t<E_Int>& pgoids)
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
    E_Int* pPGi = _ng->PHs.get_facets_ptr(i);

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

template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t, typename crd_t>
void hierarchical_mesh<ELT_t, STYPE, ngo_t, crd_t>::__conformize_next_lvl(Vector_t<E_Int> &molec, E_Int PGi, E_Int i)
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
template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t, typename crd_t>
void hierarchical_mesh<ELT_t, STYPE, ngo_t, crd_t>::filter_ngon(ngon_type& filtered_ng)
{
  filtered_ng.PGs = _ng->PGs;
        
  //E_Int sz = _PGtree.get_parent_size();
  E_Int sz = _ng->PHs.size();
        
  for (int i = 0; i < sz; i++)
  {
    if (_PHtree.is_enabled(i) == true)
    {
      E_Int* p = _ng->PHs.get_facets_ptr(i);
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
template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t, typename crd_t>
void hierarchical_mesh<ELT_t, STYPE, ngo_t, crd_t>::get_higher_level_neighbours(E_Int PHi, E_Int PGi, E_Int* neighbours, E_Int& nb_neighbours)
{
  // Warning : assume same orientation for all the descendant of a given PG

  E_Int* children = _PGtree.children(PGi);
  E_Int nb_children = _PGtree.nb_children(PGi);
  
  for (int i = 0; i < nb_children; ++i)
  {
    E_Int PH = (_F2E(0,PGi) == PHi) ? _F2E(1,children[i]) : _F2E(0,children[i]);
    neighbours[nb_neighbours++] = PH;
  }
}

///
template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t, typename crd_t>
void hierarchical_mesh<ELT_t, STYPE, ngo_t, crd_t>::get_enabled_neighbours(E_Int PHi, E_Int* neighbours, E_Int& nb_neighbours)
{
  // fill in up to 24 enabled neighbours and gives the number of enabled neighbours
  E_Int* p = _ng->PHs.get_facets_ptr(PHi);
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
/*template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t, typename crd_t>
bool hierarchical_mesh<K_MESH::Hexahedron, ngon_type, K_FLD::FloatArray>::enabled_neighbours(E_Int PHi)
{
  // if at least 1 neighbour is disabled, it returns false, otherwise it returns true
  bool enabled = true;
  E_Int* p = _ng->PHs.get_facets_ptr(PHi);
  E_Int nb_faces = _ng->PHs.stride(PHi);
    
  for (int i = 0; i < nb_faces; ++i)
  {
    E_Int PGi = p[i] - 1;
    E_Int PH = _F2E(0,PGi);

    if (PH == E_IDX_NONE) continue;

    if (PH != PHi && !_PHtree.is_enabled(PH))
      enabled = false;
    else
    {
      PH = _F2E(1,PGi);
            
      if (PH == E_IDX_NONE) continue;
          
      if (!_PHtree.is_enabled(PH))
        enabled = false;            
    }

  }
    
  return enabled;
}*/

///
template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t, typename crd_t>
void hierarchical_mesh<ELT_t, STYPE, ngo_t, crd_t>::smooth(output_type& adap_incr)
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
        E_Int* p = _PHtree.children(father);
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
template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t, typename crd_t>
bool hierarchical_mesh<ELT_t, STYPE, ngo_t, crd_t>::is_initialised()
{
  if (_initialized) return true;
  else return false;
}


///
template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t, typename crd_t>
void hierarchical_mesh<ELT_t, STYPE, ngo_t, crd_t>::get_cell_center(E_Int PHi, E_Float* center)
{    

  ELT_t::iso_barycenter(*_crd, _ng->PGs, _ng->PHs.get_facets_ptr(PHi), _ng->PHs.stride(PHi), 1, center);
}

#ifdef DEBUG_2019
///
template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t, typename crd_t>
void hierarchical_mesh<ELT_t, STYPE, ngo_t, crd_t>::quality_measure()
{
  E_Float qmean(.0), qmin(1.0), qmax(.0);
  E_Int* face;
  E_Int nb_T(0);
  E_Float Vol;
  E_Float VolT(.0);
    
  for (int i=0; i<_ng->PHs.size(); i++)
  {
    if (_PHtree.is_enabled(i) && _ng->PHs.stride(i)==4)
    {
      ++nb_T;
      face= _ng->PHs.get_facets_ptr(i);

      ELT_t e(_ng->PGs,face);
      
      E_Float q = e.quality(_crd, &Vol);

      qmean += q;
      qmin=std::min(qmin,q);
      qmax=std::max(qmax,q);
      
      VolT += Vol;
    }
  }
  
  qmean = qmean/nb_T;
   

  std::cout << "Qualité_m= " << qmean << std::endl;
  std::cout << "nb_T= " << nb_T << std::endl;
  std::cout << "Qualité_min= " << qmin << "  Qualité_max= " << qmax << std::endl; 
  std::cout.precision(12);
  std::cout << "Volume Total= " << VolT << std::endl;
}

template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t, typename crd_t>
bool hierarchical_mesh<ELT_t, STYPE, ngo_t, crd_t>::is_hybrid(const ngon_type& ng){
  E_Int s1(0), s2(0), s3(0), s4(0);  
  E_Int err = 0;
  for (E_Int i = 0; (i < ng.PHs.size()) && !err; ++i){
        if (K_MESH::Polyhedron<0>::is_HX8(ng.PGs, ng.PHs.get_facets_ptr(i), ng.PHs.stride(i)) ) ++s1;
        else if (K_MESH::Polyhedron<0>::is_TH4(ng.PGs, ng.PHs.get_facets_ptr(i), ng.PHs.stride(i)) ) ++s2;
        else if (K_MESH::Polyhedron<0>::is_PY5(ng.PGs, ng.PHs.get_facets_ptr(i), ng.PHs.stride(i)) ) ++s3;
        else if (K_MESH::Polyhedron<0>::is_PR6(ng.PGs, ng.PHs.get_facets_ptr(i), ng.PHs.stride(i)) ) ++s4;
        else{
            std::cout << "i= " << i << " ng.PHs.stride(i)= " << ng.PHs.stride(i) << std::endl;
            return false;
        }

    }
    return true;
}


template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t, typename crd_t>
void hierarchical_mesh<ELT_t, STYPE, ngo_t, crd_t>::check_vol(E_Float Vol_init, bool is_conformize){
  std::vector<E_Float> vols;
  E_Int a=ngon_t<K_FLD::IntArray>::volumes< DELAUNAY::Triangulator >(_crd, _ng, vols);
  //std::vector<E_Float>::iterator min = std::min_element(std::begin(vols), std::end(vols));
  a++;
  E_Float Vol(0);
  std::cout << "vols size= " << vols.size() << std::endl;
  std::cout << "_ng->PHs size= " << _ng->PHs.size() << std::endl;
  if (is_conformize)
  for (int i=0; i<vols.size(); i++){
//      if (_PHtree.is_enabled(i)){
          Vol += vols[i];
//      }      
  }
  else{
    for (int i=0; i<vols.size(); i++){
      if (_PHtree.is_enabled(i)){
          Vol += vols[i];
      }      
    }
  }
  std::cout << "Vol= "; 
  std::cout.precision(12);
  std::cout << Vol << std::endl;
  E_Float eps = ::fabs(Vol_init-Vol)/Vol_init;
  
  std::cout << "eps= ";
  std::cout.precision(12);
  std::cout << eps << std::endl;
  //std::cout << "min= " << *min << std::endl;
}

template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t, typename crd_t>
E_Float hierarchical_mesh<ELT_t, STYPE, ngo_t, crd_t>::vol_init(){
  std::vector<E_Float> vols;
  E_Int a=ngon_t<K_FLD::IntArray>::volumes< DELAUNAY::Triangulator >(_crd, _ng, vols);
  a++;
  E_Float Vol(0);
  std::cout << "vols size= " << vols.size() << std::endl;
  std::cout << "_ng->PHs size= " << _ng->PHs.size() << std::endl;
  for (int i=0; i<vols.size(); i++){
          Vol += vols[i];      
  }
  return Vol;
}

template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t, typename crd_t>
E_Int hierarchical_mesh<ELT_t, STYPE, ngo_t, crd_t>::control_children(E_Float tol){
    std::vector<E_Float> vols;
    E_Int a=ngon_t<K_FLD::IntArray>::volumes< DELAUNAY::Triangulator >(_crd, _ng, vols);
    a++;
    E_Int err(0);
    E_Float res(0);
    E_Int ind;
    E_Int s1(0), s2(0), s3(0);
    E_Int nb_phs = _ng->PHs.size();
    for (int i=0; i<nb_phs; i++){
        if (_PHtree.nb_children(i)!=0){
            E_Float vol_children(0);
            E_Float vol_father= vols[i];
            for (int j=0; j<_PHtree.nb_children(i); j++){
                E_Int* child= _PHtree.children(i);
                E_Int child_j= *(child+j);
                E_Float vol_child=vols[child_j];
                vol_children += vol_child;
            }
            if (::fabs(vol_father-vol_children)>tol){
                err += 1;
                if (K_MESH::Polyhedron<0>::is_HX8(_ng->PGs, _ng->PHs.get_facets_ptr(i), _ng->PHs.stride(i)) ) ++s1;
                else if (K_MESH::Polyhedron<0>::is_TH4(_ng->PGs, _ng->PHs.get_facets_ptr(i), _ng->PHs.stride(i)) ) ++s2;
                else if (K_MESH::Polyhedron<0>::is_PY5(_ng->PGs, _ng->PHs.get_facets_ptr(i), _ng->PHs.stride(i)) ) ++s3;
                if (::fabs(vol_father-vol_children)> res){
                    res=::fabs(vol_father-vol_children);
                    ind=i;
                }
                    
            }
        }
    }
    std::cout << "err tol= " << err <<std::endl;
    std::cout << "///////" << std::endl;
    std::cout << "s1= " << s1 << std::endl;
    std::cout << "s2= " << s2 << std::endl;
    std::cout << "s3= " << s3 << std::endl;
    std::cout << "///////" << std::endl;
    return ind;
}

template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t, typename crd_t>
void hierarchical_mesh<ELT_t, STYPE, ngo_t, crd_t>::control_tree(E_Int k){
    E_Float vol_father;      
    DELAUNAY::Triangulator dt;
    E_Float Gdum[3];
    //E_Float vol(0);
    E_Int err = K_MESH::Polyhedron<UNKNOWN>::metrics<DELAUNAY::Triangulator>(dt, _crd, _ng->PGs, _ng->PHs.get_facets_ptr(k), _ng->PHs.stride(k), vol_father, Gdum);
    err++;
    std::cout << "vol= " << vol_father << std::endl;
    std::cout << "k= " << k << std::endl;
    std::cout << "nb_children= " << _PHtree.nb_children(k) << std::endl;
    E_Float v(0);
    for (int j=0; j< _PHtree.nb_children(k); j++){
        E_Int* child= _PHtree.children(k);
        E_Int child_j= *(child+j);
        if (_PHtree.nb_children(child_j)!=0){
            vol_enf(child_j, v);
        }
        else{
        
        E_Float vol_child;
        E_Int err = K_MESH::Polyhedron<UNKNOWN>::metrics<DELAUNAY::Triangulator>(dt, _crd, _ng->PGs, _ng->PHs.get_facets_ptr(child_j), _ng->PHs.stride(child_j), vol_child, Gdum);
        err++;
        v += vol_child;
        }
    }
    std::cout << "vol_child= " << v << std::endl;

}


template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t, typename crd_t>
void hierarchical_mesh<ELT_t, STYPE, ngo_t, crd_t>::vol_enf(E_Int k, E_Float &v){
    E_Float vol_child;
    DELAUNAY::Triangulator dt;
    E_Float Gdum[3];
    for (int j=0; j< _PHtree.nb_children(k); j++){
        E_Int* child= _PHtree.children(k);
        E_Int child_j= *(child+j);
        if (_PHtree.nb_children(child_j)!=0){
            vol_enf(child_j, v);
        }
        else{
            E_Int err = K_MESH::Polyhedron<UNKNOWN>::metrics<DELAUNAY::Triangulator>(dt, _crd, _ng->PGs, _ng->PHs.get_facets_ptr(child_j), _ng->PHs.stride(child_j), vol_child, Gdum);
            err++;
            v += vol_child;
        }
    }
}

///
template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t, typename crd_t>
void hierarchical_mesh<ELT_t, STYPE, ngo_t, crd_t>::verif3()
{
  E_Int err(0);  
  E_Int nb_elts = _ng->PHs.size();

  for (int i=0; i< nb_elts; i++){
      E_Int level= _PHtree.get_level(i);
      E_Int nb_faces= _ng->PHs.stride(i);
      E_Int* pPGi = _ng->PHs.get_facets_ptr(i);
      if (K_MESH::Polyhedron<0>::is_PY5(_ng->PGs, _ng->PHs.get_facets_ptr(i), _ng->PHs.stride(i))){
           continue;
       }     
      for (int j=0; j< nb_faces; j++){
        E_Int PGj = *(pPGi +j) - 1;
        E_Int PHn = NEIGHBOR(i, _F2E, PGj);
        if (PHn==E_IDX_NONE || K_MESH::Polyhedron<0>::is_PY5(_ng->PGs, _ng->PHs.get_facets_ptr(PHn), _ng->PHs.stride(PHn))){
            continue;
        }
        E_Int leveln= _PHtree.get_level(PHn);
        if (::fabs(leveln-level)>1){
            err += 1;
        }
      }
  }
  std::cout << "err level= " << err << std::endl;
}

///
template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t, typename crd_t>
void hierarchical_mesh<ELT_t, STYPE, ngo_t, crd_t>::verif4()
{
  E_Int err(0);  
  E_Int nb_elts = _ng->PHs.size();

  for (int i=0; i< nb_elts; i++){
    //E_Int level= _hmesh._PHtree.get_level(i);
    //E_Int nb_faces= _hmesh._ng->PHs.stride(i);
    E_Int* pPGi = _ng->PHs.get_facets_ptr(i);
    if (K_MESH::Polyhedron<0>::is_PY5(_ng->PGs, _ng->PHs.get_facets_ptr(i), _ng->PHs.stride(i)))
      continue;
 
    if (K_MESH::Polyhedron<0>::is_HX8(_ng->PGs, _ng->PHs.get_facets_ptr(i), _ng->PHs.stride(i))){
      for (int j=0; j< 3; j++){
        E_Int PGj = *(pPGi +2*j) - 1;
        E_Int PHn = NEIGHBOR(i, _F2E, PGj);
        if (PHn==E_IDX_NONE || K_MESH::Polyhedron<0>::is_PY5(_ng->PGs, _ng->PHs.get_facets_ptr(PHn), _ng->PHs.stride(PHn)))
          continue;
        
        E_Int leveln= _PHtree.get_level(PHn);
        for (int k=0; k< 4; k++){
          E_Int PGk= *(pPGi +(j+k+2) % 6) - 1;
          E_Int PHk = NEIGHBOR(i, _F2E, PGk);
          if (PHk==E_IDX_NONE || K_MESH::Polyhedron<0>::is_PY5(_ng->PGs, _ng->PHs.get_facets_ptr(PHk), _ng->PHs.stride(PHk)))
            continue;      
 
          E_Int levelk= _PHtree.get_level(PHk);
          if (::fabs(leveln-levelk)>1){
            err += 1;
          }
        }
      }
                 
      for (int j=0; j< 3; j++){
        E_Int PGj = *(pPGi +2*j+1) - 1;
        E_Int PHn = NEIGHBOR(i, _F2E, PGj);
        if (PHn==E_IDX_NONE || K_MESH::Polyhedron<0>::is_PY5(_ng->PGs, _ng->PHs.get_facets_ptr(PHn), _ng->PHs.stride(PHn)))
          continue;
        
        E_Int leveln= _PHtree.get_level(PHn);
        for (int k=0; k< 4; k++){
          E_Int PGk= *(pPGi +(j+k+1) % 6) - 1;
          E_Int PHk = NEIGHBOR(i, _F2E, PGk);
          if (PHk==E_IDX_NONE || K_MESH::Polyhedron<0>::is_PY5(_ng->PGs, _ng->PHs.get_facets_ptr(PHk), _ng->PHs.stride(PHk)))
            continue;      
 
          E_Int levelk= _PHtree.get_level(PHk);
          if (::fabs(leveln-levelk)>1){
            err += 1;
          }
        }
      }
    }
    else if (K_MESH::Polyhedron<0>::is_TH4(_ng->PGs, _ng->PHs.get_facets_ptr(i), _ng->PHs.stride(i))){
      for (int j=0; j< 4; j++){
        E_Int PGj = *(pPGi +j) - 1;
        E_Int PHn = NEIGHBOR(i, _F2E, PGj);
        if (PHn==E_IDX_NONE || K_MESH::Polyhedron<0>::is_PY5(_ng->PGs, _ng->PHs.get_facets_ptr(PHn), _ng->PHs.stride(PHn)))
          continue;      
        
        E_Int leveln= _PHtree.get_level(PHn);
        for (int k=0; k< 3; k++){
          E_Int PGk= *(pPGi +(j+k+1) % 4) - 1;
          E_Int PHk = NEIGHBOR(i, _F2E, PGk);
          if (PHk==E_IDX_NONE || K_MESH::Polyhedron<0>::is_PY5(_ng->PGs, _ng->PHs.get_facets_ptr(PHk), _ng->PHs.stride(PHk)))
            continue;      
           
          E_Int levelk= _PHtree.get_level(PHk);
          if (::fabs(leveln-levelk)>1){
            err += 1;
          }
        }
      } 
    }
  }
          
//          if (::fabs(level-leveln)==1){

  //        E_Int leveln= _hmesh._PHtree.get_level(PHn);
  //        if (::fabs(leveln-level)>1){
  //            err += 1;
  //        }
  std::cout << "err level geom2= " << err << std::endl;
}

template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t, typename crd_t>
void hierarchical_mesh<ELT_t, STYPE, ngo_t, crd_t>::verif5()
{
    Vector_t<E_Int> lmax(_crd->cols(),0);
    Vector_t<E_Int> lmin(_crd->cols(),E_IDX_NONE);
    
    E_Int nb_phs= _ng->PHs.size();
    for (int i=0; i< nb_phs; i++){
//        E_Int lvlmin(E_IDX_NONE);
//        E_Int lvlmax(0);
        if (_PHtree.is_enabled(i)){
            E_Int level= _PHtree.get_level(i);
            E_Int nb_faces= _ng->PHs.stride(i);
            E_Int* face= _ng->PHs.get_facets_ptr(i);
            for (int j=0; j< nb_faces; j++){
                E_Int PGj= *(face+j)-1;
                E_Int nb_nodes= _ng->PGs.stride(PGj);
                E_Int* nodes= _ng->PGs.get_facets_ptr(PGj);
                for (int k=0; k< nb_nodes; k++){
                    E_Int pNk= *(nodes+k)-1;
                    lmax[pNk]= std::max(lmax[pNk], level);
                    lmin[pNk]= std::min(lmin[pNk], level);
                }
            }
        }
    }
    
    E_Int nb_nodes= _crd->cols();
    E_Int err(0);
    for (int i=0; i< nb_nodes; i++){
        if (lmax[i]-lmin[i]>1){
            err += 1;
        }
    }
    std::cout << "err nodal level= " << err << std::endl; 
}


template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t, typename crd_t>
void hierarchical_mesh<ELT_t, STYPE, ngo_t, crd_t>::debug_ph_actif(E_Int j, const ngon_type& ng, tree<arr_t> PHtree, ngon_unit& ph)
{
    if (PHtree.is_enabled(j)){
        ph.add(ng.PHs.stride(j), ng.PHs.get_facets_ptr(j));
        ph._facet[0] += 1;
    }
    else if (PHtree.nb_children(j) != 0){
        E_Int *children = PHtree.children(j);
        for (int i=0; i< PHtree.nb_children(j); i++){
            debug_ph_actif(*(children+i), ng, PHtree, ph);
        }
    }
}

// fixme : debug par encore au point
template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t, typename crd_t>
void hierarchical_mesh<ELT_t, STYPE, ngo_t, crd_t>::debug_draw(E_Int j, const ngon_type& ng, tree<arr_t> PHtree, const K_FLD::FloatArray& crd)
{
    ngon_unit ph;
    //ph._NGON[0]=1;
    if (PHtree.is_enabled(j)){
        ph.add(ng.PHs.stride(j), ng.PHs.get_facets_ptr(j));
    }
    else if (PHtree.nb_children(j) != 0){
        E_Int *children = PHtree.children(j);
        for (int i=0; i< PHtree.nb_children(j); i++){
            debug_ph_actif(*(children+i), ng, PHtree, ph);
        }
    }
    //ph._NGON[1] = ph._NGON.size() - 2;
    ph.updateFacets();
//    std::cout << "ph size= " << ph.size() << std::endl;
//    E_Int* PGi = ph.get_facets_ptr(0);
//    for (int i=0; i< ph.stride(0); i++){
//         std::cout << "PG[" << i << "]= " << *(PGi+i)-1 << std::endl;
//    }   
    ngon_type one_ph(ng.PGs, ph);
 
    NGDBG::draw_PH("draw_out.plt", crd, one_ph, 0);
  //static void draw_PH(const char* fname, const K_FLD::FloatArray& coord, const ngon_type& ng, E_Int i);
}


#endif

}

#endif

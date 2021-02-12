/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr), Alexis Gay (alexis.gay@onera.fr), Alexis Rouil (alexis.rouil@onera.fr)

#ifndef NUGA_HIERACHICAL_MESH_HXX
#define NUGA_HIERACHICAL_MESH_HXX

#if defined (DEBUG_HIERARCHICAL_MESH) || defined (OUTPUT_ITER_MESH)
#if defined (VISUAL) || defined(NETBEANSZ)
#include "Nuga/include/medit.hxx"
#else
#include "IO/io.h"
#include "Nuga/include/NGON_debug.h"
using NGDBG = NGON_debug<K_FLD::FloatArray, K_FLD::IntArray>;
#endif
#endif

#include "Nuga/include/tree.hxx"
#include "Nuga/include/subdivision.hxx"

#include "Nuga/include/IdTool.h"
#include "Nuga/include/Triangulator.h"
#include "Nuga/include/Basic.h"
#include "Nuga/include/Prism.h"
#include "Nuga/include/refiner.hxx"
#include "Nuga/include/macros.h"

#include "Nuga/include/join_sensor.hxx"
#include "Nuga/include/communicator.hxx"
#include "Nuga/include/join_t.hxx"


namespace NUGA
{

using crd_t = K_FLD::FloatArray;

///
template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t = ngon_type>
class hierarchical_mesh
{
  public:

    using elt_t = ELT_t;
    static constexpr  eSUBDIV_TYPE SUBTYPE = STYPE;

    using self_t = hierarchical_mesh<ELT_t, STYPE, ngo_t>;
    using subdiv_t = subdiv_pol<ELT_t, STYPE>;
    using pg_arr_t = typename subdiv_t::pg_arr_t;
    using ph_arr_t = typename subdiv_t::ph_arr_t;
    using output_t = typename sensor_output_data<STYPE>::type;
    using pg_tree_t = tree<pg_arr_t>; 
    using ph_tree_t = tree<ph_arr_t>;

    using bc_data_t = std::vector<std::vector<E_Int>>;

    // mutli-zone stuff
    using jsensor_t = join_sensor<self_t>;
    using join_data_t = std::vector<std::pair<E_Int, std::vector<E_Int>>>;
    using jcom_t = jsensor_com_agent<self_t, typename jsensor_t::input_t>;
    using communicator_t = NUGA::communicator<jcom_t>;

    crd_t                     _crd;             // Coordinates
    ngo_t                     _ng;              // NGON mesh
    pg_tree_t                 _PGtree;          // Polygons hierarchy
    ph_tree_t                 _PHtree;          // Polyhedra hierarchy
    K_FLD::IntArray           _F2E;             // neighborhood data structure : 2 X (nb_pgs) array : for each PG gives the left and right cells ids that share this PG
    bool                      _initialized;     // flag to avoid initialization more than once.
    refiner<ELT_t, STYPE>     _refiner;         // refinement methods must stay static, so this object is here to store _ecenter (not relevant in hmesh)

    E_Int _idx_start;

    // BCs
    bc_data_t BCptLists;

    // JOINS
    E_Int            zid;
    join_t<self_t>* join;
    jsensor_t*       jsensor;
    communicator_t*  COM;

    //hack
    mutable std::vector<E_Int> _pgnids, _phnids;

    // for fields projetion
    E_Int                     _nb_phs0;         // intial nb of PHs
    Vector_t<bool>            _enabledPHi;      // current PH enabling status

    ///
    hierarchical_mesh(crd_t& crd, ngo_t & ng):_crd(crd), _ng(ng), _PGtree(ng.PGs), _PHtree(ng.PHs), _initialized(false), zid(0), join(nullptr), jsensor(nullptr), COM(nullptr), _idx_start(0) { init(); }
    ///
    hierarchical_mesh(crd_t& crd, K_FLD::IntArray& cnt, E_Int idx_start, bc_data_t& bcptlists) :_crd(crd), _ng(cnt), _PGtree(_ng.PGs), _PHtree(_ng.PHs), _initialized(false), _idx_start(idx_start), BCptLists(bcptlists), zid(0), join(nullptr), jsensor(nullptr), COM(nullptr) { init(); }

    //multi-zone constructor
    hierarchical_mesh(E_Int id, K_FLD::FloatArray& crd, K_FLD::IntArray& cnt, const bc_data_t& bcptlists, const join_data_t& jdata, E_Int idx_start, communicator_t* com);

    ~hierarchical_mesh()
    {
      if (join != nullptr) delete join;
      if (jsensor != nullptr) delete jsensor;
      if (COM != nullptr)
        if (COM->agents[zid] != nullptr) delete COM->agents[zid];
    }

    ///
    E_Int init();
    ///
    /*E_Int relocate (crd_t& crd, ngo_t & ng) {
      _crd = &crd;
      _ng = &ng;
      _PGtree.set_entities(ng.PGs);
      _PHtree.set_entities(ng.PHs);

      if (ng.PGs.size() != _PGtree.size())
        return 1; //trying to relocate on the wrong mesh
      if (ng.PHs.size() != _PHtree.size())
        return 1; //trying to relocate on the wrong mesh
      
      return 0;
    }*/

    ///
    E_Int adapt(output_t& adap_incr, bool do_agglo);
  
    /// face-conformity
    void conformize(ngo_t& ngo, Vector_t<E_Int>& pgoids) const;
    /// Keep only enabled PHs
    void extract_enabled_phs(ngon_type& filtered_ng) const ;
    ///
    void extract_enabled_pgs_descendance(E_Int PGi, bool reverse, std::vector<E_Int>& pointlist);
    ///
    void __extract_enabled_pgs_descendance(E_Int PGi, NUGA::reordering_func F, bool reverse, std::vector<E_Int>& pointlist);
    ///
    void extract_plan(E_Int PGi, bool reverse, E_Int i0, pg_arr_t& plan) const;
    ///
    void get_cell_center(E_Int PHi, E_Float* center) const ;
    ///
    template <typename InputIterator> void get_enabled_neighbours(E_Int PHi, InputIterator neighbours, E_Int& nb_neighbours) const ;
    ///
    template <typename InputIterator> void get_higher_level_neighbours(E_Int PHi, E_Int PGi, InputIterator neighbours, E_Int& nb_neighbours) const ;
    ///
    void enable_PGs();

    void update_BCs();
    
    ///
    bool is_initialised() const ;

    ///
    E_Int project_cell_center_sol_order1(std::vector<std::vector<E_Float>>& fields);

private:
    ///
    void __init();
    ///
    template <typename PG_t>
    void __append_children_plan(E_Int PGi, bool reverse, int i0, std::map<E_Int, std::vector<E_Int>>& lvl_to_plan, E_Int lvl);
 
};

///
template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t>
hierarchical_mesh<ELT_t, STYPE, ngo_t>::hierarchical_mesh
(E_Int id, K_FLD::FloatArray& crd, K_FLD::IntArray& cnt, const bc_data_t& bcptlists, const join_data_t& jdata, E_Int idx_start, communicator_t* com) :
  _crd(crd), _ng(cnt), _PGtree(_ng.PGs), _PHtree(_ng.PHs), _initialized(false), _idx_start(idx_start), BCptLists(bcptlists), zid(id), join(nullptr), jsensor(nullptr), COM(com)
{
  //std::cout << "hierarchical_mesh : begin " << std::endl;
  init();
  //std::cout << "hierarchical_mesh : jdata sz" << jdata.size() << std::endl;
  join_data_t jmp = jdata;
  if (!jmp.empty()) // join is specified
  {
    //std::cout << "hierarchical_mesh : join specified " << std::endl;
    jsensor = new jsensor_t(*this);
    join = new join_t<self_t>(id, *this, idx_start);

    assert((E_Int)COM->agents.size() > zid); //COM agents must be resized by the caller before this call

    COM->agents[zid] = new jcom_t(zid, join, jsensor);
    //std::cout << "hierarchical_mesh : jmp sz " << jmp.size() << std::endl;
    for (auto j : jmp)
    {
      E_Int joinedZid = j.first;
      std::vector<E_Int>& ptlist = j.second;

      join->link(joinedZid, ptlist);
      COM->agents[zid]->link(joinedZid);
    }
  }
}

// default implementation for any basci element in ISO mode
template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t>
void hierarchical_mesh<ELT_t, STYPE, ngo_t>::__init()
{
  // sort the PGs of the initial NGON
  E_Int nb_phs = _ng.PHs.size();
  for (E_Int i = 0; i < nb_phs; ++i)
    ELT_t::reorder_pgs(_ng,_F2E,i);
}

template <> inline
void hierarchical_mesh<K_MESH::Polyhedron<0>, ISO_HEX, ngon_type>::__init()
{
  // nothing to do
}


/*template <> inline
void hierarchical_mesh<K_MESH::Hexahedron, DIR, ngon_type>::__init()
{
  // alexis : todo
  
  E_Int nb_phs = _ng.PHs.size();
  
  //type init
  _ng.PHs._type.clear();
  _ng.PHs._type.resize(nb_phs, (E_Int)K_MESH::Polyhedron<0>::eType::HEXA);
  
  // first reorder : by opposite pair
  E_Int generators[2], HX6opposites[6], *pg(generators), *qopp(HX6opposites);
  for (E_Int i=0; i < nb_phs; ++i)
  {
    K_MESH::Polyhedron<0>::is_prismN(_ng.PGs, _ng.PHs.get_facets_ptr(i), _ng.PHs.stride(i), pg, qopp);
    // alexis : utiliser pg et qopp pour réordonner
  }
  
  // promote eventually to layer type + 2nd reorder
  for (E_Int i=0; i < _F2E.cols(); ++i)
  {
    // alexis : si i n'est pas une frontière => continue
    
    E_Int PHi = (_F2E(0,i) != IDX_NONE) ? _F2E(0,i) : _F2E(1,i);
    
    E_Int PHcur = PHi;
//    E_Int Basecur = i;
    
    while (true) // climb over the layer
    {
      if (_ng.PHs._type[PHcur] == K_MESH::Polyhedron<0>::eType::LAYER) //already set
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
    K_MESH::Hexahedron::reorder_pgs(_ng,_F2E,i);
  
  //detect now any layer cell
//  double aniso_ratio = 0.2; //fixme : parametre a externaliser
  
  //to do : si aniso => faire  ng.PHs._type =  K_MESH::Polyhedron<0>::eType::LAYER
  
}*/

///
template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t>
E_Int hierarchical_mesh<ELT_t, STYPE, ngo_t>::init()
{
  if (_initialized) return 0;

  if (_ng.PGs.size() == 0) return 1;

  _nb_phs0 = _ng.PHs.size();
  
  E_Int err(0);
  
  // We reorient the PG of our NGON
  _ng.flag_externals(1);
  DELAUNAY::Triangulator dt;
  bool has_been_reversed;
  err = ngon_type::reorient_skins(dt, _crd, _ng, has_been_reversed);
  if (err)
    return 1;
  
  //F2E
  ngon_unit neighbors;
  _ng.build_ph_neighborhood(neighbors);
  _ng.build_F2E(neighbors, _F2E);
  
  __init(); // for pure type == reorder_pgs

  _initialized = true;
  
  return err;
}

#ifdef OUTPUT_ITER_MESH
 static E_Int iter = 1;
#endif

template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t>
E_Int hierarchical_mesh<ELT_t, STYPE, ngo_t>::adapt(output_t& adap_incr, bool do_agglo)
{
  E_Int fmax{1}, cmax{1}; //initialized wit 1 to do refinement at first iter

  // infinite loop but breaking test is done at each iteration (and terminates at some point)
  while (true)
  {
    E_Int nb_phs0 = _ng.PHs.size();
    // refine Faces : create missing children (those PGi with _PGtree.nb_children(PGi) == 0)
    if (fmax > 0 || cmax > 0) refiner<ELT_t, STYPE>::refine_Faces(adap_incr, _ng, _PGtree, _crd, _F2E, _refiner._ecenter);
    // refine Cells with missing children (those PHi with _PHtree.nb_children(PHi) == 0)
    if (cmax > 0) refiner<ELT_t, STYPE>::refine_PHs(adap_incr, _ng, _PGtree, _PHtree, _crd, _F2E);
        
    //std::cout << "update cell_adap_incr, enable the right PHs & their levels" << std::endl;
    adap_incr.cell_adap_incr.resize(_ng.PHs.size(),0);// resize to new size

    for (E_Int PHi = 0; PHi < nb_phs0; ++PHi)
    {
      E_Int& adincrPHi = adap_incr.cell_adap_incr[PHi];
      if (adincrPHi == 0) continue;

      if (!_PHtree.is_enabled(PHi))
      {
        adincrPHi = 0; // reset in cas it has something : do not allow splitting of disabled entities
        continue;
      }

      if (adincrPHi > 0) // refinement : activate the children, transfer the adapincr & set their level
      {
        E_Int nb_child = _PHtree.nb_children(PHi);
        const E_Int* children = _PHtree.children(PHi);
        for (E_Int j = 0; j < nb_child; ++j)
        {
          _PHtree.enable(*(children+j));
          adap_incr.cell_adap_incr[*(children+j)] = adincrPHi - 1;

          E_Int lvl_p1 = _PHtree.get_level(PHi) + 1;
          _PHtree.set_level(*(children+j), lvl_p1); //fixme : do it somewhere else ?
        }
        adincrPHi = 0;//reset incr
      }
      else // agglomeration : activate the father, transfer that adap incr
      {
        E_Int father = _PHtree.parent(PHi);
        _PHtree.enable(father);
        adap_incr.cell_adap_incr[father] = adincrPHi + 1;
        // reset incr on children
        E_Int nb_child = _PHtree.nb_children(PHi);
        const E_Int* children = _PHtree.children(PHi);
        for (E_Int j = 0; j < nb_child; ++j)
          adap_incr.cell_adap_incr[*(children + j)] = 0;
      }
    }

    //std::cout << "enable_PGs..." << std::endl;
    enable_PGs();
    
    _ng.PGs.updateFacets();
    _ng.PHs.updateFacets();

#ifdef DEBUG_HIERARCHICAL_MESH     
    if (! _ng.attributes_are_consistent()) return false;
#endif
    
#ifdef OUTPUT_ITER_MESH
    ngon_type filtered_ng;
    extract_enabled_phs(filtered_ng);

    std::ostringstream o;
 
#if defined (VISUAL) || defined(NETBEANSZ)
    o << "NGON_it_" << iter; // we create a file at each iteration
    medith::write(o.str().c_str(), _crd, filtered_ng);
#else
    o << "NGON_it_" << iter << ".plt"; // we create a file at each iteration
    K_FLD::IntArray cnto;
    filtered_ng.export_to_array(cnto);
    MIO::write(o.str().c_str(), _crd, cnto, "NGON");
#endif

    ++iter;
#endif

    // get the extrema values
    cmax = *std::max_element(ALL(adap_incr.cell_adap_incr));
    E_Int cmin = (!do_agglo) ? 0 : *std::min_element(ALL(adap_incr.cell_adap_incr));
    fmax = *std::max_element(ALL(adap_incr.face_adap_incr));
    //E_Int fmin = (!do_agglo) ? 0 : *std::min_element(ALL(adap_incr.face_adap_incr));

    //std::cout << "fmin/fmax/cmin/cmax : " << fmin << "/" << fmax << "/" << cmin << "/" << cmax << std::endl;

    if (cmax == 0 && cmin == 0 /*&& fmin == 0 && fmax == 0*/) return 0; // no adaptation required anymore
  }

  return 1;
}


///
template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t>
void hierarchical_mesh<ELT_t, STYPE, ngo_t>::conformize(ngo_t& ngo, Vector_t<E_Int>& pgoids) const
{
  ngon_unit new_phs;
  Vector_t<E_Int> molec, ids;

  E_Int nb_phs = _ng.PHs.size();
  for (E_Int i = 0; i < nb_phs; ++i)
  {
    if (!_PHtree.is_enabled(i)) continue;

    molec.clear();
    E_Int s = _ng.PHs.stride(i);
    const E_Int* pPGi = _ng.PHs.get_facets_ptr(i);

    for (E_Int j = 0; j < s; ++j)
    {
      E_Int PGi = *(pPGi + j) - 1;

      if (_PGtree.is_enabled(PGi))
        molec.push_back(PGi + 1);
      else // get enabled descendants
      {
        _PGtree.get_enabled_descendants(PGi, ids);

#ifdef DEBUG_HIERARCHICAL_MESH
        if (ids.empty())
        {
          E_Int pid = _PGtree.parent(PGi);
          E_Int gpid = E_IDX_NONE;
          E_Int nbc = _PGtree.nb_children(PGi);
          bool pid_is_enabled = (pid != E_IDX_NONE) ? _PGtree.is_enabled(pid) : false;

          bool gpid_is_enabled = false;
          if (pid != E_IDX_NONE){
            gpid = _PGtree.parent(pid);
            gpid_is_enabled = (gpid != E_IDX_NONE) ? _PGtree.is_enabled(gpid) : false;
          }

          std::cout << "faulty PGi : " << PGi << std::endl;
          std::cout << "PHi : " << i << std::endl;
          std::cout << "F2E : " << _F2E(0,PGi) << "/" << _F2E(1,PGi) << std::endl;
          std::cout << "parent ? : " << pid << " enabled ? : " << pid_is_enabled << std::endl;
          std::cout << "grand parent ? : " << gpid << " enabled ? : " << gpid_is_enabled<< std::endl;
          std::cout << "nbc ? : " << nbc << std::endl;
          if (nbc != 0)
          {
            const E_Int* childz = _PGtree.children(PGi);
            for (size_t i=0; i < nbc; ++i)
              std::cout << "child : " << i << " : " << childz[i] << std::endl;
          }

          medith::write("faultyPG", _crd, _ng.PGs.get_facets_ptr(PGi), _ng.PGs.stride(PGi), 1);
          medith::write("PH557", _crd, _ng, i);
          medith::write("PH6", _crd, _ng, _F2E(0,PGi));
        }
#endif
        assert (!ids.empty());

        K_CONNECT::IdTool::shift(ids, 1);
        molec.insert(molec.end(), ALL(ids));
      }
    }

    new_phs.add(molec.size(), &molec[0]);  //alexis : set _type for children ??
  }

  ngo.PGs = _ng.PGs; // on copie tous les polygones
  ngo.PHs = new_phs;
  ngo.PHs.updateFacets();

  _pgnids.clear();
  _phnids.clear();
  ngo.remove_unreferenced_pgs(_pgnids, _phnids);

  //for history (BC and Join preserving)

  pgoids.clear();
  K_CONNECT::IdTool::init_inc(pgoids, _ng.PGs.size());
  
  Vector_t<E_Int> old_pgoids;
  _PGtree.get_oids(old_pgoids);
  
  for (size_t i = 0; i < _pgnids.size(); ++i)
  {
    //old_pgoids cannot have IDX_NONE : new entities must be self referring
    if (_pgnids[i] != IDX_NONE)pgoids[_pgnids[i]] = old_pgoids[i];
  }
}

///
template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t>
void hierarchical_mesh<ELT_t, STYPE, ngo_t>::extract_enabled_phs(ngon_type& filtered_ng) const
{
  filtered_ng.PGs = _ng.PGs;
        
  //E_Int sz = _PGtree.get_parent_size();
  E_Int sz = _ng.PHs.size();
        
  for (int i = 0; i < sz; i++)
  {
    if (_PHtree.is_enabled(i) == true)
    {
      const E_Int* p = _ng.PHs.get_facets_ptr(i);
      E_Int s = _ng.PHs.stride(i);
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
void hierarchical_mesh<ELT_t, STYPE, ngo_t>::extract_enabled_pgs_descendance(E_Int PGi, bool reverse, std::vector<E_Int>& ids)
{
  ids.clear();
  if (_PGtree.is_enabled(PGi))
    return;

  reordering_func F{ nullptr };
  if (_ng.PGs.stride(PGi) == 3)
    F = subdiv_pol<K_MESH::Triangle, STYPE>::reorder_children;
  else if (_ng.PGs.stride(PGi) == 4)
    F = subdiv_pol<K_MESH::Quadrangle, STYPE>::reorder_children;
  
  __extract_enabled_pgs_descendance(PGi, F, reverse, ids);
}

///
template <>
void hierarchical_mesh<K_MESH::Polyhedron<0>, NUGA::ISO_HEX, ngon_type>::extract_enabled_pgs_descendance(E_Int PGi, bool reverse, std::vector<E_Int>& pointlist)
{
  //todo
}

///
template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t>
void hierarchical_mesh<ELT_t, STYPE, ngo_t>::__extract_enabled_pgs_descendance(E_Int PGi, reordering_func F, bool reverse, std::vector<E_Int>& ids)
{
  //
  if (_PGtree.is_enabled(PGi))
  {
    ids.push_back(PGi);
    return;
  }

  E_Int nbc = _PGtree.nb_children(PGi);
  const E_Int* pchild = _PGtree.children(PGi);

  STACK_ARRAY(E_Int, nbc, children);
  for (E_Int i = 0; i < nbc; ++i) children[i] = pchild[i];

  // to put in asked ref frame
  F(children.get(), reverse, 0);// 0 because shift_geom must have been called

  //
  for (E_Int i = 0; i < nbc; ++i)
    __extract_enabled_pgs_descendance(children[i], F, reverse, ids);
}

///
template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t>
template <typename InputIterator>
void hierarchical_mesh<ELT_t, STYPE, ngo_t>::get_higher_level_neighbours
(E_Int PHi, E_Int PGi, InputIterator neighbours, E_Int& nb_neighbours) const
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
template <typename InputIterator>
void hierarchical_mesh<ELT_t, STYPE, ngo_t>::get_enabled_neighbours
(E_Int PHi, InputIterator neighbours, E_Int& nb_neighbours) const
{
  // fill in up to 24 enabled neighbours and gives the number of enabled neighbours
  const E_Int* p = _ng.PHs.get_facets_ptr(PHi);
  E_Int nb_faces = _ng.PHs.stride(PHi);
    
#ifdef DEBUG_HIERARCHICAL_MESH
  assert (nb_neighbours == 0);
#endif
  for (int i = 0; i < nb_faces; ++i)
  {
    E_Int PGi = p[i] - 1;
      
    E_Int PH = NEIGHBOR(PHi, _F2E, PGi);

    if (PH == IDX_NONE) // returns only the enabled neighbours
      continue;
        
    if ( (_PHtree.is_enabled(PH)) || ((_PHtree.get_level(PH) > 0) && (_PHtree.is_enabled(_PHtree.parent(PH)))) )
      neighbours[nb_neighbours++] = PH; // no children : add the PH
    else
      get_higher_level_neighbours(PHi, PGi, neighbours, nb_neighbours); // add the 4 PH of higher level
  }
}

///
template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t>
bool hierarchical_mesh<ELT_t, STYPE, ngo_t>::is_initialised() const
{
  return _initialized;
}

///
template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t>
void hierarchical_mesh<ELT_t, STYPE, ngo_t>::get_cell_center(E_Int PHi, E_Float* center) const
{    
  ELT_t::iso_barycenter(_crd, _ng.PGs, _ng.PHs.get_facets_ptr(PHi), _ng.PHs.stride(PHi), 1, center);
}


///
template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t>
void hierarchical_mesh<ELT_t, STYPE, ngo_t>::extract_plan(E_Int PGi, bool reverse, E_Int i0, pg_arr_t& plan) const
{
  plan.clear();

  /*E_Int ret{ 0 };*/

  if (_ng.PGs.stride(PGi) == 3)
    /*ret = */join_plan<pg_arr_t>::extract_compact_enabled_tree(_PGtree, PGi, subdiv_pol<K_MESH::Triangle, STYPE>::reorder_children, reverse, i0, plan);
  else if (_ng.PGs.stride(PGi) == 4)
    /*ret = */join_plan<pg_arr_t>::extract_compact_enabled_tree(_PGtree, PGi, subdiv_pol<K_MESH::Quadrangle, STYPE>::reorder_children, reverse, i0, plan);
}

template <>
void hierarchical_mesh<K_MESH::Polyhedron<0>, NUGA::ISO_HEX, ngon_type>::extract_plan(E_Int PGi, bool reverse, E_Int i0, pg_arr_t& plan) const
{
  //todo
}

///
template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t>
void hierarchical_mesh<ELT_t, STYPE, ngo_t>::enable_PGs()
{
  // INFO : at the end, some join PG can be wrongly disbaled
  // (those which ensure conformal join with the other side)
  // indeed, for some join, children must be enabled instead of natural face
  // but these are fixed in join_sensor<mesh_t>::update()

  //reset
  _PGtree.reset_enabled(_ng.PGs.size(), false);

  for (E_Int PHi = 0; PHi < _ng.PHs.size(); ++PHi)
  {
    if (!_PHtree.is_enabled(PHi)) continue;

    const E_Int* faces = _ng.PHs.get_facets_ptr(PHi);
    E_Int nfaces = _ng.PHs.stride(PHi);

    for (E_Int j = 0; j < nfaces; ++j)
    {
      E_Int PGi = faces[j] - 1;
      E_Int PHn = NEIGHBOR(PHi, _F2E, PGi);

      if (PHn == IDX_NONE)
        _PGtree.enable(PGi);
      else if (_PHtree.is_enabled(PHn))
        _PGtree.enable(PGi);
    }
  }
}

//
template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t>
void hierarchical_mesh<ELT_t, STYPE, ngo_t>::update_BCs()
{
  //std::cout << "join_t<mesh_t>::update() : begin" << std::endl;
  std::vector<E_Int> ids;
  E_Int nb_pgs = _ng.PGs.size();

  // update pointlists
  for (size_t i=0; i < BCptLists.size(); ++i)
  {
    std::vector<E_Int>& ptlist = BCptLists[i];
    std::vector<E_Int> new_ptlist;

    for (size_t i = 0; i < ptlist.size(); ++i)
    {
      E_Int PGi = ptlist[i] - _idx_start;
      //std::cout << "PGi : " << PGi << std::endl;

      if (PGi < 0 || PGi >= nb_pgs)
      {
        std::cout << "update_BCs : WARNING : wrong PG id : " << PGi << std::endl;
        continue;
      }
      
      if (_PGtree.is_enabled(PGi))
        new_ptlist.push_back(PGi);
      else // look in the genealogy where are the enabled
      {
        ids.clear();
        extract_enabled_pgs_descendance(PGi, false/*reverse not required*/, ids);

        if (!ids.empty()) //refinement
          new_ptlist.insert(new_ptlist.end(), ALL(ids));
        else //agglo : get the enabled ancestor
        {
          E_Int pid{IDX_NONE};
          _PGtree.get_enabled_parent(PGi, pid);
          assert (pid != IDX_NONE);
          new_ptlist.push_back(pid);
        }
      }
    }

    //std::cout << "old BC sz vs new : " << ptlist.size() << " vs " << new_ptlist.size() << std::endl;

    ptlist = new_ptlist;
    // remove duplicates due to agglo
    K_CONNECT::IdTool::compress_unic(ptlist);
    // for the outside world
    K_CONNECT::IdTool::shift(ptlist, _idx_start);
  }
}

///
template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t>
E_Int hierarchical_mesh<ELT_t, STYPE, ngo_t>::project_cell_center_sol_order1
(std::vector<std::vector<E_Float>>& cfields)
{
  size_t nbf = cfields.size();

  if (nbf == 0) return 0;

  //check all fields have same size
  size_t fsz{0};
  for (size_t f=0; f < nbf; ++f)
  {
    if (fsz==0)fsz = cfields[f].size();
    else if (fsz != cfields[f].size())
    {
      std::cout << "project_cell_center_sol_order1 : ERROR. input fields sizes are inconsistent between them." << std::endl;
      return 1;
    }
  }

  if (_enabledPHi.empty()) // one-shot adaption or first time in adaptCellsDyn
    _enabledPHi.resize(_nb_phs0, true);
  else
  {
    /*std::cout << "next pass !!!" << std::endl;
    std::cout << "nb phs : " << _ng.PHs.size() << std::endl;
    std::cout << "_enabledPHi sz : " << _enabledPHi.size() << std::endl;
    std::cout << "fsz : " << fsz << std::endl;*/
  }

  size_t nb_enabledi = std::count(ALL(_enabledPHi), true);
  
  if (fsz != nb_enabledi)
  {
    std::cout << "project_cell_center_sol_order1 : ERROR. input fields sizes are inconsistent with hmesh." << std::endl;
    return 1;
  }

  // compute indir field to PH in _enabledPHi
  std::vector<E_Int> PHid1(fsz);
  size_t c{0};
  for (size_t i=0; i < _enabledPHi.size(); ++i)
  {
    if (_enabledPHi[i]) PHid1[c++]=i;
  }

  // compute indir field to PH and reverse in current enabled vector
  std::vector<E_Int> nPHid, nfid(_ng.PHs.size(), IDX_NONE);
  c=0;
  for (E_Int i=0; i < _ng.PHs.size(); ++i)
  {
    if (_PHtree.is_enabled(i)) 
    {
      nfid[i] = nPHid.size();
      nPHid.push_back(i);
    }
  }
  E_Int fsz2 = nPHid.size();

  // initialize new fields
  std::vector<std::vector<E_Float>> new_fields(nbf);
  for (size_t f = 0; f < nbf; ++f)
    new_fields[f].resize(fsz2, 0.);

  std::vector<E_Int> children;
  std::vector<bool> to_agglo(fsz2, false);

  // loop on old fields and apply to new
  for (size_t i=0; i < PHid1.size(); ++i)
  {
    E_Int PHj = PHid1[i];
    assert (_enabledPHi[PHj] == true); //by defininition

    if (_PHtree.is_enabled(PHj)) // still enabled, just pass the values
    {
      for (size_t f = 0; f < nbf; ++f)
        new_fields[f][nfid[PHj]] = cfields[f][i];
    }
    else // refinement :look in the genealogy where are the enabled
    {
      _PHtree.get_enabled_descendants(PHj, children);

      if (!children.empty())
      {
        for (size_t c=0; c < children.size(); ++c)
        {
          for (size_t f = 0; f < nbf; ++f)
            new_fields[f][nfid[children[c]]] = cfields[f][i];
        }
      }
      else //agglo : get the enabled ancestor
      {
        E_Int pid{IDX_NONE};
        _PHtree.get_enabled_parent(PHj, pid);
        assert (pid != IDX_NONE);

        //compute PHj volume
        E_Float v;
        K_MESH::Polyhedron<0>::volume<DELAUNAY::Triangulator>(_crd, _ng.PGs, _ng.PHs.get_facets_ptr(PHj), _ng.PHs.stride(PHj), v, true);

        for (size_t f = 0; f < nbf; ++f)
          new_fields[f][nfid[pid]] += cfields[f][i] * ::fabs(v); // accumulate mass

        to_agglo[nfid[pid]] = true;

      }
    }
  }

  // go back to density for agglomerated cells
  for (size_t i=0; i < to_agglo.size(); ++i) // i is ith field value
  {
    if (!to_agglo[i]) continue;

    E_Int PH = nPHid[i];

    E_Float v;
    K_MESH::Polyhedron<0>::volume<DELAUNAY::Triangulator>(_crd, _ng.PGs, _ng.PHs.get_facets_ptr(PH), _ng.PHs.stride(PH), v, true);

    if (::fabs(v) < ZERO_M) continue;

    for (size_t f = 0; f < nbf; ++f)
      new_fields[f][i] /= v;
  }

  cfields = new_fields;
  _enabledPHi = _PHtree.get_enabled();

  return 0;
}

}

#endif

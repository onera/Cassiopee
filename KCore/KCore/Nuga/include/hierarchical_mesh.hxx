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


#define NEIGHBOR(PHi, _F2E, PGi) ( (_F2E(0,PGi) == PHi) ? _F2E(1,PGi) : _F2E(0,PGi) )

namespace NUGA
{
  
  enum eSUBDIV_TYPE { ISO, ANISO};

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

template <> 
class hmesh_trait<K_MESH::Basic, ISO>
{
  public:
    using arr_type = K_FLD::IntArray;
    
    static const E_Int PHNBC = 8;
    static const E_Int PGNBC = 4;
};

template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t = ngon_type, typename crd_t = K_FLD::FloatArray>
class hierarchical_mesh
{
  public:

    using elt_type = ELT_t;
    using htrait = hmesh_trait<ELT_t, STYPE>;
    using arr_t = typename htrait::arr_type; //InttArray or ngon_unit
    using self_type = hierarchical_mesh<ELT_t, STYPE, ngo_t, crd_t>;
    
    crd_t&                    _crd;             // Coordinates
    ngo_t&                    _ng;              // NGON mesh
    tree<arr_t>               _PGtree, _PHtree; // Polygons/Polyhedra hierarchy
    K_FLD::IntArray           _F2E;             // neighborhood data structure : 2 X (nb_pgs) array : for each PG gives the left and right cells ids that share this PG
    bool                      _initialized;     // flag to avoid initialization more than once.

    ///
    hierarchical_mesh(crd_t& crd, ngo_t & ng):_crd(crd), _ng(ng), _PGtree(ng.PGs, htrait::PGNBC), _PHtree(ng.PHs, htrait::PHNBC), _initialized(false){}

    ///
    E_Int init();
    ///
    E_Int adapt(Vector_t<E_Int>& adap_incr, bool do_agglo);
  
    /// face-conformity
    void conformize();
    /// Keep only enabled PHs
    void filter_ngon(ngon_type& filtered_ng);
    ///
    static void refine_PGs(const Vector_t<E_Int> &PHadap, 
            ngo_t& ng, tree<arr_t> & PGtree, crd_t& crd, K_FLD::IntArray & F2E, 
            std::map<K_MESH::NO_Edge,E_Int>& ecenter);
    ///
    static void refine_PHs(const Vector_t<E_Int> &PHadap, 
                            ngo_t& ng, tree<arr_t> & PGtree, tree<arr_t> & PHtree,crd_t& crd, K_FLD::IntArray & F2E);
    ///
    static void get_nodes_PHi(E_Int* nodes, E_Int PHi, E_Int centroid_id, E_Int* BOT, E_Int* TOP, E_Int* LEFT, E_Int* RIGHT, E_Int* FRONT, E_Int* BACK,
                              ngo_t& ng, K_FLD::FloatArray& crd, K_FLD::IntArray & F2E, tree<arr_t> & PGtree);
    ///
    static void get_nodes_PHi_T(E_Int* nodes, E_Int PHi, E_Int centroid_id, E_Int* BOT, E_Int* F1, E_Int* F2, E_Int* F3,
                                ngo_t& ng, K_FLD::FloatArray& crd, K_FLD::IntArray & F2E, tree<arr_t> & PGtree);
    ///
    static void retrieve_ordered_data(E_Int PGi, E_Int i0, bool reorient, E_Int* four_childrenPG, E_Int* LNODES,
                                      ngo_t& ng, tree<arr_t> & PGtree);
    ///
    static bool need_a_reorient(E_Int PGi, E_Int PHi, bool oriented_if_R, 
                        K_FLD::IntArray & F2E);
    ///
    static E_Int get_i0(E_Int* pFace, E_Int common_node, E_Int* nodes, E_Int nb_edges_face);
    ///
    static void update_F2E(E_Int PHi, E_Int PHchildr0, E_Int* INT, E_Int* BOT, E_Int* TOP, E_Int* LEFT, E_Int* RIGHT, E_Int* FRONT, E_Int* BACK,
                           ngo_t& ng, K_FLD::IntArray & F2E, tree<arr_t> & PGtree);
    ///
    static void update_F2E_T(E_Int PHi, E_Int PHchildr0, E_Int* INT, E_Int* BOT, E_Int* F1, E_Int* F2, E_Int* F3, E_Int ndiag, 
                            ngo_t& ng, K_FLD::IntArray & F2E, tree<arr_t> & PGtree);
    ///
    static void update_children_F2E(E_Int PGi, E_Int side, tree<arr_t> & PGtree, K_FLD::IntArray & F2E);
    ///
    void get_cell_center(E_Int PHi, E_Float* center);
    ///
    void get_enabled_neighbours(E_Int PHi, E_Int* neighbours, E_Int& nb_neighbours);
    ///
    void get_higher_level_neighbours(E_Int PHi, E_Int PGi, E_Int* neighbours, E_Int& nb_neighbours);
    ///
    bool enabled_neighbours(E_Int PHi); // return true if the PHi-th PH only has enabled neighbours
    ///
    void smooth(Vector_t<E_Int>& adap_incr);
    ///
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
#endif

  private:
    std::map<K_MESH::NO_Edge,E_Int> _ecenter;

    ///
    void __conformize_next_lvl(Vector_t<E_Int>& molec, E_Int PGi, E_Int i);
        
    static void __compute_edge_center(const Vector_t<E_Int> &PGlist, 
                                      ngo_t& ng, crd_t& crd, 
                                      std::map<K_MESH::NO_Edge,E_Int> & ecenter);
    //void __compute_face_centers(K_FLD::FloatArray& crd, const typename ngo_t::unit_type& pgs, const Vector_t<E_Int> &PGlist, Vector_t<E_Int>& fcenter);
    static inline void __compute_face_center(const crd_t& crd, const E_Int* nodes, E_Int nb_nodes, E_Float* C);
    static void __compute_cell_center(const crd_t& crd, const E_Int* nodes27, E_Float* C);
  
};

///
template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t, typename crd_t>
E_Int hierarchical_mesh<ELT_t, STYPE, ngo_t, crd_t>::init()
{
  if (_initialized) return 0;
  
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
  
  // sort the PGs of the initial NGON
  for (int i = 0; i < (int)_ng.PHs.size(); ++i)
  {
    ELT_t::reorder_pgs(_ng,_F2E,i);
  }

  _initialized = true;
  
  return err;
}

#ifdef OUTPUT_ITER_MESH
 static E_Int iter = 1;
#endif

template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t, typename crd_t>
E_Int hierarchical_mesh<ELT_t, STYPE, ngo_t, crd_t>::adapt(Vector_t<E_Int>& adap_incr, bool do_agglo)
{
  E_Int err(0);

  Vector_t<E_Int> PHadap, PHagglo, PHref, PGref;


  
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
    PHadap.clear();
    PHagglo.clear();

    // PHadap : Get the list of enabled PHs for which adap_incr[PHi] != 0
    if (max > 0) // there is a need to refine
    {
      for (size_t i = 0; i < _ng.PHs.size(); ++i)
        if (adap_incr[i] > 0 && _PHtree.is_enabled(i)) PHadap.push_back(i);
    }
    
    if (min < 0) // there is a need to agglomerate
    {
      for (size_t i = 0; i < _ng.PHs.size(); ++i)
        if (adap_incr[i] == -1) PHagglo.push_back(i);
    }
    
    if (PHadap.empty() && PHagglo.empty()) break; // adapted

    // refine PGs : create missing children (those PGi with _PGtree.nb_children(PGi) == 0)
    refine_PGs(PHadap, _ng, _PGtree, _crd, _F2E, _ecenter);

//    std::cout << "is hybrid after refine_PGs ? " << is_hybrid(_ng) << std::endl;

    // refine PHs with missing children (those PHi with _PHtree.nb_children(PHi) == 0)
    refine_PHs(PHadap, _ng, _PGtree, _PHtree, _crd, _F2E);

//    std::cout << "is hybrid after refine_PHs ? " << is_hybrid(_ng) << std::endl;

    // enable the right PHs & their levels, disable subdivided PHs
    adap_incr.resize(_ng.PHs.size(),0);
    for (size_t i = 0; i < PHadap.size(); ++i)
    {
      E_Int PHi = PHadap[i];

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
    
    _ng.PGs.updateFacets();
    _ng.PHs.updateFacets();

#ifdef DEBUG_HIERARCHICAL_MESH  
//    std::cout << "is hybrid at the end of adap ? " << is_hybrid(_ng) << std::endl;
    if (! _ng.attributes_are_consistent()) return false;
#endif
    
#ifdef OUTPUT_ITER_MESH
    ngon_type filtered_ng;
    filter_ngon(filtered_ng);

    std::ostringstream o;
    o << "NGON_it_" << iter << ".plt"; // we create a file at each iteration
    K_FLD::IntArray cnto;
    filtered_ng.export_to_array(cnto);
    MIO::write(o.str().c_str(), _crd, cnto, "NGON");     
    ++iter;
#endif
    
    if (!loop) break;
  }

  return err;
}

template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t, typename crd_t>
void hierarchical_mesh<ELT_t, STYPE, ngo_t, crd_t>::conformize()
{
  ngon_unit new_phs;
  Vector_t<E_Int> molec;

  E_Int nb_phs = _ng.PHs.size();
  for (E_Int i = 0; i < nb_phs; ++i)
  {
    if (!_PHtree.is_enabled(i)) continue;
    
    molec.clear();
    E_Int s = _ng.PHs.stride(i);
    E_Int* pPGi = _ng.PHs.get_facets_ptr(i);

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
        else // append the 4 children
        {
          E_Int nbc = _PGtree.nb_children(PGi);
          for (E_Int c=0; c < nbc; ++c){
            E_Int PG_f= *(_PGtree.children(PGi)+c);
            __conformize_next_lvl(molec, PG_f, PHn);
          }
        }        
      }
    }
    
    new_phs.add(molec.size(), &molec[0]);  
  }

  _ng.PHs = new_phs;
  _ng.PHs.updateFacets();
  
  std::vector<E_Int> pgnids, phnids;
  _ng.remove_unreferenced_pgs(pgnids, phnids);
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
  filtered_ng.PGs = _ng.PGs;
        
  //E_Int sz = _PGtree.get_parent_size();
  E_Int sz = _ng.PHs.size();
        
  for (int i = 0; i < sz; i++)
  {
    if (_PHtree.is_enabled(i) == true)
    {
      E_Int* p = _ng.PHs.get_facets_ptr(i);
      E_Int s = _ng.PHs.stride(i);
      filtered_ng.PHs.add(s,p);
    }
  }
    
  filtered_ng.PGs.updateFacets();
  filtered_ng.PHs.updateFacets();
    
  Vector_t<E_Int> pgnids, phnids;
  filtered_ng.remove_unreferenced_pgs(pgnids, phnids);
}


///
template <>
void hierarchical_mesh<K_MESH::Hexahedron, eSUBDIV_TYPE::ISO, ngon_type, K_FLD::FloatArray>::__compute_cell_center
(const K_FLD::FloatArray& crd, const E_Int* nodes27, E_Float* centroid)
{
  E_Int new_bary[8];
  new_bary[0] = nodes27[0];
  new_bary[1] = nodes27[1];
  new_bary[2] = nodes27[2];
  new_bary[3] = nodes27[3];
  new_bary[4] = nodes27[4];
  new_bary[5] = nodes27[5];
  new_bary[6] = nodes27[6];
  new_bary[7] = nodes27[7];

  K_MESH::Polyhedron<STAR_SHAPED>::iso_barycenter(crd, new_bary, 8, 1,centroid);
}

///
template <>
void hierarchical_mesh<K_MESH::Hexahedron, eSUBDIV_TYPE::ISO, ngon_type, K_FLD::FloatArray>::refine_PGs
(const Vector_t<E_Int> &PHadap, 
 ngon_type& ng, tree<arr_t> & PGtree, K_FLD::FloatArray& crd, K_FLD::IntArray & F2E,         
 std::map<K_MESH::NO_Edge,E_Int>& ecenter)
{
  
  E_Int nb_phs = PHadap.size();
  
  // Gets PGs to refine
  Vector_t<E_Int> PGref;
  {
    E_Int nb_pgs(ng.PGs.size()), nb_pgs_ref(0);
    Vector_t<bool> is_PG_to_refine(nb_pgs, false);
    //
    for (E_Int i = 0; i < nb_phs; ++i)
    {
      E_Int PHi = PHadap[i];

      E_Int nb_faces = ng.PHs.stride(PHi); 
      E_Int* faces = ng.PHs.get_facets_ptr(PHi);
      
#ifdef DEBUG_HIERARCHICAL_MESH
      assert (nb_faces == 6);
#endif
      for (E_Int j = 0; j < nb_faces; ++j)
      {
        E_Int PGi = * (faces + j) - 1;
        
        if (PGtree.nb_children(PGi) == 0) // leaf PG => to refine
        {
          is_PG_to_refine[PGi] = true;
          ++nb_pgs_ref;
        }
      }
    }

    PGref.reserve(nb_pgs_ref);
    for (E_Int i = 0; i < nb_pgs; ++i)
    {
      if (is_PG_to_refine[i])
        PGref.push_back(i);
    }
  }

  E_Int nb_pgs_ref = PGref.size();
  // Compute Edges refine points
  __compute_edge_center(PGref, ng, crd, ecenter);
    
  // Reserve space for children in the tree
  PGtree.resize(PGref);
  // And in the mesh : each Q4 is split into 4 Q4
  E_Int nb_pgs0 = ng.PGs.size();
  ng.PGs.expand_n_fixed_stride(4*nb_pgs_ref, 4/*Q4 stride*/);

  F2E.resize(2,nb_pgs0+4*nb_pgs_ref,E_IDX_NONE);

#ifdef DEBUG_HIERARCHICAL_MESH
  Vector_t<E_Int> ids;
#endif
  
  // face centers  
  Vector_t<E_Int> fcenter(nb_pgs_ref, E_IDX_NONE);  
  E_Int pos = crd.cols();
  crd.resize(3, pos + nb_pgs_ref);  
  
  K_MESH::NO_Edge noE;  
  
  // Refine them
#ifndef DEBUG_HIERARCHICAL_MESH
#pragma omp parallel for private(noE)
#endif
  for (E_Int i = 0; i < nb_pgs_ref; ++i)
  {
    E_Int PGi = PGref[i];
    E_Int PGichildr[4];
    E_Int q9[9];
   
    E_Int* nodes = ng.PGs.get_facets_ptr(PGi);
    
    // Centroid calculation
    __compute_face_center(crd, ng.PGs.get_facets_ptr(PGi), ng.PGs.stride(PGi), crd.col(pos+i));
    fcenter[i] = pos+i;
    
#ifdef DEBUG_HIERARCHICAL_MESH
    assert (ng.PGs.stride(PGi) == 4);
#endif
    
    for (E_Int n=0; n < 4; ++n)
    {
        
      PGichildr[n] = nb_pgs0 + 4*i + n;
      // children have the L & R elements of the father
      F2E(0,PGichildr[n]) = F2E(0,PGi);
      F2E(1,PGichildr[n]) = F2E(1,PGi);
        
      q9[n] = *(nodes + n);
      noE.setNodes(*(nodes + n), *(nodes + (n+1)%4));
      q9[n+4] = ecenter[noE];
    }
    
    q9[8] = fcenter[i] + 1;
    
    // set them in _ng.PGs
    E_Int* q41 = ng.PGs.get_facets_ptr(PGichildr[0]);
    E_Int* q42 = ng.PGs.get_facets_ptr(PGichildr[1]);
    E_Int* q43 = ng.PGs.get_facets_ptr(PGichildr[2]);
    E_Int* q44 = ng.PGs.get_facets_ptr(PGichildr[3]);
    
    NUGA::Q9::splitQ4(crd, q9, q41, q42, q43, q44);
    
#ifdef DEBUG_HIERARCHICAL_MESH
//    ngon_unit pgs;
//    pgs.add(4,q41);
//    pgs.add(4,q42);
//    pgs.add(4,q43);
//    pgs.add(4,q44);
//    NGDBG::draw_PGT3s(_crd,pgs);
#endif   

#ifdef DEBUG_HIERARCHICAL_MESH
//    std::cout << "q41 nb nodes : " << _ng.PGs.stride(PGichildr[0]) << std::endl;
//    for (size_t i=0; i < 9; ++i)
//      std::cout << q9[i] << "/";
//    std::cout << std::endl;
//    
//    std::cout << q41[0] << "/" << q41[1] << "/" << q41[2] << "/" << q41[3] << std::endl;
#endif

    // set them in the tree
    PGtree.set_children(PGi, PGichildr, 4);
   
#ifdef DEBUG_HIERARCHICAL_MESH
    ids.insert(ids.end(), PGichildr, PGichildr+4);
#endif
    
  }
  
#ifdef DEBUG_HIERARCHICAL_MESH
  //NGDBG::draw_PGs(_crd, _ng.PGs, ids, false);
#endif

}

///
template <>
void hierarchical_mesh<K_MESH::Tetrahedron, eSUBDIV_TYPE::ISO, ngon_type, K_FLD::FloatArray>::refine_PGs
(const Vector_t<E_Int> &PHadap,
 ngon_type& ng, tree<arr_t> & PGtree, K_FLD::FloatArray& crd, K_FLD::IntArray & F2E,         
 std::map<K_MESH::NO_Edge,E_Int>& ecenter)
{
  
  E_Int nb_phs = PHadap.size();
  
  // Gets PGs to refine
  Vector_t<E_Int> PGref;
  {
    E_Int nb_pgs(ng.PGs.size()), nb_pgs_ref(0);
    Vector_t<bool> is_PG_to_refine(nb_pgs, false);
    //
    for (E_Int i = 0; i < nb_phs; ++i)
    {
      E_Int PHi = PHadap[i];

      E_Int nb_faces = ng.PHs.stride(PHi); 
      E_Int* faces = ng.PHs.get_facets_ptr(PHi);
      
#ifdef DEBUG_HIERARCHICAL_MESH
      //assert (nb_faces == 6);
      assert (nb_faces == 4);
#endif
      for (E_Int j = 0; j < nb_faces; ++j)
      {
        E_Int PGi = * (faces + j) - 1;
        
        if (PGtree.nb_children(PGi) == 0) // leaf PG => to refine
        {
          is_PG_to_refine[PGi] = true;
          ++nb_pgs_ref;
        }
      }
    }

    PGref.reserve(nb_pgs_ref);
    for (E_Int i = 0; i < nb_pgs; ++i)
    {
      if (is_PG_to_refine[i])
        PGref.push_back(i);
    }
  }

  E_Int nb_pgs_ref = PGref.size();
  // Compute Edges refine points
  __compute_edge_center(PGref, ng, crd, ecenter);
    
  // Reserve space for children in the tree
  PGtree.resize(PGref);
  // And in the mesh : each Q4 is split into 4 Q4
  E_Int nb_pgs0 = ng.PGs.size();
  //_ng.PGs.expand_n_fixed_stride(4*nb_pgs_ref, 4/*Q4 stride*/);
  ng.PGs.expand_n_fixed_stride(4*nb_pgs_ref, 3/*Q4 stride*/);

  F2E.resize(2,nb_pgs0+4*nb_pgs_ref,E_IDX_NONE);

#ifdef DEBUG_HIERARCHICAL_MESH
  Vector_t<E_Int> ids;
#endif

  K_MESH::NO_Edge noE;  
  
  // Refine them
#ifndef DEBUG_HIERARCHICAL_MESH
#pragma omp parallel for private(noE)
#endif
  for (E_Int i = 0; i < nb_pgs_ref; ++i)
  {
    E_Int PGi = PGref[i];
    E_Int PGichildr[4];
    E_Int q6[6]; // définition 6 noeuds subdivision de PGi
   
    E_Int* nodes = ng.PGs.get_facets_ptr(PGi);
    
#ifdef DEBUG_HIERARCHICAL_MESH
    assert (ng.PGs.stride(PGi) == 3);
#endif
    
    for (E_Int n=0; n < 3; ++n) // définition des PGs enfants
    {
        
      PGichildr[n] = nb_pgs0 + 4*i + n;
      // children have the L & R elements of the father
      F2E(0,PGichildr[n]) = F2E(0,PGi);
      F2E(1,PGichildr[n]) = F2E(1,PGi);
        
      q6[n] = *(nodes + n);
      noE.setNodes(*(nodes + n), *(nodes + (n+1)%3));
      q6[n+3] = ecenter[noE]; // définition du centre de chaque arête (tétra 3 arêtes)
    }
    
    PGichildr[3] = nb_pgs0 + 4*i + 3;
    F2E(0,PGichildr[3]) = F2E(0,PGi);
    F2E(1,PGichildr[3]) = F2E(1,PGi);
    
    // set them in _ng.PGs
    E_Int* q41 = ng.PGs.get_facets_ptr(PGichildr[0]);
    E_Int* q42 = ng.PGs.get_facets_ptr(PGichildr[1]);
    E_Int* q43 = ng.PGs.get_facets_ptr(PGichildr[2]);
    E_Int* q44 = ng.PGs.get_facets_ptr(PGichildr[3]);
    
    NUGA::Q9::splitQ4T(crd, q6, q41, q42, q43, q44);
    
    //std::cout << "q41 : " << q41[0] << "/" << 
    
#ifdef DEBUG_HIERARCHICAL_MESH
//    ngon_unit pgs;
//    pgs.add(3,q41);
//    pgs.add(3,q42);
//    pgs.add(3,q43);
//    pgs.add(3,q44);
//    NGDBG::draw_PGT3s(_crd,pgs);
//    ngon_unit pgs2;
//    pgs2.add(3, _ng.PGs.get_facets_ptr(PGi));
//    NGDBG::draw_PGT3s(_crd,pgs2);
#endif   

#ifdef DEBUG_HIERARCHICAL_MESH
//    std::cout << "q41 nb nodes : " << _ng.PGs.stride(PGichildr[0]) << std::endl;
//    for (size_t i=0; i < 6; ++i)
//      std::cout << q6[i] << "/";
//    std::cout << std::endl;
//    
//    std::cout << q41[0] << "/" << q41[1] << "/" << q41[2] << "/" << q41[3] << std::endl;
#endif

    // set them in the tree
    PGtree.set_children(PGi, PGichildr, 4);
   
#ifdef DEBUG_HIERARCHICAL_MESH
    ids.insert(ids.end(), PGichildr, PGichildr+4);
#endif
    
  }
  
#ifdef DEBUG_HIERARCHICAL_MESH
  NGDBG::draw_PGs(crd, ng.PGs, ids, false);
#endif

}
  

template <>
void hierarchical_mesh<K_MESH::Basic, eSUBDIV_TYPE::ISO, ngon_type, K_FLD::FloatArray>::refine_PGs
(const Vector_t<E_Int> &PHadap, 
 ngon_type& ng, tree<arr_t> & PGtree, K_FLD::FloatArray& crd, K_FLD::IntArray & F2E,                 
 std::map<K_MESH::NO_Edge,E_Int>& ecenter)
{
    Vector_t<E_Int> PHadap1, PHadap2;
    E_Int nb_phs = PHadap.size();
    for (E_Int i=0; i< nb_phs; i++){
        E_Int PHi=PHadap[i];
        E_Int s= ng.PHs.stride(PHi);
        if (K_MESH::Polyhedron<0>::is_HX8(ng.PGs, ng.PHs.get_facets_ptr(PHi), s) ==true){
            PHadap1.push_back(PHi);
        }
        else if (K_MESH::Polyhedron<0>::is_TH4(ng.PGs, ng.PHs.get_facets_ptr(PHi), s) ==true){
            PHadap2.push_back(PHi);
        }
    }

    hierarchical_mesh<K_MESH::Hexahedron, eSUBDIV_TYPE::ISO, ngon_type, K_FLD::FloatArray>::refine_PGs(PHadap1, 
                                                                                                        ng, PGtree, crd, F2E
                                                                                                        , ecenter);
    hierarchical_mesh<K_MESH::Tetrahedron, eSUBDIV_TYPE::ISO, ngon_type, K_FLD::FloatArray>::refine_PGs(PHadap2, 
                                                                                                        ng, PGtree, crd, F2E
                                                                                                        , ecenter);    
}

template <>
void hierarchical_mesh<K_MESH::Hexahedron, eSUBDIV_TYPE::ISO, ngon_type, K_FLD::FloatArray>::retrieve_ordered_data
(E_Int PGi, E_Int i0, bool reorient, E_Int* four_childrenPG, E_Int* LNODES,
ngo_t& ng, tree<arr_t> & PGtree)
{
  E_Int* pN = ng.PGs.get_facets_ptr(PGi);
  E_Int nb_edges = ng.PGs.stride(PGi);
    
  for (int i = 0; i < nb_edges; i++)
  {
    LNODES[i] = pN[i];
    four_childrenPG[i] = *(PGtree.children(PGi)+i); // the four children of PGi
  }

  E_Int* pNFils0 = ng.PGs.get_facets_ptr(four_childrenPG[0]);
  E_Int* pNFils2 = ng.PGs.get_facets_ptr(four_childrenPG[2]);
    

#ifdef DEBUG_HIERARCHICAL_MESH    
  E_Int* pNFils1 = ng.PGs.get_facets_ptr(four_childrenPG[1]);
  E_Int* pNFils3 = ng.PGs.get_facets_ptr(four_childrenPG[3]);

  assert(pNFils0[2] == pNFils2[0]);
  assert(pNFils1[0] == pNFils0[1]);
  assert(pNFils1[2] == pNFils2[1]);
  assert(pNFils2[3] == pNFils3[2]);
  assert(pNFils0[3] == pNFils3[0]);
#endif

    // we got the 4 to 8 thanks to child 0 and child 2 (never swapped)
  LNODES[4] = pNFils0[1];
  LNODES[5] = pNFils2[1];
  LNODES[6] = pNFils2[3];
  LNODES[7] = pNFils0[3];
  LNODES[8] = pNFils0[2]; // centroid
    
    
  K_CONNECT::IdTool::right_shift<4>(&LNODES[0],i0);
  K_CONNECT::IdTool::right_shift<4>(&LNODES[4],i0);     
  K_CONNECT::IdTool::right_shift<4>(&four_childrenPG[0],i0);
    
  if (reorient == true)
  {
    std::reverse(&LNODES[1], &LNODES[1] + 3 );
    std::reverse(&LNODES[4], &LNODES[4] + 4 );
    std::swap(four_childrenPG[1],four_childrenPG[3]);
  }
    
}
template <>
void hierarchical_mesh<K_MESH::Tetrahedron, eSUBDIV_TYPE::ISO, ngon_type, K_FLD::FloatArray>::retrieve_ordered_data(E_Int PGi, E_Int i0, bool reorient, E_Int* four_childrenPG, E_Int* LNODES,
ngo_t& ng, tree<arr_t> & PGtree)
{
  E_Int* pN = ng.PGs.get_facets_ptr(PGi);
  E_Int nb_edges = ng.PGs.stride(PGi);
    
  for (int i = 0; i < nb_edges; i++)
  {
    LNODES[i] = pN[i];
    four_childrenPG[i] = *(PGtree.children(PGi)+i); // the four children of PGi
  }
  four_childrenPG[3] = *(PGtree.children(PGi)+3);
  E_Int* pNFils3 = ng.PGs.get_facets_ptr(four_childrenPG[3]);
    

#ifdef DEBUG_HIERARCHICAL_MESH    
//  E_Int* pNFils1 = _ng.PGs.get_facets_ptr(four_childrenPG[1]);
//  pNFils3 = _ng.PGs.get_facets_ptr(four_childrenPG[3]);

//  //assert(pNFils0[2] == pNFils2[0]);
//  assert(pNFils1[0] == pNFils0[1]);
//  assert(pNFils1[2] == pNFils2[1]);
//  assert(pNFils2[3] == pNFils3[2]);
//  assert(pNFils0[3] == pNFils3[0]);
#endif

    // we got the 4 to 8 thanks to child 0 and child 2 (never swapped)
  LNODES[3] = pNFils3[0];
  LNODES[4] = pNFils3[1];
  LNODES[5] = pNFils3[2];

#ifdef DEBUG_2019 
//  std::cout << "avant shift avec i0= "<< i0 << std::endl;
//  for (int i=0; i<6; i++){
//      std::cout << "LNODES["<< i << "]= " << LNODES[i] << std::endl;
//  }  
//    
#endif
  K_CONNECT::IdTool::right_shift<3>(&LNODES[0],i0);
  K_CONNECT::IdTool::right_shift<3>(&LNODES[3],i0);     
  K_CONNECT::IdTool::right_shift<3>(&four_childrenPG[0],i0);
  
#ifdef DEBUG_2019 
//  std::cout << "après shift" << std::endl;
//  for (int i=0; i<6; i++){
//      std::cout << "LNODES["<< i << "]= " << LNODES[i] << std::endl;
//  }  
//  
//  
//  std::cout << "reorient= " << reorient << std::endl; 
#endif
  if (reorient == true)
  {
    std::swap(LNODES[1], LNODES[2]);
    std::swap(LNODES[3], LNODES[5]);
    std::swap(four_childrenPG[1],four_childrenPG[2]);
  }
}
///
template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t, typename crd_t>
bool hierarchical_mesh<ELT_t, STYPE, ngo_t, crd_t>::need_a_reorient(E_Int PGi, E_Int PHi, bool oriented_if_R,
                                                                    K_FLD::IntArray & F2E)
{
  if (F2E(1,PGi) == PHi && oriented_if_R == true) return false;
  else if (F2E(0,PGi) == PHi && oriented_if_R == false) return false;
  else return true;
}

///
template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t, typename crd_t>
E_Int hierarchical_mesh<ELT_t, STYPE, ngo_t, crd_t>::get_i0(E_Int* pFace, E_Int common_node, E_Int* nodes, E_Int nb_edges_face)
{
  for (int i = 0; i < nb_edges_face; i++)
    if (pFace[i] == nodes[common_node]) return i; 
  return -1;
}

///
template <>
void hierarchical_mesh<K_MESH::Hexahedron, eSUBDIV_TYPE::ISO, ngon_type, K_FLD::FloatArray>::get_nodes_PHi
(E_Int* nodes, E_Int PHi, E_Int centroidId, E_Int* BOT, E_Int* TOP, E_Int* LEFT, E_Int* RIGHT, E_Int* FRONT, E_Int* BACK,
  ngo_t& ng, K_FLD::FloatArray& crd, K_FLD::IntArray & F2E, tree<arr_t> & PGtree)
{
  E_Int* pPGi = ng.PHs.get_facets_ptr(PHi);
  E_Int PGi = pPGi[0] - 1;
  E_Int* pN = ng.PGs.get_facets_ptr(PGi);

  nodes[0] = *pN; // 0 -> PHi(0,0)

  if (F2E(1,PGi) == PHi) // for BOT, PH is the right element : well oriented
    for (int k = 1; k < 4; k++) nodes[k] =*(pN+k);

  else // otherwise : wrong orientation (swap 1 & 3)
  { 
    nodes[1] = *(pN+3);
    nodes[2] = *(pN+2);
    nodes[3] = *(pN+1);
  }
  // let's fill in nodes
  E_Int tmp[9];
  bool reorient;

  // BOT
  E_Int i0 = 0;

  reorient = need_a_reorient(pPGi[0]-1,PHi,true, F2E);
  retrieve_ordered_data(pPGi[0]-1,i0,reorient,BOT,tmp, ng, PGtree);

  nodes[8] = tmp[4];
  nodes[9] = tmp[5];
  nodes[10] = tmp[6];
  nodes[11] = tmp[7];
  nodes[12] = tmp[8];

#ifdef DEBUG_HIERARCHICAL_MESH
  assert(tmp[0] == nodes[0]);
  assert(tmp[1] == nodes[1]);
  assert(tmp[2] == nodes[2]);
  assert(tmp[3] == nodes[3]);
#endif

  // LEFT 
  E_Int* p = ng.PGs.get_facets_ptr(pPGi[2]-1);
  i0 = get_i0(p,0,nodes,4);// common point : nodes[0]

  reorient = need_a_reorient(pPGi[2]-1,PHi,true, F2E);
  retrieve_ordered_data(pPGi[2]-1,i0,reorient,LEFT,tmp, ng, PGtree);

  nodes[21] = tmp[5];
  nodes[7] = tmp[2];
  nodes[16] = tmp[6];
  nodes[4] = tmp[3];
  nodes[18] = tmp[7];
  nodes[25] = tmp[8];

#ifdef DEBUG_HIERARCHICAL_MESH
  assert(tmp[0] == nodes[0]);
  assert(tmp[1] == nodes[3]);
  assert(tmp[4] == nodes[11]);
#endif

  // RIGHT
  p = ng.PGs.get_facets_ptr(pPGi[3]-1);
  i0 = get_i0(p,1,nodes,4);

  reorient = need_a_reorient(pPGi[3]-1,PHi,false, F2E);
  retrieve_ordered_data(pPGi[3]-1,i0,reorient,RIGHT,tmp, ng, PGtree);

  nodes[20] = tmp[5];
  nodes[6] = tmp[2];
  nodes[14] = tmp[6];
  nodes[5] = tmp[3];
  nodes[19] = tmp[7];
  nodes[23] = tmp[8];

#ifdef DEBUG_HIERARCHICAL_MESH
  assert(tmp[0] == nodes[1]);
  assert(tmp[1] == nodes[2]);
  assert(tmp[4] == nodes[9]);
#endif

  // TOP
  p = ng.PGs.get_facets_ptr(pPGi[1]-1);
  i0 = get_i0(p,4,nodes,4);

  reorient = need_a_reorient(pPGi[1]-1,PHi,false, F2E);
  retrieve_ordered_data(pPGi[1]-1,i0,reorient,TOP,tmp, ng, PGtree);

  nodes[13] = tmp[4];
  nodes[15] = tmp[6];
  nodes[17] = tmp[8];

#ifdef DEBUG_HIERARCHICAL_MESH
  assert(tmp[0] == nodes[4]);
  assert(tmp[1] == nodes[5]);
  assert(tmp[2] == nodes[6]);
  assert(tmp[3] == nodes[7]);
  assert(tmp[5] == nodes[14]);
  assert(tmp[7] == nodes[16]);
#endif

  // FRONT
  p = ng.PGs.get_facets_ptr(pPGi[4]-1);
  i0 = get_i0(p,1,nodes,4);

  reorient = need_a_reorient(pPGi[4]-1,PHi,true, F2E);
  retrieve_ordered_data(pPGi[4]-1,i0,reorient,FRONT,tmp, ng, PGtree);

  nodes[22] = tmp[8];

#ifdef DEBUG_HIERARCHICAL_MESH
  assert(tmp[0] == nodes[1]);
  assert(tmp[1] == nodes[0]);
  assert(tmp[2] == nodes[4]);
  assert(tmp[3] == nodes[5]);
  assert(tmp[4] == nodes[8]);
  assert(tmp[5] == nodes[18]);
  assert(tmp[6] == nodes[13]);
  assert(tmp[7] == nodes[19]);
#endif

  // BACK
  p = ng.PGs.get_facets_ptr(pPGi[5]-1);
  i0 = get_i0(p,2,nodes,4);

  reorient = need_a_reorient(pPGi[5]-1,PHi,false, F2E);
  retrieve_ordered_data(pPGi[5]-1,i0,reorient,BACK,tmp, ng, PGtree);

  nodes[24] = tmp[8];

#ifdef DEBUG_HIERARCHICAL_MESH
  assert(tmp[0] == nodes[2]);
  assert(tmp[1] == nodes[3]);
  assert(tmp[2] == nodes[7]);
  assert(tmp[3] == nodes[6]);
  assert(tmp[4] == nodes[10]);
  assert(tmp[5] == nodes[21]);
  assert(tmp[6] == nodes[15]);
  assert(tmp[7] == nodes[20]);
#endif

  // Centroid calculation
  __compute_cell_center(crd, nodes, crd.col(centroidId));
  
  nodes[26] = centroidId+1;

#ifdef DEBUG_HIERARCHICAL_MESH
//        for (int i = 0 ; i < 27; i++){
//            for (int j = 0; j < 27; j++){
//                if (nodes[i] == nodes[j] && i < j )
//                    std::cout << std::endl << " couple " << "      (" << i << " , " << j << ")"; 
//            }
//        }
//        std::cout << std::endl;
#endif
}

template <>
void hierarchical_mesh<K_MESH::Tetrahedron, eSUBDIV_TYPE::ISO, ngon_type, K_FLD::FloatArray>::get_nodes_PHi_T
(E_Int* nodes, E_Int PHi, E_Int centroidId, E_Int* BOT, E_Int* F1, E_Int* F2, E_Int* F3,
ngo_t& ng, K_FLD::FloatArray& crd, K_FLD::IntArray & F2E, tree<arr_t> & PGtree)
{      
  E_Int* pPGi = ng.PHs.get_facets_ptr(PHi);
  E_Int PGi = pPGi[0] - 1;
  E_Int* pN = ng.PGs.get_facets_ptr(PGi);

  nodes[0] = *pN; // 0 -> PHi(0,0)

  if (F2E(1,PGi) == PHi) // for BOT, PH is the right element : well oriented
    for (int k = 1; k < 3; k++) nodes[k] =*(pN+k);

  else // otherwise : wrong orientation (swap 1 & 3)
  { 
    nodes[1] = *(pN+2);
    nodes[2] = *(pN+1);
  }
  // let's fill in nodes
  E_Int tmp[6];
  bool reorient;

  // BOT
  E_Int i0 = 0;

  reorient = need_a_reorient(pPGi[0]-1,PHi,true, F2E);
  retrieve_ordered_data(pPGi[0]-1,i0,reorient,BOT,tmp, ng, PGtree);

  nodes[4] = tmp[3];
  nodes[5] = tmp[4];
  nodes[6] = tmp[5];

#ifdef DEBUG_2019   
//  std::cout << "nodes[0]= " << nodes[0] << std::endl;
//  std::cout << "nodes[1]= " << nodes[1] << std::endl;
//  std::cout << "nodes[2]= " << nodes[2] << std::endl << std::endl;
//  for (int i=0; i<6; i++){
//    std::cout << "tmp["<< i<<"]= " << tmp[i] << std::endl;
//  }
//
#endif
  
#ifdef DEBUG_HIERARCHICAL_MESH
  assert(tmp[0] == nodes[0]);
  assert(tmp[1] == nodes[1]);
  assert(tmp[2] == nodes[2]);
  //assert(tmp[3] == nodes[3]);
#endif

  // F1
  E_Int* p = ng.PGs.get_facets_ptr(pPGi[1]-1);
  i0 = get_i0(p,0,nodes,4);// common point : nodes[0]

  reorient = need_a_reorient(pPGi[1]-1,PHi,false, F2E);
  retrieve_ordered_data(pPGi[1]-1,i0,reorient,F1,tmp, ng, PGtree);

  nodes[7] = tmp[5];
  nodes[8] = tmp[4];
  nodes[3] = tmp[2];
  

#ifdef DEBUG_HIERARCHICAL_MESH
  assert(tmp[0] == nodes[0]);
  assert(tmp[1] == nodes[1]);
  assert(tmp[3] == nodes[4]);
#endif

  // F2
  for (int i=0; i<6; i++){
      tmp[i]=0;
  }
  p = ng.PGs.get_facets_ptr(pPGi[2]-1);
  i0 = get_i0(p,1,nodes,4);

  reorient = need_a_reorient(pPGi[2]-1,PHi,false, F2E);
  retrieve_ordered_data(pPGi[2]-1,i0,reorient,F2,tmp, ng, PGtree);

  nodes[9] = tmp[4];
  

#ifdef DEBUG_HIERARCHICAL_MESH
  assert(tmp[0] == nodes[1]);
  assert(tmp[1] == nodes[2]);
  assert(tmp[2] == nodes[3]);
  assert(tmp[5] == nodes[8]);
  assert(tmp[3] == nodes[5]);
#endif

  // F3
  p = ng.PGs.get_facets_ptr(pPGi[3]-1);
  i0 = get_i0(p,2,nodes,4);

  reorient = need_a_reorient(pPGi[3]-1,PHi,false, F2E);
  retrieve_ordered_data(pPGi[3]-1,i0,reorient,F3,tmp, ng, PGtree);

#ifdef DEBUG_HIERARCHICAL_MESH
  assert(tmp[0] == nodes[2]);
  assert(tmp[1] == nodes[0]);
  assert(tmp[3] == nodes[6]);
  assert(tmp[2] == nodes[3]);
  assert(tmp[4] == nodes[7]);
  assert(tmp[5] == nodes[9]);
#endif

  // Centroid calculation
  //__compute_cell_center(_crd, nodes, _crd.col(centroidId));
  
  //nodes[26] = centroidId+1;

#ifdef DEBUG_HIERARCHICAL_MESH
//        for (int i = 0 ; i < 27; i++){
//            for (int j = 0; j < 27; j++){
//                if (nodes[i] == nodes[j] && i < j )
//                    std::cout << std::endl << " couple " << "      (" << i << " , " << j << ")"; 
//            }
//        }
//        std::cout << std::endl;
#endif
}

///
template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t, typename crd_t>
void hierarchical_mesh<ELT_t, STYPE, ngo_t, crd_t>::update_children_F2E(E_Int PGi, E_Int side, tree<arr_t> & PGtree, K_FLD::IntArray & F2E)
{
  if (PGtree.children(PGi) != nullptr)
  {
    E_Int* p = PGtree.children(PGi);
    E_Int nb_children = PGtree.nb_children(PGi);
    for (int j = 0; j < nb_children; ++j)
      F2E(side,p[j]) = F2E(side,PGi);
  }
}

///
template <>
void hierarchical_mesh<K_MESH::Hexahedron, eSUBDIV_TYPE::ISO, ngon_type, K_FLD::FloatArray>::update_F2E
(E_Int PHi, E_Int PHchildr0, E_Int* INT, E_Int* BOT, E_Int* TOP, E_Int* LEFT, E_Int* RIGHT, E_Int* FRONT, E_Int* BACK, 
ngo_t& ng, K_FLD::IntArray & F2E, tree<arr_t> & PGtree)
{
  E_Int elt0 = PHchildr0;
  E_Int* pPGi = ng.PHs.get_facets_ptr(PHi);
  // BOT and TOP
  for (int i = 0; i < 4; i++)
  {
    E_Int sid = (F2E(1,pPGi[0]-1) == PHi) ? 1 : 0;
    F2E(sid,BOT[i]) = elt0+i;
    update_children_F2E(BOT[i],sid, PGtree, F2E);

    sid = (F2E(0,pPGi[1]-1) == PHi) ? 0 : 1;
    F2E(sid,TOP[i]) = elt0+4+i;
    update_children_F2E(TOP[i],sid, PGtree, F2E);
  }
  // LEFT
  E_Int sid = (F2E(1,pPGi[2]-1) == PHi) ? 1: 0;
  F2E(sid,LEFT[0]) = elt0;
  F2E(sid,LEFT[1]) = elt0+3;
  F2E(sid,LEFT[2]) = elt0+7;
  F2E(sid,LEFT[3]) = elt0+4;
  for (int i = 0; i < 4; ++i){
    update_children_F2E(LEFT[i], sid, PGtree, F2E);
  }
  // RIGHT
  sid = (F2E(0,pPGi[3]-1) == PHi) ? 0 : 1;
  F2E(sid,RIGHT[0]) = elt0+1;
  F2E(sid,RIGHT[1]) = elt0+2;
  F2E(sid,RIGHT[2]) = elt0+6;
  F2E(sid,RIGHT[3]) = elt0+5;
  for (int i = 0; i < 4; ++i){
    update_children_F2E(RIGHT[i],sid, PGtree, F2E);
  }
  // FRONT
  sid = (F2E(1,pPGi[4]-1) == PHi) ? 1: 0;
  F2E(sid,FRONT[0]) = elt0+1;
  F2E(sid,FRONT[1]) = elt0;
  F2E(sid,FRONT[2]) = elt0+4;
  F2E(sid,FRONT[3]) = elt0+5;
  for (int i = 0; i < 4; ++i){
    update_children_F2E(FRONT[i],sid, PGtree, F2E);
  }
  // BACK
  sid = (F2E(0,pPGi[5]-1) == PHi) ? 0 : 1;
  F2E(sid,BACK[0]) = elt0+2;
  F2E(sid,BACK[1]) = elt0+3;
  F2E(sid,BACK[2]) = elt0+7;
  F2E(sid,BACK[3]) = elt0+6;
  for (int i = 0; i < 4; ++i)
    update_children_F2E(BACK[i],sid, PGtree, F2E);
  // INTERNAL faces
  F2E(0,INT[0]) = elt0+1;
  F2E(1,INT[0]) = elt0+2;

  F2E(0,INT[1]) = elt0;
  F2E(1,INT[1]) = elt0+3;

  F2E(0,INT[2]) = elt0+4;
  F2E(1,INT[2]) = elt0+7;

  F2E(0,INT[3]) = elt0+5;
  F2E(1,INT[3]) = elt0+6;

  F2E(0,INT[4]) = elt0;
  F2E(1,INT[4]) = elt0+4;

  F2E(0,INT[5]) = elt0+1;
  F2E(1,INT[5]) = elt0+5;

  F2E(0,INT[6]) = elt0+2;
  F2E(1,INT[6]) = elt0+6;

  F2E(0,INT[7]) = elt0+3;
  F2E(1,INT[7]) = elt0+7;

  F2E(0,INT[8]) = elt0;
  F2E(1,INT[8]) = elt0+1;

  F2E(0,INT[9]) = elt0+3;
  F2E(1,INT[9]) = elt0+2;

  F2E(0,INT[10]) = elt0+7;
  F2E(1,INT[10]) = elt0+6;

  F2E(0,INT[11]) = elt0+4;
  F2E(1,INT[11]) = elt0+5;

}

template <>
void hierarchical_mesh<K_MESH::Tetrahedron, eSUBDIV_TYPE::ISO, ngon_type, K_FLD::FloatArray>::update_F2E_T
(E_Int PHi, E_Int PHchildr0, E_Int* INT, E_Int* BOT, E_Int* F1, E_Int* F2, E_Int* F3, E_Int ndiag,
ngo_t& ng, K_FLD::IntArray & F2E, tree<arr_t> & PGtree)
{
  E_Int elt0 = PHchildr0;
  E_Int* pPGi = ng.PHs.get_facets_ptr(PHi);
  // BOT and TOP
  for (int i = 0; i < 3; i++)
  {
    E_Int sid = (F2E(1,pPGi[0]-1) == PHi) ? 1 : 0;
    F2E(sid,BOT[i]) = elt0+i;
    update_children_F2E(BOT[i],sid, PGtree, F2E);
  }
  E_Int sid = (F2E(1,pPGi[0]-1) == PHi) ? 1 : 0;
  
  
  if (ndiag==1){
  F2E(sid,BOT[3]) = elt0+5;
  update_children_F2E(BOT[3],sid, PGtree, F2E);
  // F1
  sid = (F2E(1,pPGi[1]-1) == PHi) ? 1: 0;
  F2E(sid,F1[0]) = elt0;
  F2E(sid,F1[1]) = elt0+1;
  F2E(sid,F1[2]) = elt0+3;
  F2E(sid,F1[3]) = elt0+4;
  for (int i = 0; i < 4; ++i)
    update_children_F2E(F1[i], sid, PGtree, F2E);
  // F2
  sid = (F2E(0,pPGi[2]-1) == PHi) ? 0 : 1;
  F2E(sid,F2[0]) = elt0+1;
  F2E(sid,F2[1]) = elt0+2;
  F2E(sid,F2[2]) = elt0+3;
  F2E(sid,F2[3]) = elt0+6;
  for (int i = 0; i < 4; ++i)
    update_children_F2E(F2[i],sid, PGtree, F2E);
  // F3
  sid = (F2E(1,pPGi[3]-1) == PHi) ? 1: 0;
  F2E(sid,F3[0]) = elt0+2;
  F2E(sid,F3[1]) = elt0;
  F2E(sid,F3[2]) = elt0+3;
  F2E(sid,F3[3]) = elt0+7;
  for (int i = 0; i < 4; ++i)
    update_children_F2E(F3[i],sid, PGtree, F2E);
  // INTERNAL faces
  F2E(0,INT[0]) = elt0;
  F2E(1,INT[0]) = elt0+4;

  F2E(0,INT[1]) = elt0+1;
  F2E(1,INT[1]) = elt0+5;

  F2E(0,INT[2]) = elt0+2;
  F2E(1,INT[2]) = elt0+6;

  F2E(0,INT[3]) = elt0+3;
  F2E(1,INT[3]) = elt0+7;

  F2E(0,INT[4]) = elt0+5;
  F2E(1,INT[4]) = elt0+6;

  F2E(0,INT[5]) = elt0+4;
  F2E(1,INT[5]) = elt0+7;

  F2E(0,INT[6]) = elt0+4;
  F2E(1,INT[6]) = elt0+5;

  F2E(0,INT[7]) = elt0+7;
  F2E(1,INT[7]) = elt0+6;
  }
  else if (ndiag==2) {
  F2E(sid,BOT[3]) = elt0+4;
  update_children_F2E(BOT[3],sid, PGtree, F2E);
  // F1
  sid = (F2E(1,pPGi[1]-1) == PHi) ? 1: 0;
  F2E(sid,F1[0]) = elt0;
  F2E(sid,F1[1]) = elt0+1;
  F2E(sid,F1[2]) = elt0+3;
  F2E(sid,F1[3]) = elt0+5;
  for (int i = 0; i < 4; ++i)
    update_children_F2E(F1[i], sid, PGtree, F2E);
  // F2
  sid = (F2E(0,pPGi[2]-1) == PHi) ? 0 : 1;
  F2E(sid,F2[0]) = elt0+1;
  F2E(sid,F2[1]) = elt0+2;
  F2E(sid,F2[2]) = elt0+3;
  F2E(sid,F2[3]) = elt0+7;
  for (int i = 0; i < 4; ++i)
    update_children_F2E(F2[i],sid, PGtree, F2E);
  // F3
  sid = (F2E(1,pPGi[3]-1) == PHi) ? 1: 0;
  F2E(sid,F3[0]) = elt0+2;
  F2E(sid,F3[1]) = elt0;
  F2E(sid,F3[2]) = elt0+3;
  F2E(sid,F3[3]) = elt0+6;
  for (int i = 0; i < 4; ++i)
    update_children_F2E(F3[i],sid, PGtree, F2E);
  // INTERNAL faces
  F2E(0,INT[0]) = elt0;
  F2E(1,INT[0]) = elt0+4;

  F2E(0,INT[1]) = elt0+1;
  F2E(1,INT[1]) = elt0+5;

  F2E(0,INT[2]) = elt0+2;
  F2E(1,INT[2]) = elt0+6;

  F2E(0,INT[3]) = elt0+3;
  F2E(1,INT[3]) = elt0+7;

  F2E(0,INT[4]) = elt0+5;
  F2E(1,INT[4]) = elt0+7;

  F2E(0,INT[5]) = elt0+4;
  F2E(1,INT[5]) = elt0+6;

  F2E(0,INT[6]) = elt0+5;
  F2E(1,INT[6]) = elt0+4;

  F2E(0,INT[7]) = elt0+7;
  F2E(1,INT[7]) = elt0+6;
  }
  else {
  F2E(sid,BOT[3]) = elt0+6;
  update_children_F2E(BOT[3],sid, PGtree, F2E);
  // F1
  sid = (F2E(1,pPGi[1]-1) == PHi) ? 1: 0;
  F2E(sid,F1[0]) = elt0;
  F2E(sid,F1[1]) = elt0+1;
  F2E(sid,F1[2]) = elt0+3;
  F2E(sid,F1[3]) = elt0+7;
  for (int i = 0; i < 4; ++i)
    update_children_F2E(F1[i], sid, PGtree, F2E);
  // F2
  sid = (F2E(0,pPGi[2]-1) == PHi) ? 0 : 1;
  F2E(sid,F2[0]) = elt0+1;
  F2E(sid,F2[1]) = elt0+2;
  F2E(sid,F2[2]) = elt0+3;
  F2E(sid,F2[3]) = elt0+5;
  for (int i = 0; i < 4; ++i)
    update_children_F2E(F2[i],sid, PGtree, F2E);
  // F3
  sid = (F2E(1,pPGi[3]-1) == PHi) ? 1: 0;
  F2E(sid,F3[0]) = elt0+2;
  F2E(sid,F3[1]) = elt0;
  F2E(sid,F3[2]) = elt0+3;
  F2E(sid,F3[3]) = elt0+4;
  for (int i = 0; i < 4; ++i)
    update_children_F2E(F3[i],sid, PGtree, F2E);
  // INTERNAL faces
  F2E(0,INT[0]) = elt0;
  F2E(1,INT[0]) = elt0+4;

  F2E(0,INT[1]) = elt0+1;
  F2E(1,INT[1]) = elt0+5;

  F2E(0,INT[2]) = elt0+2;
  F2E(1,INT[2]) = elt0+6;

  F2E(0,INT[3]) = elt0+3;
  F2E(1,INT[3]) = elt0+7;

  F2E(0,INT[4]) = elt0+6;
  F2E(1,INT[4]) = elt0+4;

  F2E(0,INT[5]) = elt0+5;
  F2E(1,INT[5]) = elt0+7;

  F2E(0,INT[6]) = elt0+6;
  F2E(1,INT[6]) = elt0+5;

  F2E(0,INT[7]) = elt0+4;
  F2E(1,INT[7]) = elt0+7;
  }
}

///
template <>
void hierarchical_mesh<K_MESH::Hexahedron, eSUBDIV_TYPE::ISO, ngon_type, K_FLD::FloatArray>::refine_PHs
(const Vector_t<E_Int> &PHadap, 
 ngo_t& ng, tree<arr_t> & PGtree, tree<arr_t> & PHtree, K_FLD::FloatArray& crd, K_FLD::IntArray & F2E)
{
  E_Int nb_phs = PHadap.size();
  // internal PGs created (12 per PH)
  // Reserve space for internal faces in the tree
  E_Int current_sz = PGtree.get_parent_size();

  PGtree.resize_hierarchy(current_sz+12*nb_phs);

  E_Int nb_pgs0 = ng.PGs.size();

  // And in the mesh : each H27 is split in 36 Pgs including 12 new internal Q4

  ng.PGs.expand_n_fixed_stride(12*nb_phs, 4/*Q4 stride*/);

  F2E.resize(2,nb_pgs0+12*nb_phs,E_IDX_NONE);

  // for the children PH
  // Reserve space for children in the tree
  PHtree.resize(PHadap);
  // And in the mesh : each H27 is split into 8 H27
  E_Int nb_phs0 = ng.PHs.size();
  ng.PHs.expand_n_fixed_stride(8*nb_phs, 6/*H27 stride*/);

  E_Int pos = crd.cols();
  crd.resize(3, pos + nb_phs);
  
  int j;
#ifndef DEBUG_HIERARCHICAL_MESH
#pragma omp parallel for private (j)
#endif
  for (E_Int i = 0; i < nb_phs; ++i)
  {
    E_Int PHi = PHadap[i];
    E_Int nodes[27]; // 0-26 to crd1 indexes
    E_Int BOT[4],TOP[4],LEFT[4],RIGHT[4],FRONT[4],BACK[4];// children lists des faces

    get_nodes_PHi(nodes, PHi, pos+i, BOT,TOP,LEFT,RIGHT,FRONT,BACK,
                  ng, crd, F2E, PGtree);

    E_Int points[27];
    K_MESH::Hexahedron::get_internal(nodes,points);

    E_Int PGichildr[4];
    E_Int INT[12];
    // INTERNAL faces
    for ( j = 0; j < 3; ++j)
    {
      PGichildr[0] = nb_pgs0 + 12*i + 4*j;
      PGichildr[1] = nb_pgs0 + 12*i + 4*j + 1;
      PGichildr[2] = nb_pgs0 + 12*i + 4*j + 2;
      PGichildr[3] = nb_pgs0 + 12*i + 4*j + 3;

      E_Int* q41 = ng.PGs.get_facets_ptr(PGichildr[0]);
      E_Int* q42 = ng.PGs.get_facets_ptr(PGichildr[1]);
      E_Int* q43 = ng.PGs.get_facets_ptr(PGichildr[2]);
      E_Int* q44 = ng.PGs.get_facets_ptr(PGichildr[3]);

      NUGA::Q9::splitQ4(crd, &points[9*j], q41, q42, q43, q44);

      INT[4*j]   = PGichildr[0];
      INT[4*j+1] = PGichildr[1];
      INT[4*j+2] = PGichildr[2];
      INT[4*j+3] = PGichildr[3];
    }

    // the 8 children of PH
    E_Int PHichildr[8];

    for (int j = 0; j < 8; ++j)
      PHichildr[j] = nb_phs0 + 8*i + j;
    
    E_Int* h271 = ng.PHs.get_facets_ptr(PHichildr[0]);
    E_Int* h272 = ng.PHs.get_facets_ptr(PHichildr[1]);
    E_Int* h273 = ng.PHs.get_facets_ptr(PHichildr[2]);
    E_Int* h274 = ng.PHs.get_facets_ptr(PHichildr[3]);        
    E_Int* h275 = ng.PHs.get_facets_ptr(PHichildr[4]);
    E_Int* h276 = ng.PHs.get_facets_ptr(PHichildr[5]);
    E_Int* h277 = ng.PHs.get_facets_ptr(PHichildr[6]);
    E_Int* h278 = ng.PHs.get_facets_ptr(PHichildr[7]);

    NUGA::H27::splitH27(crd, INT, BOT, TOP, LEFT, RIGHT, FRONT, BACK, h271, h272, h273, h274, h275, h276, h277, h278);

    // set them in the tree
    PHtree.set_children(PHi, PHichildr, 8);

    // update F2E
    update_F2E(PHi,PHichildr[0],INT,BOT,TOP,LEFT,RIGHT,FRONT,BACK, ng, F2E, PGtree);
  }  
}

///
template <>
void hierarchical_mesh<K_MESH::Tetrahedron, eSUBDIV_TYPE::ISO, ngon_type, K_FLD::FloatArray>::refine_PHs
(const Vector_t<E_Int> &PHadap, 
ngo_t& ng, tree<arr_t> & PGtree, tree<arr_t> & PHtree, K_FLD::FloatArray& crd, K_FLD::IntArray & F2E)
{
  E_Int nb_phs = PHadap.size();

  // internal PGs created (12 per PH)
  // Reserve space for internal faces in the tree
  E_Int current_sz = PGtree.get_parent_size();

  PGtree.resize_hierarchy(current_sz+8*nb_phs);

  E_Int nb_pgs0 = ng.PGs.size();

  // And in the mesh : each H27 is split in 36 Pgs including 12 new internal Q4

  //_ng.PGs.expand_n_fixed_stride(8*nb_phs, 4/*Q4 stride*/);
  ng.PGs.expand_n_fixed_stride(8*nb_phs, 3/*Q4 stride*/);

  F2E.resize(2,nb_pgs0+8*nb_phs,E_IDX_NONE);

  // for the children PH
  // Reserve space for children in the tree
  PHtree.resize(PHadap);
  // And in the mesh : each H27 is split into 8 H27
  E_Int nb_phs0 = ng.PHs.size();
  ng.PHs.expand_n_fixed_stride(8*nb_phs, 4/*TH4 stride*/);

  E_Int pos = crd.cols();
  //_crd.resize(3, pos + nb_phs);
  
  int j;
#ifndef DEBUG_HIERARCHICAL_MESH
#pragma omp parallel for private (j)
#endif
  
#ifdef DEBUG_2019 
  //std::cout << "nb_phs= " << nb_phs << std::endl;
  //E_Int n1(0), n2(0), n3(0);
#endif
  
  for (E_Int i = 0; i < nb_phs; ++i)
  {
    E_Int PHi = PHadap[i];
    E_Int nodes[10]; // 0-26 to crd1 indexes
    E_Int BOT[4], F1[4], F2[4], F3[4];// children lists

    get_nodes_PHi_T(nodes, PHi, pos+i, BOT, F1, F2, F3,
                    ng, crd, F2E, PGtree);
    
    /// Choix de la diagonale ici param nodes[10]
//    E_Float d;
       
    E_Float d1=K_FUNC::sqrDistance(crd.col(nodes[8]-1),crd.col(nodes[6]-1), 3);
    E_Float d2=K_FUNC::sqrDistance(crd.col(nodes[7]-1),crd.col(nodes[5]-1), 3);
    E_Float d3=K_FUNC::sqrDistance(crd.col(nodes[9]-1),crd.col(nodes[4]-1), 3);

    E_Int ndiag= ((d1 <= d2) && (d1 <= d3) ) ? 1 : ((d2 <= d1) && (d2 <= d3) ) ? 2 : 3;

#ifdef DEBUG_2019     
//    d=K_FUNC::sqrDistance(_crd.col(nodes[5]-1),_crd.col(nodes[9]-1), 3)
//
//    std::cout << "d1= " << d1 << std::endl;
//    std::cout << "d2= " << d2 << std::endl;
//    std::cout << "d3= " << d3 << std::endl;
//    std::cout << "d= " << d << std::endl;
    //ndiag=3;
//    if (ndiag==1){
//        n1 += 1;
//    }
//    else if (ndiag==2) {
//        n2 += 1;
//    }
//    else {
//        n3 += 1;
//    }
#endif
    
    E_Int INT[8];
    // INTERNAL faces
    
    for ( j = 0; j < 8; ++j)
    {
        INT[j]= nb_pgs0 + 8*i + j; 
    }
    
    E_Int* q41 = ng.PGs.get_facets_ptr(INT[0]);
    E_Int* q42 = ng.PGs.get_facets_ptr(INT[1]);
    E_Int* q43 = ng.PGs.get_facets_ptr(INT[2]);
    E_Int* q44 = ng.PGs.get_facets_ptr(INT[3]);
    E_Int* q45 = ng.PGs.get_facets_ptr(INT[4]);
    E_Int* q46 = ng.PGs.get_facets_ptr(INT[5]);
    E_Int* q47 = ng.PGs.get_facets_ptr(INT[6]);
    E_Int* q48 = ng.PGs.get_facets_ptr(INT[7]);

    q41[0]=nodes[4]; q41[1]=nodes[6]; q41[2]=nodes[7];
    q42[0]=nodes[5]; q42[1]=nodes[4]; q42[2]=nodes[8];
    q43[0]=nodes[6]; q43[1]=nodes[5]; q43[2]=nodes[9];
    q44[0]=nodes[8]; q44[1]=nodes[7]; q44[2]=nodes[9];
    
    if (ndiag==1){
    q45[0]=nodes[6]; q45[1]=nodes[8]; q45[2]=nodes[5];
    q46[0]=nodes[7]; q46[1]=nodes[8]; q46[2]=nodes[6];
    q47[0]=nodes[4]; q47[1]=nodes[6]; q47[2]=nodes[8];
    q48[0]=nodes[8]; q48[1]=nodes[6]; q48[2]=nodes[9];
    }
    else if (ndiag==2){
    q45[0]=nodes[5]; q45[1]=nodes[7]; q45[2]=nodes[8];
    q46[0]=nodes[6]; q46[1]=nodes[7]; q46[2]=nodes[5];
    q47[0]=nodes[4]; q47[1]=nodes[7]; q47[2]=nodes[5];
    q48[0]=nodes[5]; q48[1]=nodes[7]; q48[2]=nodes[9];
    }
    else {
    q45[0]=nodes[4]; q45[1]=nodes[9]; q45[2]=nodes[6];
    q46[0]=nodes[8]; q46[1]=nodes[9]; q46[2]=nodes[4];
    q47[0]=nodes[5]; q47[1]=nodes[9]; q47[2]=nodes[4];
    q48[0]=nodes[4]; q48[1]=nodes[9]; q48[2]=nodes[7];    
    }


    // the 8 children of PH
    E_Int PHichildr[8];

    for (int j = 0; j < 8; ++j)
      PHichildr[j] = nb_phs0 + 8*i + j;
    
    E_Int* h271 = ng.PHs.get_facets_ptr(PHichildr[0]);
    E_Int* h272 = ng.PHs.get_facets_ptr(PHichildr[1]);
    E_Int* h273 = ng.PHs.get_facets_ptr(PHichildr[2]);
    E_Int* h274 = ng.PHs.get_facets_ptr(PHichildr[3]);
    E_Int* h275 = ng.PHs.get_facets_ptr(PHichildr[4]);
    E_Int* h276 = ng.PHs.get_facets_ptr(PHichildr[5]);
    E_Int* h277 = ng.PHs.get_facets_ptr(PHichildr[6]);
    E_Int* h278 = ng.PHs.get_facets_ptr(PHichildr[7]);

    NUGA::H27::splitT10(crd, INT, BOT, F1, F2, F3, h271, h272, h273, h274, h275, h276, h277, h278, ndiag);

    // set them in the tree
    PHtree.set_children(PHi, PHichildr, 8);

    // update F2E
    update_F2E_T(PHi,PHichildr[0],INT,BOT, F1, F2, F3, ndiag, ng, F2E, PGtree);

#ifdef DEBUG_2019    
//  E_Float Vol;
//  E_Int* face;
//  //std::vector<std::pair<E_Float,E_Int> > palma;
//  for (int j=0; j<8; j++){
//      K_MESH::Tetrahedron T;
//      //std::cout << "Tetra " << i << "    "; 
//      face= _ng.PHs.get_facets_ptr(PHichildr[j]);
//      T.init(_ng.PGs,face);
//      T.Hm(_crd);
//      E_Float q=T.qualityTet(_crd, &Vol);
//      //palma.insert(std::make_pair(q,j));
//      std::cout << "q[" << j << "]= " << q << std::endl;
//  }
//  //std::sort(palma.begin(),palma.end());
//  for (int j=0; j<8; j++){
//      //std::cout << palma[j].first << "  ";
//      //std::cout << palma[j].second << std::endl;
//  }
#endif
  }
//  std::cout << "n1= " << n1 << "   n2= " << n2 << "   n3= "<< n3 << std::endl;
  //std::cout << "nb_phs= " << nb_phs << std::endl;
}

template <>
void hierarchical_mesh<K_MESH::Basic, eSUBDIV_TYPE::ISO, ngon_type, K_FLD::FloatArray>::refine_PHs
(const Vector_t<E_Int> &PHadap, 
 ngo_t& ng, tree<arr_t> & PGtree, tree<arr_t> & PHtree, K_FLD::FloatArray& crd, K_FLD::IntArray & F2E)
{
    Vector_t<E_Int> PHadap1, PHadap2;
    E_Int nb_phs = PHadap.size();
    for (E_Int i=0; i< nb_phs; i++){
        E_Int PHi=PHadap[i];
        E_Int s= ng.PHs.stride(PHi);
        if (K_MESH::Polyhedron<0>::is_HX8(ng.PGs, ng.PHs.get_facets_ptr(PHi), s)){
            PHadap1.push_back(PHi);
        }
        else if (K_MESH::Polyhedron<0>::is_TH4(ng.PGs, ng.PHs.get_facets_ptr(PHi), s)){
            PHadap2.push_back(PHi);
        }
    }
    //std::cout << "ok" << std::endl;
//    std::cout << "0)" << std::endl;
//    std::cout << "ng.PG size= " << ng.PGs.size() <<std::endl;
//    std::cout << "ng.PH size= " << ng.PHs.size() <<std::endl;
//    std::cout << "PHtree size_children= " << PHtree.size_children() <<std::endl;
//    std::cout << "PHtree size_parent= " << PHtree.size_parent() <<std::endl;
    //NUGA::array_trait<K_FLD::IntArray>::size(ng.PHs)
//    std::cout << "Hexa" << std::endl;
//    std::cout << "PHadap1 size= " << PHadap1.size() << std::endl;
    hierarchical_mesh<K_MESH::Hexahedron, eSUBDIV_TYPE::ISO, ngon_type, K_FLD::FloatArray>::refine_PHs(PHadap1, 
                                                                                                        ng, PGtree, PHtree, crd, F2E);
    
//    if (is_hybrid(ng))
//      std::cout << "hybrid" << std::endl;
//    else
//      std::cout << "false" << std::endl;
    
    //std::cout << "ok" << std::endl;
//    std::cout << "PHtree size_children= " << PHtree.size_children() <<std::endl;
//    std::cout << "PHtree size_parent= " << PHtree.size_parent() <<std::endl;
//    std::cout << "Tetra" << std::endl;
//    std::cout << "PHadap2 size= " << PHadap2.size() << std::endl;
    hierarchical_mesh<K_MESH::Tetrahedron, eSUBDIV_TYPE::ISO, ngon_type, K_FLD::FloatArray>::refine_PHs(PHadap2, 
                                                                                                        ng, PGtree, PHtree, crd, F2E);
    
//    if (is_hybrid(ng))
//      std::cout << "hybrid" << std::endl;
//    else
//      std::cout << "false" << std::endl;
    
//    std::cout << "PHtree size_children= " << PHtree.size_children() <<std::endl;
//    std::cout << "PHtree size_parent= " << PHtree.size_parent() <<std::endl;
//    std::cout << "1)" << std::endl;
//    std::cout << "ng.PG size= " << ng.PGs.size() <<std::endl;
//    std::cout << "ng.PH size= " << ng.PHs.size() <<std::endl;
}
///
template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t, typename crd_t>
void hierarchical_mesh<ELT_t, STYPE, ngo_t, crd_t>::__compute_edge_center
(const Vector_t<E_Int> &PGlist, 
ngo_t& ng, crd_t& crd,         
std::map<K_MESH::NO_Edge,E_Int> & ecenter
)
{
  for (int i = 0; i < PGlist.size(); i++)
  {
    E_Int PGi = PGlist[i];

    E_Int nb_nodes = ng.PGs.stride(PGi);
    E_Int* nodes = ng.PGs.get_facets_ptr(PGlist[i]);
                
    for (E_Int j = 0; j < nb_nodes; j++)
    {

      E_Int ind_point1 = *(nodes+j);
      E_Int ind_point2 = *(nodes+(j+1)%nb_nodes);

      K_MESH::NO_Edge no_edge(ind_point1,ind_point2); // not oriented (a,b) = (b,a)

      auto it = ecenter.find(no_edge);
      if (it == ecenter.end())
      {
        E_Float mid[3];
        K_FUNC::sum<3>(0.5, crd.col(ind_point1-1), 0.5, crd.col(ind_point2-1), mid);
        crd.pushBack(mid, mid+3);
        ecenter[no_edge] = crd.cols();
      }   
    }
  }
}

/////
//template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t, typename crd_t>
//void hierarchical_mesh<ELT_t, STYPE, ngo_t, crd_t>::__compute_face_centers
//(K_FLD::FloatArray& crd, const typename ngo_t::unit_type& PGs, const Vector_t<E_Int> &PGlist, Vector_t<E_Int>& fcenter)
//{
//  E_Int nb_pgs = PGlist.size();
//  
//  fcenter.resize(nb_pgs, E_IDX_NONE);
//  
//  E_Int pos = crd.cols();
//  
//  crd.resize(3, pos + nb_pgs);
//  
//  for (size_t i = 0; i < nb_pgs; ++i)
//  {
//    E_Int PGi = PGlist[i];
//    // Centroid calculation
//    __compute_face_center(crd, PGs.get_facets_ptr(PGi), PGs.stride(PGi), crd.col(pos));
//    
//    fcenter[i] = pos++;
//  }
//}

///
template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t, typename crd_t>
void hierarchical_mesh<ELT_t, STYPE, ngo_t, crd_t>::__compute_face_center
(const crd_t& crd, const E_Int* nodes, E_Int nb_nodes, E_Float* C)
{   
  K_MESH::Polygon::iso_barycenter<crd_t,3>(crd, nodes, nb_nodes, 1, C);
}

///
template <typename ELT_t, eSUBDIV_TYPE STYPE, typename ngo_t, typename crd_t>
void hierarchical_mesh<ELT_t, STYPE, ngo_t, crd_t>::get_cell_center(E_Int PHi, E_Float* center)
{    

  ELT_t::iso_barycenter(_crd, _ng.PGs, _ng.PHs.get_facets_ptr(PHi), _ng.PHs.stride(PHi), 1, center);
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
  E_Int* p = _ng.PHs.get_facets_ptr(PHi);
  E_Int nb_faces = _ng.PHs.stride(PHi);
    
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
  E_Int* p = _ng.PHs.get_facets_ptr(PHi);
  E_Int nb_faces = _ng.PHs.stride(PHi);
    
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
void hierarchical_mesh<ELT_t, STYPE, ngo_t, crd_t>::smooth(std::vector<E_Int>& adap_incr)
{
  std::stack<E_Int> stck;

  for (E_Int i = 0; i < _ng.PHs.size();  i++){
    if (adap_incr[i] != 0){
      stck.push(i);
    }
  }
  
  while (!stck.empty()){

    E_Int ind_PHi = stck.top(); // index of ith PH
    stck.pop();
    E_Int s = _ng.PHs.stride(ind_PHi);

    
    E_Int* neighbours= new E_Int[4*s];
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
    
  for (int i=0; i<_ng.PHs.size(); i++)
  {
    if (_PHtree.is_enabled(i) && _ng.PHs.stride(i)==4)
    {
      ++nb_T;
      face= _ng.PHs.get_facets_ptr(i);

      ELT_t e(_ng.PGs,face);
      
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
  std::cout << "_ng.PHs size= " << _ng.PHs.size() << std::endl;
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
  std::cout << "_ng.PHs size= " << _ng.PHs.size() << std::endl;
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
    E_Int nb_phs = _ng.PHs.size();
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
                if (K_MESH::Polyhedron<0>::is_HX8(_ng.PGs, _ng.PHs.get_facets_ptr(i), _ng.PHs.stride(i)) ) ++s1;
                else if (K_MESH::Polyhedron<0>::is_TH4(_ng.PGs, _ng.PHs.get_facets_ptr(i), _ng.PHs.stride(i)) ) ++s2;
                else if (K_MESH::Polyhedron<0>::is_PY5(_ng.PGs, _ng.PHs.get_facets_ptr(i), _ng.PHs.stride(i)) ) ++s3;
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
    E_Int err = K_MESH::Polyhedron<UNKNOWN>::metrics<DELAUNAY::Triangulator>(dt, _crd, _ng.PGs, _ng.PHs.get_facets_ptr(k), _ng.PHs.stride(k), vol_father, Gdum);
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
        E_Int err = K_MESH::Polyhedron<UNKNOWN>::metrics<DELAUNAY::Triangulator>(dt, _crd, _ng.PGs, _ng.PHs.get_facets_ptr(child_j), _ng.PHs.stride(child_j), vol_child, Gdum);
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
            E_Int err = K_MESH::Polyhedron<UNKNOWN>::metrics<DELAUNAY::Triangulator>(dt, _crd, _ng.PGs, _ng.PHs.get_facets_ptr(child_j), _ng.PHs.stride(child_j), vol_child, Gdum);
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
  E_Int nb_elts = _ng.PHs.size();

  for (int i=0; i< nb_elts; i++){
      E_Int level= _PHtree.get_level(i);
      E_Int nb_faces= _ng.PHs.stride(i);
      E_Int* pPGi = _ng.PHs.get_facets_ptr(i);
      if (K_MESH::Polyhedron<0>::is_PY5(_ng.PGs, _ng.PHs.get_facets_ptr(i), _ng.PHs.stride(i))){
           continue;
       }     
      for (int j=0; j< nb_faces; j++){
        E_Int PGj = *(pPGi +j) - 1;
        E_Int PHn = NEIGHBOR(i, _F2E, PGj);
        if (PHn==E_IDX_NONE || K_MESH::Polyhedron<0>::is_PY5(_ng.PGs, _ng.PHs.get_facets_ptr(PHn), _ng.PHs.stride(PHn))){
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
  E_Int nb_elts = _ng.PHs.size();

  for (int i=0; i< nb_elts; i++){
    //E_Int level= _hmesh._PHtree.get_level(i);
    //E_Int nb_faces= _hmesh._ng.PHs.stride(i);
    E_Int* pPGi = _ng.PHs.get_facets_ptr(i);
    if (K_MESH::Polyhedron<0>::is_PY5(_ng.PGs, _ng.PHs.get_facets_ptr(i), _ng.PHs.stride(i)))
      continue;
 
    if (K_MESH::Polyhedron<0>::is_HX8(_ng.PGs, _ng.PHs.get_facets_ptr(i), _ng.PHs.stride(i))){
      for (int j=0; j< 3; j++){
        E_Int PGj = *(pPGi +2*j) - 1;
        E_Int PHn = NEIGHBOR(i, _F2E, PGj);
        if (PHn==E_IDX_NONE || K_MESH::Polyhedron<0>::is_PY5(_ng.PGs, _ng.PHs.get_facets_ptr(PHn), _ng.PHs.stride(PHn)))
          continue;
        
        E_Int leveln= _PHtree.get_level(PHn);
        for (int k=0; k< 4; k++){
          E_Int PGk= *(pPGi +(j+k+2) % 6) - 1;
          E_Int PHk = NEIGHBOR(i, _F2E, PGk);
          if (PHk==E_IDX_NONE || K_MESH::Polyhedron<0>::is_PY5(_ng.PGs, _ng.PHs.get_facets_ptr(PHk), _ng.PHs.stride(PHk)))
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
        if (PHn==E_IDX_NONE || K_MESH::Polyhedron<0>::is_PY5(_ng.PGs, _ng.PHs.get_facets_ptr(PHn), _ng.PHs.stride(PHn)))
          continue;
        
        E_Int leveln= _PHtree.get_level(PHn);
        for (int k=0; k< 4; k++){
          E_Int PGk= *(pPGi +(j+k+1) % 6) - 1;
          E_Int PHk = NEIGHBOR(i, _F2E, PGk);
          if (PHk==E_IDX_NONE || K_MESH::Polyhedron<0>::is_PY5(_ng.PGs, _ng.PHs.get_facets_ptr(PHk), _ng.PHs.stride(PHk)))
            continue;      
 
          E_Int levelk= _PHtree.get_level(PHk);
          if (::fabs(leveln-levelk)>1){
            err += 1;
          }
        }
      }
    }
    else if (K_MESH::Polyhedron<0>::is_TH4(_ng.PGs, _ng.PHs.get_facets_ptr(i), _ng.PHs.stride(i))){
      for (int j=0; j< 4; j++){
        E_Int PGj = *(pPGi +j) - 1;
        E_Int PHn = NEIGHBOR(i, _F2E, PGj);
        if (PHn==E_IDX_NONE || K_MESH::Polyhedron<0>::is_PY5(_ng.PGs, _ng.PHs.get_facets_ptr(PHn), _ng.PHs.stride(PHn)))
          continue;      
        
        E_Int leveln= _PHtree.get_level(PHn);
        for (int k=0; k< 3; k++){
          E_Int PGk= *(pPGi +(j+k+1) % 4) - 1;
          E_Int PHk = NEIGHBOR(i, _F2E, PGk);
          if (PHk==E_IDX_NONE || K_MESH::Polyhedron<0>::is_PY5(_ng.PGs, _ng.PHs.get_facets_ptr(PHk), _ng.PHs.stride(PHk)))
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
#endif

}

#endif

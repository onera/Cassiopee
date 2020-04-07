/*
 
 
 
              NUGA 
 
 
 
 */
//Authors : Sâm Landier (sam.landier@onera.fr)

#ifndef REFINER_HXX
#define REFINER_HXX

#include "MeshElement/Edge.h"
#include "MeshElement/Basic.h"
#include "T6.hxx"
#include "Q9.hxx"
#include "h27.hxx"
#include "sensor.hxx"

namespace NUGA
{

// general case : nb of varying children
template <typename ELT_t, eSUBDIV_TYPE STYPE>
struct NBC
{
  enum { nbc = -1 };
  using arr_t = ngon_unit;
};

//
template <>
struct NBC<K_MESH::Hexahedron, ISO>
{
  enum { nbc = 8 };
  using arr_t = K_FLD::IntArray;
  static constexpr E_Int nbi = 12; // nb of inner faces
};

//
template <>
struct NBC<K_MESH::Tetrahedron, ISO>
{
  enum { nbc = 8 };
  using arr_t = K_FLD::IntArray;
};

// isotropic HEXA subdivision => 4 Quadrangles children => fixed stride array
template <>
struct NBC<K_MESH::Quadrangle, ISO>
{
  enum { nbc = 4 };
  using arr_t = K_FLD::IntArray;
};

// isotropic Triangle subdivision => 4 Trianges children => fixed stride array
template <>
struct NBC<K_MESH::Triangle, ISO>
{
  enum { nbc = 4 };
  using arr_t = K_FLD::IntArray;
};

///
template <typename ELT_t>
class refine_point_computer
{
  public:
    static void compute_center(const K_FLD::FloatArray& crd, const E_Int* nodes, E_Int nb_nodes, E_Int idx_start, E_Float* C);
};
  
///
template <>
class refine_point_computer<K_MESH::NO_Edge>
{
public:
  static void compute_centers(const Vector_t<E_Int> &PGlist, 
                              ngon_type& ng, K_FLD::FloatArray& crd, 
                              std::map<K_MESH::NO_Edge,E_Int> & ecenter)
  {
    for (size_t i = 0; i < PGlist.size(); i++)
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
};

///
template <typename ELT_t>
inline
void refine_point_computer<ELT_t>::compute_center
(const K_FLD::FloatArray& crd, const E_Int* nodes, E_Int nb_nodes, E_Int idx_start, E_Float* C)
{   
  ELT_t::iso_barycenter(crd, nodes, nb_nodes, idx_start, C);
}
  
 ///
template <typename ELT_t, eSUBDIV_TYPE STYPE>
class refiner
{
  public : 
    using sensor_output_t = typename sensor_output_data<STYPE>::type;

  public:
    std::map<K_MESH::NO_Edge, E_Int> _ecenter;

  public:
    
    // Generic wrapper for cells (calls relevant(s) refine_PGs)
    template <typename arr_t>
    static void refine_Faces(const sensor_output_t &adap_incr, 
                           ngon_type& ng, tree<arr_t> & PGtree, K_FLD::FloatArray& crd, K_FLD::IntArray & F2E,         
                           std::map<K_MESH::NO_Edge,E_Int>& ecenter);

    /// Generic wrapper for cells
    template <typename arr_t>
    static void refine_PHs(const sensor_output_t &adap_incr, ngon_type& ng, tree<arr_t> & PGtree, tree<arr_t> & PHtree,
      K_FLD::FloatArray& crd, K_FLD::IntArray & F2E);

    // T3/Q4/PG
    template <typename arr_t>
    static void refine_PGs(const sensor_output_t &adap_incr,
      ngon_type& ng, tree<arr_t> & PGtree, K_FLD::FloatArray& crd, K_FLD::IntArray & F2E,
      std::map<K_MESH::NO_Edge, E_Int>& ecenter);

    static void refine_PG(K_FLD::FloatArray& crd, ngon_unit& PGs, E_Int PGi, E_Int posC, E_Int posChild,
                          std::map<K_MESH::NO_Edge, E_Int>& ecenter);

    
      
    //private:
    
    ///
    template <typename arr_t>
    static void get_PGs_to_refine(const ngon_type& ng, const tree<arr_t> & PGtree, const sensor_output_t &adap_incr, Vector_t<E_Int> & PGlist);
    ///
    ///
    template <typename arr_t>
    static void get_PHs_to_refine(const ngon_type& ng, const tree<arr_t> & PHtree, const sensor_output_t &adap_incr, Vector_t<E_Int> & PHlist);
    ///
    template <typename arr_t>
    static void reserve_mem_PGs(K_FLD::FloatArray& crd, ngon_unit& PGs, const Vector_t<E_Int>& PGlist, tree<arr_t> & PGtree, K_FLD::IntArray& F2E);
    ///
    template <typename arr_t>
    static void reserve_mem_PHs(K_FLD::FloatArray& crd, ngon_type& ng, const Vector_t<E_Int>& PHlist, tree<arr_t> & PGtree, tree<arr_t> & PHtree, K_FLD::IntArray& F2E);
    ///    
    template <typename arr_t>
    static void retrieve_ordered_data(E_Int PGi, E_Int i0, bool reorient, E_Int* four_childrenPG, E_Int* LNODES, ngon_type& ng, tree<arr_t> & PGtree);
        
    ///
    template <typename arr_t>
    static void get_nodes_PHi (E_Int* nodes, E_Int PHi, E_Int centroidId, E_Int** FACES,
                          ngon_type& ng, K_FLD::FloatArray& crd, K_FLD::IntArray & F2E, tree<arr_t> & PGtree);
  
    ///
    template <typename arr_t>
    static void update_F2E (E_Int PHi, E_Int *children, E_Int* INT, E_Int** FACES, E_Int ndiag, ngon_type& ng, K_FLD::IntArray & F2E, tree<arr_t> & PGtree);
        
    ///
    template <typename arr_t>
    static void update_children_F2E(E_Int PGi, E_Int side, tree<arr_t> & PGtree, K_FLD::IntArray & F2E);

};

///
template <typename ELT_t, eSUBDIV_TYPE STYPE>
template <typename arr_t>
void refiner<ELT_t, STYPE>::get_PGs_to_refine
(const ngon_type& ng, const tree<arr_t> & PGtree, const sensor_output_t &adap_incr, Vector_t<E_Int> & PGlist)
{
  E_Int nb_phs = ng.PHs.size();
  // Gets PGs to refine
  E_Int nb_pgs(ng.PGs.size()), nb_pgs_ref(0);
  Vector_t<bool> is_PG_to_refine(nb_pgs, false);
  //
  for (E_Int i = 0; i < nb_phs; ++i)
  {
    if (adap_incr[i] <= 0) continue;

    E_Int nb_faces = ng.PHs.stride(i);
    const E_Int* faces = ng.PHs.get_facets_ptr(i);

    for (E_Int j = 0; j < nb_faces; ++j)
    {
      E_Int PGi = *(faces + j) - 1;

      if (PGtree.nb_children(PGi) == 0) // leaf PG => to refine
      {
        if (ng.PGs.stride(PGi) == ELT_t::NB_NODES)
        {
          ++nb_pgs_ref;
          is_PG_to_refine[PGi] = true;
        }
      }
    }
  }

  PGlist.reserve(nb_pgs_ref);
  for (E_Int i = 0; i < nb_pgs; ++i)
  {
    if (is_PG_to_refine[i]) {
      PGlist.push_back(i);
    }
  }
}


///
template <typename ELT_t, eSUBDIV_TYPE STYPE>
template <typename arr_t>
void refiner<ELT_t, STYPE>::get_PHs_to_refine
(const ngon_type& ng, const tree<arr_t> & PHtree, const sensor_output_t &adap_incr, Vector_t<E_Int> & PHlist)
{
  for (size_t i = 0; i < adap_incr.size(); ++i)
    if (adap_incr[i] > 0 && PHtree.is_enabled(i) && (ELT_t::isOfType(ng.PGs, ng.PHs.get_facets_ptr(i), ng.PHs.stride(i))))
      PHlist.push_back(i);
}


/// Default impl for 
template <typename ELT_t, eSUBDIV_TYPE STYPE>
template <typename arr_t>
void refiner<ELT_t, STYPE>::reserve_mem_PGs
(K_FLD::FloatArray& crd, ngon_unit& PGs, const Vector_t<E_Int>& PGlist, tree<arr_t> & PGtree, K_FLD::IntArray& F2E)
{
  const E_Int &nbc = NBC<ELT_t, STYPE>::nbc;
  E_Int nb_pgs_ref = PGlist.size();

  PGtree.resize(PGlist, nbc);
  
  E_Int nb_pgs0 = PGs.size();

  // And in the mesh : each ELT(T3/) is split into nbc
  PGs.expand_n_fixed_stride(nbc * nb_pgs_ref, ELT_t::NB_NODES);

  F2E.resize(2, nb_pgs0 + nbc * nb_pgs_ref, E_IDX_NONE);
}

/// Impl for PG with centroid point as refinement
template <>
template <typename arr_t>
void refiner<K_MESH::Quadrangle, ISO>::reserve_mem_PGs
(K_FLD::FloatArray& crd, ngon_unit& PGs, const Vector_t<E_Int>& PGlist, tree<arr_t> & PGtree, K_FLD::IntArray& F2E)
{
  const E_Int &nbc = NBC<K_MESH::Quadrangle, ISO>::nbc;
  E_Int nb_pgs_ref = PGlist.size();

  PGtree.resize(PGlist, nbc);
  
  E_Int nb_pgs0 = PGs.size();

  // And in the mesh : each ELT(T3/) is split into nbc
  PGs.expand_n_fixed_stride(nbc * nb_pgs_ref, K_MESH::Quadrangle::NB_NODES);

  F2E.resize(2, nb_pgs0 + nbc * nb_pgs_ref, E_IDX_NONE); 

  // reserve space for face centers
  crd.resize(3, crd.cols() + nb_pgs_ref);
}

// Default impl for mono-bound-type element (TETRA, HEXA)
template <typename ELT_t, eSUBDIV_TYPE STYPE>
template <typename arr_t>
void refiner<ELT_t, STYPE>::reserve_mem_PHs
(K_FLD::FloatArray& crd, ngon_type& ng, const Vector_t<E_Int>& PHadap, tree<arr_t> & PGtree, tree<arr_t> & PHtree, K_FLD::IntArray& F2E)
{
  const E_Int &phnbc = NBC<ELT_t, STYPE>::nbc;
  const E_Int& nbi = NBC<ELT_t, STYPE>::nbi;

  using face_t = typename ELT_t::boundary_type;

  E_Int nb_phs_ref = PHadap.size();

  // internal PGs created (12 per PH)
  // Reserve space for internal faces in the tree
  E_Int current_sz = PGtree.size();
  assert(current_sz != -1);

  PGtree.resize_hierarchy(current_sz + nbi * nb_phs_ref);

  E_Int nb_pgs0 = ng.PGs.size();

  // each Hex/ISO is split in 36 PGs including 12 new internal Q4
  ng.PGs.expand_n_fixed_stride(nbi * nb_phs_ref, face_t::NB_NODES);

  F2E.resize(2, nb_pgs0 + nbi * nb_phs_ref, E_IDX_NONE);

  // for the children PH
  // Reserve space for children in the tree
  PHtree.resize(PHadap, phnbc);
  // each Hex/ISO is split into 8 Hex
  ng.PHs.expand_n_fixed_stride(phnbc * nb_phs_ref, ELT_t::NB_BOUNDS/*H27 stride*/);

  E_Int pos = crd.cols();
  crd.resize(3, pos + nb_phs_ref);
}

/// default implementation : T3/Q4 ISO case
template <typename ELT_t, eSUBDIV_TYPE STYPE>
template <typename arr_t>
void refiner<ELT_t, STYPE>::refine_PGs
(const sensor_output_t& adap_incr,
 ngon_type& ng, tree<arr_t>& PGtree, K_FLD::FloatArray& crd, K_FLD::IntArray& F2E,
 std::map<K_MESH::NO_Edge, E_Int>& ecenter)
{
  std::vector<E_Int> PGref;
  get_PGs_to_refine(ng, PGtree, adap_incr, PGref);
  E_Int nb_pgs_ref = PGref.size();
  if (nb_pgs_ref == 0) return;

  // Compute Edges refine points
  refine_point_computer<K_MESH::NO_Edge>::compute_centers(PGref, ng, crd, ecenter);

  E_Int pos = crd.cols(); // first face center in crd (if relevant)
  E_Int nb_pgs0 = ng.PGs.size();
  // Reserve space for children in the tree
  reserve_mem_PGs(crd, ng.PGs, PGref, PGtree, F2E);

  const E_Int &nbc = NBC<ELT_t, STYPE>::nbc;

  // Refine them
#ifndef DEBUG_HIERARCHICAL_MESH
#pragma omp parallel for
#endif
  for (E_Int i = 0; i < nb_pgs_ref; ++i)
  {
    E_Int PGi = PGref[i];
    E_Int firstChild = nb_pgs0 + nbc * i;

    refine_PG(crd, ng.PGs, PGi, pos + i, firstChild, ecenter);

    for (E_Int n = 0; n < nbc; ++n)
    {
      // children have the L & R elements of the father
      F2E(0, firstChild + n) = F2E(0, PGi);
      F2E(1, firstChild + n) = F2E(1, PGi);
    }

    // set them in the tree
    PGtree.set_children(PGi, firstChild, nbc);
  }
}

template<>
void refiner<K_MESH::Quadrangle, ISO>::refine_PG
(K_FLD::FloatArray& crd, ngon_unit& PGs, E_Int PGi, E_Int posC, E_Int posChild, std::map<K_MESH::NO_Edge, E_Int>& ecenter)
{
  E_Int* nodes = PGs.get_facets_ptr(PGi);
  const E_Int &nbc = NBC<K_MESH::Quadrangle, ISO>::nbc;

  // Centroid calculation
  NUGA::refine_point_computer<K_MESH::Polygon>::compute_center(crd, nodes, K_MESH::Quadrangle::NB_NODES, 1, crd.col(posC));

  E_Int q9[9];
  K_MESH::NO_Edge noE;
  for (E_Int n = 0; n < K_MESH::Quadrangle::NB_NODES; ++n)
  {
    q9[n] = *(nodes + n);
    noE.setNodes(*(nodes + n), *(nodes + (n + 1) % 4));
    q9[n + K_MESH::Quadrangle::NB_NODES] = ecenter[noE];
  }
  q9[8] = posC + 1;

  NUGA::Q9::split<eSUBDIV_TYPE::ISO>(PGs, q9, posChild);
}

///
template<>
void refiner<K_MESH::Triangle, ISO>::refine_PG
(K_FLD::FloatArray& crd, ngon_unit& PGs, E_Int PGi, E_Int posC, E_Int posChild, std::map<K_MESH::NO_Edge, E_Int>& ecenter)
{
  E_Int* nodes = PGs.get_facets_ptr(PGi);
  const E_Int &nbc = NBC<K_MESH::Triangle, ISO>::nbc;

  E_Int T6[6];
  K_MESH::NO_Edge noE;
  for (E_Int n = 0; n < K_MESH::Triangle::NB_NODES; ++n)
  {
    T6[n] = *(nodes + n);
    noE.setNodes(*(nodes + n), *(nodes + (n + 1) % K_MESH::Triangle::NB_NODES));
    T6[n + K_MESH::Triangle::NB_NODES] = ecenter[noE];
  }

  NUGA::T6::split<eSUBDIV_TYPE::ISO>(PGs, T6, posChild);
}

///
template <>
template <typename arr_t>
void refiner<K_MESH::Quadrangle, eSUBDIV_TYPE::DIR>::refine_PGs
(const sensor_output_t& adap_incr, ngon_type& ng, tree<arr_t>& PGtree, K_FLD::FloatArray& crd, 
 K_FLD::IntArray& F2E, std::map<K_MESH::NO_Edge,E_Int>& ecenter)
{
  //alexis : todo:
  
  //extraire d'abord la partii ISO de adap_incr et creer un vector<E_Int> avec cette partie, appler la version ISO de refine_Faces dessus
  
  // faire le reste en XY
}

///
template<>
template <typename arr_t>
void refiner<K_MESH::Triangle, DIR>::refine_PGs
(const dir_incr_type &adap_incr,
 ngon_type& ng, tree<arr_t> & PGtree, K_FLD::FloatArray& crd, K_FLD::IntArray & F2E,         
 std::map<K_MESH::NO_Edge,E_Int>& ecenter)
{
  // do not handle DIR for Triangle yet
  refiner<K_MESH::Triangle, ISO>::refine_PGs<arr_t>(adap_incr._adap_incr, ng, PGtree, crd, F2E, ecenter);
}

// default implementation : ISO case , for all basic element types (Basic, Tetra, Pyra, Penta, Hexa)
template <typename ELT_t, eSUBDIV_TYPE STYPE>
template <typename arr_t>
void refiner<ELT_t, STYPE>::refine_Faces
(const sensor_output_t &adap_incr, 
 ngon_type& ng, tree<arr_t> & PGtree, K_FLD::FloatArray& crd, K_FLD::IntArray & F2E,                 
 std::map<K_MESH::NO_Edge,E_Int>& ecenter)
{ 
  refiner<K_MESH::Quadrangle, STYPE>::refine_PGs(adap_incr, ng, PGtree, crd, F2E, ecenter);
  refiner<K_MESH::Triangle, STYPE>::refine_PGs(adap_incr, ng, PGtree, crd, F2E, ecenter);  
}

///
template <>
template <typename arr_t>
void refiner<K_MESH::Quadrangle, eSUBDIV_TYPE::ISO>::retrieve_ordered_data
(E_Int PGi, E_Int i0, bool reorient, E_Int* four_childrenPG, E_Int* LNODES,
 ngon_type& ng, tree<arr_t> & PGtree)
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
template <typename arr_t>
void refiner<K_MESH::Triangle, eSUBDIV_TYPE::ISO>::retrieve_ordered_data
(E_Int PGi, E_Int i0, bool reorient, E_Int* four_childrenPG, E_Int* LNODES,
ngon_type& ng, tree<arr_t> & PGtree)
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
  E_Int* pNFils1 = ng.PGs.get_facets_ptr(four_childrenPG[1]);
  E_Int* pNFils2 = ng.PGs.get_facets_ptr(four_childrenPG[2]);
  E_Int* pNFils0 = ng.PGs.get_facets_ptr(four_childrenPG[0]);

  assert(pNFils0[2] == pNFils2[0]);
  assert(pNFils1[0] == pNFils0[1]);
  assert(pNFils1[2] == pNFils2[1]);
  assert(pNFils2[0] == pNFils3[2]);
  assert(pNFils0[2] == pNFils3[2]);
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
inline bool need_a_reorient(E_Int PGi, E_Int PHi, bool oriented_if_R, K_FLD::IntArray & F2E)
{
  if (F2E(1,PGi) == PHi && oriented_if_R == true) return false;
  else if (F2E(0,PGi) == PHi && oriented_if_R == false) return false;
  else return true;
}

///
inline E_Int get_i0(E_Int* pFace, E_Int common_node, E_Int* nodes, E_Int nb_edges_face)
{
  for (int i = 0; i < nb_edges_face; i++)
    if (pFace[i] == nodes[common_node]) return i; 
#ifdef DEBUG_HIERARCHICAL_MESH
  assert(false);
#endif
  return -1;
}

///
template <>
template <typename arr_t>
void refiner<K_MESH::Hexahedron, eSUBDIV_TYPE::ISO>::get_nodes_PHi
(E_Int* nodes, E_Int PHi, E_Int centroidId, E_Int** FACES,
  ngon_type& ng, K_FLD::FloatArray& crd, K_FLD::IntArray & F2E, tree<arr_t> & PGtree)
{
  E_Int* BOT = &FACES[0][0];  
  E_Int* TOP = &FACES[1][0]; 
  E_Int* LEFT = &FACES[2][0];  
  E_Int* RIGHT = &FACES[3][0];  
  E_Int* FRONT = &FACES[4][0];  
  E_Int* BACK = &FACES[5][0];   

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
  refiner<K_MESH::Quadrangle, eSUBDIV_TYPE::ISO>::retrieve_ordered_data(pPGi[0]-1,i0,reorient,BOT,tmp, ng, PGtree);

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
  refiner<K_MESH::Quadrangle, eSUBDIV_TYPE::ISO>::retrieve_ordered_data(pPGi[2]-1,i0,reorient,LEFT,tmp, ng, PGtree);

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
  refiner<K_MESH::Quadrangle, eSUBDIV_TYPE::ISO>::retrieve_ordered_data(pPGi[3]-1,i0,reorient,RIGHT,tmp, ng, PGtree);

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
  refiner<K_MESH::Quadrangle, eSUBDIV_TYPE::ISO>::retrieve_ordered_data(pPGi[1]-1,i0,reorient,TOP,tmp, ng, PGtree);

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
  refiner<K_MESH::Quadrangle, eSUBDIV_TYPE::ISO>::retrieve_ordered_data(pPGi[4]-1,i0,reorient,FRONT,tmp, ng, PGtree);

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
  refiner<K_MESH::Quadrangle, eSUBDIV_TYPE::ISO>::retrieve_ordered_data(pPGi[5]-1,i0,reorient,BACK,tmp, ng, PGtree);

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
  NUGA::refine_point_computer<K_MESH::Hexahedron>::compute_center(crd, nodes, 8, 1, crd.col(centroidId));
  
  nodes[26] = centroidId+1;
}

///
template <>
template <typename arr_t>
void refiner<K_MESH::Tetrahedron, eSUBDIV_TYPE::ISO>::get_nodes_PHi
(E_Int* nodes, E_Int PHi, E_Int centroidId, E_Int** FACES,
ngon_type& ng, K_FLD::FloatArray& crd, K_FLD::IntArray & F2E, tree<arr_t> & PGtree)
{      
  E_Int* BOT = &FACES[0][0];  
  E_Int* F1 = &FACES[1][0]; 
  E_Int* F2 = &FACES[2][0];  
  E_Int* F3 = &FACES[3][0];  
    
  E_Int* pPGi = ng.PHs.get_facets_ptr(PHi);
  E_Int PGi = pPGi[0] - 1;
  E_Int* pN = ng.PGs.get_facets_ptr(PGi);

  nodes[0] = *pN; // 0 -> PHi(0,0)

//  std::cout << "orient= " << F2E(1,PGi) << std::endl;
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
  refiner<K_MESH::Triangle, eSUBDIV_TYPE::ISO>::retrieve_ordered_data(pPGi[0]-1,i0,reorient,BOT,tmp, ng, PGtree);

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
  refiner<K_MESH::Triangle, eSUBDIV_TYPE::ISO>::retrieve_ordered_data(pPGi[1]-1,i0,reorient,F1,tmp, ng, PGtree);

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
  refiner<K_MESH::Triangle, eSUBDIV_TYPE::ISO>::retrieve_ordered_data(pPGi[2]-1,i0,reorient,F2,tmp, ng, PGtree);

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
  refiner<K_MESH::Triangle, eSUBDIV_TYPE::ISO>::retrieve_ordered_data(pPGi[3]-1,i0,reorient,F3,tmp, ng, PGtree);

#ifdef DEBUG_HIERARCHICAL_MESH
  assert(tmp[0] == nodes[2]);
  assert(tmp[1] == nodes[0]);
  assert(tmp[3] == nodes[6]);
  assert(tmp[2] == nodes[3]);
  assert(tmp[4] == nodes[7]);
  assert(tmp[5] == nodes[9]);
#endif
}

template <>
template <typename arr_t>
void refiner<K_MESH::Pyramid, eSUBDIV_TYPE::ISO>::get_nodes_PHi
(E_Int* nodes, E_Int PHi, E_Int centroidId, E_Int** FACES,
ngon_type& ng, K_FLD::FloatArray& crd, K_FLD::IntArray & F2E, tree<arr_t> & PGtree)
{      
  E_Int* BOT = &FACES[0][0];  
  E_Int* F1 = &FACES[1][0]; 
  E_Int* F2 = &FACES[2][0];  
  E_Int* F3 = &FACES[3][0];  
  E_Int* F4 = &FACES[4][0];  
    
  E_Int* pPGi = ng.PHs.get_facets_ptr(PHi);
  E_Int PGi = pPGi[0] - 1;
  E_Int* pN = ng.PGs.get_facets_ptr(PGi);

//  for (int i=0; i< ng.PGs.stride(PGi); i++){
//      std::cout << "pN= " << *(pN+i) << std::endl;
//  }
  nodes[0] = *pN; // 0 -> PHi(0,0)

  if (F2E(1,PGi) == PHi){ // for BOT, PH is the right element : well oriented
    for (int k = 1; k < 4; k++) nodes[k] =*(pN+k);
  }
  else // otherwise : wrong orientation (swap 1 & 3)
  { 
    nodes[1] = *(pN+3);
    nodes[2] = *(pN+2);
    nodes[3] = *(pN+1);
  }
//  for (int i=0; i< ng.PGs.stride(PGi); i++){
//      std::cout << "nodes[" << i << "]= " << nodes[i] << std::endl;
//  }
  
  // let's fill in nodes
  E_Int tmp[9];
  bool reorient;

  // BOT
  E_Int i0 = 0;

  reorient = need_a_reorient(pPGi[0]-1,PHi,true, F2E);
  refiner<K_MESH::Quadrangle, eSUBDIV_TYPE::ISO>::retrieve_ordered_data(pPGi[0]-1,i0,reorient,BOT,tmp, ng, PGtree);

  nodes[5] = tmp[4];
  nodes[6] = tmp[5];
  nodes[7] = tmp[6];
  nodes[8] = tmp[7];
  nodes[9] = tmp[8];

#ifdef DEBUG_2019   
//  std::cout << "nodes[0]= " << nodes[0] << std::endl;
//  std::cout << "nodes[1]= " << nodes[1] << std::endl;
//  std::cout << "nodes[2]= " << nodes[2] << std::endl << std::endl;
//  std::cout << "face BOT stride= " << ng.PGs.stride(PGi) << std::endl;
//  for (int i=0; i<4; i++){
//    std::cout << "nodes["<< i<<"]= " << nodes[i] << std::endl;
//  }
//  for (int i=0; i<9; i++){
//    std::cout << "tmp["<< i<<"]= " << tmp[i] << std::endl;
//  }
//
#endif
  
#ifdef DEBUG_HIERARCHICAL_MESH
  assert(tmp[0] == nodes[0]);
  assert(tmp[1] == nodes[1]);
  assert(tmp[2] == nodes[2]);
  assert(tmp[3] == nodes[3]);
#endif
  E_Int tmp2[6];

  // F1
  E_Int* p = ng.PGs.get_facets_ptr(pPGi[1]-1);
  i0 = get_i0(p,0,nodes,4);// common point : nodes[0]

  reorient = need_a_reorient(pPGi[1]-1,PHi,false, F2E);
  refiner<K_MESH::Triangle, eSUBDIV_TYPE::ISO>::retrieve_ordered_data(pPGi[1]-1,i0,reorient,F1,tmp2, ng, PGtree);

  nodes[4] = tmp2[2];
  nodes[10] = tmp2[4];
  nodes[13] = tmp2[5];
  

#ifdef DEBUG_HIERARCHICAL_MESH
  assert(tmp2[0] == nodes[0]);
  assert(tmp2[1] == nodes[1]);
  assert(tmp2[3] == nodes[5]);
#endif

  // F2
  for (int i=0; i<6; i++){
      tmp2[i]=0;
  }
  p = ng.PGs.get_facets_ptr(pPGi[2]-1);
  i0 = get_i0(p,1,nodes,4);

  reorient = need_a_reorient(pPGi[2]-1,PHi,false, F2E);
  refiner<K_MESH::Triangle, eSUBDIV_TYPE::ISO>::retrieve_ordered_data(pPGi[2]-1,i0,reorient,F2,tmp2, ng, PGtree);  

#ifdef DEBUG_HIERARCHICAL_MESH
  assert(tmp2[0] == nodes[1]);
  assert(tmp2[1] == nodes[2]);
  assert(tmp2[2] == nodes[4]);
  assert(tmp2[3] == nodes[6]);
  assert(tmp2[5] == nodes[10]);
#endif

  // F3
  p = ng.PGs.get_facets_ptr(pPGi[3]-1);
  i0 = get_i0(p,2,nodes,4);

  reorient = need_a_reorient(pPGi[3]-1,PHi,false, F2E);
  refiner<K_MESH::Triangle, eSUBDIV_TYPE::ISO>::retrieve_ordered_data(pPGi[3]-1,i0,reorient,F3,tmp2, ng, PGtree);
  
  nodes[11] = tmp2[5];
  nodes[12] = tmp2[4];
  
#ifdef DEBUG_HIERARCHICAL_MESH
  assert(tmp2[0] == nodes[2]);
  assert(tmp2[1] == nodes[3]);
  assert(tmp2[2] == nodes[4]);
  assert(tmp2[3] == nodes[7]);
#endif
  
  // F4
  p = ng.PGs.get_facets_ptr(pPGi[4]-1);
  i0 = get_i0(p,3,nodes,4);

  reorient = need_a_reorient(pPGi[4]-1,PHi,false, F2E);
  refiner<K_MESH::Triangle, eSUBDIV_TYPE::ISO>::retrieve_ordered_data(pPGi[4]-1,i0,reorient,F4,tmp2, ng, PGtree);

#ifdef DEBUG_HIERARCHICAL_MESH
  assert(tmp2[0] == nodes[3]);
  assert(tmp2[1] == nodes[0]);
  assert(tmp2[2] == nodes[4]);
  assert(tmp2[3] == nodes[8]);
  assert(tmp2[4] == nodes[13]);
  assert(tmp2[5] == nodes[12]);
#endif

  // Centroid calculation
  //__compute_cell_center(_crd, nodes, _crd.col(centroidId));
  
  //nodes[26] = centroidId+1;
}

template <>
template <typename arr_t>
void refiner<K_MESH::Prism, eSUBDIV_TYPE::ISO>::get_nodes_PHi
(E_Int* nodes, E_Int PHi, E_Int centroidId, E_Int** FACES,
ngon_type& ng, K_FLD::FloatArray& crd, K_FLD::IntArray & F2E, tree<arr_t> & PGtree)
{      
  E_Int* BOT = &FACES[0][0];  
  E_Int* F1 = &FACES[1][0]; 
  E_Int* F2 = &FACES[2][0];  
  E_Int* F3 = &FACES[3][0];  
  E_Int* TOP = &FACES[4][0];      
    
  E_Int* pPGi = ng.PHs.get_facets_ptr(PHi);
  E_Int PGi = pPGi[0] - 1;
  E_Int* pN = ng.PGs.get_facets_ptr(PGi);

  for (size_t i=0; i < 18; ++i)
  	nodes[i] = E_IDX_NONE;

  //sam
  assert (ng.PGs.stride(PGi) == 3);

#ifdef DEBUG_2019  
//  std::cout << "PHi= " << PHi << std::endl;
//  for (int j=0; j< ng.PHs.stride(PHi); j++){
//      E_Int PGj= pPGi[j]-1;
//    for (int i=0; i< ng.PGs.stride(PGj); i++){
//        E_Int* pNj= ng.PGs.get_facets_ptr(PGj);
//        std::cout << "pN= " << *(pNj+i)-1 << std::endl;
//    }
//      std::cout << "/////////" << std::endl;
//  }
#endif
  
  nodes[0] = *pN; // 0 -> PHi(0,0)

  if (F2E(1,PGi) == PHi){ // for BOT, PH is the right element : well oriented
    for (int k = 1; k < 3; k++) nodes[k] =*(pN+k);
  }
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
  refiner<K_MESH::Triangle, eSUBDIV_TYPE::ISO>::retrieve_ordered_data(pPGi[0]-1,i0,reorient,BOT,tmp, ng, PGtree);

  nodes[6] = tmp[3];
  nodes[7] = tmp[4];
  nodes[8] = tmp[5];


#ifdef DEBUG_2019   
//  std::cout << "nodes[0]= " << nodes[0] << std::endl;
//  std::cout << "nodes[1]= " << nodes[1] << std::endl;
//  std::cout << "nodes[2]= " << nodes[2] << std::endl << std::endl;
//  std::cout << "face BOT stride= " << ng.PGs.stride(PGi) << std::endl;
//  for (int i=0; i<6; i++){
//    std::cout << "nodes["<< i<<"]= " << nodes[i] << std::endl;
//  }
//  for (int i=0; i<6; i++){
//    std::cout << "tmp["<< i<<"]= " << tmp[i] << std::endl;
//  }
//
#endif
  
#ifdef DEBUG_HIERARCHICAL_MESH
  assert(tmp[0] == nodes[0]);
  assert(tmp[1] == nodes[1]);
  assert(tmp[2] == nodes[2]);
#endif
  
#ifdef DEBUG_2019  
//    for (int i=0; i<4; i++){
//      E_Int* pNj= ng.PGs.get_facets_ptr(BOT[i]);
//      for (int j=0; j< ng.PGs.stride(BOT[i]); j++){
//          std::cout << "BOT[" << i <<"]= " << *(pNj+j)-1 << std::endl;
//      }
//      std::cout << "/////" << std::endl;
//  } 
//
//  for (int i=0; i<3; i++){
//      std::cout << "nodes[" << i << "]= " << nodes[i]-1 << std::endl;
//  }
//  for (int i=0; i<3; i++){
//      std::cout << "nodes[" << i+6 << "]= " << nodes[i+6]-1 << std::endl;
//  }
//  std::cout << "////////" << std::endl;
#endif
  
  E_Int tmp2[9];
  

  // F1
  E_Int* p = ng.PGs.get_facets_ptr(pPGi[1]-1);
  assert (ng.PGs.stride(pPGi[1]-1) == 4);
  i0 = get_i0(p,0,nodes,4);// common point : nodes[0]

  reorient = need_a_reorient(pPGi[1]-1,PHi,false, F2E);
  refiner<K_MESH::Quadrangle, eSUBDIV_TYPE::ISO>::retrieve_ordered_data(pPGi[1]-1,i0,reorient,F1,tmp2, ng, PGtree);

  nodes[15] = tmp2[7];
  nodes[12] = tmp2[8];
  nodes[16] = tmp2[5];
  nodes[3] = tmp2[3];
  nodes[9] = tmp2[6];
  nodes[4] = tmp2[2];
  
#ifdef DEBUG_HIERARCHICAL_MESH
  assert(tmp2[0] == nodes[0]);
  assert(tmp2[4] == nodes[6]);
  assert(tmp2[1] == nodes[1]);
#endif

  // F2
  for (int i=0; i<9; i++){
      tmp2[i]=E_IDX_NONE;
  }
  p = ng.PGs.get_facets_ptr(pPGi[2]-1);
  assert (ng.PGs.stride(pPGi[2]-1) == 4);//sam
  i0 = get_i0(p,1,nodes,4);

  reorient = need_a_reorient(pPGi[2]-1,PHi,false, F2E);
  refiner<K_MESH::Quadrangle, eSUBDIV_TYPE::ISO>::retrieve_ordered_data(pPGi[2]-1,i0,reorient,F2,tmp2, ng, PGtree);  

  nodes[13] = tmp2[8];
  nodes[17] = tmp2[5];
  nodes[10] = tmp2[6];
  nodes[5] = tmp2[2];
  
#ifdef DEBUG_HIERARCHICAL_MESH
  assert(tmp2[0] == nodes[1]);
  assert(tmp2[4] == nodes[7]);
  assert(tmp2[1] == nodes[2]);
  assert(tmp2[7] == nodes[16]);
  assert(tmp2[3] == nodes[4]);
#endif

  // F3
  for (int i=0; i<9; i++){
      tmp2[i]=E_IDX_NONE;
  }
  p = ng.PGs.get_facets_ptr(pPGi[3]-1);
  i0 = get_i0(p,2,nodes,4);

  reorient = need_a_reorient(pPGi[3]-1,PHi,false, F2E);
  refiner<K_MESH::Quadrangle, eSUBDIV_TYPE::ISO>::retrieve_ordered_data(pPGi[3]-1,i0,reorient,F3,tmp2, ng, PGtree);
  
  nodes[14] = tmp2[8];
  nodes[11] = tmp2[6];
  
#ifdef DEBUG_HIERARCHICAL_MESH
  assert(tmp2[0] == nodes[2]);
  assert(tmp2[4] == nodes[8]);
  assert(tmp2[1] == nodes[0]);
  assert(tmp2[7] == nodes[17]);
  assert(tmp2[5] == nodes[15]);
  assert(tmp2[3] == nodes[5]);
  assert(tmp2[2] == nodes[3]);
#endif
  
  // TOP
  p = ng.PGs.get_facets_ptr(pPGi[4]-1);
  assert (ng.PGs.stride(pPGi[4]-1) == 3);//sam
  i0 = get_i0(p,3,nodes,3);

  reorient = need_a_reorient(pPGi[4]-1,PHi,false, F2E);
  refiner<K_MESH::Triangle, eSUBDIV_TYPE::ISO>::retrieve_ordered_data(pPGi[4]-1,i0,reorient,TOP,tmp, ng, PGtree);

#ifdef DEBUG_HIERARCHICAL_MESH
  assert(tmp[0] == nodes[3]);
  assert(tmp[1] == nodes[4]);
  assert(tmp[2] == nodes[5]);
  assert(tmp[3] == nodes[9]);
  assert(tmp[4] == nodes[10]);
  assert(tmp[5] == nodes[11]);

  for (size_t i=0; i < 18; ++i)
  	assert(nodes[i] != E_IDX_NONE);
#endif
}

///
template <typename ELT_t, eSUBDIV_TYPE STYPE>
template <typename arr_t>
void refiner<ELT_t, STYPE>::update_children_F2E(E_Int PGi, E_Int side, tree<arr_t> & PGtree, K_FLD::IntArray & F2E)
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
template <typename arr_t>
void refiner<K_MESH::Hexahedron, eSUBDIV_TYPE::ISO>::update_F2E
(E_Int PHi, E_Int *children, E_Int* INT, E_Int** FACES, E_Int ndiag,
ngon_type& ng, K_FLD::IntArray & F2E, tree<arr_t>& PGtree)
{
  E_Int* BOT = &FACES[0][0];  
  E_Int* TOP = &FACES[1][0]; 
  E_Int* LEFT = &FACES[2][0];  
  E_Int* RIGHT = &FACES[3][0];  
  E_Int* FRONT = &FACES[4][0];  
  E_Int* BACK = &FACES[5][0];      
    
  E_Int* pPGi = ng.PHs.get_facets_ptr(PHi);
  // BOT and TOP
  for (int i = 0; i < 4; i++)
  {
    E_Int sid = (F2E(1,pPGi[0]-1) == PHi) ? 1 : 0;
    F2E(sid,BOT[i]) = *(children+i);
    update_children_F2E(BOT[i],sid, PGtree, F2E);

    sid = (F2E(0,pPGi[1]-1) == PHi) ? 0 : 1;
    F2E(sid,TOP[i]) = *(children+4+i);
    update_children_F2E(TOP[i],sid, PGtree, F2E);
  }
  // LEFT
  E_Int sid = (F2E(1,pPGi[2]-1) == PHi) ? 1: 0;
  F2E(sid,LEFT[0]) = *(children);
  F2E(sid,LEFT[1]) = *(children+3);
  F2E(sid,LEFT[2]) = *(children+7);
  F2E(sid,LEFT[3]) = *(children+4);
  for (int i = 0; i < 4; ++i){
    update_children_F2E(LEFT[i], sid, PGtree, F2E);
  }
  // RIGHT
  sid = (F2E(0,pPGi[3]-1) == PHi) ? 0 : 1;
  F2E(sid,RIGHT[0]) = *(children+1);
  F2E(sid,RIGHT[1]) = *(children+2);
  F2E(sid,RIGHT[2]) = *(children+6);
  F2E(sid,RIGHT[3]) = *(children+5);
  for (int i = 0; i < 4; ++i){
    update_children_F2E(RIGHT[i],sid, PGtree, F2E);
  }
  // FRONT
  sid = (F2E(1,pPGi[4]-1) == PHi) ? 1: 0;
  F2E(sid,FRONT[0]) = *(children+1);
  F2E(sid,FRONT[1]) = *(children);
  F2E(sid,FRONT[2]) = *(children+4);
  F2E(sid,FRONT[3]) = *(children+5);
  for (int i = 0; i < 4; ++i){
    update_children_F2E(FRONT[i],sid, PGtree, F2E);
  }
  // BACK
  sid = (F2E(0,pPGi[5]-1) == PHi) ? 0 : 1;
  F2E(sid,BACK[0]) = *(children+2);
  F2E(sid,BACK[1]) = *(children+3);
  F2E(sid,BACK[2]) = *(children+7);
  F2E(sid,BACK[3]) = *(children+6);
  for (int i = 0; i < 4; ++i)
    update_children_F2E(BACK[i],sid, PGtree, F2E);
  // INTERNAL faces
  F2E(0,INT[0]) = *(children+1);
  F2E(1,INT[0]) = *(children+2);

  F2E(0,INT[1]) = *(children);
  F2E(1,INT[1]) = *(children+3);

  F2E(0,INT[2]) = *(children+4);
  F2E(1,INT[2]) = *(children+7);

  F2E(0,INT[3]) = *(children+5);
  F2E(1,INT[3]) = *(children+6);

  F2E(0,INT[4]) = *(children);
  F2E(1,INT[4]) = *(children+4);

  F2E(0,INT[5]) = *(children+1);
  F2E(1,INT[5]) = *(children+5);

  F2E(0,INT[6]) = *(children+2);
  F2E(1,INT[6]) = *(children+6);

  F2E(0,INT[7]) = *(children+3);
  F2E(1,INT[7]) = *(children+7);

  F2E(0,INT[8]) = *(children);
  F2E(1,INT[8]) = *(children+1);

  F2E(0,INT[9]) = *(children+3);
  F2E(1,INT[9]) = *(children+2);

  F2E(0,INT[10]) = *(children+7);
  F2E(1,INT[10]) = *(children+6);

  F2E(0,INT[11]) = *(children+4);
  F2E(1,INT[11]) = *(children+5);

}

template <>
template <typename arr_t>
void refiner<K_MESH::Tetrahedron, eSUBDIV_TYPE::ISO>::update_F2E
(E_Int PHi, E_Int *children, E_Int* INT, E_Int** FACES, E_Int ndiag,
ngon_type& ng, K_FLD::IntArray & F2E, tree<arr_t>& PGtree)
{
  E_Int* BOT = &FACES[0][0];  
  E_Int* F1 = &FACES[1][0]; 
  E_Int* F2 = &FACES[2][0];  
  E_Int* F3 = &FACES[3][0];      
    
    
  E_Int* pPGi = ng.PHs.get_facets_ptr(PHi);
  // BOT and TOP
  for (int i = 0; i < 3; i++)
  {
    E_Int sid = (F2E(1,pPGi[0]-1) == PHi) ? 1 : 0;
    F2E(sid,BOT[i]) = *(children+i);
    update_children_F2E(BOT[i],sid, PGtree, F2E);
  }
  E_Int sid = (F2E(1,pPGi[0]-1) == PHi) ? 1 : 0;
  
  
  if (ndiag==1){
  F2E(sid,BOT[3]) = *(children+5);
  update_children_F2E(BOT[3],sid, PGtree, F2E);
  // F1
  sid = (F2E(1,pPGi[1]-1) == PHi) ? 1: 0;
  F2E(sid,F1[0]) = *(children);
  F2E(sid,F1[1]) = *(children+1);
  F2E(sid,F1[2]) = *(children+3);
  F2E(sid,F1[3]) = *(children+4);
  for (int i = 0; i < 4; ++i)
    update_children_F2E(F1[i], sid, PGtree, F2E);
  // F2
  sid = (F2E(0,pPGi[2]-1) == PHi) ? 0 : 1;
  F2E(sid,F2[0]) = *(children+1);
  F2E(sid,F2[1]) = *(children+2);
  F2E(sid,F2[2]) = *(children+3);
  F2E(sid,F2[3]) = *(children+6);
  for (int i = 0; i < 4; ++i)
    update_children_F2E(F2[i],sid, PGtree, F2E);
  // F3
  sid = (F2E(1,pPGi[3]-1) == PHi) ? 1: 0;
  F2E(sid,F3[0]) = *(children+2);
  F2E(sid,F3[1]) = *(children);
  F2E(sid,F3[2]) = *(children+3);
  F2E(sid,F3[3]) = *(children+7);
  for (int i = 0; i < 4; ++i)
    update_children_F2E(F3[i],sid, PGtree, F2E);
  // INTERNAL faces
  F2E(0,INT[0]) = *(children);
  F2E(1,INT[0]) = *(children+4);

  F2E(0,INT[1]) = *(children+1);
  F2E(1,INT[1]) = *(children+5);

  F2E(0,INT[2]) = *(children+2);
  F2E(1,INT[2]) = *(children+6);

  F2E(0,INT[3]) = *(children+3);
  F2E(1,INT[3]) = *(children+7);

  F2E(0,INT[4]) = *(children+5);
  F2E(1,INT[4]) = *(children+6);

  F2E(0,INT[5]) = *(children+4);
  F2E(1,INT[5]) = *(children+7);

  F2E(0,INT[6]) = *(children+4);
  F2E(1,INT[6]) = *(children+5);

  F2E(0,INT[7]) = *(children+7);
  F2E(1,INT[7]) = *(children+6);
  }
  else if (ndiag==2) {
  F2E(sid,BOT[3]) = *(children+4);
  update_children_F2E(BOT[3],sid, PGtree, F2E);
  // F1
  sid = (F2E(1,pPGi[1]-1) == PHi) ? 1: 0;
  F2E(sid,F1[0]) = *(children);
  F2E(sid,F1[1]) = *(children+1);
  F2E(sid,F1[2]) = *(children+3);
  F2E(sid,F1[3]) = *(children+5);
  for (int i = 0; i < 4; ++i)
    update_children_F2E(F1[i], sid, PGtree, F2E);
  // F2
  sid = (F2E(0,pPGi[2]-1) == PHi) ? 0 : 1;
  F2E(sid,F2[0]) = *(children+1);
  F2E(sid,F2[1]) = *(children+2);
  F2E(sid,F2[2]) = *(children+3);
  F2E(sid,F2[3]) = *(children+7);
  for (int i = 0; i < 4; ++i)
    update_children_F2E(F2[i],sid, PGtree, F2E);
  // F3
  sid = (F2E(1,pPGi[3]-1) == PHi) ? 1: 0;
  F2E(sid,F3[0]) = *(children+2);
  F2E(sid,F3[1]) = *(children);
  F2E(sid,F3[2]) = *(children+3);
  F2E(sid,F3[3]) = *(children+6);
  for (int i = 0; i < 4; ++i)
    update_children_F2E(F3[i],sid, PGtree, F2E);
  // INTERNAL faces
  F2E(0,INT[0]) = *(children);
  F2E(1,INT[0]) = *(children+4);

  F2E(0,INT[1]) = *(children+1);
  F2E(1,INT[1]) = *(children+5);

  F2E(0,INT[2]) = *(children+2);
  F2E(1,INT[2]) = *(children+6);

  F2E(0,INT[3]) = *(children+3);
  F2E(1,INT[3]) = *(children+7);

  F2E(0,INT[4]) = *(children+5);
  F2E(1,INT[4]) = *(children+7);

  F2E(0,INT[5]) = *(children+4);
  F2E(1,INT[5]) = *(children+6);

  F2E(0,INT[6]) = *(children+5);
  F2E(1,INT[6]) = *(children+4);

  F2E(0,INT[7]) = *(children+7);
  F2E(1,INT[7]) = *(children+6);
  }
  else {
  F2E(sid,BOT[3]) = *(children+6);
  update_children_F2E(BOT[3],sid, PGtree, F2E);
  // F1
  sid = (F2E(1,pPGi[1]-1) == PHi) ? 1: 0;
  F2E(sid,F1[0]) = *(children);
  F2E(sid,F1[1]) = *(children+1);
  F2E(sid,F1[2]) = *(children+3);
  F2E(sid,F1[3]) = *(children+7);
  for (int i = 0; i < 4; ++i)
    update_children_F2E(F1[i], sid, PGtree, F2E);
  // F2
  sid = (F2E(0,pPGi[2]-1) == PHi) ? 0 : 1;
  F2E(sid,F2[0]) = *(children+1);
  F2E(sid,F2[1]) = *(children+2);
  F2E(sid,F2[2]) = *(children+3);
  F2E(sid,F2[3]) = *(children+5);
  for (int i = 0; i < 4; ++i)
    update_children_F2E(F2[i],sid, PGtree, F2E);
  // F3
  sid = (F2E(1,pPGi[3]-1) == PHi) ? 1: 0;
  F2E(sid,F3[0]) = *(children+2);
  F2E(sid,F3[1]) = *(children);
  F2E(sid,F3[2]) = *(children+3);
  F2E(sid,F3[3]) = *(children+4);
  for (int i = 0; i < 4; ++i)
    update_children_F2E(F3[i],sid, PGtree, F2E);
  // INTERNAL faces
  F2E(0,INT[0]) = *(children);
  F2E(1,INT[0]) = *(children+4);

  F2E(0,INT[1]) = *(children+1);
  F2E(1,INT[1]) = *(children+5);

  F2E(0,INT[2]) = *(children+2);
  F2E(1,INT[2]) = *(children+6);

  F2E(0,INT[3]) = *(children+3);
  F2E(1,INT[3]) = *(children+7);

  F2E(0,INT[4]) = *(children+6);
  F2E(1,INT[4]) = *(children+4);

  F2E(0,INT[5]) = *(children+5);
  F2E(1,INT[5]) = *(children+7);

  F2E(0,INT[6]) = *(children+6);
  F2E(1,INT[6]) = *(children+5);

  F2E(0,INT[7]) = *(children+4);
  F2E(1,INT[7]) = *(children+7);
  }
}

template <>
template <typename arr_t>
void refiner<K_MESH::Pyramid, eSUBDIV_TYPE::ISO>::update_F2E
(E_Int PHi, E_Int* children, E_Int* INT, E_Int** FACES, E_Int ndiag,
ngon_type& ng, K_FLD::IntArray & F2E, tree<arr_t>& PGtree)
{
  E_Int* BOT = &FACES[0][0];  
  E_Int* F1 = &FACES[1][0]; 
  E_Int* F2 = &FACES[2][0];  
  E_Int* F3 = &FACES[3][0];  
  E_Int* F4 = &FACES[4][0];      
    
  E_Int* pPGi = ng.PHs.get_facets_ptr(PHi);
  // BOT
  for (int i = 0; i < 4; i++)
  {
    E_Int sid = (F2E(1,pPGi[0]-1) == PHi) ? 1 : 0;
    F2E(sid,BOT[i]) = *(children+i);
    update_children_F2E(BOT[i],sid, PGtree, F2E);
  }
  // F1
  E_Int sid = (F2E(1,pPGi[1]-1) == PHi) ? 1: 0;
  F2E(sid,F1[0]) = *children;
  F2E(sid,F1[1]) = *(children+1);
  F2E(sid,F1[2]) = *(children+4);
  F2E(sid,F1[3]) = *(children+6);
  for (int i = 0; i < 4; ++i){
    update_children_F2E(F1[i], sid, PGtree, F2E);
  }
  // F2
  sid = (F2E(0,pPGi[2]-1) == PHi) ? 0 : 1;
  F2E(sid,F2[0]) = *(children+1);
  F2E(sid,F2[1]) = *(children+2);
  F2E(sid,F2[2]) = *(children+4);
  F2E(sid,F2[3]) = *(children+7);
  for (int i = 0; i < 4; ++i){
    update_children_F2E(F2[i],sid, PGtree, F2E);
  }
  // F3
  sid = (F2E(1,pPGi[3]-1) == PHi) ? 1: 0;
  F2E(sid,F3[0]) = *(children+2);
  F2E(sid,F3[1]) = *(children+3);
  F2E(sid,F3[2]) = *(children+4);
  F2E(sid,F3[3]) = *(children+8);
  for (int i = 0; i < 4; ++i){
    update_children_F2E(F3[i],sid, PGtree, F2E);
  }
  // F4
  sid = (F2E(0,pPGi[4]-1) == PHi) ? 0 : 1;
  F2E(sid,F4[0]) = *(children+3);
  F2E(sid,F4[1]) = *(children);
  F2E(sid,F4[2]) = *(children+4);
  F2E(sid,F4[3]) = *(children+9);
  for (int i = 0; i < 4; ++i){
    update_children_F2E(F4[i],sid, PGtree, F2E);
  }
  // INTERNAL faces
//  for (int i=0; i<12; i++){
//      std::cout << "INT["<< i << "]= " << INT[i] << std::endl;
//  }
//  std::cout << "PG size= " << ng.PGs.size() << std::endl;
//  std::cout << "F2E size= " << F2E.cols() << std::endl;
  F2E(0,INT[0]) = *(children);
  F2E(1,INT[0]) = *(children+6);

  F2E(0,INT[1]) = *(children+1);
  F2E(1,INT[1]) = *(children+6);

  F2E(0,INT[2]) = *(children+1);
  F2E(1,INT[2]) = *(children+7);

  F2E(0,INT[3]) = *(children+2);
  F2E(1,INT[3]) = *(children+7);

  F2E(0,INT[4]) = *(children+2);
  F2E(1,INT[4]) = *(children+8);

  F2E(0,INT[5]) = *(children+3);
  F2E(1,INT[5]) = *(children+8);

  F2E(0,INT[6]) = *(children+3);
  F2E(1,INT[6]) = *(children+9);

  F2E(0,INT[7]) = *(children);
  F2E(1,INT[7]) = *(children+9);

  F2E(0,INT[8]) = *(children+5);
  F2E(1,INT[8]) = *(children+6);

  F2E(0,INT[9]) = *(children+5);
  F2E(1,INT[9]) = *(children+7);

  F2E(0,INT[10]) = *(children+5);
  F2E(1,INT[10]) = *(children+8);

  F2E(0,INT[11]) = *(children+5);
  F2E(1,INT[11]) = *(children+9);
  
  F2E(0,INT[12]) = *(children+5);
  F2E(1,INT[12]) = *(children+4);
}

template <>
template<typename arr_t>
void refiner<K_MESH::Prism, eSUBDIV_TYPE::ISO>::update_F2E
(E_Int PHi, E_Int* children, E_Int* INT, E_Int** FACES, E_Int ndiag,
ngon_type& ng, K_FLD::IntArray & F2E, tree<arr_t>& PGtree)
{
  E_Int* BOT = &FACES[0][0];  
  E_Int* F1 = &FACES[1][0]; 
  E_Int* F2 = &FACES[2][0];  
  E_Int* F3 = &FACES[3][0];  
  E_Int* TOP = &FACES[4][0];     

  E_Int* pPGi = ng.PHs.get_facets_ptr(PHi);
  // BOT
  E_Int sid = (F2E(1,pPGi[0]-1) == PHi) ? 1 : 0;
  for (int i = 0; i < 4; i++)
  {
    F2E(sid,BOT[i]) = *(children+i);
    update_children_F2E(BOT[i],sid, PGtree, F2E);
  }
  
  
  // F1
  sid = (F2E(1,pPGi[1]-1) == PHi) ? 1: 0;
  F2E(sid,F1[0]) = *children;
  F2E(sid,F1[1]) = *(children+1);
  F2E(sid,F1[2]) = *(children+5);
  F2E(sid,F1[3]) = *(children+4);
  for (int i = 0; i < 4; ++i){
    update_children_F2E(F1[i], sid, PGtree, F2E);
  }
  // F2
  sid = (F2E(0,pPGi[2]-1) == PHi) ? 0 : 1;
  F2E(sid,F2[0]) = *(children+1);
  F2E(sid,F2[1]) = *(children+2);
  F2E(sid,F2[2]) = *(children+6);
  F2E(sid,F2[3]) = *(children+5);
  for (int i = 0; i < 4; ++i){
    update_children_F2E(F2[i],sid, PGtree, F2E);
  }
  // F3
  sid = (F2E(1,pPGi[3]-1) == PHi) ? 1: 0;
  F2E(sid,F3[0]) = *(children+2);
  F2E(sid,F3[1]) = *(children);
  F2E(sid,F3[2]) = *(children+4);
  F2E(sid,F3[3]) = *(children+6);
  for (int i = 0; i < 4; ++i){
    update_children_F2E(F3[i],sid, PGtree, F2E);
  }
  // TOP
  sid = (F2E(0,pPGi[4]-1) == PHi) ? 0 : 1;
  F2E(sid,TOP[0]) = *(children+4);
  F2E(sid,TOP[1]) = *(children+5);
  F2E(sid,TOP[2]) = *(children+6);
  F2E(sid,TOP[3]) = *(children+7);
  for (int i = 0; i < 4; ++i){
    update_children_F2E(TOP[i],sid, PGtree, F2E);
  }
  // INTERNAL faces

  F2E(0,INT[0]) = *children;
  F2E(1,INT[0]) = *(children+4);

  F2E(0,INT[1]) = *(children+1);
  F2E(1,INT[1]) = *(children+5);

  F2E(0,INT[2]) = *(children+2);
  F2E(1,INT[2]) = *(children+6);

  F2E(0,INT[3]) = *(children+3);
  F2E(1,INT[3]) = *(children+7);

  F2E(0,INT[4]) = *(children);
  F2E(1,INT[4]) = *(children+3);

  F2E(0,INT[5]) = *(children+1);
  F2E(1,INT[5]) = *(children+3);

  F2E(0,INT[6]) = *(children+2);
  F2E(1,INT[6]) = *(children+3);

  F2E(0,INT[7]) = *(children+4);
  F2E(1,INT[7]) = *(children+7);

  F2E(0,INT[8]) = *(children+5);
  F2E(1,INT[8]) = *(children+7);

  F2E(0,INT[9]) = *(children+6);
  F2E(1,INT[9]) = *(children+7);
}

///
template <>
template <typename arr_t>
void refiner<K_MESH::Hexahedron, eSUBDIV_TYPE::ISO>::refine_PHs
(const Vector_t<E_Int> &adap_incr, 
 ngon_type& ng, tree<arr_t> & PGtree, tree<arr_t> & PHtree, K_FLD::FloatArray& crd, K_FLD::IntArray & F2E)
{
  assert(adap_incr.size() <= ng.PHs.size());
  //
  std::vector<E_Int> PHadap;
  get_PHs_to_refine(ng, PHtree, adap_incr, PHadap);
  E_Int nb_phs = PHadap.size();
  if (nb_phs == 0) return;

  E_Int pos = crd.cols(); // first face center in crd (if relevant)
  E_Int nb_pgs0 = ng.PGs.size();
  E_Int nb_phs0 = ng.PHs.size();
  // Reserve space for children in the tree
  reserve_mem_PHs(crd, ng, PHadap, PGtree, PHtree, F2E);

  static constexpr E_Int nbc = NBC<K_MESH::Hexahedron, ISO>::nbc;

#ifndef DEBUG_HIERARCHICAL_MESH
#pragma omp parallel for
#endif
  for (E_Int i = 0; i < nb_phs; ++i)
  {
    E_Int FACES[K_MESH::Hexahedron::NB_BOUNDS][NBC<K_MESH::Hexahedron, ISO>::nbc];
    E_Int PHi = PHadap[i];
    E_Int nodes[27]; // 0-26 to crd1 indexes
    ////
    E_Int* BOT = &FACES[0][0];
    E_Int* TOP = &FACES[1][0];
    E_Int* LEFT = &FACES[2][0];
    E_Int* RIGHT = &FACES[3][0];
    E_Int* FRONT = &FACES[4][0];
    E_Int* BACK = &FACES[5][0];
            
    E_Int *FACESindex[6] = {FACES[0], FACES[1], FACES[2], FACES[3], FACES[4], FACES[5]}; // convert int[][] to int** ---> see stackoverflow
    get_nodes_PHi(nodes, PHi, pos+i, FACESindex, ng, crd, F2E, PGtree);

    E_Int points[27];
    K_MESH::Hexahedron::get_internal(nodes,points);

    E_Int INT[12];
    // INTERNAL faces
    for ( size_t j = 0; j < 3; ++j)
    {
      E_Int firstPGChild = nb_pgs0 + 12 * i + 4 * j;
      
      NUGA::Q9::split<eSUBDIV_TYPE::ISO>(ng.PGs, &points[9*j], firstPGChild);

      INT[4 * j]     = firstPGChild;
      INT[4 * j + 1] = firstPGChild + 1;
      INT[4 * j + 2] = firstPGChild + 2;
      INT[4 * j + 3] = firstPGChild + 3;
    }

    // the 8 children of PH
    E_Int firstPHchild = nb_phs0 + 8 * i;
    NUGA::H27::splitH27(crd, ng.PHs, INT, BOT, TOP, LEFT, RIGHT, FRONT, BACK, firstPHchild);

    // set them in the tree
    PHtree.set_children(PHi, firstPHchild, 8);

    // update F2E
    E_Int ndiag(0); // spécialisation argument ndiag pour tétra
    E_Int PHchildren[8];
    for (size_t k = 0; k < 8; ++k) PHchildren[k] = nb_phs0 + 8 * i + k;
    update_F2E(PHi, PHchildren, INT, FACESindex, ndiag, ng, F2E, PGtree);
  }  
}

///
template <>
template <typename arr_t>
void refiner<K_MESH::Tetrahedron, eSUBDIV_TYPE::ISO>::refine_PHs
(const Vector_t<E_Int> &adap_incr, 
 ngon_type& ng, tree<arr_t> & PGtree, tree<arr_t> & PHtree, K_FLD::FloatArray& crd, K_FLD::IntArray & F2E)
{
  assert(adap_incr.size() <= ng.PHs.size());

  //alexis : fixme : remplacer l'appel ais_HX8 pour l'emploi de ng.PHs.type
  std::vector<E_Int> PHadap;
  for (size_t i = 0; i < adap_incr.size(); ++i)
    if (adap_incr[i] > 0 && PHtree.is_enabled(i) && (K_MESH::Polyhedron<0>::is_TH4(ng.PGs, ng.PHs.get_facets_ptr(i), ng.PHs.stride(i))))
      PHadap.push_back(i);
  
  E_Int nb_phs = PHadap.size();
  if (nb_phs == 0) return;

  // internal PGs created (12 per PH)
  // Reserve space for internal faces in the tree
  E_Int current_sz = PGtree.size();
  assert(current_sz != -1);

  PGtree.resize_hierarchy(current_sz+8*nb_phs);

  E_Int nb_pgs0 = ng.PGs.size();

  // And in the mesh : each H27 is split in 36 Pgs including 12 new internal Q4

  //_ng.PGs.expand_n_fixed_stride(8*nb_phs, 4/*Q4 stride*/);
  ng.PGs.expand_n_fixed_stride(8*nb_phs, 3/*Q4 stride*/);

  F2E.resize(2,nb_pgs0+8*nb_phs,E_IDX_NONE);

  // for the children PH
  // Reserve space for children in the tree
  PHtree.resize(PHadap, 8);
  // And in the mesh : each H27 is split into 8 H27
  E_Int nb_phs0 = ng.PHs.size();
  ng.PHs.expand_n_fixed_stride(8*nb_phs, 4/*TH4 stride*/);

  E_Int pos = crd.cols();
  
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
    E_Int FACES[4][4];     
    E_Int PHi = PHadap[i];
    E_Int nodes[10]; // 0-26 to crd1 indexes
 
    E_Int* BOT = &FACES[0][0];
    E_Int* F1 = &FACES[1][0];
    E_Int* F2 = &FACES[2][0];
    E_Int* F3 = &FACES[3][0];
            
    E_Int *FACESindex[4] = {FACES[0], FACES[1], FACES[2], FACES[3]}; // convert int[][] to int** ---> see stackoverflow
    get_nodes_PHi(nodes, PHi, pos+i, FACESindex, ng, crd, F2E, PGtree);
    
    /// Choix de la diagonale ici param nodes[10]
//    E_Float d;
       
    E_Float d1=K_FUNC::sqrDistance(crd.col(nodes[8]-1),crd.col(nodes[6]-1), 3);
    E_Float d2=K_FUNC::sqrDistance(crd.col(nodes[7]-1),crd.col(nodes[5]-1), 3);
    E_Float d3=K_FUNC::sqrDistance(crd.col(nodes[9]-1),crd.col(nodes[4]-1), 3);

    E_Int ndiag= ((d1 <= d2) && (d1 <= d3) ) ? 1 : ((d2 <= d1) && (d2 <= d3) ) ? 2 : 3;
    
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
    E_Int PHchildren[8];
    for (size_t k = 0; k < 8; ++k) PHchildren[k] = nb_phs0 + 8 * i + k;
    
    E_Int* h271 = ng.PHs.get_facets_ptr(PHchildren[0]);
    E_Int* h272 = ng.PHs.get_facets_ptr(PHchildren[1]);
    E_Int* h273 = ng.PHs.get_facets_ptr(PHchildren[2]);
    E_Int* h274 = ng.PHs.get_facets_ptr(PHchildren[3]);
    E_Int* h275 = ng.PHs.get_facets_ptr(PHchildren[4]);
    E_Int* h276 = ng.PHs.get_facets_ptr(PHchildren[5]);
    E_Int* h277 = ng.PHs.get_facets_ptr(PHchildren[6]);
    E_Int* h278 = ng.PHs.get_facets_ptr(PHchildren[7]);

    NUGA::H27::splitT10(crd, INT, BOT, F1, F2, F3, h271, h272, h273, h274, h275, h276, h277, h278, ndiag);

    // set them in the tree
    PHtree.set_children(PHi, PHchildren, 8);

    // update F2E
    update_F2E(PHi, PHchildren, INT, FACESindex, ndiag, ng, F2E, PGtree);
  }
}

///
template <>
template <typename arr_t>
void refiner<K_MESH::Pyramid, eSUBDIV_TYPE::ISO>::refine_PHs
(const Vector_t<E_Int> &adap_incr, 
 ngon_type& ng, tree<arr_t> & PGtree, tree<arr_t> & PHtree, K_FLD::FloatArray& crd, K_FLD::IntArray & F2E)
{
  assert(adap_incr.size() <= ng.PHs.size());

  //alexis : fixme : remplacer l'appel is_PY5 pour l'emploi de ng.PHs.type
  std::vector<E_Int> PHadap;
  for (size_t i = 0; i < adap_incr.size(); ++i)
    if (adap_incr[i] > 0 && PHtree.is_enabled(i) && (K_MESH::Polyhedron<0>::is_PY5(ng.PGs, ng.PHs.get_facets_ptr(i), ng.PHs.stride(i))))
      PHadap.push_back(i);
  
  E_Int nb_phs = PHadap.size();
  if (nb_phs == 0) return;
  // internal PGs created (12 per PH)
  // Reserve space for internal faces in the tree
  E_Int current_sz = PGtree.size();
  assert(current_sz != -1);

  PGtree.resize_hierarchy(current_sz+13*nb_phs); // pyra faces internes 13

  E_Int nb_pgs0 = ng.PGs.size();

  // And in the mesh : each H27 is split in 36 Pgs including 12 new internal Q4
  ng.PGs.expand_n_fixed_stride(12*nb_phs, 3/*Q4 stride*/); // 12 triangle stride 3 
  ng.PGs.expand_n_fixed_stride(1*nb_phs, 4/*Q4 stride*/);  // 1 quad stride 4  

  F2E.resize(2,nb_pgs0+13*nb_phs,E_IDX_NONE);
//  std::cout << "ng.PG size= " << ng.PGs.size() <<std::endl;
  // for the children PH
  // Reserve space for children in the tree
  PHtree.resize(PHadap, 10);

  // And in the mesh : each H27 is split into 8 H27
  E_Int nb_phs0 = ng.PHs.size();
  ng.PHs.expand_n_fixed_stride(6*nb_phs, 5/*H27 stride*/); // 6 pyra stride 5
  ng.PHs.expand_n_fixed_stride(4*nb_phs, 4/*H27 stride*/); // 4 tetra stride 4

  E_Int pos = crd.cols();
  
  int j;

#ifndef DEBUG_HIERARCHICAL_MESH
#pragma omp parallel for private (j)
#endif
  for (E_Int i = 0; i < nb_phs; ++i)
  {  
    E_Int FACES[5][4];   
    //std::cout << "k= " << i << " sur " << nb_phs << std::endl;
    E_Int PHi = PHadap[i];
    E_Int nodes[14]; // 0-26 to crd1 indexes
    //E_Int BOT[4], F1[4], F2[4], F3[4], F4[4];// children lists des faces
    
    E_Int* BOT = &FACES[0][0];
    E_Int* F1 = &FACES[1][0];
    E_Int* F2 = &FACES[2][0];
    E_Int* F3 = &FACES[3][0];
    E_Int* F4 = &FACES[4][0];
            
    E_Int *FACESindex[5] = {FACES[0], FACES[1], FACES[2], FACES[3], FACES[4]}; // convert int[][] to int** ---> see stackoverflow
    get_nodes_PHi(nodes, PHi, pos+i, FACESindex, ng, crd, F2E, PGtree);

    E_Int INT[13];
    // INTERNAL faces  
    for ( j = 0; j < 12; ++j){
      INT[j]= nb_pgs0 + 12*i + j;
    }

    INT[12]= nb_pgs0 + 12*nb_phs + i;
    
    E_Int* q41 = ng.PGs.get_facets_ptr(INT[0]);
    E_Int* q42 = ng.PGs.get_facets_ptr(INT[1]);
    E_Int* q43 = ng.PGs.get_facets_ptr(INT[2]);
    E_Int* q44 = ng.PGs.get_facets_ptr(INT[3]);
    E_Int* q45 = ng.PGs.get_facets_ptr(INT[4]);
    E_Int* q46 = ng.PGs.get_facets_ptr(INT[5]);
    E_Int* q47 = ng.PGs.get_facets_ptr(INT[6]);
    E_Int* q48 = ng.PGs.get_facets_ptr(INT[7]);
    E_Int* q49 = ng.PGs.get_facets_ptr(INT[8]);
    E_Int* q410 = ng.PGs.get_facets_ptr(INT[9]);
    E_Int* q411 = ng.PGs.get_facets_ptr(INT[10]);
    E_Int* q412 = ng.PGs.get_facets_ptr(INT[11]);
    E_Int* q413 = ng.PGs.get_facets_ptr(INT[12]);

    q41[0]=nodes[5]; q41[1]=nodes[9]; q41[2]=nodes[13];
    q42[0]=nodes[9]; q42[1]=nodes[5]; q42[2]=nodes[10];
    q43[0]=nodes[6]; q43[1]=nodes[9]; q43[2]=nodes[10];
    q44[0]=nodes[9]; q44[1]=nodes[6]; q44[2]=nodes[11];    
    q45[0]=nodes[7]; q45[1]=nodes[9]; q45[2]=nodes[11];
    q46[0]=nodes[9]; q46[1]=nodes[7]; q46[2]=nodes[12];
    q47[0]=nodes[8]; q47[1]=nodes[9]; q47[2]=nodes[12];
    q48[0]=nodes[9]; q48[1]=nodes[8]; q48[2]=nodes[13];
    
    q49[0]=nodes[10]; q49[1]=nodes[13]; q49[2]=nodes[9];
    q410[0]=nodes[11]; q410[1]=nodes[10]; q410[2]=nodes[9];
    q411[0]=nodes[12]; q411[1]=nodes[11]; q411[2]=nodes[9];
    q412[0]=nodes[13]; q412[1]=nodes[12]; q412[2]=nodes[9];
    
    q413[0]=nodes[13]; q413[1]=nodes[10]; q413[2]=nodes[11]; q413[3]=nodes[12];

    // the 8 children of PH
    E_Int PHichildr[10];

    for (int j = 0; j < 6; ++j){
      PHichildr[j] = nb_phs0 + 6*i + j;
    }
    for (int j = 0; j < 4; ++j){
      PHichildr[6+j] = nb_phs0 + 6*nb_phs + 4*i +j;
    }
    
    E_Int* h271 = ng.PHs.get_facets_ptr(PHichildr[0]);
    E_Int* h272 = ng.PHs.get_facets_ptr(PHichildr[1]);
    E_Int* h273 = ng.PHs.get_facets_ptr(PHichildr[2]);
    E_Int* h274 = ng.PHs.get_facets_ptr(PHichildr[3]);        
    E_Int* h275 = ng.PHs.get_facets_ptr(PHichildr[4]);
    E_Int* h276 = ng.PHs.get_facets_ptr(PHichildr[5]);
    E_Int* h277 = ng.PHs.get_facets_ptr(PHichildr[6]);
    E_Int* h278 = ng.PHs.get_facets_ptr(PHichildr[7]);
    E_Int* h279 = ng.PHs.get_facets_ptr(PHichildr[8]);
    E_Int* h2710 = ng.PHs.get_facets_ptr(PHichildr[9]);


    NUGA::H27::splitP13(crd, INT, BOT, F1, F2, F3, F4, h271, h272, h273, h274, h275, h276, h277, h278, h279, h2710);

    // set them in the tree
    PHtree.set_children(PHi, PHichildr, 10);

    // update F2E
    E_Int ndiag(0);
    update_F2E(PHi,PHichildr,INT, FACESindex, ndiag,ng, F2E, PGtree);
    //std::cout << "is hybrid after refine pyra i=  "<< i << "? " << (is_hybrid(ng) ? "oui" : "non") << std::endl;

  }  
  
}

///
template <>
template <typename arr_t>
void refiner<K_MESH::Prism, eSUBDIV_TYPE::ISO>::refine_PHs
(const Vector_t<E_Int> &adap_incr, 
 ngon_type& ng, tree<arr_t> & PGtree, tree<arr_t> & PHtree, K_FLD::FloatArray& crd, K_FLD::IntArray & F2E)
{
  assert(adap_incr.size() <= ng.PHs.size());

  //alexis : fixme : remplacer l'appel is_PR6 pour l'emploi de ng.PHs.type
  std::vector<E_Int> PHadap;
  for (size_t i = 0; i < adap_incr.size(); ++i)
    if (adap_incr[i] > 0 && PHtree.is_enabled(i) && (K_MESH::Polyhedron<0>::is_PR6(ng.PGs, ng.PHs.get_facets_ptr(i), ng.PHs.stride(i))))
      PHadap.push_back(i);
  
  E_Int nb_phs = PHadap.size();
  if (nb_phs == 0) return;
  // internal PGs created (12 per PH)
  // Reserve space for internal faces in the tree
  E_Int current_sz = PGtree.size();
  assert(current_sz != -1);

  PGtree.resize_hierarchy(current_sz+10*nb_phs); // prism faces internes 10  (4)

  E_Int nb_pgs0 = ng.PGs.size();
  // faces internes : 4 TRI + 6 QUAD
  ng.PGs.expand_n_fixed_stride(4*nb_phs, 3/*T3 stride*/);
  ng.PGs.expand_n_fixed_stride(6*nb_phs, 4/*Q4 stride*/);
  
  F2E.resize(2,nb_pgs0+10*nb_phs,E_IDX_NONE);
  PHtree.resize(PHadap, 8);

  E_Int nb_phs0 = ng.PHs.size();
  ng.PHs.expand_n_fixed_stride(8*nb_phs, 5); // 8 prism
  E_Int pos = crd.cols();
  
  int j;

#ifndef DEBUG_HIERARCHICAL_MESH
#pragma omp parallel for private (j)
#endif
  for (E_Int i = 0; i < nb_phs; ++i)
  {  
    E_Int FACES[5][4];  
    E_Int PHi = PHadap[i];
    E_Int nodes[18]; // 0-26 to crd1 indexes
    //E_Int BOT[4], F1[4], F2[4], F3[4], TOP[4];// children lists des faces
    E_Int* BOT = &FACES[0][0];
    E_Int* F1 = &FACES[1][0];
    E_Int* F2 = &FACES[2][0];
    E_Int* F3 = &FACES[3][0];
    E_Int* TOP = &FACES[4][0];
            
    E_Int *FACESindex[5] = {FACES[0], FACES[1], FACES[2], FACES[3], FACES[4]}; // convert int[][] to int** ---> see stackoverflow  
    get_nodes_PHi(nodes, PHi, pos+i, FACESindex,
                  ng, crd, F2E, PGtree);

    E_Int INT[10];  
    for ( j = 0; j < 4; ++j){
      INT[j]= nb_pgs0 + 4*i + j;
    }

    for ( j = 0; j < 6; ++j){
      INT[4+j]= nb_pgs0 + 4*nb_phs + 6*i +j;
    }
    
    E_Int* q41 = ng.PGs.get_facets_ptr(INT[0]);     // triangle
    E_Int* q42 = ng.PGs.get_facets_ptr(INT[1]);
    E_Int* q43 = ng.PGs.get_facets_ptr(INT[2]);
    E_Int* q44 = ng.PGs.get_facets_ptr(INT[3]);
    
    E_Int* q45 = ng.PGs.get_facets_ptr(INT[4]);     // quad
    E_Int* q46 = ng.PGs.get_facets_ptr(INT[5]);
    E_Int* q47 = ng.PGs.get_facets_ptr(INT[6]);
    E_Int* q48 = ng.PGs.get_facets_ptr(INT[7]);
    E_Int* q49 = ng.PGs.get_facets_ptr(INT[8]);
    E_Int* q410 = ng.PGs.get_facets_ptr(INT[9]);
    
    q41[0]=nodes[15]; q41[1]=nodes[12]; q41[2]=nodes[14];
    q42[0]=nodes[12]; q42[1]=nodes[16]; q42[2]=nodes[13];
    q43[0]=nodes[14]; q43[1]=nodes[13]; q43[2]=nodes[17];
    q44[0]=nodes[12]; q44[1]=nodes[13]; q44[2]=nodes[14];
    
    q45[0]=nodes[6]; q45[1]=nodes[8]; q45[2]=nodes[14]; q45[3]=nodes[12];  // étage 1
    q46[0]=nodes[7]; q46[1]=nodes[6]; q46[2]=nodes[12]; q46[3]=nodes[13];
    q47[0]=nodes[8]; q47[1]=nodes[7]; q47[2]=nodes[13]; q47[3]=nodes[14];
    
    q48[0]=nodes[12]; q48[1]=nodes[14]; q48[2]=nodes[11]; q48[3]=nodes[9];  // étage 2    
    q49[0]=nodes[13]; q49[1]=nodes[12]; q49[2]=nodes[9]; q49[3]=nodes[10];
    q410[0]=nodes[14]; q410[1]=nodes[13]; q410[2]=nodes[10]; q410[3]=nodes[11]; 
    
    // the 8 children of PH
    E_Int PHichildr[8];
    for (int j = 0; j < 8; ++j){
      PHichildr[j] = nb_phs0 + 8*i + j;
    }
        
    E_Int* h271 = ng.PHs.get_facets_ptr(PHichildr[0]);
    E_Int* h272 = ng.PHs.get_facets_ptr(PHichildr[1]);
    E_Int* h273 = ng.PHs.get_facets_ptr(PHichildr[2]);
    E_Int* h274 = ng.PHs.get_facets_ptr(PHichildr[3]);        
    E_Int* h275 = ng.PHs.get_facets_ptr(PHichildr[4]);
    E_Int* h276 = ng.PHs.get_facets_ptr(PHichildr[5]);
    E_Int* h277 = ng.PHs.get_facets_ptr(PHichildr[6]);
    E_Int* h278 = ng.PHs.get_facets_ptr(PHichildr[7]);

    NUGA::H27::splitPr18(crd, INT, BOT, F1, F2, F3, TOP, h271, h272, h273, h274, h275, h276, h277, h278);

    // set them in the tree
    PHtree.set_children(PHi, PHichildr, 8);

    // update F2E
    E_Int ndiag(0);
    update_F2E(PHi,PHichildr,INT, FACESindex, ndiag,ng, F2E, PGtree);

  }  
  
}

///
template <>
template <typename arr_t>
void refiner<K_MESH::Basic, eSUBDIV_TYPE::ISO>::refine_PHs
(const sensor_output_t &adap_incr, 
 ngon_type& ng, tree<arr_t> & PGtree, tree<arr_t> & PHtree, K_FLD::FloatArray& crd, K_FLD::IntArray & F2E)
{  
  assert(adap_incr.size() <= ng.PHs.size());

  refiner<K_MESH::Hexahedron, eSUBDIV_TYPE::ISO>::refine_PHs(adap_incr, ng, PGtree, PHtree, crd, F2E);

  refiner<K_MESH::Tetrahedron, eSUBDIV_TYPE::ISO>::refine_PHs(adap_incr, ng, PGtree, PHtree, crd, F2E);

  refiner<K_MESH::Pyramid, eSUBDIV_TYPE::ISO>::refine_PHs(adap_incr, ng, PGtree, PHtree, crd, F2E);

  refiner<K_MESH::Prism, eSUBDIV_TYPE::ISO>::refine_PHs(adap_incr, ng, PGtree, PHtree, crd, F2E);
    
}


// DIRECTIONNEL
//default impl : DIR IS ISO for all elements but HEXA "Layer"
template <>
template <typename arr_t>
void refiner<K_MESH::Tetrahedron, DIR>::refine_PHs
(const sensor_output_t &adap_incr, 
 ngon_type& ng, tree<arr_t> & PGtree, tree<arr_t> & PHtree, K_FLD::FloatArray& crd, K_FLD::IntArray & F2E)
{
  refiner<K_MESH::Tetrahedron, ISO>::refine_PHs(adap_incr, ng, PGtree, PHtree, crd, F2E);
}

///
template <>
template <typename arr_t>
void refiner<K_MESH::Pyramid, DIR>::refine_PHs
(const sensor_output_t &adap_incr,
 ngon_type& ng, tree<arr_t> & PGtree, tree<arr_t> & PHtree, K_FLD::FloatArray& crd, K_FLD::IntArray & F2E)
{
  refiner<K_MESH::Pyramid, ISO>::refine_PHs(adap_incr, ng, PGtree, PHtree, crd, F2E);
}

///
template <>
template <typename arr_t>
void refiner<K_MESH::Prism, DIR>::refine_PHs
(const sensor_output_t &adap_incr,
 ngon_type& ng, tree<arr_t> & PGtree, tree<arr_t> & PHtree, K_FLD::FloatArray& crd, K_FLD::IntArray & F2E)
{
  refiner<K_MESH::Prism, ISO>::refine_PHs(adap_incr, ng, PGtree, PHtree, crd, F2E);
}

///
template <>
template <typename arr_t>
void refiner<K_MESH::Hexahedron, DIR>::refine_PHs
(const sensor_output_t &adap_incr,
 ngon_type& ng, tree<arr_t> & PGtree, tree<arr_t> & PHtree, K_FLD::FloatArray& crd, K_FLD::IntArray & F2E)
{
  //alexis : todo
}

///
template <>
template <typename arr_t>
void refiner<K_MESH::Basic, DIR>::refine_PHs
(const sensor_output_t &adap_incr, 
 ngon_type& ng, tree<arr_t> & PGtree, tree<arr_t> & PHtree, K_FLD::FloatArray& crd, K_FLD::IntArray & F2E)
{
  //alexis : todo
}

} // NUGA


#endif /* REFINER_HXX */


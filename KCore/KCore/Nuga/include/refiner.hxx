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
#include "HX27.hxx"
#include "TH10.hxx"
#include "PR18.hxx"
#include "PY13.hxx"
#include "PHQ4.hxx"

#include "sensor.hxx"

namespace NUGA
{

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
    static void reserve_mem_PGs(K_FLD::FloatArray& crd, ngon_unit& PGs, const Vector_t<E_Int>& PGlist, tree<arr_t> & PGtree, K_FLD::IntArray& F2E, std::vector<E_Int>& childpos);
    ///
    template <typename arr_t>
    static void reserve_mem_PHs(K_FLD::FloatArray& crd, ngon_type& ng, const Vector_t<E_Int>& PHlist, tree<arr_t> & PGtree, tree<arr_t> & PHtree,
                                K_FLD::IntArray& F2E, std::vector<E_Int>& intpos, std::vector<E_Int>& childpos);

    ///
    template <typename arr_t>
    static void __reserve_mem_single_bound_type_PHs
    (K_FLD::FloatArray& crd, ngon_type& ng, const Vector_t<E_Int>& PHlist, tree<arr_t> & PGtree, tree<arr_t> & PHtree,
      K_FLD::IntArray& F2E, std::vector<E_Int>& intpos, std::vector<E_Int>& childpos);
    
};

/// Impl for basic elements
template <typename ELT_t, eSUBDIV_TYPE STYPE>
template <typename arr_t>
void refiner<ELT_t, STYPE>::get_PGs_to_refine
(const ngon_type& ng, const tree<arr_t> & PGtree, const sensor_output_t &adap_incr, Vector_t<E_Int> & PGlist)
{
  E_Int nb_phs = ng.PHs.size();
  // Gets PGs to refine
  E_Int nb_pgs(ng.PGs.size());
  Vector_t<E_Int> is_PG_to_refine = adap_incr.face_adap_incr;
  is_PG_to_refine.resize(nb_pgs, 0);

  // disable non-admissible elements (face_adap_incr may contain info fo Q4 & T3)
  for (E_Int i = 0; i < nb_pgs; ++i)
  {
    if (ng.PGs.stride(i) != ELT_t::NB_NODES)
      is_PG_to_refine[i] = 0;
  }

  //
  for (E_Int i = 0; i < nb_phs; ++i)
  {
    if (adap_incr.cell_adap_incr[i] <= 0) continue;

    E_Int nb_faces = ng.PHs.stride(i);
    const E_Int* faces = ng.PHs.get_facets_ptr(i);

    for (E_Int j = 0; j < nb_faces; ++j)
    {
      E_Int PGi = *(faces + j) - 1;
      if (is_PG_to_refine[PGi] != 0) continue;

      if (PGtree.nb_children(PGi) == 0) // leaf PG => to refine
      {
        if (ng.PGs.stride(PGi) == ELT_t::NB_NODES)
          is_PG_to_refine[PGi] = 1;
      }
    }
  }

  E_Int nb_pgs_ref = std::count(is_PG_to_refine.begin(), is_PG_to_refine.end(), 1);
  PGlist.reserve(nb_pgs_ref);

  for (E_Int i = 0; i < nb_pgs; ++i)
  {
    if (is_PG_to_refine[i]) {
      PGlist.push_back(i);
    }
  }
}

///
template <>
template <typename arr_t>
void refiner<K_MESH::Polygon, ISO_HEX>::get_PGs_to_refine
(const ngon_type& ng, const tree<arr_t> & PGtree, const sensor_output_t &adap_incr, Vector_t<E_Int> & PGlist)
{
  //todo JP
}


///
template <typename ELT_t, eSUBDIV_TYPE STYPE>
template <typename arr_t>
void refiner<ELT_t, STYPE>::get_PHs_to_refine
(const ngon_type& ng, const tree<arr_t> & PHtree, const sensor_output_t &adap_incr, Vector_t<E_Int> & PHlist)
{
  for (size_t i = 0; i < adap_incr.cell_adap_incr.size(); ++i)
    if (adap_incr.cell_adap_incr[i] > 0 && (PHtree.nb_children(i)==0)  && (ELT_t::is_of_type(ng.PGs, ng.PHs.get_facets_ptr(i), ng.PHs.stride(i))))
      PHlist.push_back(i);
}

///
template <>
template <typename arr_t>
void refiner<K_MESH::Polyhedron<0>, ISO_HEX>::get_PHs_to_refine
(const ngon_type& ng, const tree<arr_t> & PHtree, const sensor_output_t &adap_incr, Vector_t<E_Int> & PHlist)
{
  //todo JP
}


/// Impl for  T3/IS
template <>
template <typename arr_t>
void refiner<K_MESH::Triangle, ISO>::reserve_mem_PGs
(K_FLD::FloatArray& crd, ngon_unit& PGs, const Vector_t<E_Int>& PGlist, tree<arr_t> & PGtree, K_FLD::IntArray& F2E, std::vector<E_Int>&childpos)
{
  const E_Int &nbc = subdiv_pol<K_MESH::Triangle, ISO>::NBC;
  E_Int nb_pgs_ref = PGlist.size();

  PGtree.resize(PGlist, nbc);
  
  E_Int nb_pgs0 = PGs.size();

  // And in the mesh : each ELT(T3) is split into nbc
  PGs.expand_n_fixed_stride(nbc * nb_pgs_ref, K_MESH::Triangle::NB_NODES);

  childpos.resize(nb_pgs_ref + 1);//one-pass-the-end size
  E_Int s(nb_pgs0);
  for (size_t i = 0; i < nb_pgs_ref + 1; ++i)
  {
    childpos[i] = s;
    s += nbc;
  }

  F2E.resize(2, nb_pgs0 + nbc * nb_pgs_ref, E_IDX_NONE);
}

/// Impl for Q4/ISO
template <>
template <typename arr_t>
void refiner<K_MESH::Quadrangle, ISO>::reserve_mem_PGs
(K_FLD::FloatArray& crd, ngon_unit& PGs, const Vector_t<E_Int>& PGlist, tree<arr_t> & PGtree, K_FLD::IntArray& F2E, std::vector<E_Int>&childpos)
{
  const E_Int &nbc = subdiv_pol<K_MESH::Quadrangle, ISO>::NBC;
  E_Int nb_pgs_ref = PGlist.size();

  PGtree.resize(PGlist, nbc);
  
  E_Int nb_pgs0 = PGs.size();

  // And in the mesh : each ELT(T3/) is split into nbc
  PGs.expand_n_fixed_stride(nbc * nb_pgs_ref, K_MESH::Quadrangle::NB_NODES);

  childpos.resize(nb_pgs_ref + 1);//one-pass-the-end size
  E_Int s(nb_pgs0);
  for (size_t i = 0; i < nb_pgs_ref+1; ++i)
  {
    childpos[i] = s;
    s += nbc;
  }

  F2E.resize(2, nb_pgs0 + nbc * nb_pgs_ref, E_IDX_NONE); 

  // reserve space for face centers
  crd.resize(3, crd.cols() + nb_pgs_ref);
}


/// Impl for PG/ISO_HEX
template <>
template <typename arr_t>
void refiner<K_MESH::Polygon, ISO_HEX>::reserve_mem_PGs
(K_FLD::FloatArray& crd, ngon_unit& PGs, const Vector_t<E_Int>& PGlist, tree<arr_t> & PGtree, K_FLD::IntArray& F2E, std::vector<E_Int>& childpos)
{
  // info : nbc is variable and equal to number of nodes

  E_Int nb_pgs_ref = PGlist.size();
  E_Int nb_pgs0 = PGs.size();

  std::vector<E_Int> pregnant;
  subdiv_pol<K_MESH::Polygon, ISO_HEX>::nbc_list(PGs, PGlist, pregnant);

  PGtree.resize(PGlist, pregnant);

  // expand ng.PGs for children with total nb of QUADS
  // todo JP

  childpos.resize(nb_pgs_ref + 1);//one-pass-the-end size
  E_Int s(nb_pgs0);
  for (size_t i = 0; i < nb_pgs_ref; ++i)
  {
    childpos[i] = s;
    s += pregnant[i];
  }
  childpos[nb_pgs_ref] = s;

  F2E.resize(2, s, E_IDX_NONE);
  
  // reserve space for face centers
  crd.resize(3, crd.cols() + nb_pgs_ref);
}

// used for HEXA or TETRA
template <typename ELT_t, eSUBDIV_TYPE STYPE>
template <typename arr_t>
void refiner<ELT_t, STYPE>::__reserve_mem_single_bound_type_PHs
(K_FLD::FloatArray& crd, ngon_type& ng, const Vector_t<E_Int>& PHlist, tree<arr_t> & PGtree, tree<arr_t> & PHtree,
  K_FLD::IntArray& F2E, std::vector<E_Int>& intpos, std::vector<E_Int>& childpos)
{
  const E_Int &nbc = subdiv_pol<ELT_t, ISO>::PHNBC;
  const E_Int& nbi = subdiv_pol<ELT_t, ISO>::NBI;

  using face_t = typename ELT_t::boundary_type;

  E_Int nb_phs_ref = PHlist.size();

  // internal PGs created (12 per PH)
  // Reserve space for internal faces in the tree
  E_Int current_sz = PGtree.size();
  assert(current_sz != -1);

  //just sync the tree attributes (internals have no parents so do not a call to PGtree.resize)
  PGtree.resize_hierarchy(current_sz + nbi * nb_phs_ref); //fixme : in hmesh.init ?

  E_Int nb_pgs0 = ng.PGs.size();

  // each Hex/ISO is split in 36 PGs including 12 new internal Q4
  ng.PGs.expand_n_fixed_stride(nbi * nb_phs_ref, face_t::NB_NODES);

  intpos.resize(nb_phs_ref + 1);//one-pass-the-end size
  E_Int s(nb_pgs0);
  for (E_Int i = 0; i < nb_phs_ref + 1; ++i)
  {
    intpos[i] = s;
    s += nbi;
  }

  F2E.resize(2, nb_pgs0 + nbi * nb_phs_ref, E_IDX_NONE);

  // for the children PH
  // Reserve space for children in the tree
  PHtree.resize(PHlist, nbc);
  // each Hex/ISO is split into 8 Hex
  E_Int nb_phs0(ng.PHs.size());
  ng.PHs.expand_n_fixed_stride(nbc * nb_phs_ref, ELT_t::NB_BOUNDS);

  childpos.resize(nb_phs_ref + 1);//one-pass-the-end size
  s = nb_phs0;
  for (size_t i = 0; i < nb_phs_ref + 1; ++i)
  {
    childpos[i] = s;
    s += nbc;
  }
}


// Impl for TETRA/ISO
template <>
template <typename arr_t>
void refiner<K_MESH::Tetrahedron, ISO>::reserve_mem_PHs
(K_FLD::FloatArray& crd, ngon_type& ng, const Vector_t<E_Int>& PHadap, tree<arr_t> & PGtree, tree<arr_t> & PHtree, K_FLD::IntArray& F2E,
  std::vector<E_Int>& intpos, std::vector<E_Int>& childpos)
{
  refiner<K_MESH::Tetrahedron, ISO>::__reserve_mem_single_bound_type_PHs(crd, ng, PHadap, PGtree, PHtree, F2E, intpos, childpos);
}

// Impl for HEXA
template <>
template <typename arr_t>
void refiner<K_MESH::Hexahedron, ISO>::reserve_mem_PHs
(K_FLD::FloatArray& crd, ngon_type& ng, const Vector_t<E_Int>& PHadap, tree<arr_t> & PGtree, tree<arr_t> & PHtree, K_FLD::IntArray& F2E,
  std::vector<E_Int>& intpos, std::vector<E_Int>& childpos)
{
  refiner<K_MESH::Hexahedron, ISO>::__reserve_mem_single_bound_type_PHs(crd, ng, PHadap, PGtree, PHtree, F2E, intpos, childpos);
  crd.resize(3, crd.cols() + PHadap.size()); // room for centroids
}

// Impl for Prism/ISO
template <>
template <typename arr_t>
void refiner<K_MESH::Prism, ISO>::reserve_mem_PHs
(K_FLD::FloatArray& crd, ngon_type& ng, const Vector_t<E_Int>& PHadap, tree<arr_t> & PGtree, tree<arr_t> & PHtree, K_FLD::IntArray& F2E,
  std::vector<E_Int>& intpos, std::vector<E_Int>& childpos)
{
  using ELT_t = K_MESH::Prism;

  const E_Int &nbc = subdiv_pol<ELT_t, ISO>::PHNBC;
  const E_Int& nbi = subdiv_pol<ELT_t, ISO>::NBI;

  E_Int nb_phs_ref = PHadap.size();
  // internal PGs created (12 per PH)
  // Reserve space for internal faces in the tree
  E_Int current_sz = PGtree.size();
  assert(current_sz != -1);

  //just sync the tree attributes (internals have no parents so do not a call to PGtree.resize)
  PGtree.resize_hierarchy(current_sz + nbi * nb_phs_ref); // //fixme : in hmesh.init ?

  E_Int nb_pgs0 = ng.PGs.size();
  // faces internes : 4 TRI + 6 QUAD
  ng.PGs.expand_n_ab_fixed_stride(nb_phs_ref, 4, 3, 6, 4);

  intpos.resize(nb_phs_ref + 1);//one-pass-the-end size
  E_Int s(nb_pgs0);
  for (size_t i = 0; i < nb_phs_ref + 1; ++i)
  {
    intpos[i] = s;
    s += nbi;
  }

  F2E.resize(2, nb_pgs0 + nbi * nb_phs_ref, E_IDX_NONE);
  PHtree.resize(PHadap, nbc);

  //
  E_Int nb_phs0(ng.PHs.size());
  ng.PHs.expand_n_fixed_stride(nbc * nb_phs_ref, ELT_t::NB_BOUNDS);

  childpos.resize(nb_phs_ref + 1);//one-pass-the-end size
  s = nb_phs0;
  for (size_t i = 0; i < nb_phs_ref + 1; ++i)
  {
    childpos[i] = s;
    s += nbc;
  }
}

// Impl for Prism/ISO
template <>
template <typename arr_t>
void refiner<K_MESH::Pyramid, ISO>::reserve_mem_PHs
(K_FLD::FloatArray& crd, ngon_type& ng, const Vector_t<E_Int>& PHadap, tree<arr_t> & PGtree, tree<arr_t> & PHtree, K_FLD::IntArray& F2E,
  std::vector<E_Int>& intpos, std::vector<E_Int>& childpos)
{
  using ELT_t = K_MESH::Pyramid;

  const E_Int &nbc = subdiv_pol<ELT_t, ISO>::PHNBC;
  const E_Int& nbi = subdiv_pol<ELT_t, ISO>::NBI;

  E_Int nb_phs_ref = PHadap.size();
  // internal PGs created (12 per PH)
  // Reserve space for internal faces in the tree
  E_Int current_sz = PGtree.size();
  assert(current_sz != -1);

  //just sync the tree attributes (internals have no parents so do not a call to PGtree.resize)
  PGtree.resize_hierarchy(current_sz + nbi * nb_phs_ref); //fixme : in hmesh.init ?

  E_Int nb_pgs0 = ng.PGs.size();

  // 13 faces internes : 12 TRI + 1 QUAD
  ng.PGs.expand_n_ab_fixed_stride(nb_phs_ref, 12, 3, 1, 4);

  intpos.resize(nb_phs_ref + 1);//one-pass-the-end size
  E_Int s(nb_pgs0);
  for (size_t i = 0; i < nb_phs_ref + 1; ++i)
  {
    intpos[i] = s;
    s += nbi;
  }

  F2E.resize(2, nb_pgs0 + nbi * nb_phs_ref, E_IDX_NONE);
  //  std::cout << "ng.PG size= " << ng.PGs.size() <<std::endl;
  // for the children PH
  // Reserve space for children in the tree
  PHtree.resize(PHadap, nbc);

  // 10 children : 6 PYRA + 4 TETRA
  E_Int nb_phs0(ng.PHs.size());
  ng.PHs.expand_n_ab_fixed_stride(nb_phs_ref, 6, 5, 4, 4);

  childpos.resize(nb_phs_ref + 1);//one-pass-the-end size
  s = nb_phs0;
  for (size_t i = 0; i < nb_phs_ref + 1; ++i)
  {
    childpos[i] = s;
    s += nbc;
  }
}


// Impl for PH/ISO_HEX
template <>
template <typename arr_t>
void refiner<K_MESH::Polyhedron<0>, ISO_HEX>::reserve_mem_PHs
(K_FLD::FloatArray& crd, ngon_type& ng, const Vector_t<E_Int>& PHadap, tree<arr_t> & PGtree, tree<arr_t> & PHtree, K_FLD::IntArray& F2E,
  std::vector<E_Int>& intpos, std::vector<E_Int>& childpos)
{
  //
  E_Int nb_phs_ref = PHadap.size();

  E_Int nbi_tot = subdiv_pol<K_MESH::Polyhedron<0>, ISO_HEX>::nbi_sum(ng.PHs, PHadap);

  // reserve space for internal faces in the tree
  E_Int current_sz = PGtree.size();
  assert(current_sz != -1);

  // just sync the tree attributes (internals have no parents so do not a call to PGtree.resize)
  PGtree.resize_hierarchy(current_sz + nbi_tot); //fixme : in hmesh.init ?

  E_Int nb_pgs0 = ng.PGs.size();

  // expand PGs for internals QUADS
  ng.PGs.expand_n_fixed_stride(nbi_tot, 4);

  //
  F2E.resize(2, nb_pgs0 + nbi_tot, E_IDX_NONE);

  // reserve space for children in the tree
  std::vector<E_Int> pregnant;
  subdiv_pol<K_MESH::Polyhedron<0>, ISO_HEX>::nbc_list(ng.PHs, PHadap, pregnant);

  PHtree.resize(PHadap, pregnant);

  // expand ng.PHs for children with total nb of HEXAS
  // todo JP

  // allocate space in coordinates for centroids
  crd.resize(3, crd.cols() + nb_phs_ref);
}

/// default implementation : T3/Q4-ISO, PG-ISO_HEX cases
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

  // Reserve space for children in the tree
  std::vector<E_Int> childpos;
  reserve_mem_PGs(crd, ng.PGs, PGref, PGtree, F2E, childpos);

  // Refine them
#ifndef DEBUG_HIERARCHICAL_MESH
//#pragma omp parallel for
#endif
  for (E_Int i = 0; i < nb_pgs_ref; ++i)
  {
    E_Int PGi = PGref[i];
    E_Int firstChild = childpos[i];
    
    refine_PG(crd, ng.PGs, PGi, pos + i, firstChild, ecenter);

    E_Int nbc = childpos[i+1]- firstChild;

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

///
template <>
template <typename arr_t>
void refiner<K_MESH::Quadrangle, eSUBDIV_TYPE::DIR>::refine_PGs
(const sensor_output_t& adap_incr, ngon_type& ng, tree<arr_t>& PGtree, K_FLD::FloatArray& crd,
  K_FLD::IntArray& F2E, std::map<K_MESH::NO_Edge, E_Int>& ecenter)
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
  std::map<K_MESH::NO_Edge, E_Int>& ecenter)
{
  // do not handle DIR for Triangle yet
  //refiner<K_MESH::Triangle, ISO>::refine_PGs<arr_t>(adap_incr._adap_incr, ng, PGtree, crd, F2E, ecenter);
}

///
template<>
void refiner<K_MESH::Quadrangle, ISO>::refine_PG
(K_FLD::FloatArray& crd, ngon_unit& PGs, E_Int PGi, E_Int posC, E_Int posChild, std::map<K_MESH::NO_Edge, E_Int>& ecenter)
{
  E_Int* nodes = PGs.get_facets_ptr(PGi);
  const E_Int &nbc = subdiv_pol<K_MESH::Quadrangle, ISO>::NBC;

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
  const E_Int &nbc = subdiv_pol<K_MESH::Triangle, ISO>::NBC;

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
template<>
void refiner<K_MESH::Polygon, ISO_HEX>::refine_PG
(K_FLD::FloatArray& crd, ngon_unit& PGs, E_Int PGi, E_Int posC, E_Int posChild, std::map<K_MESH::NO_Edge, E_Int>& ecenter)
{
  E_Int* nodes = PGs.get_facets_ptr(PGi);
  E_Int nnodes = PGs.stride(PGi);

  // Centroid calculation
  NUGA::refine_point_computer<K_MESH::Polygon>::compute_center(crd, nodes, K_MESH::Quadrangle::NB_NODES, 1, crd.col(posC));

  //refining face structure : nodes + mid points + centroid
  E_Int nrefnodes = nnodes * 2 + 1;
  STACK_ARRAY(E_Int, nrefnodes, refE); // emulate E_Int refE[refnodes]

  //todo JP : remplir refE. A toi de choisir une convention d'ordre de stockage coherente avec l'appel qui suit à K_MESH::Polygon::split,
  // fonction que tu dois definir egalement.
  // si le dernier noeud est le centroid, faire : refE[nrefnodes-1]= posC + 1 

  K_MESH::Polygon::split<ISO_HEX>(refE.get(), nrefnodes, PGs, posChild);
}

/// default implementation : ISO case , for all basic element types (Basic, Tetra, Pyra, Penta, Hexa)
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

// Polyhedron/ISO_HEX impl
template <>
template <typename arr_t>
void refiner<K_MESH::Polyhedron<0>, ISO_HEX>::refine_Faces
(const sensor_output_t &adap_incr,
  ngon_type& ng, tree<arr_t> & PGtree, K_FLD::FloatArray& crd, K_FLD::IntArray & F2E,
  std::map<K_MESH::NO_Edge, E_Int>& ecenter)
{
  refiner<K_MESH::Polygon, ISO_HEX>::refine_PGs(adap_incr, ng, PGtree, crd, F2E, ecenter);
}

/// HEXA/TETRA/PYRA/PRISM ISO, PH-ISO_HEX impl
template <typename ELT_t, eSUBDIV_TYPE STYPE>
template <typename arr_t>
void refiner<ELT_t, STYPE>::refine_PHs
(const sensor_output_t &adap_incr, ngon_type& ng, tree<arr_t> & PGtree, tree<arr_t> & PHtree,
  K_FLD::FloatArray& crd, K_FLD::IntArray & F2E)
{
  assert(adap_incr.cell_adap_incr.size() <= ng.PHs.size());
  //
  std::vector<E_Int> PHadap;
  get_PHs_to_refine(ng, PHtree, adap_incr, PHadap);
  E_Int nb_phs_ref = PHadap.size();
  if (nb_phs_ref == 0) return;

  E_Int pos = crd.cols(); // first cell center in crd (if relevant)
  E_Int nb_pgs0 = ng.PGs.size();
  E_Int nb_phs0 = ng.PHs.size();
  // Reserve space for children in the tree
  std::vector<E_Int> intpos, childpos;
  reserve_mem_PHs(crd, ng, PHadap, PGtree, PHtree, F2E, intpos, childpos);

  using spliting_t = splitting_t<ELT_t, STYPE, 1>; // HX27, TH10, PY13 or PR18

#ifndef DEBUG_HIERARCHICAL_MESH
//#pragma omp parallel for
#endif
  for (E_Int i = 0; i < nb_phs_ref; ++i)
  {
    E_Int PHi = PHadap[i];
   
    spliting_t elt(crd, ng, PHi, pos+i, F2E, PGtree);
 
    elt.split(ng, PHi, PHtree, PGtree, F2E, intpos[i], childpos[i]);
  }  
}

///
template <>
template <typename arr_t>
void refiner<K_MESH::Basic, eSUBDIV_TYPE::ISO>::refine_PHs
(const sensor_output_t &adap_incr, 
 ngon_type& ng, tree<arr_t> & PGtree, tree<arr_t> & PHtree, K_FLD::FloatArray& crd, K_FLD::IntArray & F2E)
{  
  assert(adap_incr.cell_adap_incr.size() <= ng.PHs.size());

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
  //refiner<K_MESH::Tetrahedron, ISO>::refine_PHs(adap_incr.cell_adap_incr, ng, PGtree, PHtree, crd, F2E);
}

///
template <>
template <typename arr_t>
void refiner<K_MESH::Pyramid, DIR>::refine_PHs
(const sensor_output_t &adap_incr,
 ngon_type& ng, tree<arr_t> & PGtree, tree<arr_t> & PHtree, K_FLD::FloatArray& crd, K_FLD::IntArray & F2E)
{
  //refiner<K_MESH::Pyramid, ISO>::refine_PHs(adap_incr.cell_adap_incr, ng, PGtree, PHtree, crd, F2E);
}

///
template <>
template <typename arr_t>
void refiner<K_MESH::Prism, DIR>::refine_PHs
(const sensor_output_t &adap_incr,
 ngon_type& ng, tree<arr_t> & PGtree, tree<arr_t> & PHtree, K_FLD::FloatArray& crd, K_FLD::IntArray & F2E)
{
  //refiner<K_MESH::Prism, ISO>::refine_PHs(adap_incr.cell_adap_incr, ng, PGtree, PHtree, crd, F2E);
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

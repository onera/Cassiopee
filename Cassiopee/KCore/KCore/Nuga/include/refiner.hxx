/*    
    Copyright 2013-2025 Onera.

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

#ifndef REFINER_HXX
#define REFINER_HXX

#include "Nuga/include/Edge.h"
#include "Nuga/include/Basic.h"
#include "T6.hxx"
#include "Q9.hxx"
#include "Q6.hxx"
#include "HX27.hxx"
#include "TH10.hxx"
#include "PR18.hxx"
#include "PY13.hxx"
#include "PHQ4.hxx"
#include "HX18.hxx"
#include "HX12.hxx"

#include "sensor.hxx"

namespace NUGA
{

  ///
  class refine_point_computer
  {
  private:
    static void refine_edge(std::map<K_MESH::NO_Edge, E_Int>& ecenter, E_Int n1, E_Int n2,
      K_MESH::NO_Edge e, K_FLD::FloatArray& crd)
    {
      e.setNodes(n1, n2);
      auto it = ecenter.find(e);
      if (it == ecenter.end()) {
        E_Float mid[3];
        NUGA::sum<3>(0.5, crd.col(n1 - 1), 0.5, crd.col(n2 - 1), mid);
        crd.pushBack(mid, mid + 3);
        ecenter[e] = crd.cols();
      }
    }

  public:
        static void compute_centers(const Vector_t<E_Int> &PGlist,
            const Vector_t<NUGA::eDIR> &PG_directive,
            ngon_type& ng, K_FLD::FloatArray& crd,
            std::map<K_MESH::NO_Edge, E_Int> & ecenter)
        {
            K_MESH::NO_Edge e;
            E_Int stride, n1, n2;
            E_Int PGi;

            for (size_t i = 0; i < PGlist.size(); i++) {
                PGi = PGlist[i];
                stride = ng.PGs.stride(PGi);
                const E_Int* pN = ng.PGs.get_facets_ptr(PGlist[i]);
                if (PG_directive[i] == XY) {
                    for (E_Int j = 0; j < stride; j++) {
                        n1 = *(pN + j);
                        n2 = *(pN + (j + 1) % stride);
                        refine_edge(ecenter, n1, n2, e, crd);
                    }
                } else if (PG_directive[i] == Xd) {
                    n1 = pN[0]; n2 = pN[1];
                    refine_edge(ecenter, n1, n2, e, crd);
                    n1 = pN[2]; n2 = pN[3];
                    refine_edge(ecenter, n1, n2, e, crd);
                } else if (PG_directive[i] == Y) {
                    n1 = pN[1]; n2 = pN[2];
                    refine_edge(ecenter, n1, n2, e, crd);
                    n1 = pN[3]; n2 = pN[0];
                    refine_edge(ecenter, n1, n2, e, crd);
                }
            }
        }
    };

   ///
  template <typename ELT_t, eSUBDIV_TYPE STYPE>
  class refiner
  {
  public:
    using output_t = adap_incr_type<STYPE>;

  public:
    std::map<K_MESH::NO_Edge, E_Int> _ecenter;

  public:

    // Generic wrapper for cells (calls relevant(s) refine_PGs)
    template <typename arr_t>
    static void refine_Faces(const output_t &adap_incr,
                             ngon_type& ng, tree<arr_t> & PGtree, K_FLD::FloatArray& crd, K_FLD::IntArray & F2E,
                             std::map<K_MESH::NO_Edge, E_Int>& ecenter);

    /// Generic wrapper for cells
    template <typename pg_arr_t, typename ph_arr_t>
    static void refine_PHs(output_t &adap_incr, ngon_type& ng, tree<pg_arr_t> & PGtree, tree<ph_arr_t> & PHtree,
                           K_FLD::FloatArray& crd, K_FLD::IntArray & F2E);

    // T3/Q4/PG
    template <typename arr_t>
    static void refine_PGs(Vector_t<E_Int>& PG_to_ref, Vector_t<NUGA::eDIR> & PG_directive,
                           ngon_type& ng, tree<arr_t> & PGtree, K_FLD::FloatArray& crd, K_FLD::IntArray & F2E,
                           std::map<K_MESH::NO_Edge, E_Int>& ecenter);

    static void refine_PG(K_FLD::FloatArray& crd, ngon_unit& PGs, E_Int PGi, E_Int posC, E_Int posChild, eDIR d,
                          const std::map<K_MESH::NO_Edge, E_Int>& ecenter);



    //private:
    template <typename pg_arr_t, typename ph_arr_t>
    static void __refine_PHS_hexa_dir (output_t &adap_incr,
                                ngon_type& ng, tree<pg_arr_t> & PGtree, tree<ph_arr_t> & PHtree, K_FLD::FloatArray& crd, K_FLD::IntArray & F2E);

    ///
    template <typename arr_t>
    static E_Int extract_PGs_to_refine(K_FLD::FloatArray& crd/*hack for CLEF*/, const ngon_type& ng, const tree<arr_t> & PGtree, const output_t &adap_incr, Vector_t<E_Int>& PG_to_ref, Vector_t<NUGA::eDIR> & PG_directive);
    ///
    ///
    template <typename arr_t>
    static void extract_PHs_to_refine(const ngon_type& ng, const tree<arr_t> & PHtree, const output_t &adap_incr, Vector_t<E_Int> & PH_to_ref, Vector_t<NUGA::eDIR> & PH_directive);
    ///
    template <typename arr_t>
    static void reserve_mem_PGs(K_FLD::FloatArray& crd, ngon_unit& PGs, const Vector_t<E_Int>& PG_to_ref, const Vector_t<NUGA::eDIR> & PG_directive, tree<arr_t> & PGtree, K_FLD::IntArray& F2E, std::vector<E_Int>& childpos, std::vector<E_Int>& crdpos);
    ///
    template <typename pg_arr_t, typename ph_arr_t>
    static void reserve_mem_PHs(K_FLD::FloatArray& crd, ngon_type& ng, const Vector_t<E_Int> & PH_to_ref, const Vector_t<NUGA::eDIR> & PH_directive, tree<pg_arr_t> & PGtree, tree<ph_arr_t> & PHtree,
      K_FLD::IntArray& F2E, std::vector<E_Int>& intpos, std::vector<E_Int>& childpos);

    ///
    template <typename pg_arr_t, typename ph_arr_t>
    static void __reserve_mem_single_bound_type_PHs
    (K_FLD::FloatArray& crd, ngon_type& ng, const Vector_t<E_Int> & PH_to_ref, const Vector_t<NUGA::eDIR> & PH_directive, tree<pg_arr_t> & PGtree, tree<ph_arr_t> & PHtree,
      K_FLD::IntArray& F2E, std::vector<E_Int>& intpos, std::vector<E_Int>& childpos);

  };

  /// Impl for basic elements
  template <typename ELT_t, eSUBDIV_TYPE STYPE>
  template <typename arr_t>
  E_Int refiner<ELT_t, STYPE>::extract_PGs_to_refine
  (K_FLD::FloatArray& crd/*hack for CLEF*/, const ngon_type& ng, const tree<arr_t> & PGtree, const output_t &adap_incr, Vector_t<E_Int>& PG_to_ref, Vector_t<NUGA::eDIR> & PG_directive)
  {
    PG_to_ref.clear();
    PG_directive.clear();

    E_Int nb_phs = ng.PHs.size();
    // Gets PGs to refine
    E_Int nb_pgs(ng.PGs.size());
    Vector_t<bool> is_PG_to_refine(nb_pgs, false);

    // disable non-admissible elements (face_adap_incr may contain info fo Q4 & T3)
    // IMPORTANT : need to be over all PGS (not only those belonging to cells to refine) 
    for (E_Int i = 0; i < nb_pgs; ++i)
    {
      if (ng.PGs.stride(i) != ELT_t::NB_NODES) continue;

      ASSERT_IN_VECRANGE(adap_incr.face_adap_incr, i)

      if (adap_incr.face_adap_incr[i] != 0) // IMPORTANT : CAN HAVE SOMETHING UPON ENTRY (e.g. WHEN SYNCING JOINS)
        is_PG_to_refine[i] = true;
    }

    // select cell-to-refine faces, but discard those having already children and not being of input type (Tri or Quad)
    for (E_Int i = 0; i < nb_phs; ++i)
    {
      ASSERT_IN_VECRANGE(adap_incr.cell_adap_incr, i)
      
      if (!(adap_incr.cell_adap_incr[i] > 0)) continue;

      E_Int nb_faces = ng.PHs.stride(i);
      const E_Int* faces = ng.PHs.get_facets_ptr(i);

      for (E_Int j = 0; j < nb_faces; ++j)
      {
        E_Int PGi = *(faces + j) - 1;

        if (ng.PGs.stride(PGi) != ELT_t::NB_NODES) continue;
        if (PGtree.nb_children(PGi) != 0) continue;

        if (STYPE != DIR) is_PG_to_refine[PGi] = true;
      }
    }

    E_Int nb_pgs_ref = std::count(is_PG_to_refine.begin(), is_PG_to_refine.end(), true);

    if (nb_pgs_ref == 0) return 0;

    PG_to_ref.reserve(nb_pgs_ref);
    PG_directive.reserve(nb_pgs_ref);

    for (E_Int i = 0; i < nb_pgs; ++i)
    {
      // 1. adap_incr.get_face_dir(i)  must return a value
      if (is_PG_to_refine[i])
      {
        PG_to_ref.push_back(i);
        
        auto fdir = XY;
        
        if (STYPE == DIR_PROTO)
        {
          E_Int nnodes = ng.PGs.stride(i);
          const E_Int* nodes = ng.PGs.get_facets_ptr(i);

          fdir = get_dir(crd, nodes, nnodes);
        }
        else if (STYPE == DIR) {
          fdir = adap_incr.get_face_dir(i);
        }

        assert(fdir != NONE);

        PG_directive.push_back(fdir);
      }
    }

    return PG_to_ref.size();
  }

  ///
  template <>
  template <typename arr_t>
  E_Int refiner<K_MESH::Polygon, ISO_HEX>::extract_PGs_to_refine
  (K_FLD::FloatArray& dummy/*hack for CLEF*/, const ngon_type& ng, const tree<arr_t> & PGtree, const output_t &adap_incr, Vector_t<E_Int>& PG_to_ref, Vector_t<NUGA::eDIR> & PG_directive)
  {
    PG_to_ref.clear();
    PG_directive.clear();

    //
    E_Int nb_pgs(ng.PGs.size());
    auto face_adap_incr = adap_incr.face_adap_incr; // copy : input is preserved
    face_adap_incr.resize(nb_pgs, 0);

    //
    E_Int nb_phs = ng.PHs.size();
    for (E_Int i = 0; i < nb_phs; ++i)
    {
      if (adap_incr.cell_adap_incr[i] <= 0) continue;

      E_Int nb_faces = ng.PHs.stride(i);             // this i-th PH is to be refined
      const E_Int* faces = ng.PHs.get_facets_ptr(i);

      for (E_Int j = 0; j < nb_faces; ++j)           // browse its faces
      {
        E_Int PGi = *(faces + j) - 1;
        face_adap_incr[PGi] = (PGtree.nb_children(PGi) == 0);
      }
    }

    for (size_t k = 0; k < face_adap_incr.size(); ++k)
      if (face_adap_incr[k] != 0)
      {
        PG_to_ref.push_back(k);
        PG_directive.push_back(XY);
      }
    return PG_to_ref.size();
  }


  ///
  template <typename ELT_t, eSUBDIV_TYPE STYPE>
  template <typename arr_t>
  void refiner<ELT_t, STYPE>::extract_PHs_to_refine
  (const ngon_type& ng, const tree<arr_t> & PHtree, const output_t &adap_incr, Vector_t<E_Int> & PH_to_ref, Vector_t<NUGA::eDIR> & PH_directive)
  {
    for (size_t i = 0; i < adap_incr.cell_adap_incr.size(); ++i)
    {
      if (adap_incr.cell_adap_incr[i] > 0 && (PHtree.nb_children(i) == 0) && (ELT_t::is_of_type(ng.PGs, ng.PHs.get_facets_ptr(i), ng.PHs.stride(i))))
      {
        PH_to_ref.push_back(i);
        PH_directive.push_back(XYZ/*adap_incr.cell_adap_incr[i]*/);
      }
    }
  }

  ///
  template <>
  template <typename arr_t>
  void refiner<K_MESH::Hexahedron, DIR>::extract_PHs_to_refine
  (const ngon_type& ng, const tree<arr_t> & PHtree, const output_t &adap_incr, Vector_t<E_Int> & PH_to_ref, Vector_t<NUGA::eDIR> & PH_directive)
  {
    for (size_t i = 0; i < adap_incr.cell_adap_incr.size(); ++i)
      if (adap_incr.cell_adap_incr[i] > 0 && (PHtree.nb_children(i) == 0) && (K_MESH::Hexahedron::is_of_type(ng.PGs, ng.PHs.get_facets_ptr(i), ng.PHs.stride(i))))
      {
        PH_to_ref.push_back(i);
        
        eDIR d = NONE;
        auto v = adap_incr.cell_adap_incr[i];

        int count = 0;
        if (v.n[0]) ++count;
        if (v.n[1]) ++count;
        if (v.n[2]) ++count;

        assert(count != 0);

        if (count == 3) d = XYZ;
        else if (count == 2) d = XY;
        else if (count == 1) d = Xd;
 
        PH_directive.push_back(d);
      }
  }

  ///
  template <>
  template <typename arr_t>
  void refiner<K_MESH::Hexahedron, DIR_PROTO>::extract_PHs_to_refine
  (const ngon_type& ng, const tree<arr_t> & PHtree, const output_t &adap_incr, Vector_t<E_Int> & PH_to_ref, Vector_t<NUGA::eDIR> & PH_directive)
  {
    refiner<K_MESH::Hexahedron, DIR>::extract_PHs_to_refine(ng, PHtree, adap_incr, PH_to_ref, PH_directive);
  }

  ///
  template <>
  template <typename arr_t>
  void refiner<K_MESH::Polyhedron<0>, ISO_HEX>::extract_PHs_to_refine
  (const ngon_type& ng, const tree<arr_t> & PHtree, const output_t &adap_incr, Vector_t<E_Int> & PH_to_ref, Vector_t<NUGA::eDIR> & PH_directive)
  {
    for (size_t i = 0; i < adap_incr.cell_adap_incr.size(); ++i)
    {
      if (adap_incr.cell_adap_incr[i] > 0 && PHtree.is_enabled(i))
      {
        PH_to_ref.push_back(i);
      }
    }
  }


  /// Impl for  T3/IS
  template <>
  template <typename arr_t>
  void refiner<K_MESH::Triangle, ISO>::reserve_mem_PGs
  (K_FLD::FloatArray& crd, ngon_unit& PGs, const Vector_t<E_Int>& PGlist, const Vector_t<NUGA::eDIR> & PG_directive, tree<arr_t> & PGtree, K_FLD::IntArray& F2E, std::vector<E_Int>&childpos, std::vector<E_Int>& dummy)
  {
    const E_Int &nbc = subdiv_pol<K_MESH::Triangle, ISO>::NBC;
    E_Int nb_pgs_ref = PGlist.size();

    PGtree.resize(PGlist, nbc);

    E_Int nb_pgs0 = PGs.size();

    // And in the mesh : each ELT(T3) is split into nbc
    PGs.expand_n_fixed_stride(nbc * nb_pgs_ref, K_MESH::Triangle::NB_NODES);

    childpos.resize(nb_pgs_ref + 1);//one-pass-the-end size
    E_Int s(nb_pgs0);
    for (E_Int i = 0; i < nb_pgs_ref + 1; ++i)
    {
      childpos[i] = s;
      s += nbc;
    }

    F2E.resize(2, nb_pgs0 + nbc * nb_pgs_ref, IDX_NONE);
  }

  /// Impl for Q4/ISO
  template <>
  template <typename arr_t>
  void refiner<K_MESH::Quadrangle, ISO>::reserve_mem_PGs
  (K_FLD::FloatArray& crd, ngon_unit& PGs, const Vector_t<E_Int>& PGlist, const Vector_t<NUGA::eDIR> & PG_directive, tree<arr_t> & PGtree, K_FLD::IntArray& F2E, std::vector<E_Int>&childpos, std::vector<E_Int>& crdpos)
  {
    const E_Int &nbc = subdiv_pol<K_MESH::Quadrangle, ISO>::NBC;
    E_Int nb_pgs_ref = PGlist.size();

    PGtree.resize(PGlist, nbc);

    E_Int nb_pgs0 = PGs.size();

    // And in the mesh : each ELT(T3/) is split into nbc
    PGs.expand_n_fixed_stride(nbc * nb_pgs_ref, K_MESH::Quadrangle::NB_NODES);

    childpos.resize(nb_pgs_ref + 1);//one-pass-the-end size
    E_Int s(nb_pgs0);
    for (E_Int i = 0; i < nb_pgs_ref + 1; ++i)
    {
      childpos[i] = s;
      s += nbc;
    }

    F2E.resize(2, nb_pgs0 + nbc * nb_pgs_ref, IDX_NONE);

    // reserve space for face centers
    E_Int pos0 = crd.cols();
    crdpos.resize(nb_pgs_ref, IDX_NONE);
    crd.resize(3, crd.cols() + nb_pgs_ref);
    for (size_t i = 0; i < crdpos.size(); ++i)
      crdpos[i] = pos0++;
  }


  /// Impl for PG/ISO_HEX
  template <>
  template <typename arr_t>
  void refiner<K_MESH::Polygon, ISO_HEX>::reserve_mem_PGs
  (K_FLD::FloatArray& crd, ngon_unit& PGs, const Vector_t<E_Int>& PGlist, const Vector_t<NUGA::eDIR> & PG_directive, tree<arr_t> & PGtree, K_FLD::IntArray& F2E, std::vector<E_Int>& childpos, std::vector<E_Int>& crdpos)
  {
    // info : nbc is variable and equal to number of nodes

    //E_Int nb_pts_ref = crd.cols();
    E_Int nb_pgs_ref = PGlist.size();
    E_Int nb_pgs0 = PGs.size();

    std::vector<E_Int> pregnant;
    subdiv_pol<K_MESH::Polygon, ISO_HEX>::nbc_list(PGs, PGlist, PG_directive, pregnant); // pregnant[i] is number of children on PGlist[i]

    PGtree.resize(PGlist, pregnant);

    // expand ng.PGs for children with total nb of QUADS
    E_Int nb_tot_child = 0;
    for (E_Int i = 0; i < nb_pgs_ref; ++i)
      nb_tot_child += pregnant[i];
    PGs.expand_n_fixed_stride(nb_tot_child, 4);

    childpos.resize(nb_pgs_ref + 1);//one-pass-the-end size
    E_Int s(nb_pgs0);
    for (E_Int i = 0; i < nb_pgs_ref; ++i)
    {
      childpos[i] = s; // numero du polygone ou commence les enfants
      s += pregnant[i];
    }
    childpos[nb_pgs_ref] = s; // s is the total number of faces including new 
    F2E.resize(2, s, IDX_NONE);

    // reserve space for face centers
    E_Int pos0 = crd.cols();
    crdpos.resize(nb_pgs_ref, IDX_NONE);
    crd.resize(3, crd.cols() + nb_pgs_ref); // reserve space for external face centers
    for (size_t i = 0; i < crdpos.size(); ++i)
      crdpos[i] = pos0++;
  }

  /// Impl for QUAD/DIR
  template <>
  template <typename arr_t>
  void refiner<K_MESH::Quadrangle, DIR>::reserve_mem_PGs
  (K_FLD::FloatArray& crd, ngon_unit& PGs, const Vector_t<E_Int>& PGlist, const Vector_t<NUGA::eDIR> & PG_directive, tree<arr_t> & PGtree, K_FLD::IntArray& F2E, std::vector<E_Int>& childpos, std::vector<E_Int>& crdpos)
  {
    // info : nbc is variable for directionnal : either 2 or 4

    E_Int nb_pgs_ref = PGlist.size();
    E_Int nb_pgs0 = PGs.size();

    std::vector<E_Int> pregnant;
    subdiv_pol<K_MESH::Quadrangle, DIR>::nbc_list(PGs, PGlist, PG_directive, pregnant);

    PGtree.resize(PGlist, pregnant);

    // expand ng.PGs for children with total nb of QUADS
    E_Int nb_tot_child = 0;
    for (E_Int i = 0; i < nb_pgs_ref; ++i)
      nb_tot_child += pregnant[i];
    PGs.expand_n_fixed_stride(nb_tot_child, 4);

    childpos.resize(nb_pgs_ref + 1);//one-pass-the-end size
    E_Int s(nb_pgs0);
    for (E_Int i = 0; i < nb_pgs_ref; ++i)
    {
      childpos[i] = s;
      s += pregnant[i];
    }
    childpos[nb_pgs_ref] = s;

    F2E.resize(2, s, IDX_NONE);

    // add centers for ISO
    crdpos.resize(nb_pgs_ref, IDX_NONE);
    E_Int pos0 = crd.cols();
    for (E_Int i = 0; i < nb_pgs_ref; ++i)
    {
      if (pregnant[i] == 4) // ISO : need a center
        crdpos[i] = pos0++;
    }

    // reserve space for face centers
    crd.resize(3, pos0);
  }

  /// Impl for TRI/DIR
  template <>
  template <typename arr_t>
  void refiner<K_MESH::Triangle, DIR>::reserve_mem_PGs
  (K_FLD::FloatArray& crd, ngon_unit& PGs, const Vector_t<E_Int>& PGlist, const Vector_t<NUGA::eDIR> & PG_directive, tree<arr_t> & PGtree, K_FLD::IntArray& F2E, std::vector<E_Int>& childpos, std::vector<E_Int>& crdpos)
  {
    //todo
  }
  template <>
  template <typename arr_t>
  void refiner<K_MESH::Triangle, DIR_PROTO>::reserve_mem_PGs
  (K_FLD::FloatArray& crd, ngon_unit& PGs, const Vector_t<E_Int>& PGlist, const Vector_t<NUGA::eDIR> & PG_directive, tree<arr_t> & PGtree, K_FLD::IntArray& F2E, std::vector<E_Int>& childpos, std::vector<E_Int>& crdpos)
  {
    //todo
  }

  template <>
  template <typename arr_t>
  void refiner<K_MESH::Quadrangle, DIR_PROTO>::reserve_mem_PGs
  (K_FLD::FloatArray& crd, ngon_unit& PGs, const Vector_t<E_Int>& PGlist, const Vector_t<NUGA::eDIR> & PG_directive, tree<arr_t> & PGtree, K_FLD::IntArray& F2E, std::vector<E_Int>& childpos, std::vector<E_Int>& crdpos)
  {
    refiner<K_MESH::Quadrangle, DIR>::reserve_mem_PGs (crd, PGs, PGlist, PG_directive, PGtree, F2E, childpos, crdpos);
  }

  // used for HEXA or TETRA
  template <typename ELT_t, eSUBDIV_TYPE STYPE>
  template <typename pg_arr_t, typename ph_arr_t>
  void refiner<ELT_t, STYPE>::__reserve_mem_single_bound_type_PHs
  (K_FLD::FloatArray& crd, ngon_type& ng, const Vector_t<E_Int> & PH_to_ref, const Vector_t<NUGA::eDIR> & PH_directive, tree<pg_arr_t> & PGtree, tree<ph_arr_t> & PHtree,
    K_FLD::IntArray& F2E, std::vector<E_Int>& intpos, std::vector<E_Int>& childpos)
  {
    const E_Int &nbc = subdiv_pol<ELT_t, ISO>::PHNBC;
    const E_Int& nbi = subdiv_pol<ELT_t, ISO>::NBI;

    using face_t = typename ELT_t::boundary_type;

    E_Int nb_phs_ref = PH_to_ref.size();

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

    F2E.resize(2, nb_pgs0 + nbi * nb_phs_ref, IDX_NONE);

    // for the children PH
    // Reserve space for children in the tree
    PHtree.resize(PH_to_ref, nbc);
    // each Hex/ISO is split into 8 Hex
    E_Int nb_phs0(ng.PHs.size());
    ng.PHs.expand_n_fixed_stride(nbc * nb_phs_ref, ELT_t::NB_BOUNDS);

    childpos.resize(nb_phs_ref + 1);//one-pass-the-end size
    s = nb_phs0;
    for (E_Int i = 0; i < nb_phs_ref + 1; ++i)
    {
      childpos[i] = s;
      s += nbc;
    }
  }


  // Impl for TETRA/ISO
  template <>
  template <typename pg_arr_t, typename ph_arr_t>
  void refiner<K_MESH::Tetrahedron, ISO>::reserve_mem_PHs
  (K_FLD::FloatArray& crd, ngon_type& ng, const Vector_t<E_Int> & PH_to_ref, const Vector_t<NUGA::eDIR> & PH_directive, tree<pg_arr_t> & PGtree, tree<ph_arr_t> & PHtree, K_FLD::IntArray& F2E,
    std::vector<E_Int>& intpos, std::vector<E_Int>& childpos)
  {
    refiner<K_MESH::Tetrahedron, ISO>::__reserve_mem_single_bound_type_PHs(crd, ng, PH_to_ref, PH_directive, PGtree, PHtree, F2E, intpos, childpos);
  }

  // Impl for HEXA
  template <>
  template <typename pg_arr_t, typename ph_arr_t>
  void refiner<K_MESH::Hexahedron, ISO>::reserve_mem_PHs
  (K_FLD::FloatArray& crd, ngon_type& ng, const Vector_t<E_Int> & PH_to_ref, const Vector_t<NUGA::eDIR> & PH_directive, tree<pg_arr_t> & PGtree, tree<ph_arr_t> & PHtree, K_FLD::IntArray& F2E,
    std::vector<E_Int>& intpos, std::vector<E_Int>& childpos)
  {
    refiner<K_MESH::Hexahedron, ISO>::__reserve_mem_single_bound_type_PHs(crd, ng, PH_to_ref, PH_directive, PGtree, PHtree, F2E, intpos, childpos);
    crd.resize(3, crd.cols() + PH_to_ref.size()); // room for centroids
  }

  // Impl for Prism/ISO
  template <>
  template <typename pg_arr_t, typename ph_arr_t>
  void refiner<K_MESH::Prism, ISO>::reserve_mem_PHs
  (K_FLD::FloatArray& crd, ngon_type& ng, const Vector_t<E_Int> & PH_to_ref, const Vector_t<NUGA::eDIR> & PH_directive,
    tree<pg_arr_t> & PGtree, tree<ph_arr_t> & PHtree, K_FLD::IntArray& F2E,
    std::vector<E_Int>& intpos, std::vector<E_Int>& childpos)
  {
    using ELT_t = K_MESH::Prism;

    const E_Int &nbc = subdiv_pol<ELT_t, ISO>::PHNBC;
    const E_Int& nbi = subdiv_pol<ELT_t, ISO>::NBI;

    E_Int nb_phs_ref = PH_to_ref.size();
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
    for (E_Int i = 0; i < nb_phs_ref + 1; ++i)
    {
      intpos[i] = s;
      s += nbi;
    }

    F2E.resize(2, nb_pgs0 + nbi * nb_phs_ref, IDX_NONE);
    PHtree.resize(PH_to_ref, nbc);

    //
    E_Int nb_phs0(ng.PHs.size());
    ng.PHs.expand_n_fixed_stride(nbc * nb_phs_ref, ELT_t::NB_BOUNDS);

    childpos.resize(nb_phs_ref + 1);//one-pass-the-end size
    s = nb_phs0;
    for (E_Int i = 0; i < nb_phs_ref + 1; ++i)
    {
      childpos[i] = s;
      s += nbc;
    }
  }

  // Impl for Prism/ISO
  template <>
  template <typename pg_arr_t, typename ph_arr_t>
  void refiner<K_MESH::Pyramid, ISO>::reserve_mem_PHs
  (K_FLD::FloatArray& crd, ngon_type& ng, const Vector_t<E_Int> & PH_to_ref, const Vector_t<NUGA::eDIR> & PH_directive,
    tree<pg_arr_t> & PGtree, tree<ph_arr_t> & PHtree, K_FLD::IntArray& F2E,
    std::vector<E_Int>& intpos, std::vector<E_Int>& childpos)
  {
    using ELT_t = K_MESH::Pyramid;

    const E_Int &nbc = subdiv_pol<ELT_t, ISO>::PHNBC;
    const E_Int& nbi = subdiv_pol<ELT_t, ISO>::NBI;

    E_Int nb_phs_ref = PH_to_ref.size();
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
    for (E_Int i = 0; i < nb_phs_ref + 1; ++i)
    {
      intpos[i] = s;
      s += nbi;
    }

    F2E.resize(2, nb_pgs0 + nbi * nb_phs_ref, IDX_NONE);
    //  std::cout << "ng.PG size= " << ng.PGs.size() <<std::endl;
    // for the children PH
    // Reserve space for children in the tree
    PHtree.resize(PH_to_ref, nbc);

    // 10 children : 6 PYRA + 4 TETRA
    E_Int nb_phs0(ng.PHs.size());
    ng.PHs.expand_n_ab_fixed_stride(nb_phs_ref, 6, 5, 4, 4);

    childpos.resize(nb_phs_ref + 1);//one-pass-the-end size
    s = nb_phs0;
    for (E_Int i = 0; i < nb_phs_ref + 1; ++i)
    {
      childpos[i] = s;
      s += nbc;
    }
  }


  // Impl for PH/ISO_HEX
  template <>
  template <typename pg_arr_t, typename ph_arr_t>
  void refiner<K_MESH::Polyhedron<0>, ISO_HEX>::reserve_mem_PHs
  (K_FLD::FloatArray& crd, ngon_type& ng, const Vector_t<E_Int> & PH_to_ref, const Vector_t<NUGA::eDIR> & PH_directive,
    tree<pg_arr_t> & PGtree, tree<ph_arr_t> & PHtree, K_FLD::IntArray& F2E,
    std::vector<E_Int>& intpos, std::vector<E_Int>& childpos)
  {
    //
    E_Int nb_phs_ref = PH_to_ref.size();

    E_Int nbi_tot = subdiv_pol<K_MESH::Polyhedron<0>, ISO_HEX>::nbi_sum(ng, PH_to_ref);

    // reserve space for internal faces in the tree
    E_Int current_sz = PGtree.size();
    assert(current_sz != -1);

    // just sync the tree attributes (internals have no parents so do not a call to PGtree.resize)
    PGtree.resize_hierarchy(current_sz + nbi_tot); //fixme : in hmesh.init ?

    E_Int nb_pgs0 = ng.PGs.size();

    // expand PGs for internals QUADS
    ng.PGs.expand_n_fixed_stride(nbi_tot, 4);

    //
    F2E.resize(2, nb_pgs0 + nbi_tot, IDX_NONE);

    // reserve space for children in the tree
    std::vector<E_Int> pregnant;
    subdiv_pol<K_MESH::Polyhedron<0>, ISO_HEX>::nbc_list(ng, PH_to_ref, pregnant);

    PHtree.resize(PH_to_ref, pregnant);

    // expand ng.PHs for children with total nb of HEXAS

    E_Int nb_phs0 = ng.PHs.size();
    E_Int s(nb_phs0);

    childpos.resize(nb_phs_ref + 1);//one-pass-the-end size
    for (E_Int i = 0; i < nb_phs_ref; ++i)
    {
      childpos[i] = s;  // numero du polygone ou commence les enfants
      s += pregnant[i];
    }

    childpos[nb_phs_ref] = s;
    std::vector<E_Int> unodes;
    E_Int nbFineCells = 0;
    for (size_t j = 0; j < PH_to_ref.size(); ++j)
    {
      E_Int* first_pg = ng.PHs.get_facets_ptr(PH_to_ref[j]);
      E_Int  nb_pgs = ng.PHs.stride(PH_to_ref[j]);
      unodes.clear();
      K_MESH::Polyhedron<UNKNOWN>::unique_nodes(ng.PGs, first_pg, nb_pgs, unodes);

      nbFineCells += unodes.size();
    }

    // l'ordre est donne par le retour de unique_nodes pour les cellules a diviser 
    //
    std::vector<E_Int> arity, pregnant_c(nbFineCells);
    E_Int ind = 0;

    for (size_t j = 0; j < PH_to_ref.size(); ++j)
    {
      E_Int* first_pg = ng.PHs.get_facets_ptr(PH_to_ref[j]);
      E_Int  nb_pgs = ng.PHs.stride(PH_to_ref[j]);

      unodes.clear();
      arity.clear();
      K_MESH::Polyhedron<UNKNOWN>::unique_nodes_arity(ng.PGs, first_pg, nb_pgs, unodes, arity);

      for (size_t k = 0; k < unodes.size(); k++)
      {
        // number of faces of the future fine PHc based on the k-th node arity of the j-th coarse-cell to adapt
        pregnant_c[ind] = 2 * arity[k];
        ++ind;
      }
    }

    ng.PHs.expand_variable_stride(nbFineCells, &pregnant_c[0]); // expand PGs for internals QUADS

    intpos.resize(nb_phs_ref + 1);  //one-pass-the-end size
    s = nb_pgs0;
    for (E_Int i = 0; i < nb_phs_ref; ++i)
    {
      intpos[i] = s;
      E_Int PHi = PH_to_ref[i];
      s += subdiv_pol<K_MESH::Polyhedron<0>, ISO_HEX>::nbi(ng.PGs, ng.PHs.get_facets_ptr(PHi), ng.PHs.stride(PHi));
    }
    intpos[nb_phs_ref] = s;

    // allocate space in coordinates for centroids
    crd.resize(3, crd.cols() + nb_phs_ref);
  }

  // Impl for HEXA/DIR
  template <>
  template <typename pg_arr_t, typename ph_arr_t>
  void refiner<K_MESH::Hexahedron, DIR>::reserve_mem_PHs
  (K_FLD::FloatArray& crd, ngon_type& ng, const Vector_t<E_Int> & PH_to_ref, const Vector_t<NUGA::eDIR> & PH_directive, tree<pg_arr_t> & PGtree, tree<ph_arr_t> & PHtree, K_FLD::IntArray& F2E,
    std::vector<E_Int>& intpos, std::vector<E_Int>& childpos)
  {
    //
    E_Int nb_phs_ref = PH_to_ref.size();

    E_Int nbi_tot = subdiv_pol<K_MESH::Hexahedron, DIR>::nbi_sum(ng, PH_to_ref, PH_directive);

    // reserve space for internal faces in the tree
    E_Int current_sz = PGtree.size();
    assert(current_sz != -1);

    // just sync the tree attributes (internals have no parents so do not a call to PGtree.resize)
    PGtree.resize_hierarchy(current_sz + nbi_tot); //fixme : in hmesh.init ?

    E_Int nb_pgs0 = ng.PGs.size();

    // expand PGs for internals QUADS
    ng.PGs.expand_n_fixed_stride(nbi_tot, 4);

    //
    F2E.resize(2, nb_pgs0 + nbi_tot, IDX_NONE);

    // reserve space for children in the tree
    std::vector<E_Int> pregnant;
    subdiv_pol<K_MESH::Hexahedron, DIR>::nbc_list(ng, PH_to_ref, PH_directive, pregnant);

    PHtree.resize(PH_to_ref, pregnant);
    E_Int nb_phs0{ ng.PHs.size() };
    E_Int nb_tot_child = 0;
    for (E_Int i = 0; i < nb_phs_ref; ++i)
      nb_tot_child += pregnant[i];
    ng.PHs.expand_n_fixed_stride(nb_tot_child, K_MESH::Hexahedron::NB_BOUNDS);

    intpos.resize(nb_phs_ref + 1);//one-pass-the-end size
    E_Int s(nb_pgs0), nb_xyz(0);
    for (E_Int i = 0; i < nb_phs_ref; ++i)
    {
      intpos[i] = s;
      if (PH_directive[i] == XYZ) {
        s += 12;
        ++nb_xyz;
      }
      else if (PH_directive[i] == XY || PH_directive[i] == XZ || PH_directive[i] == YZ)
        s += 4;
      else if (PH_directive[i] == Xd || PH_directive[i] == Y || PH_directive[i] == Z)
        s += 1;
    }
    intpos[nb_phs_ref] = s;

    childpos.resize(nb_phs_ref + 1);//one-pass-the-end size
    s = nb_phs0;
    for (E_Int i = 0; i < nb_phs_ref; ++i)
    {
      childpos[i] = s;
      s += pregnant[i];
    }
    childpos[nb_phs_ref] = s;

    // allocate space in coordinates for XYZ centroids
    crd.resize(3, crd.cols() + nb_xyz);
  }

  ///
  template <>
  template <typename pg_arr_t, typename ph_arr_t>
  void refiner<K_MESH::Hexahedron, DIR_PROTO>::reserve_mem_PHs
  (K_FLD::FloatArray& crd, ngon_type& ng, const Vector_t<E_Int> & PH_to_ref, const Vector_t<NUGA::eDIR> & PH_directive, tree<pg_arr_t> & PGtree, tree<ph_arr_t> & PHtree, K_FLD::IntArray& F2E,
    std::vector<E_Int>& intpos, std::vector<E_Int>& childpos)
  {
    refiner<K_MESH::Hexahedron, DIR>::reserve_mem_PHs (crd, ng, PH_to_ref, PH_directive, PGtree, PHtree, F2E, intpos, childpos);
  }

  ///
  template<> inline
    void refiner<K_MESH::Quadrangle, ISO>::refine_PG
    (K_FLD::FloatArray& crd, ngon_unit& PGs, E_Int PGi, E_Int posC, E_Int posChild, eDIR d, const std::map<K_MESH::NO_Edge, E_Int>& ecenter)
  {
    E_Int* nodes = PGs.get_facets_ptr(PGi);
    assert(PGs.stride(PGi) == 4);
    //const E_Int &nbc = subdiv_pol<K_MESH::Quadrangle, ISO>::NBC;

    // Centroid calculation
    K_MESH::Polygon::iso_barycenter(crd, nodes, 4, 1, crd.col(posC));

    E_Int q9[9];
    K_MESH::NO_Edge noE;
    for (E_Int n = 0; n < K_MESH::Quadrangle::NB_NODES; ++n)
    {
      q9[n] = *(nodes + n);
      noE.setNodes(*(nodes + n), *(nodes + (n + 1) % 4));
      auto it = ecenter.find(noE);
      assert(it != ecenter.end());
      q9[n + K_MESH::Quadrangle::NB_NODES] = it->second;
    }
    q9[8] = posC + 1;

    NUGA::Q9::split<eSUBDIV_TYPE::ISO>(PGs, q9, posChild);
  }

  ///
  template<> inline
    void refiner<K_MESH::Quadrangle, DIR>::refine_PG
    (K_FLD::FloatArray& crd, ngon_unit& PGs, E_Int PGi, E_Int posC, E_Int posChild, eDIR d, const std::map<K_MESH::NO_Edge, E_Int>& ecenter)
  {
    E_Int* nodes = PGs.get_facets_ptr(PGi);
    //const E_Int &nbc = subdiv_pol<K_MESH::Quadrangle, ISO>::NBC;

    E_Int q6[6];
    K_MESH::NO_Edge noE;

    E_Int m1, m2;

    if (d == Xd)
    {
      noE.setNodes(*(nodes), *(nodes + 1));
      auto it = ecenter.find(noE);
      assert(it != ecenter.end());
      m1 = it->second;
      noE.setNodes(*(nodes + 2), *(nodes + 3));
      it = ecenter.find(noE);
      assert(it != ecenter.end());
      m2 = it->second;
    }
    else //d == Y
    {
      noE.setNodes(*(nodes + 1), *(nodes + 2));
      auto it = ecenter.find(noE);
      assert(it != ecenter.end());
      m1 = it->second;
      noE.setNodes(*(nodes), *(nodes + 3));
      it = ecenter.find(noE);
      assert(it != ecenter.end());
      m2 = it->second;
    }

    for (E_Int n = 0; n < K_MESH::Quadrangle::NB_NODES; ++n)
      q6[n] = *(nodes + n);
    q6[4] = m1;
    q6[5] = m2;

    NUGA::Q6::split(PGs, q6, d, posChild);
  }

  template<> inline
    void refiner<K_MESH::Quadrangle, DIR_PROTO>::refine_PG
    (K_FLD::FloatArray& crd, ngon_unit& PGs, E_Int PGi, E_Int posC, E_Int posChild, eDIR d, const std::map<K_MESH::NO_Edge, E_Int>& ecenter)
  {
    refiner<K_MESH::Quadrangle, DIR>::refine_PG(crd, PGs, PGi, posC, posChild, d, ecenter);
  }

  ///
  template<> inline
    void refiner<K_MESH::Triangle, ISO>::refine_PG
    (K_FLD::FloatArray& crd, ngon_unit& PGs, E_Int PGi, E_Int posC, E_Int posChild, eDIR d, const std::map<K_MESH::NO_Edge, E_Int>& ecenter)
  {
    E_Int* nodes = PGs.get_facets_ptr(PGi);
    //const E_Int &nbc = subdiv_pol<K_MESH::Triangle, ISO>::NBC;

    E_Int T6[6];
    K_MESH::NO_Edge noE;
    for (E_Int n = 0; n < K_MESH::Triangle::NB_NODES; ++n)
    {
      T6[n] = *(nodes + n);
      noE.setNodes(*(nodes + n), *(nodes + (n + 1) % K_MESH::Triangle::NB_NODES));
      auto it = ecenter.find(noE);
      assert(it != ecenter.end());
      T6[n + K_MESH::Triangle::NB_NODES] = it->second;
    }

    NUGA::T6::split<eSUBDIV_TYPE::ISO>(PGs, T6, posChild);
  }

  ///
  template<> inline
    void refiner<K_MESH::Polygon, ISO_HEX>::refine_PG
    (K_FLD::FloatArray& crd, ngon_unit& PGs, E_Int PGi, E_Int posC, E_Int posChild, eDIR d, const std::map<K_MESH::NO_Edge, E_Int>& ecenter)
  {
    //   remplir refE. A toi de choisir une convention d'ordre de stockage coherente avec l'appel qui suit
    //   (K_MESH::Polygon::split): si le dernier noeud est le centroid, faire : refE[nrefnodes-1]= posC + 1
    E_Int* nodes = PGs.get_facets_ptr(PGi);
    E_Int nnodes = PGs.stride(PGi);

    // Centroid calculation
    K_MESH::Polygon::iso_barycenter(crd, nodes, nnodes, 1, crd.col(posC));

    //refining face structure : nodes + mid points + centroid
    E_Int nrefnodes = nnodes * 2 + 1;
    STACK_ARRAY(E_Int, nrefnodes, refE); // emulate E_Int refE[refnodes]

    //todo JP : remplir refE. A toi de choisir une convention d'ordre de stockage coherente avec l'appel qui suit à K_MESH::Polygon::split,
    // fonction que tu dois definir egalement.
    // si le dernier noeud est le centroid, faire : refE[nrefnodes-1]= posC + 1 

    K_MESH::NO_Edge noE;
    for (E_Int n = 0; n < nnodes; ++n)
    {
      refE[2 * n] = *(nodes + n);
      noE.setNodes(*(nodes + n), *(nodes + (n + 1) % nnodes));
      auto it = ecenter.find(noE);
      assert(it != ecenter.end());
      refE[2 * n + 1] = it->second;
    }
    refE[2 * nnodes] = posC + 1;

    K_MESH::Polygon::split<ISO_HEX>(refE.get(), nrefnodes, PGs, posChild);
  }

  ///
  template<> inline
    void refiner<K_MESH::Polygon, ISO>::refine_PG
    (K_FLD::FloatArray& crd, ngon_unit& PGs, E_Int PGi, E_Int posC, E_Int posChild, eDIR d, const std::map<K_MESH::NO_Edge, E_Int>& ecenter)
  {
    assert(false);
  }

  ///
  template<> inline
    void refiner<K_MESH::Polygon, DIR_PROTO>::refine_PG
    (K_FLD::FloatArray& crd, ngon_unit& PGs, E_Int PGi, E_Int posC, E_Int posChild, eDIR d, const std::map<K_MESH::NO_Edge, E_Int>& ecenter)
  {
    assert(false);
  }

  ///
  template<> inline
    void refiner<K_MESH::Triangle, DIR_PROTO>::refine_PG
    (K_FLD::FloatArray& crd, ngon_unit& PGs, E_Int PGi, E_Int posC, E_Int posChild, eDIR d, const std::map<K_MESH::NO_Edge, E_Int>& ecenter)
  {
    //todo
  }

  template<> inline
    void refiner<K_MESH::Triangle, ISO_HEX>::refine_PG
    (K_FLD::FloatArray& crd, ngon_unit& PGs, E_Int PGi, E_Int posC, E_Int posChild, eDIR d, const std::map<K_MESH::NO_Edge, E_Int>& ecenter)
  {
    refiner<K_MESH::Polygon, ISO_HEX>::refine_PG(crd, PGs, PGi, posC, posChild, d, ecenter);
  }

  template<> inline
    void refiner<K_MESH::Quadrangle, ISO_HEX>::refine_PG
    (K_FLD::FloatArray& crd, ngon_unit& PGs, E_Int PGi, E_Int posC, E_Int posChild, eDIR d, const std::map<K_MESH::NO_Edge, E_Int>& ecenter)
  {
    refiner<K_MESH::Polygon, ISO_HEX>::refine_PG(crd, PGs, PGi, posC, posChild, d, ecenter);
  }

  /// default implementation : T3/Q4-ISO, PG-ISO_HEX cases
  template <typename ELT_t, eSUBDIV_TYPE STYPE>
  template <typename arr_t>
  void refiner<ELT_t, STYPE>::refine_PGs
  (Vector_t<E_Int>& PG_to_ref, Vector_t<NUGA::eDIR> & PG_directive,
    ngon_type& ng, tree<arr_t>& PGtree, K_FLD::FloatArray& crd, K_FLD::IntArray& F2E,
    std::map<K_MESH::NO_Edge, E_Int>& ecenter)
  {
    // Compute Edges refine points
    refine_point_computer::compute_centers(PG_to_ref, PG_directive, ng, crd, ecenter);

    //E_Int pos = crd.cols(); // first face center in crd (if relevant)

    // Reserve space for children in the tree
    std::vector<E_Int> childpos, crdpos;
    reserve_mem_PGs(crd, ng.PGs, PG_to_ref, PG_directive, PGtree, F2E, childpos, crdpos);

    // Refine them
#ifndef DEBUG_HIERARCHICAL_MESH
//#pragma omp parallel for
#endif
    for (size_t i = 0; i < PG_to_ref.size(); ++i)
    {
      E_Int PGi = PG_to_ref[i];
      E_Int firstChild = childpos[i];

      E_Int nbc = childpos[i + 1] - firstChild;

      E_Int C = !crdpos.empty() ? crdpos[i] : IDX_NONE;

      if (nbc == 2 && STYPE == DIR)
        refiner<ELT_t, DIR_PROTO>::refine_PG(crd, ng.PGs, PGi, C, firstChild, PG_directive[i], ecenter);
      else if (nbc == 2 && STYPE == DIR_PROTO)
        refiner<ELT_t, DIR_PROTO>::refine_PG(crd, ng.PGs, PGi, C, firstChild, PG_directive[i], ecenter);
      else if (STYPE == ISO || STYPE == DIR_PROTO || STYPE ==DIR)
      {
        assert(nbc == 4);
        refiner<ELT_t, ISO>::refine_PG(crd, ng.PGs, PGi, C, firstChild, XY, ecenter);
      }
      else if (STYPE == ISO_HEX)
        refiner<ELT_t, ISO_HEX>::refine_PG(crd, ng.PGs, PGi, C, firstChild, XY, ecenter);
      else
        assert(false);

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

  /// default implementation : ISO case , for all basic element types (Basic, Tetra, Pyra, Penta, Hexa)
  template <typename ELT_t, eSUBDIV_TYPE STYPE>
  template <typename arr_t>
  void refiner<ELT_t, STYPE>::refine_Faces
  (const output_t &adap_incr,
    ngon_type& ng, tree<arr_t> & PGtree, K_FLD::FloatArray& crd, K_FLD::IntArray & F2E,
    std::map<K_MESH::NO_Edge, E_Int>& ecenter)
  {
    // deal with quads (ISO/DIR)
    Vector_t<E_Int> PG_to_ref;
    Vector_t<NUGA::eDIR> PG_directive;
    E_Int nb_pgs_ref = refiner<K_MESH::Quadrangle, STYPE>::extract_PGs_to_refine(crd, ng, PGtree, adap_incr, PG_to_ref, PG_directive);
    if (nb_pgs_ref != 0)
      refiner<K_MESH::Quadrangle, STYPE>::refine_PGs(PG_to_ref, PG_directive, ng, PGtree, crd, F2E, ecenter);

    // deal with tris (ISO/DIR)
    PG_to_ref.clear();
    PG_directive.clear();
    nb_pgs_ref = refiner<K_MESH::Triangle, STYPE>::extract_PGs_to_refine(crd, ng, PGtree, adap_incr, PG_to_ref, PG_directive);
    if (nb_pgs_ref != 0)
      refiner<K_MESH::Triangle, STYPE>::refine_PGs(PG_to_ref, PG_directive, ng, PGtree, crd, F2E, ecenter);
  }

  // Polyhedron/ISO_HEX impl
  template <>
  template <typename arr_t>
  void refiner<K_MESH::Polyhedron<0>, ISO_HEX>::refine_Faces
  (const output_t &adap_incr,
    ngon_type& ng, tree<arr_t> & PGtree, K_FLD::FloatArray& crd, K_FLD::IntArray & F2E,
    std::map<K_MESH::NO_Edge, E_Int>& ecenter)
  {
    Vector_t<E_Int> PG_to_ref;
    Vector_t<NUGA::eDIR> dummy;
    E_Int nb_pgs_ref = refiner<K_MESH::Polygon, ISO_HEX>::extract_PGs_to_refine(crd/*hack for CLEF*/, ng, PGtree, adap_incr, PG_to_ref, dummy);
    if (nb_pgs_ref != 0)
      refiner<K_MESH::Polygon, ISO_HEX>::refine_PGs(PG_to_ref, dummy, ng, PGtree, crd, F2E, ecenter);
  }

  /// HEXA/TETRA/PYRA/PRISM ISO, PH-ISO_HEX impl
  template <typename ELT_t, eSUBDIV_TYPE STYPE>
  template <typename pg_arr_t, typename ph_arr_t>
  void refiner<ELT_t, STYPE>::refine_PHs
  (output_t &adap_incr, ngon_type& ng, tree<pg_arr_t> & PGtree, tree<ph_arr_t> & PHtree,
    K_FLD::FloatArray& crd, K_FLD::IntArray & F2E)
  {

    assert((E_Int)adap_incr.cell_adap_incr.size() <= ng.PHs.size());
    //
    Vector_t<E_Int> PH_to_ref;
    Vector_t<NUGA::eDIR>  PH_directive;
    extract_PHs_to_refine(ng, PHtree, adap_incr, PH_to_ref, PH_directive);
    E_Int nb_phs_ref = PH_to_ref.size();
    if (nb_phs_ref == 0) return;

    E_Int pos = crd.cols(); // first cell center in crd (if relevant)
    //E_Int nb_pgs0 = ng.PGs.size();
    //E_Int nb_phs0 = ng.PHs.size();
    // Reserve space for children in the tree
    std::vector<E_Int> intpos, childpos;
    reserve_mem_PHs(crd, ng, PH_to_ref, PH_directive, PGtree, PHtree, F2E, intpos, childpos);

    using spliting_t = splitting_t<ELT_t, XYZ, 1>; // HX27, TH10, PY13 or PR18, HX18

#ifndef DEBUG_HIERARCHICAL_MESH
//#pragma omp parallel for
#endif
    for (E_Int i = 0; i < nb_phs_ref; ++i)
    {
      E_Int PHi = PH_to_ref[i];

      spliting_t elt(crd, ng, PHi, pos + i, F2E, PGtree);

      elt.split(ng, PHi, PHtree, PGtree, F2E, intpos[i], childpos[i]);
    }
  }

  ///
  template <>
  template <typename pg_arr_t, typename ph_arr_t>
  void refiner<K_MESH::Basic, eSUBDIV_TYPE::ISO>::refine_PHs
  (output_t &adap_incr,
    ngon_type& ng, tree<pg_arr_t> & PGtree, tree<ph_arr_t> & PHtree, K_FLD::FloatArray& crd, K_FLD::IntArray & F2E)
  {
    assert((E_Int)adap_incr.cell_adap_incr.size() <= ng.PHs.size());

    refiner<K_MESH::Hexahedron, eSUBDIV_TYPE::ISO>::refine_PHs(adap_incr, ng, PGtree, PHtree, crd, F2E);

    refiner<K_MESH::Tetrahedron, eSUBDIV_TYPE::ISO>::refine_PHs(adap_incr, ng, PGtree, PHtree, crd, F2E);

    refiner<K_MESH::Pyramid, eSUBDIV_TYPE::ISO>::refine_PHs(adap_incr, ng, PGtree, PHtree, crd, F2E);

    refiner<K_MESH::Prism, eSUBDIV_TYPE::ISO>::refine_PHs(adap_incr, ng, PGtree, PHtree, crd, F2E);

  }

  // DIRECTIONNEL
  
  template <typename ELT_t, eSUBDIV_TYPE STYPE>
  template <typename pg_arr_t, typename ph_arr_t>
  void refiner<ELT_t, STYPE>::__refine_PHS_hexa_dir
  (output_t &adap_incr,
    ngon_type& ng, tree<pg_arr_t> & PGtree, tree<ph_arr_t> & PHtree, K_FLD::FloatArray& crd, K_FLD::IntArray & F2E)
  {
    assert((E_Int)adap_incr.cell_adap_incr.size() <= ng.PHs.size());
    //
    Vector_t<E_Int> PH_to_ref;
    Vector_t<NUGA::eDIR>  PH_directive;
    extract_PHs_to_refine(ng, PHtree, adap_incr, PH_to_ref, PH_directive);
    E_Int nb_phs_ref = PH_to_ref.size();
    if (nb_phs_ref == 0) return;

    E_Int pos = crd.cols(); // first cell center in crd (if relevant)
                            //E_Int nb_pgs0 = ng.PGs.size();
                            //E_Int nb_phs0 = ng.PHs.size();
                            // Reserve space for children in the tree
    std::vector<E_Int> intpos, childpos;
    reserve_mem_PHs(crd, ng, PH_to_ref, PH_directive, PGtree, PHtree, F2E, intpos, childpos);

#ifndef DEBUG_HIERARCHICAL_MESH
//#pragma omp parallel for
#endif
    int xyzG = 0;
    for (E_Int i = 0; i < nb_phs_ref; ++i)
    {
      E_Int PHi = PH_to_ref[i];

      if (PH_directive[i] == XYZ)
      {
        HX27 elt(crd, ng, PHi, pos + (xyzG++), F2E, PGtree);

        elt.split(ng, PHi, PHtree, PGtree, F2E, intpos[i], childpos[i]);
      }
      else if (PH_directive[i] == XY)
      {
        bool do_reorder = (STYPE == DIR_PROTO); // do reoder (reorder_as_XY call) only if DIR_PROTO
        HX18 elt(crd, ng, PHi, pos + i, F2E, PGtree, do_reorder);

        if (elt.ok_for_split)
          elt.split(ng, PHi, PHtree, PGtree, F2E, intpos[i], childpos[i]);
        else
          adap_incr.cell_adap_incr[PHi].set(0);
      }
      else // X
      {
        HX12 elt(crd, ng, PHi, pos + i, F2E, PGtree);

        if (elt.ok_for_split)
          elt.split(ng, PHi, PHtree, PGtree, F2E, intpos[i], childpos[i]);
        else
          adap_incr.cell_adap_incr[PHi].set(0);
      }
    }
 }

  template <>
  template <typename pg_arr_t, typename ph_arr_t>
  void refiner<K_MESH::Hexahedron, DIR>::refine_PHs
  (output_t &adap_incr,
    ngon_type& ng, tree<pg_arr_t> & PGtree, tree<ph_arr_t> & PHtree, K_FLD::FloatArray& crd, K_FLD::IntArray & F2E)
  {
    refiner<K_MESH::Hexahedron, DIR>::__refine_PHS_hexa_dir <pg_arr_t, ph_arr_t> (adap_incr, ng, PGtree, PHtree, crd, F2E);
  }

  template <>
  template <typename pg_arr_t, typename ph_arr_t>
  void refiner<K_MESH::Hexahedron, DIR_PROTO>::refine_PHs
  (output_t &adap_incr,
    ngon_type& ng, tree<pg_arr_t> & PGtree, tree<ph_arr_t> & PHtree, K_FLD::FloatArray& crd, K_FLD::IntArray & F2E)
  {
    refiner<K_MESH::Hexahedron, DIR_PROTO>::__refine_PHS_hexa_dir <pg_arr_t, ph_arr_t> (adap_incr, ng, PGtree, PHtree, crd, F2E);
  }


  //default impl : DIR IS ISO for all elements but HEXA "Layer"
  /*template <>
  template <typename pg_arr_t, typename ph_arr_t>
  void refiner<K_MESH::Tetrahedron, DIR_PROTO>::refine_PHs
  (const output_t &adap_incr,
    ngon_type& ng, tree<pg_arr_t> & PGtree, tree<ph_arr_t> & PHtree, K_FLD::FloatArray& crd, K_FLD::IntArray & F2E)
  {
    //refiner<K_MESH::Tetrahedron, ISO>::refine_PHs(adap_incr.cell_adap_incr, ng, PGtree, PHtree, crd, F2E);
  }*/

  ///
  template <>
  template <typename pg_arr_t, typename ph_arr_t>
  void refiner<K_MESH::Pyramid, DIR_PROTO>::refine_PHs
  (output_t &adap_incr,
    ngon_type& ng, tree<pg_arr_t> & PGtree, tree<ph_arr_t> & PHtree, K_FLD::FloatArray& crd, K_FLD::IntArray & F2E)
  {
    //refiner<K_MESH::Pyramid, ISO>::refine_PHs(adap_incr.cell_adap_incr, ng, PGtree, PHtree, crd, F2E);
  }

  ///
  template <>
  template <typename pg_arr_t, typename ph_arr_t>
  void refiner<K_MESH::Prism, DIR_PROTO>::refine_PHs
  (output_t &adap_incr,
    ngon_type& ng, tree<pg_arr_t> & PGtree, tree<ph_arr_t> & PHtree, K_FLD::FloatArray& crd, K_FLD::IntArray & F2E)
  {
    //refiner<K_MESH::Prism, ISO>::refine_PHs(adap_incr.cell_adap_incr, ng, PGtree, PHtree, crd, F2E);
  }

  ///
  /*template <>
  template <typename pg_arr_t, typename ph_arr_t>
  void refiner<K_MESH::Hexahedron, DIR_PROTO>::refine_PHs
  (const output_t &adap_incr,
   ngon_type& ng, tree<pg_arr_t> & PGtree, tree<ph_arr_t> & PHtree, K_FLD::FloatArray& crd, K_FLD::IntArray & F2E)
  {
    //alexis : todo
  }

  ///
  template <>
  template <typename pg_arr_t, typename ph_arr_t>
  void refiner<K_MESH::Basic, DIR_PROTO>::refine_PHs
  (const output_t &adap_incr,
   ngon_type& ng, tree<pg_arr_t> & PGtree, tree<ph_arr_t> & PHtree, K_FLD::FloatArray& crd, K_FLD::IntArray & F2E)
  {
    //alexis : todo
  }*/

} // NUGA


#endif /* REFINER_HXX */

/*    
    Copyright 2013-2026 ONERA.

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

#ifndef NUGA_CLASSIFYER_HXX
#define NUGA_CLASSIFYER_HXX

//#define CLASSIFYER_DBG

#include "Nuga/include/macros.h"
#include "Nuga/include/selector.hxx"
#include "Nuga/include/collider.hxx"
#include "Nuga/include/displacement.hxx"
#include "Nuga/include/random.h"

#ifdef CLASSIFYER_DBG
#include "Nuga/include/medit.hxx"
#endif

#define TEMPLATE_TYPES template<eXPolicy POLICY, typename zmesh_t, typename bound_mesh_t>
#define TEMPLATE_CLASS classifyer<POLICY, zmesh_t, bound_mesh_t>

namespace NUGA
{
  enum eXPolicy { COLLISION, XCELLN_VAL, XCELLN_OUT, MESH };
  enum eClassify { AMBIGUOUS=-1, IN/*IN_2*/=0, IN_1=1, X=2, OUT=3, UPPER_COL=4 };

  template <eXPolicy POLICY, typename zmesh_t> // implemented for COLLISION policy. Second template arg is only there for XCELLN_OUT/MESH modes
  struct data_trait
  {
    using wdata_t = std::vector<E_Float>;
    using outdata_t = std::vector<E_Float>;
  
    static void mark_cell_w_mask(wdata_t& data, E_Int i, E_Int im, E_Int val) { assert ((size_t)i < data.size()); data[i] = val; } // minus to mark as new X for __flag_hidden_subzones
  };

  struct color_t
  {
    color_t(E_Int col) :val(col) {}
    color_t()
    {
      //todo SL
      val = OUT;  //where this is needed?
    }

    color_t& operator=(color_t& col) { val = col.val; masks = col.masks; return *this;}

    color_t& operator=(E_Int col) { val = col; return *this; } // do not touch to masks
    bool operator!=(E_Int col) const = delete; // { return (col != val); }
    bool operator==(E_Int col) const { return (col == val); }

    operator double() const { return val; }

    double val;
    std::vector<E_Int> masks;
  };

  template <typename zmesh_t>
  struct data_trait<XCELLN_VAL, zmesh_t>
  {
    using wdata_t = std::vector<color_t>;
    using outdata_t = std::vector<double>;

    static void mark_cell_w_mask(wdata_t & data, E_Int i, E_Int im, E_Int val)
    {
      // minus to mark as new X for __flag_hidden_subzones
      assert ((size_t)i < data.size());
      data[i].val = val;
      data[i].masks.push_back(im);
    }
  };

  template <typename zmesh_t>
  struct data_trait<XCELLN_OUT, zmesh_t> : public data_trait<XCELLN_VAL, zmesh_t>
  {
    struct outdata_t 
    {
      zmesh_t mesh;
      bool full_out;
      outdata_t() :mesh(), full_out(false) {};
      outdata_t(outdata_t&& d):mesh(std::move(d.mesh)), full_out(d.full_out){}
      outdata_t& operator=(outdata_t&& d)
      {
        mesh = std::move(d.mesh);
        full_out = d.full_out;
        return *this;
      }
    };
  };

  ///
  template<eXPolicy POLICY, typename zmesh_type, typename bound_mesh_type = typename NUGA::boundary_t<zmesh_type>>
  class classifyer
  {
  public:
    using zmesh_t = zmesh_type;
    using bmesh_t = bound_mesh_type;
    using wdata_t   = typename data_trait<POLICY, zmesh_t>::wdata_t;
    using outdata_t = typename data_trait<POLICY, zmesh_t>::outdata_t;

  public:
    classifyer(double RTOL):_RTOL(RTOL) {}

    void prepare(zmesh_t & z_mesh,
                 const std::vector<K_FLD::FloatArray> &mask_crds,
                 const std::vector<K_FLD::IntArray>& mask_cnts,
                 std::vector< std::vector<E_Int> > &mask_wall_ids,
                 const std::vector<K_SEARCH::BBox3D>& comp_boxes,
                 const std::vector<E_Int>& z_priorities, E_Int rank_wnp,
                 std::vector< bmesh_t*> & mask_bits, bmesh_t& WP, bmesh_t& WNP);

    E_Int compute(zmesh_t const & z_mesh, std::vector< bmesh_t*> const & mask_bits, bmesh_t& WP, bmesh_t& WNP, outdata_t& outdata);

    void finalize(zmesh_t const & z_mesh, outdata_t& outdata);

  protected:

    virtual outdata_t __process_X_cells(zmesh_t const & z_mesh, std::vector< bmesh_t*> const & mask_bits, wdata_t & wdata) = 0;

    void __compact_to_box(zmesh_t const & z_mesh, std::vector< bmesh_t*> & mask_bits,
                          bmesh_t& WP, bmesh_t& WNP);

    void __build_mask_bits(const std::vector<K_FLD::FloatArray> &mask_crds,
                           const std::vector<K_FLD::IntArray>& mask_cnts,
                           std::vector< std::vector<E_Int> > &mask_wall_ids,
                           const std::vector<E_Int>& z_priorities, E_Int rank_wnp,
                           std::vector< bmesh_t*> & mask_bits);

    void __build_mask_bits2(zmesh_t& z_mesh, double ARTOL, eMetricType mtype, double PS_MIN, // these parameters to deal with overlapping
                            const std::vector<K_FLD::FloatArray> &mask_crds,
                            const std::vector<K_FLD::IntArray>& mask_cnts,
                            std::vector< std::vector<E_Int> > &mask_wall_ids,
                            const std::vector<E_Int>& z_priorities, E_Int rank_wnp,
                            std::vector< bmesh_t*> & mask_bits,
                            bmesh_t& WP, bmesh_t& WNP);

    void __process_overlapping_boundaries(zmesh_t & z_mesh, std::vector< bmesh_t*> & mask_bits, E_Int rank_wnp, E_Float RTOL);

    bool __flag_colliding_cells(zmesh_t const & z_mesh, std::vector< bmesh_t*> const & mask_bits, E_Int im, wdata_t& wdata, E_Int XCOL);

    E_Int __flag_hidden_subzones(zmesh_t const & z_mesh, bmesh_t const & mask_bit, wdata_t& wdata, E_Int XCOL, E_Int nsubmin);

  protected:
    E_Float _RTOL;
    std::vector<K_SEARCH::BBox3D> _comp_boxes;
    std::vector<E_Int> _z_priorities;

  };

  namespace CLASSIFY // some "global" functions : they do not need all the template params of the classifyer class
  {
    template <typename T1, typename T2>
    static inline eClassify classify(T1 const& t1, T2 const& t2, bool deep);

    template <>
    inline eClassify classify(NUGA::aPolygon const& ae1, edge_mesh_t const& front, bool deep)
    {
      const E_Float* norm = ae1.get_normal();
      const E_Float* pt = ae1.get_centroid();

      E_Int sign(0);
      NUGA::random rando;

      for (E_Int i = 0; i < front.ncells(); ++i)
      {
        E_Float v1[3], v2[3], v[3];
        NUGA::diff<3>(front.crd.col(front.cnt(0, i)), pt, v1);
        NUGA::diff<3>(front.crd.col(front.cnt(1, i)), pt, v2);

        NUGA::crossProduct<3>(v1, v2, v);
        NUGA::normalize<3>(v);

        E_Float psi = NUGA::dot<3>(norm, v);
        E_Int sigpsi = zSIGN(psi, EPSILON);

        if (sigpsi == 0) //AMBIGUOUS
        {
          sign = 0;
          break;
        }

        if (sign == 0) sign = sigpsi;
        else if (sign*sigpsi < 0) // AMBIGUOUS
        {
          sign = 0;
          break;
        }
      }

      if (sign == 0) // AMBIGUOUS => deeper test to valuate sign based on visibility
      {
        if (!deep) return AMBIGUOUS;

        // pick randomly an edge, for the ray(GC) from its mid point to the centroid of ae1
        unsigned int k = rando.rand() % front.ncells();
        E_Float C[3], G[3], nrm[3];
        NUGA::sum<3>(0.5, front.crd.col(front.cnt(0, k)), 0.5, front.crd.col(front.cnt(1, k)), C);
        ae1.centroid<3>(G);
        ae1.normal<3>(nrm);

        // get the visible edge : since we are 3D, we use plane containing the edge and norm
        E_Float lambda_min(NUGA::FLOAT_MAX), R[3];
        for (E_Int j = 0; j < front.ncells(); ++j)
        {
          const E_Float* P = front.crd.col(front.cnt(0, j));
          const E_Float* Q = front.crd.col(front.cnt(1, j));
          NUGA::sum<3>(Q, nrm, R);

          E_Float lambda, UV[2], min_d;
          E_Bool parallel, coincident;
          K_MESH::Triangle::planeLineMinDistance<3>(P, Q, R, G, C, EPSILON, true, lambda, UV, parallel, coincident, min_d, true/*strict*/);

          if (lambda < lambda_min)
          {
            // is G above or under the plane ?
            E_Float PQ[3], pnorm[3], ray[3];
            NUGA::diff<3>(Q, P, PQ);
            NUGA::crossProduct<3>(PQ, nrm, pnorm);
            NUGA::diff<3>(C, G, ray);
            E_Float ps = NUGA::dot<3>(ray, pnorm);
            sign = zSIGN(ps, EPSILON); // >0 means under
            lambda_min = lambda;
          }
        }
      }

      assert(sign);
      return (sign < 0) ? OUT : IN;
    }

    template <>
    inline eClassify classify(NUGA::aPolyhedron<0> const& ae1, pg_smesh_t const& front, bool deep)
    {
      const E_Float* ae1G = ae1.get_centroid();
      
      E_Int sign(0);
      NUGA::random rando;

      for (E_Int i = 0; i < front.ncells(); ++i)
      {
        E_Float fni[3], ci[3];
        //K_MESH::Polygon pgf(front, i);
        K_MESH::Polygon::normal<K_FLD::FloatArray, 3>(front.crd, front.cnt.get_facets_ptr(i), front.cnt.stride(i), front.index_start, fni);
        K_MESH::Polygon::centroid<3>(front.crd, front.cnt.get_facets_ptr(i), front.cnt.stride(i), front.index_start, ci);

#ifdef CLASSIFYER_DBG
        E_Float l2 = sqrt(fni[0] * fni[0] + fni[1] * fni[1] + fni[2] * fni[2]);
        assert(::fabs(l2 - 1.) < EPSILON); // NOT DEGEN
#endif

        E_Float ray[3];
        NUGA::diff<3>(ci, ae1G, ray);

        E_Float psi = NUGA::dot<3>(fni, ray);
        E_Int sigpsi = zSIGN(psi, EPSILON);

        if (sigpsi == 0) //AMBIGUOUS
        {
          sign = 0;
          break;
        }

        if (sign == 0) sign = sigpsi;
        else if (sign*sigpsi < 0) // AMBIGUOUS
        {
          sign = 0;
          break;
        }
      }

      if (sign == 0) // AMBIGUOUS => deeper test to valuate sign based on visibility
      {
        if (!deep) return AMBIGUOUS;

        // pick a front face, in a REGULAR position for the ray(GC) : from the centroid G of ae1 to face centroid C
        E_Float ray[3], C[3];
        E_Int k{0};
        E_Float n[3];
        for (; (k < front.ncells()) && (sign == 0.); ++k)
        {
          K_MESH::Polygon PGk(front.cnt, k);

          K_MESH::Polygon::centroid<3>(front.crd, PGk.begin(), PGk.nb_nodes(), front.index_start, C);
          K_MESH::Polygon::normal<K_FLD::FloatArray, 3>(front.crd, PGk.begin(), PGk.nb_nodes(), front.index_start, n);
          
          NUGA::diff<3>(C, ae1G, ray);
          E_Float ps = NUGA::dot<3>(ray, n);

          sign = zSIGN(ps, EPSILON); // >0 means under
        }

        assert (k != IDX_NONE); 
        E_Float lambda_min(1.);

        // seek for closest-to-ae1G PG crossing the ray
        for (E_Int j = 0; j < front.ncells(); ++j)
        {
          if (j == k) continue;

          K_MESH::Polygon PGj(front.cnt, j);

          E_Float lambda(NUGA::FLOAT_MAX), u1;
          E_Bool overlap;
          bool isx = PGj.intersect<DELAUNAY::Triangulator>(front.crd, ae1G, C, EPSILON, true, lambda, u1, overlap);

          if (isx && lambda < lambda_min)
          {
            // is G above or under the plane ?
            E_Float nj[3];
            PGj.normal<K_FLD::FloatArray, 3>(front.crd, nj);

#ifdef CLASSIFYER_DBG
            E_Float l2 = sqrt(nj[0] * nj[0] + nj[1] * nj[1] + nj[2] * nj[2]);
            assert(::fabs(l2 - 1.) < EPSILON); // NOT DEGEN
#endif
            NUGA::diff<3>(C, ae1G, ray);
            E_Float ps = NUGA::dot<3>(ray, nj);
            sign = zSIGN(ps, EPSILON); // >0 means under
            lambda_min = lambda;
          }
        }
      }

      assert(sign);
      return (sign < 0) ? OUT : IN; //assume outward orientation
    }
    

    template <>
    inline eClassify classify(NUGA::aPolyhedron<0> const& ae1, NUGA::aPolyhedron<0> const& ae2, bool deep)
    {
      // WARNING : assume ae1 and ae2 are closed surface (polyhedra)
      // WARNING : assume no collision situation => 3 posibilities : 1 in 2, 2 in 1 or separated

      assert(ae1.m_oriented != 0);
      assert(ae2.m_oriented != 0);
      
      // 1 in 2 ?
      const E_Float* G1 = ae1.get_centroid();
      DELAUNAY::Triangulator dt;
      ae2.triangulate(dt, ae2.m_crd);
      E_Float omega = 0.;
      for (E_Int i = 0; i < ae2.nb_tris(); ++i)
      {
        E_Int T[3];
        ae2.triangle(i, T);
        omega += K_MESH::Triangle::oriented_trihedral_angle(G1, ae2.m_crd.col(T[0]), ae2.m_crd.col(T[1]), ae2.m_crd.col(T[2]));
      }

      omega = ::fabs(::fabs(omega) - 4. *NUGA::PI);

      if (omega < 1.e-13) return IN;

      // 2 in 1 ?
      const E_Float* G2 = ae2.get_centroid();
      ae1.triangulate(dt, ae1.m_crd);
      omega = 0.;
      E_Int T[3];
      for (E_Int i = 0; i < ae1.nb_tris(); ++i)
      {
        ae1.triangle(i, T);
        omega += K_MESH::Triangle::oriented_trihedral_angle(G2, ae1.m_crd.col(T[0]), ae1.m_crd.col(T[1]), ae1.m_crd.col(T[2]));
      }

      omega = ::fabs(::fabs(omega) - 4. *NUGA::PI);

      if (omega < 1.e-13) return IN_1;

      return OUT;
    }
    

    static eClassify classify(K_SEARCH::BBox3D const& t1, K_SEARCH::BBox3D const& t2)
    {
      bool is_in = t1.is_included(t2);
      return is_in ? IN : OUT;
    }

    static eClassify classify(K_SEARCH::BBox3D const& t1, std::vector<K_SEARCH::BBox3D> const& ts, std::vector<E_Int> const & ids)
    {
      eClassify ret = OUT;
      for (size_t i = 0; i < ids.size() && (ret == OUT); ++i) 
      {
        K_SEARCH::BBox3D bx = ts[ids[i]];
        bx.enlarge(0.01); //hack for 2D : eg. collar double wall can fall out of fuselage box
        ret = classify(t1, bx);
      }
      //std::cout << "nb prior masks : " << ids.size() << std::endl;
      //std::cout << "box classif : OUT ? " << (ret == OUT) << std::endl;
      return ret;
    }

    
    static inline eClassify classify2D(NUGA::aPolygon & ae1_2D, NUGA::aPolygon& ae2_2D, double ABSTOL)
    {
      DELAUNAY::Triangulator dt;

      bool is_in{ false };
      for (E_Int k = 0; k < ae1_2D.m_crd.cols(); ++k)
      {
        const double* P = ae1_2D.m_crd.col(k);
        // if P is in ae1, ae0 is a piece of ae1
        #ifdef NDEBUG
        ae2_2D.fast_is_in_pred<DELAUNAY::Triangulator, 2>(dt, ae2_2D.m_crd, P, is_in, ABSTOL);
        #else
        E_Int err = ae2_2D.fast_is_in_pred<DELAUNAY::Triangulator, 2>(dt, ae2_2D.m_crd, P, is_in, ABSTOL);
        assert(!err);
        #endif
        
        if (!is_in) break;
      }

      if (is_in) return IN;

      is_in = false;
      for (E_Int k = 0; k < ae2_2D.m_crd.cols(); ++k)
      {
        const E_Float* P = ae2_2D.m_crd.col(k);
        // if P is in ae1, ae0 is a piece of ae1
        ae1_2D.fast_is_in_pred<DELAUNAY::Triangulator, 3>(dt, ae1_2D.m_crd, P, is_in, ABSTOL);
        if (!is_in) break;
      }

      if (is_in) return IN_1;

      return OUT;
    }
    
  }

  ///
  TEMPLATE_TYPES
  void TEMPLATE_CLASS::prepare
  (zmesh_t & z_mesh, const std::vector<K_FLD::FloatArray> &mask_crds,
   const std::vector<K_FLD::IntArray>& mask_cnts,
   std::vector< std::vector<E_Int> > &mask_wall_ids,
   const std::vector<K_SEARCH::BBox3D>& comp_boxes,
   const std::vector<E_Int>& z_priorities, E_Int rank_wnp,
   std::vector< bound_mesh_t*> & mask_bits, bound_mesh_t& WP, bound_mesh_t& WNP)
  {
#ifdef CLASSIFYER_DBG
    static E_Int znb = 0;
    std::cout << "PREP_build_structures_and_reduce_to_zone : __build_mask_bits : " << znb << std::endl;
#endif

    // data to classify no-collision zones
    _z_priorities = z_priorities;
    _z_priorities.resize(rank_wnp); //truncate WNP : not required with bbox logic
    _comp_boxes = comp_boxes;
    
    // build mask data structures (mesh object) : WP are discarded. Putting first decreasing OP, then remaining WNP
    //__build_mask_bits(mask_crds, mask_cnts, mask_wall_ids, z_priorities, rank_wnp, mask_bits);
    __build_mask_bits2(z_mesh, -0.1, NUGA::ISO_MAX, 2.*NUGA::PI*15./360, mask_crds, mask_cnts, mask_wall_ids, z_priorities, rank_wnp, mask_bits, WP, WNP);

#ifdef CLASSIFYER_DBG
    for (size_t m = 0; m < mask_bits.size(); ++m) 
    {
      std::ostringstream o;
      o << "mask_init_z_" << znb << "_m_" << m;
      if (mask_bits[m] != nullptr) medith::write(o.str().c_str(), mask_bits[m]->crd, mask_bits[m]->cnt);
    }
    std::cout << "PREP_build_structures_and_reduce_to_zone : __compact_to_box : " << znb << std::endl;
#endif

    // reduce masks to pieces in zone box : TO REMOVE IF DONE BEFORE IN THE PYTHON
    __compact_to_box(z_mesh, mask_bits, WP, WNP);

#ifdef CLASSIFYER_DBG
    for (size_t m = 0; m < mask_bits.size(); ++m) {
      std::ostringstream o;
      o << "mask_inbox_z_" << znb << "_m_" << m;
      if (mask_bits[m] != nullptr) medith::write(o.str().c_str(), mask_bits[m]->crd, mask_bits[m]->cnt);
    }
    medith::write("WP", WP.crd, WP.cnt);
    medith::write("WNP", WNP.crd, WNP.cnt);
    std::cout << "PREP_build_structures_and_reduce_to_zone : __process_overlapping_boundaries : " << znb << std::endl;
#endif    

    //// 1. detecting overlaps, 
    //// 2. marking on these overlap regions zmesh' border elts as IN, 
    //// 3. discarding overlap boundaries of mask (not required anymore, neither to blank, nor to clip)
    //if (typeid(zmesh_t) == typeid(ph_mesh_t))
    //  __process_overlapping_boundaries(z_mesh, mask_bits, rank_wnp, _RTOL);

#ifdef CLASSIFYER_DBG
    for (size_t m = 0; m < mask_bits.size(); ++m) {
      std::ostringstream o;
      o << "mask_no_ovlp_z_" << znb << "_m_" << m;
      if (mask_bits[m] != nullptr) medith::write(o.str().c_str(), mask_bits[m]->crd, mask_bits[m]->cnt);
    }

    ++znb;
#endif

    // now we have reduced the meshes to the useful part, create & add localizers
    WP.build_localizer();
    WNP.build_localizer();
    for (size_t m = 0; m < mask_bits.size(); ++m)
      if (mask_bits[m] != nullptr) mask_bits[m]->build_localizer();
  }

  ///
  TEMPLATE_TYPES
  E_Int TEMPLATE_CLASS::compute
  (zmesh_t const & z_mesh, std::vector< bound_mesh_t*> const & mask_bits, bound_mesh_t& WP, bound_mesh_t& WNP, outdata_t& outdata)
  {
    E_Int ncells(z_mesh.ncells());

    // initialization of the inner data 
    wdata_t wdata(ncells, OUT);

    size_t nbits = mask_bits.size();

    std::vector< bound_mesh_t*> cpy_masks(mask_bits);
    if (WP.ncells() != 0)
    {
      cpy_masks.push_back(&WP);
      bool has_X = __flag_colliding_cells(z_mesh, cpy_masks, cpy_masks.size() - 1/*last*/, wdata, -IN);
      if (has_X)
        __flag_hidden_subzones(z_mesh, WP, wdata, IN, 2);// if there are at least 2 zones => mark as IN the hidden ones.
      cpy_masks.pop_back(); // WP appends only INs so remove it from list
    }
    if (WNP.ncells() != 0)
    {
      cpy_masks.push_back(&WNP);
      bool has_X = __flag_colliding_cells(z_mesh, cpy_masks, cpy_masks.size() - 1/*last*/, wdata, -X);
      if (has_X)
        __flag_hidden_subzones(z_mesh, WNP, wdata, X, 2);// if there are at least 2 zones => mark as IN the hidden ones.
    }

    // rearrange by decreasing mask box size (because a bigger mask potentially hides more cells) 
    STACK_ARRAY(int, nbits, maskidx);
    std::vector<std::pair<double, int>> palma;

    for (size_t i = 0; i < nbits; ++i)
    {
      if (mask_bits[i] == nullptr) continue;
      K_SEARCH::BBox3D mbx;
      mask_bits[i]->bbox(mbx);
      E_Float v = (mbx.maxB[0] - mbx.minB[0])*(mbx.maxB[1] - mbx.minB[1])*(mbx.maxB[2] - mbx.minB[2]);
      palma.push_back(std::make_pair(v, i));
    }

    std::sort(palma.rbegin(), palma.rend());
    
    wdata_t wdata0 = wdata;
    
    // nbits is unchanged => not treating again WP
    for (size_t ii = 0; ii < palma.size(); ++ii)
    {
      int i = palma[ii].second;

      wdata_t virgin_wdata= wdata0;

      // append z_xcelln with X
      bool has_X = __flag_colliding_cells(z_mesh, mask_bits, i, virgin_wdata, -X);
      if (!has_X) continue;

      // if there are at least 2 zones => mark as IN the hidden ones.
      __flag_hidden_subzones(z_mesh, *(mask_bits[i]), virgin_wdata, X, 2);

      // update wdata
      for (E_Int u = 0; u < ncells; ++u)
      {
        if (wdata[u] == (E_Int)IN) continue;
        if (virgin_wdata[u] == (E_Int)OUT) continue;

        E_Int val = virgin_wdata[u];
        data_trait<POLICY, zmesh_t>::mark_cell_w_mask(wdata, u, i, val);
      }
    }

    // set col_X color (collision mode) or polyclip (xcelln mode)
    outdata = this->__process_X_cells(z_mesh, cpy_masks, wdata);

    return 0;
  }
  
  ///
  TEMPLATE_TYPES
  void TEMPLATE_CLASS::finalize(zmesh_t const & z_mesh, outdata_t& outdata)
  {
    for (size_t i = 0; i < outdata.size(); ++i)
      outdata[i] = (outdata[i] == IN) ? 0. : (outdata[i] == OUT) ? 1. : outdata[i];

    //
    bool no_x = (*std::min_element(ALL(outdata)) == 1.);

    if (no_x)
    {
      K_SEARCH::BBox3D zbx;
      z_mesh.bbox(zbx);
      eClassify loc = NUGA::CLASSIFY::classify(zbx, _comp_boxes, _z_priorities);

      if (loc == IN)
      {
        E_Int sz = outdata.size();
        outdata.clear();
        outdata.resize(sz, 0.);
      }
    }
  }
  
  ///
  TEMPLATE_TYPES
  void TEMPLATE_CLASS::__compact_to_box
  (zmesh_t const & z_mesh, std::vector< bound_mesh_t*> & mask_bits, bound_mesh_t& WP, bound_mesh_t& WNP)
  {
    // zone reduction
    K_SEARCH::BBox3D z_box;
    z_mesh.bbox(z_box);
    if (typeid(zmesh_t) == typeid(pg_smesh_t))
      z_box.enlarge(0.01); //hack for 2D : eg. collar double wall can fall out of fuselage box

    int nmasks = mask_bits.size();
    for (E_Int m = 0; m < nmasks; ++m)
    {
      if (mask_bits[m] == nullptr) continue;
      // first coarse filtering based on brute force : localizers are not available yet
      // because we want to build them on a reduced set
      NUGA::selector::reduce_to_box(*mask_bits[m], z_box, true/*brute force*/);

      if (mask_bits[m]->ncells() == 0)
      {
        delete mask_bits[m]; mask_bits[m] = nullptr;
      }
    }

    NUGA::selector::reduce_to_box(WP, z_box, true/*brute force*/);
    NUGA::selector::reduce_to_box(WNP, z_box, true/*brute force*/);
  }

  ///
  TEMPLATE_TYPES
  void TEMPLATE_CLASS::__build_mask_bits
  (const std::vector<K_FLD::FloatArray> &mask_crds, const std::vector<K_FLD::IntArray>& mask_cnts,
   std::vector< std::vector<E_Int> > &mask_wall_ids,
   const std::vector<E_Int>& z_priorities, E_Int rank_wnp, std::vector< bound_mesh_t*> & mask_bits)
  {

  	mask_bits.clear();
  	//int nb_comps = mask_crds.size();
    
    // grabbing OP (WP are discarded) and WNP
  	for (size_t i=0; i <z_priorities.size(); ++i)
  	{
  	  int compi = z_priorities[i];
      
      //std::cout << "__build_mask_bits : comp " << compi << std::endl;
  	  
      bool z_is_prior_over_compi = ((E_Int)i < rank_wnp);

      //std::cout << "is a prioritary over this comp ? " << z_is_prior_over_compi << " rank/rank_wnp " << i << "/" << rank_wnp << std::endl;

  	  if (!z_is_prior_over_compi && mask_wall_ids[compi].empty()) continue; // no WNP
      
      //std::cout << "__build_mask_bits : nb walls : " << mask_wall_ids[compi].size() << std::endl;

  	  bound_mesh_t* bit = new bound_mesh_t(mask_crds[compi], mask_cnts[compi], 1/* ASSUME DIRECT UPEN ENTRY*/);
      
  	  E_Int nbcells = bit->ncells();
      
      //std::cout << "__build_mask_bits : 3 : nb of cells in mask : " << nbcells  << std::endl;

      // empty or only walls in it (WP are not required)
  	  bool discard_bit = ( (nbcells == 0) || (z_is_prior_over_compi && (nbcells == (E_Int)mask_wall_ids[compi].size()) ) );

  	  if (discard_bit) 
  	  {
#ifdef CLASSIFYER_DBG
        std::cout << "mask bit " << i << " is discarded" << std::endl;
#endif
  	  	delete bit; continue;
  	  }

      //nbcells = bit->ncells();
      //std::cout << "__build_mask_bits : 4 : nb of cells in mask : " << nbcells << std::endl;

#ifdef CLASSIFYER_DBG
    {
      std::ostringstream o;
      o << "mask_0_" << i;
      //medith::write<>(o.str().c_str(), bit->crd, bit->cnt);
    }
#endif

      //nbcells = bit->ncells();
      //std::cout << "__build_mask_bits : 5 : nb of cells in mask : " << nbcells << std::endl;

      // reduce to OP (discard WP) when prior comp, keep only WNP otherwise
  	  std::vector<bool> keep(nbcells, z_is_prior_over_compi);
  	  for (size_t u=0; u<mask_wall_ids[compi].size(); ++u )
        keep[mask_wall_ids[compi][u]]=!z_is_prior_over_compi;
      
  	  bit->compress(keep);
  	  if (bit->ncells() == 0) // completely gone
  	  {
  	  	delete bit; continue;
  	  }
      
      //nbcells = bit->ncells();
      //std::cout << "__build_mask_bits : 6 : nb of cells in mask : " << nbcells << std::endl;
      
      if (!z_is_prior_over_compi) // reverse WNPs
      {
        //std::cout << "reversing  : WNPs" << i << std::endl;
        if (bit->oriented == 1) bit->reverse_orient();
      }

  	  mask_bits.push_back(bit);
      
      //std::cout << "mask rank : " << mask_bits.size() -1 << std::endl;

#ifdef CLASSIFYER_DBG
    {
      //std::cout << "ouput mask_1_ " << i << std::endl;
      std::ostringstream o;
      o << "mask_1_" << i;
      //medith::write<>(o.str().c_str(), bit->crd, bit->cnt);
    }
#endif
  	}
    //std::cout << "__build_mask_bits : exit" << std::endl;
  }

  ///
  TEMPLATE_TYPES
  void TEMPLATE_CLASS::__build_mask_bits2
  (zmesh_t & z_mesh, double ARTOL, eMetricType mtype, double AMAX, // these parameters to deal with overlapping
   const std::vector<K_FLD::FloatArray> &mask_crds, const std::vector<K_FLD::IntArray>& mask_cnts,
   std::vector< std::vector<E_Int> > &mask_wall_ids,
   const std::vector<E_Int>& z_priorities, E_Int rank_wnp, std::vector< bound_mesh_t*> & mask_bits, bound_mesh_t& WP, bound_mesh_t& WNP)
  {
    mask_bits.clear();
    //int nb_comps = mask_crds.size();

    bound_mesh_t zbound;
    z_mesh.get_boundary(zbound);

    // separating WP/WNP and OVLP
    for (size_t i = 0; i <z_priorities.size(); ++i)
    {
      int compi = z_priorities[i];

      //std::cout << "__build_mask_bits : comp " << compi << std::endl;

      bool is_a_prior_mask = ((E_Int)i < rank_wnp);

      //std::cout << "is a prioritary mask ? " << is_a_prior_mask << " rank/rank_wnp " << i << "/" << rank_wnp << std::endl;

      if (!is_a_prior_mask && mask_wall_ids[compi].empty()) continue; // no WNP

      //std::cout << "__build_mask_bits : nb walls : " << mask_wall_ids[compi].size() << std::endl;

      bound_mesh_t* bit = new bound_mesh_t(mask_crds[compi], mask_cnts[compi], 1/* ASSUME DIRECT UPEN ENTRY*/);

      E_Int nbcells = bit->ncells();

      //std::cout << "__build_mask_bits : 3 : nb of cells in mask : " << nbcells  << std::endl;

      if (nbcells == 0)
      {
#ifdef CLASSIFYER_DBG
        std::cout << "mask bit " << i << " is discarded" << std::endl;
#endif
        delete bit; continue;
      }

      //nbcells = bit->ncells();
      //std::cout << "__build_mask_bits : 4 : nb of cells in mask : " << nbcells << std::endl;

#ifdef CLASSIFYER_DBG
      {
        std::ostringstream o;
        o << "mask_0_" << i;
        medith::write<>(o.str().c_str(), bit->crd, bit->cnt);
      }
#endif
      std::vector<bool> is_dw;
      NUGA::move_double_walls(bit, zbound, ARTOL, mtype, AMAX, mask_wall_ids[compi], is_dw);

#ifdef CLASSIFYER_DBG
      {
        //std::cout << "ouput mask_1_ " << i << std::endl;
        std::ostringstream o;
        o << "mask_1_" << i;
        medith::write<>(o.str().c_str(), bit->crd, bit->cnt);
      }
#endif

      //
      bound_mesh_t OVLP, WALL(*bit);
      
      std::vector<bool> keep;

      if (is_a_prior_mask)
      {
        OVLP = *bit;
        // remove WALL from OVLP
        keep.clear();
        keep.resize(nbcells, true);
        for (size_t u = 0; u < mask_wall_ids[compi].size(); ++u)
          keep[mask_wall_ids[compi][u]] = false;

        OVLP.compress(keep);
      }

      // remove OVLP and DW from WALL
      keep.clear();
      keep.resize(nbcells, false);
      for (size_t u = 0; u < mask_wall_ids[compi].size(); ++u)
      {
        E_Int wid = mask_wall_ids[compi][u];
        keep[wid] = (!is_dw.empty()) ? !is_dw[wid] : true; // is_dw is empty either in surface mode or no walls in volume mode for the current zone.
      }

      WALL.compress(keep);

      if (WALL.ncells() != 0)
      {
        if (!is_a_prior_mask)
        {
          WALL.reverse_orient();
          WNP.append(WALL);
          //bound_mesh_t* w = new bound_mesh_t(WALL);
          //mask_bits.push_back(w);

#ifdef CLASSIFYER_DBG
          {
            //std::cout << "ouput mask_1_ " << i << std::endl;
            std::ostringstream o;
            o << "wnp_" << i;
            //medith::write<>(o.str().c_str(), WALL.crd, WALL.cnt);
          }
#endif
        }
        else
        {
          WALL.reverse_orient();
          WP.append(WALL);
#ifdef CLASSIFYER_DBG
          {
            //std::cout << "ouput mask_1_ " << i << std::endl;
            std::ostringstream o;
            o << "wp_" << i;
            //medith::write<>(o.str().c_str(), WALL.crd, WALL.cnt);
          }
#endif
        }
      }

      *bit = OVLP;

      if (bit->ncells() == 0) // completely gone
      {
        delete bit; continue;
      }

      //nbcells = bit->ncells();
      //std::cout << "__build_mask_bits : 5 : nb of cells in mask : " << nbcells << std::endl;

      mask_bits.push_back(bit);

      //std::cout << "mask rank : " << mask_bits.size() -1 << std::endl;

#ifdef CLASSIFYER_DBG
      {
        //std::cout << "ouput mask_1_ " << i << std::endl;
        std::ostringstream o;
        o << "mask_2_" << i;
        medith::write<>(o.str().c_str(), bit->crd, bit->cnt);
      }
#endif
    }

#ifdef CLASSIFYER_DBG
    medith::write<>("WNP", WNP.crd, WNP.cnt);
    medith::write<>("WP", WP.crd, WP.cnt);
#endif

    //std::cout << "__build_mask_bits : exit" << std::endl;
  }

  ///
  TEMPLATE_TYPES
  void TEMPLATE_CLASS::__process_overlapping_boundaries
  (zmesh_t & z_mesh, std::vector< bound_mesh_t*> & mask_bits, E_Int rank_wnp, E_Float RTOL)
  {
  	// WARNING : assume consistent orientation

    // get zmesh boundary (not using the constructor to do so because we need the ancestor info)
    bound_mesh_t zbound;

    std::vector<E_Int> ancestor;
  	z_mesh.get_boundary(zbound, ancestor);

#ifdef CLASSIFYER_DBG
    {
      std::ostringstream o;
      o << "m";
      medith::write(o.str().c_str(), z_mesh.crd, z_mesh.cnt);
    }
    {
      std::ostringstream o;
      o << "bound";
      medith::write(o.str().c_str(), zbound.crd, zbound.cnt);
    }
#endif

  	std::vector<E_Int> is_x1, is_x2;

#ifdef CLASSIFYER_DBG
    bool has_abut_ovlp = false;
#endif
    for (size_t m=0; m < mask_bits.size(); ++m)
    {
      is_x1.clear();
      is_x2.clear();

      if (mask_bits[m] == nullptr) continue;

#ifdef CLASSIFYER_DBG
    {
      std::ostringstream o;
      o << "mask_" << m;
      medith::write<>(o.str().c_str(), mask_bits[m]->crd, mask_bits[m]->cnt);
    }
#endif
      
      bound_mesh_t& maski = *mask_bits[m];
      using elt_t = typename bound_mesh_t::elt_t;
      using loc_t = typename bound_mesh_t::loc_t;

      COLLIDE::compute_overlap<elt_t, elt_t, loc_t>(zbound.crd, zbound.cnt, 
      	                                            maski.crd, maski.cnt, *(maski.get_localizer()), 
                                                    is_x1, is_x2,
                                                    RTOL);

      // flag the attached PG of z_mesh as IN
      std::vector<E_Int> ids;
      for (size_t u=0; u< is_x1.size(); ++u) 
      {
        E_Int maskid = is_x1[u];
        if (maskid == IDX_NONE) continue;
        
        bool is_abutting = (maskid < 0);
        maskid = ::abs(maskid) - 1; //go 0-based positive

      	bool appendit = (is_abutting && (E_Int)m >=  rank_wnp); // test on rank : ie. is a WNP
      	appendit     |= (!is_abutting  && (E_Int)m <   rank_wnp); // test on rank : ie. is not a WNP, it is as OVLP
      	if (appendit) ids.push_back(ancestor[u]);
      }
      if (!ids.empty()) 
      {
        //z_mesh.set_type(IN, cids);
#ifdef CLASSIFYER_DBG
        has_abut_ovlp = true;
#endif
      }
      
      // discard overlap elts in masks
      E_Int nbcells = mask_bits[m]->ncells();
      assert (is_x2.size() == nbcells);
      //std::cout << "nbcells/is_x2 size : " << nbcells << "/" << is_x2.size() << std::endl;
      std::vector<bool> keep(nbcells, true);
  	  for (E_Int u=0; u<nbcells; ++u )
  	   	if (is_x2[u] != IDX_NONE) keep[u] = false;

#ifdef CLASSIFYER_DBG
      {
        std::ostringstream o;
        o << "mask_" << m << "_colored";
        medith::write<>(o.str().c_str(), mask_bits[m]->crd, mask_bits[m]->cnt, nullptr, 0, &keep);
      }
#endif
  	  
  	  mask_bits[m]->compress(keep);
  	  
#ifdef CLASSIFYER_DBG
      {
        std::ostringstream o;
        o << "mask_" << m << "_compressed";
        medith::write<>(o.str().c_str(), mask_bits[m]->crd, mask_bits[m]->cnt);
      }
#endif
      
      if (mask_bits[m]->ncells() == 0) // completely gone
  	  {
  	  	delete mask_bits[m]; mask_bits[m]=nullptr;
  	  }
    }

    // compress (so keep the sorting) the masks to the remaining ones
    std::vector< bound_mesh_t*> remaining_masks;
    for (size_t m=0; m < mask_bits.size(); ++m)
      if (mask_bits[m] != nullptr)
      	remaining_masks.push_back(mask_bits[m]);
    
   mask_bits = remaining_masks;

#ifdef CLASSIFYER_DBG
    if (has_abut_ovlp)
    {
      std::cout << "has ABUTTING/OVERSET" << std::endl;
      std::ostringstream o;
      o << "m_in";
      medith::write(o.str().c_str(), z_mesh.crd, z_mesh.cnt,nullptr, 0, z_mesh.e_type.empty() ? nullptr : &z_mesh.e_type);
    }
#endif
  }

  ///
  TEMPLATE_TYPES
  bool TEMPLATE_CLASS::__flag_colliding_cells
  (zmesh_t const & z_mesh, std::vector< bound_mesh_t*> const & mask_bits, E_Int im, wdata_t& data, E_Int XVAL)
  {
    assert ((size_t)im < mask_bits.size());
    const bound_mesh_t* pmask = mask_bits[im];
    if (pmask == nullptr) return false;
    if (pmask->ncells() == 0) return false;

    bool has_X = false;
    bound_mesh_t const & mask_bit = *pmask;
    
#ifdef CLASSIFYER_DBG
      medith::write<>("currentz", z_mesh.crd, z_mesh.cnt);
      medith::write<>("currentm", mask_bit.crd, mask_bit.cnt);
#endif
  
    E_Int nbcells = z_mesh.ncells();

    using loc_t = typename zmesh_t::loc_t;
    const loc_t* l= mask_bit.get_localizer();
    assert (l != nullptr);
    const loc_t& mask_loc = *l;

    mask_bit.get_nodal_metric2();
    z_mesh.get_nodal_metric2();
  
    std::vector<E_Int> cands;
    
    for (size_t i = 0; i < (size_t)nbcells; ++i)
    {
      //std::cout << i << " over " << nbcells << std::endl;
      assert (i < data.size());
      // in any POLICY, X must be re-processed to avoid missing element in current X-front being computed
      if (data[i] == (E_Int)IN) continue;

      // autonomous element
      auto ae1 = z_mesh.aelement(i);

      cands.clear();
      mask_loc.get_candidates(ae1, ae1.m_crd, cands, 1, _RTOL); //return as 1-based
      if (cands.empty()) continue;

      std::sort(ALL(cands));//fixme : here because trait::compress_ use extract_by_predicate which is not a compressing function

#ifdef CLASSIFYER_DBG
      //medith::write<ngon_type>("subj", z_mesh.crd, z_mesh.cnt, i);
      //medith::write("cands", mask_bit.crd, mask_bit.cnt, &cands, 1);
#endif

      bool is_x = NUGA::COLLIDE::get_colliding(ae1, mask_bit, cands, 1, _RTOL, true/*returns at first found*/);

#ifdef CLASSIFYER_DBG
      //medith::write("first_x", mask_bit.crd, mask_bit.cnt, &cands, 1);
#endif
      
      if (is_x) // cands[0] contains the first found
      {
#ifdef CLASSIFYER_DBG
        /*medith::write<ngon_type>("ae1", z_mesh.crd, z_mesh.cnt, i);
        medith::write<ngon_type>("mcan", mask_bit.crd, mask_bit.cnt, cands[0]-1);
        auto cand = mask_bit.aelement(cands[0] - 1);
        K_SEARCH::BBox3D bc(cand.m_crd);
        bc.enlarge(_RTOL);
        ngon_type ng2;
        K_FLD::FloatArray c2;
        bc.convert2NG(c2, ng2);
        medith::write("bc", c2, ng2);*/
#endif

        data_trait<POLICY, zmesh_t>::mark_cell_w_mask(data, i, im, XVAL);
        //z_mesh.set_flag(i, cands[0]);
      }
      
      has_X |= is_x;

    }
    
#ifdef CLASSIFYER_DBG
    std::vector<E_Int> xs, dat(data.size());
    
    for (size_t i=0; i < data.size(); ++i)
    {
      if (data[i] == XVAL) xs.push_back(i);
      dat[i] = ::fabs(data[i]);//for medit
    }
    if (xs.empty())
      std::cout << "NO COLLISIONS WITH CURRENT MASK" << std::endl;
    else
    {
      medith::write("colliding_set", z_mesh.crd, z_mesh.cnt, &xs);
      medith::write("flag_collided_zone_cells", z_mesh.crd, z_mesh.cnt, nullptr, 0, &dat);
    }
#endif

    return has_X;
  }

  ///
  TEMPLATE_TYPES
  E_Int TEMPLATE_CLASS::__flag_hidden_subzones
  (zmesh_t const & z_mesh, bound_mesh_t const & mask_bit, wdata_t& z_xcelln, E_Int XCOL/*X or IN*/, E_Int nsubmin)
  {
    // 1. compute neighborhood (lazy mode)
    auto neighborz = z_mesh.get_neighbors();
    
  	// 2. incremental coloring, starting from UPPER_COL
    
    // 2.1 : current field with only IN cells and new Xs (-X)
    std::vector<E_Float> cur_xcelln(z_xcelln.size(), OUT);

    for (size_t i=0; i < z_xcelln.size(); ++i)
    {
      if (z_xcelln[i] == (E_Int)IN) cur_xcelln[i] = IN;
      else if (z_xcelln[i] == -XCOL) // collision colors : negative to say 'new'
      {
        z_xcelln[i] = XCOL; //reset on z_xcelln
        cur_xcelln[i] = XCOL;
      }
    }
    
#ifdef CLASSIFYER_DBG
    medith::write("initial_field_coloring", z_mesh.crd, z_mesh.cnt, nullptr, 0, &cur_xcelln);
#endif

    // 2.2 : incremental coloring => new INs, some X
    assert (UPPER_COL > OUT);
    NUGA::EltAlgo<typename zmesh_t::elt_t>::coloring(*neighborz, cur_xcelln, (E_Float)OUT, (E_Float)UPPER_COL);
    
#ifdef CLASSIFYER_DBG
    medith::write("colored_field_coloring", z_mesh.crd, z_mesh.cnt, nullptr, 0, &cur_xcelln);
#endif

    // 2.3 : update z_xcelln with IN & X: RULE : IN > X > OUT
    for (size_t i=0; i < z_xcelln.size(); ++i)z_xcelln[i] = std::min(double(z_xcelln[i]), cur_xcelln[i]);
    // 2.4 : check nb of subzones. STOP if < 2
    E_Int nb_subzones = *std::max_element(ALL(cur_xcelln)) - UPPER_COL + 1;

#ifdef CLASSIFYER_DBG
    // std::set<int> all_cols(ALL(cur_xcelln));
    // for (auto c : all_cols)
    //   std::cout << "colors : " << c << std::endl;
    // std::cout << "nb of supposed sub zones: " <<nb_subzones << std::endl;
#endif

    if (nb_subzones < nsubmin) // no new subzones, just return
      return 0;
    
    // 2.5 : classify subzones
    STACK_ARRAY(eClassify, nb_subzones, z_color);
    for (E_Int i=0; i < nb_subzones; ++i)z_color[i] = AMBIGUOUS;
    
#ifdef CLASSIFYER_DBG
    std::cout << "NB SUZONES : " << nb_subzones << std::endl;
#endif
    
    E_Int ncell = z_mesh.ncells();
    E_Int missing_col=nb_subzones;

    using loc_t = typename zmesh_t::loc_t;
    const loc_t& mask_loc = *(mask_bit.get_localizer());

    std::vector<E_Int> cands;
    short iter(0);
    
    while (iter++ < 2 && missing_col)
    {
      bool deep = (iter == 1) ? false : true;
      for (E_Int i = 0; (i < ncell) && missing_col; ++i)
      {
        if (cur_xcelln[i] != XCOL) continue;

        E_Int nneighs = neighborz->stride(i);
        const E_Int* pneighs = neighborz->begin(i);

        for (E_Int j = 0; (j < nneighs); ++j)
        {
          if (pneighs[j] == IDX_NONE) continue;

          E_Int subid = cur_xcelln[pneighs[j]];

          // discard element that is already colored
          if (subid == IDX_NONE) continue;
          if (subid == X)        continue;
          if (subid == IN)       continue;

          subid -= UPPER_COL; // to be 0-based
          assert(subid > -1 && subid < nb_subzones);

          //std::cout << "subid " <<  subid << ": z_color[subid] : " << z_color[subid] << std::endl;

          if (z_color[subid] != AMBIGUOUS) continue;

          // subid is a unvisited subzone

          // GET MOLECULE of i, check position of pneighs[j] regarding that molecule

          // a. get molecule of i
          auto ae1 = z_mesh.aelement(i);
          cands.clear();
          mask_loc.get_candidates(ae1, ae1.m_crd, cands, 1, _RTOL); //return as 0-based (fixme for volumic, was 1-based)
          if (cands.empty()) continue;

#ifdef CLASSIFYER_DBG
          medith::write("cands", mask_bit.crd, mask_bit.cnt, &cands, 1);
#endif

          /*bool is_x = */NUGA::COLLIDE::get_colliding(ae1, mask_bit, cands, 1, _RTOL, false/*i.e reduce cands to true collidings*/);
          if (cands.empty()) continue;

          auto aen = z_mesh.aelement(pneighs[j]);

#ifdef CLASSIFYER_DBG
          //medith::write("colliding", ae1);
          medith::write("xmolecule", mask_bit.crd, mask_bit.cnt, &cands, 1);
          medith::write("subj", aen);
#endif

          // autonomous cutter front
          bound_mesh_t acut_front(mask_bit, cands, 1);

          // b. classify pneigh[j] with that molecule
          
          z_color[subid] = CLASSIFY::classify(aen, acut_front, deep);
          if (z_color[subid] != AMBIGUOUS)
          {
            --missing_col;

#ifdef CLASSIFYER_DBG
            std::cout << "subid : " << subid << " is in ? " << (z_color[subid] == IN) << std::endl;
#endif
          }
        }
      }
    }

    
    //if (missing_col != 0)
    //{
      // happening for unclosed portion (previously colored X) that have nothing to do with current mask
      // NOTHING TO DO
    //}
    
    // now update z_xcelln
    for (size_t i=0; i < z_xcelln.size(); ++i)
    {
      E_Int subid = cur_xcelln[i];
      
      if (subid == IDX_NONE)  continue;
      if (subid == X)         continue;
      if (subid == IN)        continue;
      
      subid -= UPPER_COL; // to be 0-based
      
      if (z_color[subid] == IN) z_xcelln[i] = IN;
    }

#ifdef CLASSIFYER_DBG
    static int counter = 0;
    std::ostringstream o;
    o << "flag_hidden_subzones_" << counter++;
    medith::write(o.str().c_str(), z_mesh.crd, z_mesh.cnt, nullptr, 0, &z_xcelln);
#endif

    return 0;
  }
}

#endif

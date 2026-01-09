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

#ifndef NUGA_SUBDIVISION_HXX
#define NUGA_SUBDIVISION_HXX

#include<vector>
#include "Nuga/include/macros.h"
#include "Nuga/include/subdiv_defs.h"

#include "Nuga/include/Tetrahedron.h"
#include "Nuga/include/Hexahedron.h"
#include "Nuga/include/Prism.h"
#include "Nuga/include/Pyramid.h"
#include "Nuga/include/Basic.h"
#include "Nuga/include/Polyhedron.h"

#include "Nuga/include/Triangle.h"
#include "Nuga/include/Quadrangle.h"
#include "Nuga/include/Polygon.h"


namespace NUGA
{

  //
  template <typename ELT_t, eSUBDIV_TYPE STYPE>
  struct subdiv_pol;

  //
  template <>
  struct subdiv_pol<K_MESH::Hexahedron, ISO>
  {
    enum { PGNBC = 4, PHNBC = 8, NBI = 12 };

    using ph_arr_t = K_FLD::IntArray;
    using pg_arr_t = K_FLD::IntArray;
  };

  //
  template <>
  struct subdiv_pol<K_MESH::Tetrahedron, ISO>
  {
    enum { PGNBC = 4, PHNBC = 8, NBI = 8 };

    using ph_arr_t = K_FLD::IntArray;
    using pg_arr_t = K_FLD::IntArray;
  };

  //
  template <>
  struct subdiv_pol<K_MESH::Prism, ISO>
  {
    enum { PGNBC = 4, PHNBC = 8, NBI = 10 }; // NBI : 4 T3 + 6 Q4

    using ph_arr_t = K_FLD::IntArray;
    using pg_arr_t = K_FLD::IntArray;
  };

  //
  template <>
  struct subdiv_pol<K_MESH::Pyramid, ISO>
  {
    enum { PGNBC = 4, PHNBC = 10, NBI = 13 }; // NBI : 12 T3 + 1 Q4

    using ph_arr_t = K_FLD::IntArray;
    using pg_arr_t = K_FLD::IntArray;

  };

  // Basic - ISO
  template <>
  struct subdiv_pol<K_MESH::Basic, ISO>
  {
    enum { PGNBC = -1, PHNBC = -1/*, NBI = -1 */ };

    using ph_arr_t = ngon_unit;
    using pg_arr_t = K_FLD::IntArray;
  };

  // ISO_HEX Poyhedron subdivision => N HEXA children , with N is the nb of nodes
  template <>
  struct subdiv_pol<K_MESH::Polyhedron<0>, ISO_HEX>
  {
    enum { PGNBC = -1, PHNBC = -1/*, NBI = -1 */ };

    using ph_arr_t = ngon_unit;
    using pg_arr_t = ngon_unit;

    static E_Int nbc(const ngon_unit& PGS, const E_Int* first_pg, E_Int nb_pgs)
    {
      std::vector<E_Int> unodes;
      K_MESH::Polyhedron<UNKNOWN>::unique_nodes(PGS, first_pg, nb_pgs, unodes);
      return (E_Int)unodes.size();
    }

    static void nbc_list(const ngon_type& ng, const std::vector<E_Int>& PHlist, std::vector<E_Int>& pregnant)
    {
      E_Int nbPHtoadapt = PHlist.size();
      pregnant.resize(nbPHtoadapt);

      E_Int PHl, nb_pgs;
      const E_Int* first_pg;
      for (E_Int l = 0; l < nbPHtoadapt; l++)
      {
        PHl = PHlist[l];
        nb_pgs = ng.PHs.stride(PHl);
        first_pg = ng.PHs.get_facets_ptr(PHl);
        pregnant[l] = nbc(ng.PGs, first_pg, nb_pgs);
      }
    }

    static E_Int nbi(const ngon_unit& PGS, const E_Int* first_pg, E_Int nb_pgs)
    {
      return (K_MESH::Polyhedron<UNKNOWN>::cumulated_arity(PGS, first_pg, nb_pgs) / 2);
    }

    static E_Int nbi_sum(const ngon_type& ng, const std::vector<E_Int>& PHlist)
    {
      E_Int sum(0);
      E_Int nbPHtoadapt = PHlist.size();
      E_Int PHl, nb_pgs;
      const E_Int* first_pg;
      for (E_Int l = 0; l < nbPHtoadapt; l++)
      {
        PHl = PHlist[l];
        nb_pgs = ng.PHs.stride(PHl);
        first_pg = ng.PHs.get_facets_ptr(PHl);

        sum = sum + nbi(ng.PGs, first_pg, nb_pgs);
      }
      return sum;
    }
  };

  // DIR for HEXA
  template <>
  struct subdiv_pol < K_MESH::Hexahedron, DIR >
  {
    enum { PGNBC = -1, PHNBC = -1/*, NBI = -1 */ };

    using ph_arr_t = ngon_unit;
    using pg_arr_t = ngon_unit;

    static E_Int nbi_sum(const ngon_type& ng, const std::vector<E_Int>& PHlist, const std::vector<eDIR>& PH_directive)
    {
      E_Int sum(0);
      E_Int nbPHtoadapt = PHlist.size();
      for (E_Int l = 0; l < nbPHtoadapt; l++)
      {
        //PHl = PHlist[l];
        const auto& dir = PH_directive[l];
        if (dir == XYZ) sum += 12;
        else if (dir == XY || dir == XZ || dir == YZ) sum += 4;
        else if (dir == Xd || dir == Y || dir == Z) sum += 1;
      }
      return sum;
    }

    static void nbc_list(const ngon_type& ng, const std::vector<E_Int>& PHlist, const std::vector<eDIR>& PH_directive, std::vector<E_Int>& pregnant)
    {
      E_Int nbPHtoadapt = PHlist.size();

      pregnant.clear();
      pregnant.resize(nbPHtoadapt);

      //E_Int PHl, nb_pgs;
      for (E_Int l = 0; l < nbPHtoadapt; l++)
      {
        //PHl = PHlist[l];
        const auto& dir = PH_directive[l];

        if (dir == XYZ) pregnant[l] = 8;
        else if (dir == XY || dir == XZ || dir == YZ) pregnant[l] = 4;
        else if (dir == Xd || dir == Y || dir == Z) pregnant[l] = 2;
      }
    }
  };

  // DIR_PROTO for HEXA
  /// same impl for both so define by inheritance
  template <>
  struct subdiv_pol < K_MESH::Hexahedron, DIR_PROTO > : public subdiv_pol < K_MESH::Hexahedron, DIR>
  {
  };

  // isotropic HEXA subdivision => 4 Quadrangles children => fixed stride array
  template <>
  struct subdiv_pol<K_MESH::Quadrangle, ISO>
  {
    enum { NBC = 4 };
    using pg_arr_t = K_FLD::IntArray;

    static void reorder_children(E_Int* child, E_Int nchildren/*dummy*/, bool reverse, E_Int i0)
    {
      K_CONNECT::IdTool::right_shift<4>(&child[0], i0);
      if (reverse)
        std::swap(child[1], child[3]);
    }

  };

  // isotropic Triangle subdivision => 4 Trianges children => fixed stride array
  template <>
  struct subdiv_pol<K_MESH::Triangle, ISO>
  {
    enum { NBC = 4 };
    using pg_arr_t = K_FLD::IntArray;

    static void reorder_children(E_Int* child, E_Int nchildren/*dummy*/, bool reverse, E_Int i0)
    {
      K_CONNECT::IdTool::right_shift<3>(&child[0], i0);
      if (reverse)
        std::swap(child[1], child[2]);
    }
  };

  // directional QUAD subdivision => 2 Quadrangles children => fixed stride array
  template <>
  struct subdiv_pol<K_MESH::Quadrangle, DIR>
  {
    using pg_arr_t = ngon_unit;

    static void reorder_children(E_Int* child, E_Int nchildren/*dummy*/, bool reverse, E_Int i0)
    {
      // has nothing to do : with 2 children, the sorting must be kept
    }

    static E_Int nbc_list(const ngon_unit& PGs, const std::vector<E_Int>& PGlist, const std::vector<eDIR>& PG_directive, std::vector<E_Int>& pregnant)
    {
      pregnant.clear();
      pregnant.resize(PGlist.size());
      for (size_t i = 0; i < PG_directive.size(); ++i)
      {
        if (PG_directive[i] == XY) pregnant[i] = 4;
        else pregnant[i] = 2;
      }

      return 0;
    }

  };

  /// same impl for both so define by inheritance
  template <>
  struct subdiv_pol<K_MESH::Quadrangle, DIR_PROTO> : public subdiv_pol<K_MESH::Quadrangle, DIR>
  {};

  // directional Triangle subdivision => 2 Trianges children => fixed stride array
  template <>
  struct subdiv_pol<K_MESH::Triangle, DIR>
  {
    using pg_arr_t = ngon_unit;

    static void reorder_children(E_Int* child, E_Int nchildren/*dummy*/, bool reverse, E_Int i0)
    {
      K_CONNECT::IdTool::right_shift<2>(&child[0], i0);
      if (reverse)
        std::swap(child[0], child[1]);
    }
  };

  /// same impl for both so define by inheritance
  template <>
  struct subdiv_pol<K_MESH::Triangle, DIR_PROTO> : public subdiv_pol<K_MESH::Triangle, DIR>
  {};

  // ISO_HEX Polygon subdivision => N quad children , with N is the nb of nodes
  template <>
  struct subdiv_pol<K_MESH::Polygon, ISO_HEX>
  {
    enum { PGNBC = -1, PHNBC = -1/*, NBI = -1 */ };
    using pg_arr_t = ngon_unit;

    static void reorder_children(E_Int* child, E_Int nchildren, bool reverse, E_Int i0)
    {
      K_CONNECT::IdTool::right_shift(&child[0], nchildren, i0);
      if (reverse)
        std::reverse(child, child+nchildren);
    }

    static E_Int nbc_list(const ngon_unit& PGs, const std::vector<E_Int>& PGlist, const std::vector<eDIR>& PG_directive, std::vector<E_Int>& pregnant)
    {
      pregnant.clear();
      pregnant.resize(PGlist.size());

      for (size_t i = 0; i < PGlist.size(); ++i)
      {
        E_Int PGi = PGlist[i];
        pregnant[i] = PGs.stride(PGi);
      }
      return 0;
    }

  };

  template <short DIM>
  struct int_tuple
  {
    E_Int n[DIM];

    void set(E_Int val) 
    {
      for (E_Int i = 0; i < DIM; i++)
        n[i] = val;
    }

    int_tuple() { n[0] = n[1] = 0; if (DIM == 3) n[2] = 0; }
    explicit int_tuple(E_Int val) { n[0] = n[1] = val; if (DIM == 3) n[2] = val; }
    ~int_tuple() {};
    
    int_tuple& operator=(E_Int val) { n[0] = n[1] = val; if (DIM == 3) n[2] = val; return *this; }
    int_tuple& operator=(const int_tuple& d) { n[0] = d.n[0];  n[1] = d.n[1]; if (DIM == 3) n[2] = d.n[2]; return *this; }
    //E_Int operator+(E_Int v) const { return max() + v; }
    int_tuple& operator+(E_Int val) { n[0] += val;  n[1] += val; if (DIM == 3) n[2] += val; return *this; }
    int_tuple& operator--()
    {
        n[0] = std::max((E_Int)0, n[0]-1);
        n[1] = std::max((E_Int)0, n[1]-1);
        if (DIM == 3)
            n[2] = std::max((E_Int)0, n[2]-1);
        return *this;
    }

    bool operator>=(E_Int v) const { return (max() >= v); }
    bool operator<=(E_Int v) const { return (max() <= v); }
    bool operator>(E_Int v) const { return (max() > v); }
    bool operator==(E_Int v) const { return (max() == v) && (min() == v); }
    bool operator!=(E_Int val) const { // WARNING : weird logic unless val is 0
      if (n[0] != val) return true;
      if (n[1] != val) return true;
      if ((DIM == 3) && (n[2] != val)) return true;
      return false;
    }
	bool operator==(const int_tuple& d) const
    {
      for (E_Int i = 0; i < DIM; i++)
        if (n[i] != d.n[i]) return false;
      return true;
    }
    bool operator!=(const int_tuple& d) const
    {
      return !(*this == d);
    }

    //DANGEROUS because weird logic
    bool operator>(const int_tuple& d) {
      if (n[0] > d.n[0]) return true;
      if (n[1] > d.n[1]) return true;
      if ((DIM == 3) && (n[2] > d.n[2])) return true;
      return false;
    }
    int_tuple& operator+=(E_Int val) { n[0] = std::max(n[0] + val, (E_Int)0); n[1] = std::max(n[1] + val, (E_Int)0);; if (DIM == 3) n[2] = std::max(n[2] + val, (E_Int)0); return *this; }
    int_tuple& operator+=(const int_tuple& d) { n[0] += d.n[0];  n[1] += d.n[1]; if (DIM == 3)  n[2] += d.n[2]; return *this; }

    int_tuple& operator/=(E_Int val) { n[0] /= val; n[1] /= val; if (DIM == 3) n[2] /= val; return *this; }
  
    int_tuple operator-(const int_tuple& d) { int_tuple res(0);  res.n[0] = n[0] - d.n[0]; res.n[1] = n[1] - d.n[1]; if (DIM == 3) res.n[2] = n[2] - d.n[2]; return res; }
    int_tuple operator-(E_Int d) const
    {
      int_tuple res(0);
      res.n[0] = std::max(n[0] - d, (E_Int)0);
      res.n[1] = std::max(n[1] - d, (E_Int)0);
      if (DIM == 3)
        res.n[2] = std::max(n[2] - d, (E_Int)0);
      return res;
    }

    E_Int max() const {
      if (DIM == 3) return std::max(n[0], std::max(n[1], n[2]));
      else return std::max(n[0], n[1]);
    }

    E_Int min() const {
      if (DIM == 3) return std::min(n[0], std::min(n[1], n[2]));
      else return std::min(n[0], n[1]);
    }
  };

  template <short DIM> inline int_tuple<DIM> max(int_tuple<DIM>&d, E_Int v) { int_tuple<DIM> res(0); for (size_t k=0; k < DIM; ++k) res.n[k] = std::max(d.n[k], v); return res; }//hack fr CLEF : l.362(hmesh.xhh)
  inline int_tuple<3> abs(int_tuple<3> d) { int_tuple<3> res(0);  res.n[0] = ::abs(d.n[0]); res.n[1] = ::abs(d.n[1]); res.n[2] = ::abs(d.n[2]); return res; }
  inline int_tuple<3> max(int_tuple<3> a, int_tuple<3> b) { int_tuple<3> res(0); res.n[0] = std::max(a.n[0], b.n[0]); res.n[1] = std::max(a.n[1], b.n[1]); res.n[2] = std::max(a.n[2], b.n[2]); return res; }
  inline int_tuple<3> min(int_tuple<3> a, int_tuple<3> b) { int_tuple<3> res(0); res.n[0] = std::min(a.n[0], b.n[0]); res.n[1] = std::min(a.n[1], b.n[1]); res.n[2] = std::min(a.n[2], b.n[2]); return res; }


  template <short DIM> inline std::ostream &operator<<(std::ostream& out, const int_tuple<DIM>& d)
  {
    out << d.n[0] << "/" << d.n[1];

    if (DIM == 3)
      out << "/" << d.n[2];
    out << std::endl;

    return out;
  }

 
  ///
  template <eSUBDIV_TYPE STYPE> // ISO impl
  struct adap_incr_type
  {
    using cell_incr_t = E_Int;
    using face_incr_t = E_Int;
    Vector_t<E_Int> cell_adap_incr;
    Vector_t<E_Int> face_adap_incr;

    E_Int cmin(E_Int k) { ASSERT_IN_VECRANGE(cell_adap_incr, k); return cell_adap_incr[k]; }
    E_Int cmax(E_Int k) { ASSERT_IN_VECRANGE(cell_adap_incr, k); return cell_adap_incr[k]; }
    E_Int fmax(E_Int k) { ASSERT_IN_VECRANGE(face_adap_incr, k); return face_adap_incr[k]; }

    NUGA::eDIR get_face_dir(E_Int k) const
    {
      ASSERT_IN_VECRANGE(face_adap_incr, k);

      if (face_adap_incr[k] != 0) return XY;
      return NONE;
    }

  };

  ///
  template<> // DIR 
  struct adap_incr_type<DIR>
  {
    using cell_incr_t = int_tuple<3>;
    using face_incr_t = int_tuple<2>;
    Vector_t<cell_incr_t> cell_adap_incr;
    Vector_t<face_incr_t> face_adap_incr;

    E_Int cmin(E_Int k) { ASSERT_IN_VECRANGE(cell_adap_incr, k); return cell_adap_incr[k].min(); }
    E_Int cmax(E_Int k) { ASSERT_IN_VECRANGE(cell_adap_incr, k); return cell_adap_incr[k].max(); }
    E_Int fmax(E_Int k) { ASSERT_IN_VECRANGE(face_adap_incr, k); return face_adap_incr[k].max(); }

    NUGA::eDIR get_face_dir(E_Int k) const
    {
      ASSERT_IN_VECRANGE(face_adap_incr, k);
      if ((face_adap_incr[k].n[0] != 0) && (face_adap_incr[k].n[1] != 0)) return XY;
      if (face_adap_incr[k].n[0] != 0) return Xd;
      if (face_adap_incr[k].n[1] != 0) return Y;
      return NONE;
    }

  };

  /// same impl for both so define by inheritance
  template<> 
  struct adap_incr_type<DIR_PROTO> : public adap_incr_type<DIR>
  {
  };



} // NUGA

#endif

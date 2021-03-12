/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

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
  enum {PGNBC=4, PHNBC=8, NBI=12};

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
  enum { PGNBC=4, PHNBC=10, NBI=13 }; // NBI : 12 T3 + 1 Q4

  using ph_arr_t = K_FLD::IntArray;
  using pg_arr_t = K_FLD::IntArray;

};

// Basic - ISO
template <>
struct subdiv_pol<K_MESH::Basic, ISO>
{
  enum { PGNBC = -1, PHNBC = -1/*, NBI = -1 */};

  using ph_arr_t = ngon_unit;
  using pg_arr_t = K_FLD::IntArray;
};

// ISO_HEX Poyhedron subdivision => N HEXA children , with N is the nb of nodes
template <>
struct subdiv_pol<K_MESH::Polyhedron<0>, ISO_HEX>
{
  enum { PGNBC = -1, PHNBC = -1/*, NBI = -1 */};

  using ph_arr_t = ngon_unit;
  using pg_arr_t = ngon_unit;

  static E_Int nbc(const ngon_unit& PGS, const E_Int* first_pg, E_Int nb_pgs)
  {
    //todo JP 
    // <=> nb of nodes : method uniques_nodes in K_MESH::Polyhdron
    return 0;
  }

  static void nbc_list(const ngon_type& ng, const std::vector<E_Int>& PHlist, std::vector<E_Int>& pregnant)
  {
    //todo JP 
    // loop on Polyhedra given in PHlist and call above nbc function and set pregnant upon exit
  }

  static E_Int nbi(const ngon_unit& PGS, const E_Int* first_pg, E_Int nb_pgs)
  {
    //todo JP
    // <=> nb of internal faces : half of the cumulated node arity
    // method : cumulated_arity in K_MESH::Polyhedron
    return 0;
  }

  static E_Int nbi_sum(const ngon_type& ng, const std::vector<E_Int>& PHlist)
  {
    E_Int sum(0);
    //todo JP
    // loop on Polyhedra given in PHlist and call above nbi function and return accumulated value
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
};

// isotropic HEXA subdivision => 4 Quadrangles children => fixed stride array
template <>
struct subdiv_pol<K_MESH::Quadrangle, ISO>
{
  enum { NBC = 4};
  using pg_arr_t = K_FLD::IntArray;

  static void reorder_children(E_Int* child, bool reverse, E_Int i0)
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

  static void reorder_children(E_Int* child, bool reverse, E_Int i0)
  {
    K_CONNECT::IdTool::right_shift<3>(&child[0], i0);
    if (reverse)
      std::swap(child[1], child[2]);
  }
};

// ISO_HEX Polygon subdivision => N quad children , with N is the nb of nodes
template <>
struct subdiv_pol<K_MESH::Polygon, ISO_HEX>
{
  enum { NBC = -1};
  using pg_arr_t = ngon_unit;

  static E_Int nbc_list(const ngon_unit& PGs, const std::vector<E_Int>& PGlist, std::vector<E_Int>& pregnant)
  {
    //todo JP 
    // loop on polygons given in PGlist and set its stride in pregnant: pregnant[i] = PGs.stride(PGlist[i])
    return 0;
  }

};
  
struct incr_type
{
  using cell_output_type = Vector_t<E_Int>;
  Vector_t<E_Int> face_adap_incr, cell_adap_incr;
};

// Directional subdivision data
struct dir_incr_type //fixme : separate also here face & cell adap incr
{
  Vector_t<E_Int>      _adap_incr;
  Vector_t<NUGA::eDIR> _ph_dir, _pg_dir;
  
  E_Int& operator[](E_Int i) {return _adap_incr[i]; };
  
  void clear() {_adap_incr.clear();}
  void resize(E_Int sz, E_Int val) {_adap_incr.resize(sz, val);}
  Vector_t<E_Int>::iterator begin() {return _adap_incr.begin();}
  Vector_t<E_Int>::iterator end() {return _adap_incr.end();}
};

  // Isotropic subdivision data : Vector_t<E_Int>
  template <eSUBDIV_TYPE STYPE>
  struct sensor_output_data;

  template<>
  struct sensor_output_data<ISO>
  {
    using type = incr_type;
  };

  template<>
  struct sensor_output_data<ISO_HEX>
  {
    using type = incr_type;
  };

  template<>
  struct sensor_output_data<DIR>
  {
    using type = dir_incr_type;
  };

  //template<>
  //struct sensor_output_data<ISO2>
  //{
  //  using type = Vector_t<E_Float>;
  //};

  } // NUGA

#endif

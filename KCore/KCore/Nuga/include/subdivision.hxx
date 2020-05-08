/*
 
 
 
              NUGA 
 
 
 
 */
//Authors : Sâm Landier (sam.landier@onera.fr)

#ifndef NUGA_SUBDIVISION_HXX
#define NUGA_SUBDIVISION_HXX

#include<vector>
#include "Nuga/include/macros.h"
#include "Nuga/include/subdiv_defs.h"

#include "MeshElement/Tetrahedron.h"
#include "MeshElement/Hexahedron.h"
#include "MeshElement/Prism.h"
#include "MeshElement/Pyramid.h"
#include "MeshElement/Basic.h"
#include "MeshElement/Polyhedron.h"

#include "MeshElement/Triangle.h"
#include "MeshElement/Quadrangle.h"
#include "MeshElement/Polygon.h"


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
  using arr_t = K_FLD::IntArray;
};

//
template <>
struct subdiv_pol<K_MESH::Tetrahedron, ISO>
{
  enum { PGNBC = 4, PHNBC = 8, NBI = 8 };
  using arr_t = K_FLD::IntArray;
};

//
template <>
struct subdiv_pol<K_MESH::Prism, ISO>
{
  enum { PGNBC = 4, PHNBC = 8, NBI = 10 }; // NBI : 4 T3 + 6 Q4
  using arr_t = K_FLD::IntArray;
};

//
template <>
struct subdiv_pol<K_MESH::Pyramid, ISO>
{
  enum { PGNBC=4, PHNBC=10, NBI=13 }; // NBI : 12 T3 + 1 Q4
  using arr_t = K_FLD::IntArray;
};

// Basic - ISO
template <>
struct subdiv_pol<K_MESH::Basic, ISO>
{
  enum { PGNBC = -1, PHNBC = -1/*, NBI = -1 */};
  using arr_t = ngon_unit;
};

// ISO_HEX Poyhedron subdivision => N HEXA children , with N is the nb of nodes
template <>
struct subdiv_pol<K_MESH::Polyhedron<0>, ISO_HEX>
{
  enum { PGNBC = -1, PHNBC = -1/*, NBI = -1 */};
  using arr_t = ngon_unit;

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
  using arr_t = ngon_unit;
};

// isotropic HEXA subdivision => 4 Quadrangles children => fixed stride array
template <>
struct subdiv_pol<K_MESH::Quadrangle, ISO>
{
  enum { NBC = 4};
  using arr_t = K_FLD::IntArray;
};

// isotropic Triangle subdivision => 4 Trianges children => fixed stride array
template <>
struct subdiv_pol<K_MESH::Triangle, ISO>
{
  enum { NBC = 4 };
  using arr_t = K_FLD::IntArray;
};

// ISO_HEX Polygon subdivision => N quad children , with N is the nb of nodes
template <>
struct subdiv_pol<K_MESH::Polygon, ISO_HEX>
{
  enum { PGNBC = -1, PHNBC = -1/*, NBI = -1 */ };
  using arr_t = ngon_unit;

  static E_Int nbc_list(const ngon_unit& PGs, const std::vector<E_Int>& PGlist, std::vector<E_Int>& pregnant)
  {
    //todo JP 
    // loop on polygons given in PGlist and set its stride in pregnant: pregnant[i] = PGs.stride(PGlist[i])
    return 0;
  }

};
  

// Directional subdivision data
  struct dir_incr_type
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
    using type = Vector_t<E_Int>;
  };

  template<>
  struct sensor_output_data<ISO_HEX>
  {
    using type = Vector_t<E_Int>;
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

  };
#endif
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

#ifndef NUGA_SENSOR_HXX
#define NUGA_SENSOR_HXX

#include<vector>
#include "Nuga/include/macros.h"
#include "subdivision.hxx"
#include "smoother.hxx"

namespace NUGA
{

/// Geometric sensor
template <typename mesh_t, typename sensor_input_t>
class sensor
{
  public:
    using output_t = adap_incr_type<mesh_t::SUBTYPE>;

  public:
    sensor(mesh_t& mesh, smoother<mesh_t>* smoother) :_hmesh(mesh), _smoother(smoother) {};
    
    virtual E_Int assign_data(const sensor_input_t& data) { _data = data; return 0; }
    
    bool compute(output_t& adap_incr, bool do_agglo);

    const smoother<mesh_t>* get_smoother() { return _smoother; }

    const mesh_t& get_hmesh() { return _hmesh; }

    virtual bool fill_adap_incr(output_t& adap_incr, bool do_agglo) = 0;
  
    virtual bool update() { return false; };

    virtual bool stop() { return false; }

    virtual ~sensor()
    {
      if (_smoother != nullptr) { delete _smoother; _smoother = nullptr; }
    }

    void append_adap_incr_w_over_connected(output_t& adap_incr);

  public:
    mesh_t &          _hmesh;
    sensor_input_t     _data;
    smoother<mesh_t>* _smoother;
};

//// fix adap incr methods : outisde the class to benefit from overloading ////////////////////////////////////////////////
//
template <typename mesh_t, typename adap_incr_t>
static void discard_disabledand_unhandled(mesh_t& hmesh, adap_incr_t& adap_incr)
{
  E_Int nb_phs = adap_incr.cell_adap_incr.size();// hmesh._ng.PHs.size();
  // prevent to adapt on unhandled elements
  for (E_Int i = 0; i < nb_phs; ++i)
  {
    if (adap_incr.cell_adap_incr[i] == 0) continue;

    E_Int nb_faces = hmesh._ng.PHs.stride(i);
    const E_Int* faces = hmesh._ng.PHs.get_facets_ptr(i);
    bool admissible_elt = K_MESH::Polyhedron<0>::is_basic(hmesh._ng.PGs, faces, nb_faces);

    if (!admissible_elt || !hmesh._PHtree.is_enabled(i)) // fixme : why the latter can happen ? due to smoothing neighbor traversal ?
      adap_incr.cell_adap_incr[i] = 0;
  }
}

///
/// ISO, ISO_HEX impl.
template <typename mesh_t, eSUBDIV_TYPE STYPE>
static void fix_adap_incr(mesh_t& hmesh, adap_incr_type<STYPE>& adap_incr)
{
  discard_disabledand_unhandled(hmesh, adap_incr);
}

/// proto DIR impl.
template <typename mesh_t>
static void fix_adap_incr(mesh_t& hmesh, adap_incr_type<DIR_PROTO>& adap_incr)
{
  // hack to make work the DIR prototype
  for (size_t i = 0; i < adap_incr.cell_adap_incr.size(); ++i)
    adap_incr.cell_adap_incr[i].n[2] = 0; // force to be XY
}

/// DIR impl.
template <typename mesh_t>
static void fix_adap_incr(mesh_t& hmesh, adap_incr_type<DIR>& adap_incr)
{
}

inline NUGA::eDIR get_dir(const K_FLD::FloatArray& crd, const E_Int* nodes, E_Int nnodes)
{
  E_Float Z[] = { 0., 0., 1. };
  E_Float n[3];

  K_MESH::Polygon::normal<K_FLD::FloatArray, 3>(crd, nodes, nnodes, 1, n);

  E_Float ps = ::fabs(NUGA::dot<3>(Z, n));

  if (ps > 0.1) return XY; // face is ortho to Z => XY

  // X or Y ?

  E_Float E[3];
  NUGA::diff<3>(crd.col(nodes[1] - 1), crd.col(nodes[0] - 1), E);
  NUGA::normalize<3>(E);

  ps = ::fabs(NUGA::dot<3>(Z, E));

  if (ps > 0.1) return Y;

  return Xd;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// 
template <typename mesh_t, typename sensor_input_t>
bool sensor<mesh_t, sensor_input_t>::compute(output_t& adap_incr, bool do_agglo)
{
  if (this->stop()) return false; // e.g. itermax is reached for geom_sensor

  //
  /* bool something_filled = */ this->fill_adap_incr(adap_incr, do_agglo);
  //if (!something_filled) return false;

  // expand adap_incr with over-connected cells : more than twice the nb of natural neighbors
  append_adap_incr_w_over_connected(adap_incr);

  //apply the 2:1 rule
  if (_smoother != nullptr) _smoother->smooth(_hmesh, adap_incr);

  // fix inconsistencies
  fix_adap_incr(_hmesh, adap_incr);

  /*std::cout << std::endl;
  std::cout << "adap incr size : " << adap_incr.cell_adap_incr.size() << std::endl;
  std::cout << "adap incr nb 0 : " << std::count_if(ALL(adap_incr.cell_adap_incr), [](int i) {return i == 0; }) << std::endl;
  std::cout << "adap incr nb 1 : " << std::count_if(ALL(adap_incr.cell_adap_incr), [](int i) {return i == 1; }) << std::endl;
  std::cout << std::endl;*/
  //detect if at least one modification is required
  //bool carry_on(false);

  if (adap_incr.cell_adap_incr.size() != 0)
  {
    for (size_t i = 0; i < adap_incr.cell_adap_incr.size(); ++i)
      if (adap_incr.cell_adap_incr[i] != 0) return true;
  }
  if (adap_incr.face_adap_incr.size() != 0)
  {
    for (size_t i = 0; i < adap_incr.face_adap_incr.size(); ++i)
      if (adap_incr.face_adap_incr[i] != 0) return true;
  }
  
  return false;
}

template <typename mesh_t, typename sensor_input_t>
void sensor<mesh_t, sensor_input_t>::append_adap_incr_w_over_connected(output_t& adap_incr)
{
  //
  for (size_t i = 0; i < adap_incr.cell_adap_incr.size(); ++i)
  {
    if (!_hmesh._PHtree.is_enabled(i)) continue;
    if (adap_incr.cell_adap_incr[i] != 0) continue;

    E_Int nbc = _hmesh._PHtree.nb_children(i);
    if (nbc != 0) continue; // deal only with leaf elts : enabling info is not up to date at this stage

    E_Int nbf = _hmesh._ng.PHs.stride(i);
    
    STACK_ARRAY(E_Int, 4 * nbf, neighbours);//fixme 4s
    E_Int nb_neighbours{ 0 };
    _hmesh.get_enabled_neighbours(i, neighbours.get(), nb_neighbours);

    if (nb_neighbours > 2 * nbf)
      adap_incr.cell_adap_incr[i] = 1;

  }
}

}
#endif

/*
 
 
 
              NUGA 
 
 
 
 */
//Authors : Sâm Landier (sam.landier@onera.fr)

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
    using output_t = typename sensor_output_data<mesh_t::SUBTYPE>::type;

  public:
    sensor(const mesh_t& mesh, smoother<mesh_t>* smoother) :_hmesh(mesh), _smoother(smoother) {};
    
    virtual E_Int assign_data(sensor_input_t& data) { _data = data; return 0; }
    
    bool compute(output_t& adap_incr, bool do_agglo);

    const smoother<mesh_t>* get_smoother() { return _smoother; }

    const mesh_t& get_hmesh() { return _hmesh; }

    virtual void fill_adap_incr(output_t& adap_incr, bool do_agglo) = 0;
  
    virtual bool update() { return false; };

    virtual bool stop() { return false; }

    virtual ~sensor()
    {
      if (_smoother != nullptr) delete _smoother; _smoother = nullptr;
    }

  protected:
    mesh_t const &    _hmesh;
    sensor_input_t     _data;
    smoother<mesh_t>* _smoother;
};

//// fix adap incr methods : outisde the class to benefit from overloading ////////////////////////////////////////////////
//
template <typename mesh_t, typename adap_incr_t>
static void discard_disabledand_unhandled(mesh_t& hmesh, adap_incr_t& adap_incr)
{
  E_Int nb_phs = hmesh._ng.PHs.size();
  // prevent to adapt on unhandled elements
  for (E_Int i = 0; i < nb_phs; ++i)
  {
    if (adap_incr[i] == 0) continue;

    E_Int nb_faces = hmesh._ng.PHs.stride(i);
    const E_Int* faces = hmesh._ng.PHs.get_facets_ptr(i);
    bool admissible_elt = K_MESH::Polyhedron<0>::is_basic(hmesh._ng.PGs, faces, nb_faces);

    if (!admissible_elt || !hmesh._PHtree.is_enabled(i)) // fixme : why the latter can happen ? due to smoothing neighbor traversal ?
      adap_incr[i] = 0;
  }
}

//
template <typename mesh_t>
static void fix_adap_incr(mesh_t& hmesh, Vector_t<E_Int>& adap_incr)
{
  discard_disabledand_unhandled(hmesh, adap_incr);
}

//
template <typename mesh_t>
static void fix_adap_incr(mesh_t& hmesh, dir_incr_type& adap_incr)
{
  discard_disabledand_unhandled(hmesh, adap_incr);

  //solve inconsistencies
  //alexis : todo

  // premiere version : desactiver l'adaptation dans les cellules XYZ (iso) qui sont connectées à une cellule "layer" par un QUAD lateral
  adap_incr._ph_dir.clear();
  E_Int nb_phs = hmesh._ng.PHs.size();
  adap_incr._ph_dir.resize(nb_phs, XYZ);

  for (E_Int i = 0; i < nb_phs; ++i)
  {
    // if type is layer => adap_incr._ph_dir[i] = XY
  }

  E_Int nb_pgs = hmesh._ng.PGs.size();
  adap_incr._pg_dir.clear();
  adap_incr._pg_dir.resize(nb_pgs, NONE);

  // boucle sur les layer

  // appel à get_local

  // remplissage approprié de adap_incr._pg_dir pour les 6 faces : X, Y, ou XY

  // boucle sur le reste : 

  // appel à get_local

  // remplissage de adap_incr._pg_dir ou desactivation de adap_incr._ph_dir

  // si NONE => remplissage
  // sinon, si valeur différente => desactivation
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// 
template <typename mesh_t, typename sensor_input_t>
bool sensor<mesh_t, sensor_input_t>::compute(output_t& adap_incr, bool do_agglo)
{
  if (this->stop()) return false; // e.g. itermax is reached for geom_sensor

  //fill in adap incr thanks to points to cell
  this->fill_adap_incr(adap_incr, do_agglo);

  //apply the 2:1 rule
  if (_smoother != nullptr) _smoother->smooth(_hmesh, adap_incr.cell_adap_incr);

  // fix inconsistencies
  fix_adap_incr(_hmesh, adap_incr.cell_adap_incr);

  //detect if at least one modification is required
  //bool carry_on(false);
  E_Int nb_elts = _hmesh._ng.PHs.size();
  for (int i = 0; i < nb_elts; ++i)
    if (adap_incr.cell_adap_incr[i] != 0) return true;
  E_Int nb_faces = _hmesh._ng.PGs.size();
  for (int i = 0; i < nb_faces; ++i)
    if (adap_incr.face_adap_incr[i] != 0) return true;
  return false;
}

}
#endif

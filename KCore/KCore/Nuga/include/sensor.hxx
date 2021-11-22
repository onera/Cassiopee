/*



--------- NUGA v1.0



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
    using output_t = incr_type<mesh_t::SUBTYPE>;

  public:
    sensor(mesh_t& mesh, smoother<mesh_t>* smoother) :_hmesh(mesh), _smoother(smoother) {};
    
    virtual E_Int assign_data(const sensor_input_t& data) { _data = data; return 0; }
    
    bool compute(output_t& adap_incr, bool do_agglo);

    const smoother<mesh_t>* get_smoother() { return _smoother; }

    const mesh_t& get_hmesh() { return _hmesh; }

    virtual bool fill_adap_incr(output_t& adap_incr, bool do_agglo) = 0;

    void append_adap_incr_w_over_connected(output_t& adap_incr);
  
    virtual bool update() { return false; };

    virtual bool stop() { return false; }

    virtual ~sensor()
    {
      if (_smoother != nullptr) { delete _smoother; _smoother = nullptr; }
    }

  protected:
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

//
template <typename mesh_t>
static void fix_adap_incr(mesh_t& hmesh, incr_type<ISO>& adap_incr)
{
  discard_disabledand_unhandled(hmesh, adap_incr);
}

template <typename mesh_t>
static void fix_adap_incr(mesh_t& hmesh, incr_type<ISO_HEX>& adap_incr)
{
  discard_disabledand_unhandled(hmesh, adap_incr);
}

inline NUGA::eDIR get_dir(const K_FLD::FloatArray& crd, const E_Int* nodes, E_Int nnodes)
{
  E_Float Z[] = { 0., 0., 1. };
  E_Float n[3];

  K_MESH::Polygon::normal<K_FLD::FloatArray, 3>(crd, nodes, nnodes, 1, n);

  double ps = ::fabs(NUGA::dot<3>(Z, n));

  if (ps > 0.1) return XY; // face is ortho to Z => XY

  // X or Y ?

  double E[3];
  NUGA::diff<3>(crd.col(nodes[1] - 1), crd.col(nodes[0] - 1), E);

  ps = ::fabs(NUGA::dot<3>(Z, E));

  if (ps > 0.1) return Y;

  return Xd;
}

//
template <typename mesh_t>
static void fix_adap_incr(mesh_t& hmesh, incr_type<DIR>& adap_incr)
{
  discard_disabledand_unhandled(hmesh, adap_incr);

  E_Int nb_phs = adap_incr.cell_adap_incr.size();// hmesh._ng.PHs.size();
  //E_Int nb_pgs = hmesh._ng.PGs.size();

  for (E_Int i = 0; i < nb_phs; ++i)
  {
    if (adap_incr.cell_adap_incr[i] == 0) continue;

    adap_incr.cell_adap_incr[i].n[2] = 0; //hack for CLEF : force to be XY
    
    E_Int nfaces = hmesh._ng.PHs.stride(i);
    const E_Int* faces = hmesh._ng.PHs.get_facets_ptr(i);

    for (size_t f = 0; f < nfaces; ++f)
    {
      E_Int PGi = faces[f] - 1;
      E_Int nnodes = hmesh._ng.PGs.stride(PGi);
      const E_Int* nodes = hmesh._ng.PGs.get_facets_ptr(PGi);

      NUGA::eDIR d = get_dir(hmesh._crd, nodes, nnodes);

      auto& ad = adap_incr.face_adap_incr.vec[PGi];
      //E_Int val = adap_incr.cell_adap_incr[i];

      if (d == XY || d == Xd)
        ad.n[0] = 1;// val;
      if (d == XY || d == Y)
        ad.n[1] = 1;// val;
    }

  }
}

//template <typename mesh_t>
//static void fix_adap_incr(mesh_t& hmesh, incr_type<dir_type>& adap_incr)
//{
//  discard_disabledand_unhandled(hmesh, adap_incr);
//
//  //solve inconsistencies
//  //alexis : todo
//
//  // premiere version : desactiver l'adaptation dans les cellules XYZ (iso) qui sont connectées à une cellule "layer" par un QUAD lateral
//  //adap_incr._ph_dir.clear();
//  E_Int nb_phs = hmesh._ng.PHs.size();
//  //adap_incr._ph_dir.resize(nb_phs, XYZ);
//
//  for (E_Int i = 0; i < nb_phs; ++i)
//  {
//    // if type is layer => adap_incr._ph_dir[i] = XY
//  }
//
//  E_Int nb_pgs = hmesh._ng.PGs.size();
//  //adap_incr._pg_dir.clear();
//  //adap_incr._pg_dir.resize(nb_pgs, NONE);
//
//  // boucle sur les layer
//
//  // appel à get_local
//
//  // remplissage approprié de adap_incr._pg_dir pour les 6 faces : X, Y, ou XY
//
//  // boucle sur le reste : 
//
//  // appel à get_local
//
//  // remplissage de adap_incr._pg_dir ou desactivation de adap_incr._ph_dir
//
//  // si NONE => remplissage
//  // sinon, si valeur différente => desactivation
//}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/// 
template <typename mesh_t, typename sensor_input_t>
bool sensor<mesh_t, sensor_input_t>::compute(output_t& adap_incr, bool do_agglo)
{
  if (this->stop()) return false; // e.g. itermax is reached for geom_sensor

  //
  bool something_filled = this->fill_adap_incr(adap_incr, do_agglo);
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
    for (int i = 0; i < adap_incr.cell_adap_incr.size(); ++i)
      if (adap_incr.cell_adap_incr[i] != 0) return true;
  }
  if (adap_incr.face_adap_incr.size() != 0)
  {
    for (int i = 0; i < adap_incr.face_adap_incr.size(); ++i)
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

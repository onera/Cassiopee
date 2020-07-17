/*
 
 
 
              NUGA 
 
 
 
 */
//Authors : Sâm Landier (sam.landier@onera.fr)

#ifndef NUGA_SHELL_SMOOTHER_HXX
#define NUGA_SHELL_SMOOTHER_HXX

#include "smoother.hxx"

namespace NUGA
{

template <typename mesh_t>
struct shell_smoother : public smoother<mesh_t>
{
  using output_t = typename sensor_output_data<mesh_t::SUBTYPE>::type;
  using cell_adap_incr_t = typename output_t::cell_output_type;

  Vector_t<E_Int> _Ln;

  shell_smoother() = default;

  void smooth(const mesh_t& hmesh, cell_adap_incr_t& adap_incr);

  void __update_nodal_levels(const mesh_t& hmesh, E_Int PHi);

  ~shell_smoother() = default;
};

template <typename mesh_t>
void shell_smoother<mesh_t>::smooth(const mesh_t& hmesh, cell_adap_incr_t& adap_incr)
{
  E_Int n_nodes = hmesh._crd.cols();
  _Ln.resize(n_nodes);

  for (size_t i = 0; i<adap_incr.size(); i++)
  {
    if (adap_incr[i] == 1) __update_nodal_levels(hmesh, i);
  }

  // NODAL SMOOTHING LOOP
  E_Int Lc;
  E_Int n(0);
  bool carry_on = true;
  while (carry_on)
  {
    n = n + 1;
    carry_on = false;
    for (int i = 0; i<hmesh._ng.PHs.size(); i++) {
      E_Int Mnodes(0);
      if (adap_incr[i]) continue;
      if (hmesh._PHtree.is_enabled(i)) {
        Lc = hmesh._PHtree.get_level(i);
        const E_Int* pPHi = hmesh._ng.PHs.get_facets_ptr(i);

        // Mnodes definition
        for (int j = 0; j<hmesh._ng.PHs.stride(i); j++) {
          E_Int PGi = *(pPHi + j) - 1;
          E_Int PHn = NEIGHBOR(i, hmesh._F2E, PGi);
          if (!(PHn == E_IDX_NONE) && !(hmesh._PHtree.is_enabled(PHn)) && !(hmesh._PHtree.get_level(PHn) == 0) && !(hmesh._PHtree.is_enabled(hmesh._PHtree.parent(PHn)))) {
            E_Int nb_ch = hmesh._PGtree.nb_children(PGi);
            for (int k = 0; k< nb_ch; k++) {
              const E_Int* enf = hmesh._PGtree.children(PGi);
              E_Int stride = hmesh._ng.PGs.stride(*(enf + k));
              const E_Int* ind = hmesh._ng.PGs.get_facets_ptr(*(enf + k));
              for (int l = 0; l<stride; l++) {
                Mnodes = std::max(_Ln[*(ind + l) - 1], Mnodes);
              }
            }
          }
          else {
            E_Int stride = hmesh._ng.PGs.stride(PGi);
            const E_Int *ind = hmesh._ng.PGs.get_facets_ptr(PGi);
            for (int l = 0; l< stride; l++) {
              Mnodes = std::max(_Ln[*(ind + l) - 1], Mnodes);
            }
          }
        }
        // fin Mnodes definition
        if (Lc < Mnodes - 1)
        {
          const E_Int* faces = hmesh._ng.PHs.get_facets_ptr(i);
          E_Int nb_faces = hmesh._ng.PHs.stride(i);
          bool admissible_elt = K_MESH::Polyhedron<0>::is_basic(hmesh._ng.PGs, faces, nb_faces);

          if (!admissible_elt)
            continue;

          adap_incr[i] = 1;
          __update_nodal_levels(hmesh, i);
          carry_on = true;
        }
      }
    }
  }
}

template <typename mesh_t>
void shell_smoother<mesh_t>::__update_nodal_levels(const mesh_t& hmesh, E_Int PHi)
{
  E_Int Lc = hmesh._PHtree.get_level(PHi);
  E_Int nb_faces = hmesh._ng.PHs.stride(PHi);
  for (int j = 0; j< nb_faces; j++) {
    const E_Int* PGj = hmesh._ng.PHs.get_facets_ptr(PHi);
    E_Int PG = PGj[j] - 1;
    const E_Int* pN = hmesh._ng.PGs.get_facets_ptr(PG);
    E_Int nb_nodes = hmesh._ng.PGs.stride(PG);
    for (int l = 0; l< nb_nodes; l++) {
      _Ln[*(pN + l) - 1] = std::max(_Ln[*(pN + l) - 1], Lc + 1);
    }
  }
}

}
#endif


#include "Nuga/include/adapt_cells.h"
#include "Nuga/include/Basic.h"
#include "Nuga/include/hierarchical_mesh.hxx"
#include "Nuga/include/geom_sensor.hxx"

using namespace NUGA;


int adapt_cells(c_phmesh_t& m, const c_crd3D_t& src_pts)
{
  using ELT_t = K_MESH::Basic;
  using hmesh_t = hierarchical_mesh<ELT_t, ISO, ngon_type>;
  using sensor_t = NUGA::geom_sensor<hmesh_t, c_phmesh_t::crd_t>;

  K_FLD::FloatArray crd;// (m.crd);      //TODO
  ngon_type ng;// (m.pgs, m.phs);    // TODO
  
  hmesh_t hmesh(crd, ng);

  sensor_t gsensor(hmesh, eSmoother::V1_NEIGH, 1/*max_pts_per_cell*/, 10/*itermax*/);

  hmesh.init();

  int err = gsensor.assign_data(src_pts);

  adaptor<hmesh_t, sensor_t>::run(hmesh, gsensor);

  std::vector<E_Int> oids;
  ngon_type ngo;
  
  hmesh.conformize(ngo, oids);

  //m.crd.release();  //TODO
  int dim{ 3 };
  bool calloc;
  hmesh._crd.relay_mem(m.crd.p, dim, m.crd.n, calloc);
  assert(calloc == m.crd.CALLOC);
  ngo.PHs.relay_mem(m.phs);
  ngo.PGs.relay_mem(m.pgs);

  return 0;
}
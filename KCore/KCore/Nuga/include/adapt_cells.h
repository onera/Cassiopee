/*



--------- NUGA v1.0



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

#ifndef NUGA_ADAPT_CELLS_H
#define NUGA_ADAPT_CELLS_H

#include "Nuga/include/cdatastruct.hxx"
//#include <vector>

namespace NUGA
{
  // adapts m with respect to input source points
  // Returns a polyhedral conformal mesh for polygons and polyhedra.
  // Just the leaves : all intermediary levels are removed
  int adapt_cells(c_phmesh_t& m,             // i/o mesh
                  const c_crd3D_t& src_pts);  // source points

  //// adapts m with respect to input source points
  //// Returns a polyhedral conformal mesh with ascendent history for polygons and polyhedra.
  //// Just the leaves : all intermediary levels are removed
  //int adapt_cells(c_phmesh_t& m,             // i/o mesh
  //  const c_crd3D_t& src_pts,                // source points
  //  std::vector<int>& pgoids,                // polygon original id in m before adaptation
  //  std::vector<int>& phoids);               // polyhedron original id in m before adaptation

  // adapts m with respect to input source points
  // PERSISTANT : returns the entire genealogy
  //int adapt_cells(const phmesh_t& m, const crd3D_t& src_pts, hmesh& hm);
}

#endif
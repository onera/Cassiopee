/*



NUGA



*/
//Authors : Sâm Landier (sam.landier@onera.fr)

#ifndef NUGA_ADAPT_CELLS_H
#define NUGA_ADAPT_CELLS_H

#include "Nuga/include/cdatastruct.hxx"

// adapts m with respect to input source points
// ONE SHOT : returns a polyhedral conformal mesh with ascendent history. Just the leaves : all intermediary levels removed
int adapt_cells(c_phmesh_t& m, const c_crd3D_t& src_pts);


// adapts m with respect to input source points
// PERSISTANT : returns the entire genealogy
//int adapt_cells(const phmesh_t& m, const crd3D_t& src_pts, hmesh& hm);

#endif
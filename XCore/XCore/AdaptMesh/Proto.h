#ifndef _PROTO_H
#define _PROTO_H

#include "Struct.h"
#include "../common/mem.h"

// Metric
void compute_ref_data(AMesh *M, E_Float **fields, E_Int nfields);
void Right_shift(E_Int *pn, E_Int pos, E_Int size);
E_Int Get_pos(E_Int e, E_Int *pn, E_Int size);

// Comm
void exchange_proc_data_d(AMesh *M, E_Float *data, E_Float ***rbuf);

// Hessian
E_Float *compute_hessian(AMesh *M, E_Float *field);

// Mesh
E_Int *get_face(E_Int face, E_Int &np, E_Int *ngon, E_Int *indPG);
E_Int *get_cell(E_Int cell, E_Int &nf, E_Int *nface, E_Int *indPH);
void make_cell_centers(AMesh *M);
E_Int get_stride(E_Int i, E_Int *indir);
E_Int *get_facets(E_Int i, E_Int *cn, E_Int *indir);
E_Int set_faces_type(AMesh *M);
E_Int set_cells_type(AMesh *M);
E_Int get_reorient(E_Int face, E_Int cell, E_Int normalIn, AMesh *M);
void reorder_cells(AMesh *M);
void compute_face_center_and_area(E_Int id, E_Int stride,
  E_Int *pn, E_Float *x, E_Float *y, E_Float *z, E_Float *fc, E_Float *fa);
void Order_tri(E_Int *local, E_Int *pn, E_Int reorient, E_Int i0);
void Order_quad(E_Int *local, E_Int *pn, E_Int reorient, E_Int i0);
void check_canon_tetra(E_Int cell, AMesh *M);
void check_canon_hexa(E_Int cell, AMesh *M);
void check_canon_cells(AMesh *M);

// Refine
void get_ref_cells_and_faces(AMesh *, std::vector<E_Int> &ref_cells,
  std::vector<E_Int> &ref_faces);
void resize_data_for_refinement(AMesh *M, E_Int nref_cells, E_Int nref_faces);
void refine_faces(const std::vector<E_Int> &ref_faces, AMesh *M);
void refine_cells(const std::vector<E_Int> &ref_cells, AMesh *M);


#endif

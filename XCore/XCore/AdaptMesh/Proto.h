#ifndef _PROTO_H
#define _PROTO_H

#include "Struct.h"

// Metric
void compute_ref_data(AMesh *M, E_Float **fields, E_Int nfields);

// Comm
void exchange_proc_data_d(AMesh *M, E_Float *data, E_Float ***rbuf);

// Hessian
E_Float *compute_hessian(AMesh *M, E_Float *field);

// Mesh
E_Int *get_face(E_Int face, E_Int &np, E_Int *ngon, E_Int *indPG);
E_Int *get_cell(E_Int cell, E_Int &nf, E_Int *nface, E_Int *indPH);
void make_cell_centers(AMesh *M);
inline E_Int get_stride(E_Int i, E_Int *indir);
inline E_Int *get_ptr(E_Int i, E_Int *cn, E_Int *indir);
E_Int set_faces_type(AMesh *M);
E_Int set_cells_type(AMesh *M);
E_Int get_reorient(E_Int face, E_Int cell, E_Int normalIn, AMesh *M);
void reorder_cells(AMesh *M);

#endif

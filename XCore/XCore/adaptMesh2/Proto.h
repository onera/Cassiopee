#ifndef _PROTO_H
#define _PROTO_H

#include "Struct.h"
#include "../common/mem.h"

// Metric
void compute_ref_data(AMesh *M, E_Float **fields, E_Int nfields);
void reconstruct_parent_quad(E_Int face, AMesh *M, E_Int pn[4]);
void smooth_ref_data_parallel(AMesh *M);
void smooth_ref_data(AMesh *M);

// Comm
void Comm_interface_data_i(AMesh *M, E_Int *data, E_Int stride, E_Int **rbuf);
void Comm_interface_data_f(AMesh *M, E_Float *data, E_Int stride, E_Float **rbuf);
AMesh *Redistribute_mesh(AMesh *M);
void Comm_waitall(AMesh *M);

// Hessian
E_Float *compute_hessian(AMesh *M, E_Float *field);

// Mesh
E_Int master_cell(E_Int cell, AMesh *M);
E_Int master_face(E_Int face, AMesh *M);
E_Int get_neighbour(E_Int cell, E_Int face, AMesh *M);
E_Int is_internal_face(E_Int face, AMesh *M);
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
E_Int check_canon_tetra(E_Int cell, AMesh *M);
E_Int check_canon_hexa(E_Int cell, AMesh *M);
E_Int check_canon_penta(E_Int cell, AMesh *M);
E_Int check_canon_pyra(E_Int cell, AMesh *M);
E_Int check_canon_cells(AMesh *M);
AMesh *init_mesh(K_FLD::FldArrayI& cn, E_Float *px, E_Float *py, E_Float *pz,
  E_Int npts);
void get_full_cell(E_Int cell, AMesh *M, E_Int &nf, E_Int *pf);
void update_boundary_faces(AMesh *M);
void ngon_print(AMesh *M);
void nface_print(AMesh *M);
const char *cell_type_to_string(E_Int type);
const char *face_type_to_string(E_Int type);
void renumber_mesh(AMesh *M, const std::vector<E_Int> &new_cells,
  const std::vector<E_Int> &new_faces, E_Int nc, E_Int nf,
  E_Int sizeNFace, E_Int sizeNGon);
void compress_mesh(AMesh *M, const std::vector<E_Int> &new_cells,
  const std::vector<E_Int> &new_faces, E_Int nc, E_Int nf,
  E_Int sizeNFace, E_Int sizeNGon);
void mesh_drop(AMesh *M);
void patch_drop(Patch *P);

// Refine
void get_ref_faces_and_cells(AMesh *, std::vector<E_Int> &ref_faces,
  std::vector<E_Int> &ref_cells);
void resize_data_for_refinement(AMesh *M, size_t nref_cells,
  size_t nref_faces);
void refine_faces(const std::vector<E_Int> &ref_faces,
  const std::vector<E_Int> &ref_patterns, size_t start, size_t stop, AMesh *M);
void refine_cells(const std::vector<E_Int> &ref_cells, size_t start,
  size_t stop, AMesh *M);
void update_external_pe_after_ref(E_Int cell, AMesh *M);
void update_external_pe_after_unref(E_Int cell, AMesh *M);
void refine_mesh(AMesh *, std::vector<E_Int> &ref_faces,
  std::vector<E_Int> &ref_cells);
void update_global_cells_after_ref(AMesh *M);
void resize_data_for_synchronisation(AMesh *M);
void update_patch_faces_after_ref(AMesh *M);
void update_global_faces_after_ref(AMesh *M);
void update_global_points_after_ref(AMesh *M);
void synchronise_patches_after_ref(AMesh *M);

// Unrefine
void get_unref_faces_and_cells(AMesh *, std::vector<E_Int> &unref_faces,
  std::vector<E_Int> &unref_cells);
void unrefine_faces(const std::vector<E_Int> &unref_faces,
  const std::vector<E_Int> &patterns, size_t start, size_t stop, AMesh *M);
void unrefine_cells(const std::vector<E_Int> &unref_cells, size_t start,
  size_t stop, AMesh *M);
void assign_pattern_to_adapt_faces(const std::vector<E_Int> &adapt_faces,
  std::vector<E_Int> &patterns, AMesh *M);
void unrefine_tetra(E_Int cell, AMesh *M);
void unrefine_pyra(E_Int cell, AMesh *M);
void unrefine_mesh(AMesh *M, std::vector<E_Int> &unref_faces,
  std::vector<E_Int> &unref_cells);

// RenumberMesh
//std::vector<E_Int> renumber_cells(AMesh *M);
//std::vector<E_Int> renumber_faces(AMesh *M,
//  const std::vector<E_Int> &new_cells);
//void renumber_mesh(AMesh *M);
//void renumber_indPH(AMesh *M, const std::vector<E_Int> &new_cells,
//  E_Int *new_indPH);
//void renumber_nface(AMesh *M, const std::vector<E_Int> &new_cells,
//  E_Int *new_indPH, E_Int *new_nface);
void renumber_indPG(AMesh *M, const std::vector<E_Int> &new_faces,
  E_Int *new_indPG);
void renumber_ngon(AMesh *M, const std::vector<E_Int> &new_faces,
  E_Int *new_indPG, E_Int *new_ngon);

void renumber_nface_shallow(AMesh *M, const std::vector<E_Int> &new_faces);
void renumber_boundary(AMesh *M, const std::vector<E_Int> &new_faces);
//void renumber_comm_patches(AMesh *M, const std::vector<E_Int> &new_faces);
void init_mesh_numbering(AMesh *M);

// Topo
E_Int Orient_boundary(AMesh *M);
E_Int Build_own_nei(AMesh *M);
void close_mesh(AMesh *M);
void Right_shift(E_Int *pn, E_Int pos, E_Int size);
E_Int Get_pos(E_Int e, E_Int *pn, E_Int size);
E_Int check_closed_cell(E_Int cell, E_Int *pf, E_Int nf, AMesh *M);
void set_cells_for_2D(AMesh *M);
E_Int check_face_aligned_with_vector(E_Int face, E_Float vec[3], AMesh *M);

// Hexa
void H18_refine(E_Int cell, AMesh *M);
void H18_unrefine(E_Int cell, AMesh *M);
void H27_refine(E_Int cell, AMesh *M);
void H27_unrefine(E_Int cell, AMesh *M);

// Penta
void Pe12_refine(E_Int cell, AMesh *M);
void Pe12_unrefine(E_Int cell, AMesh *M);
void Pe18_refine(E_Int cell, AMesh *M);
void Pe18_unrefine(E_Int cell, AMesh *M);


// Quad
void Q6_get_ordered_data(AMesh *M, E_Int NODE, E_Int reorient,
  E_Int *children, E_Int *local);
void Q6_refine(E_Int face, AMesh *M);
void Q6_unrefine(E_Int face, AMesh *M);

void Q9_get_ordered_data(AMesh *M, E_Int NODE, E_Int reorient,
  E_Int *children, E_Int *local);
void Q9_refine(E_Int face, AMesh *M);
void Q9_unrefine(E_Int face, AMesh *M);

// Tri
void T6_get_ordered_data(AMesh *M, E_Int NODE, E_Int reorient,
  E_Int *children, E_Int *local);
void T6_refine(E_Int face, AMesh *);
void T6_unrefine(E_Int face, AMesh *);

#endif

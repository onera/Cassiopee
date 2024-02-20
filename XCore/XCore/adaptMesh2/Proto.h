#ifndef _PROTO_H
#define _PROTO_H

#include "Struct.h"
#include "../common/mem.h"

// Mesh
AMesh *init_mesh(K_FLD::FldArrayI& cn, E_Float *px, E_Float *py, E_Float *pz,
  E_Int npts);
E_Int master_cell(E_Int cell, AMesh *M);
E_Int master_face(E_Int face, AMesh *M);
E_Int get_neighbour(E_Int cell, E_Int face, AMesh *M);
E_Int is_internal_face(E_Int face, AMesh *M);
E_Int *get_face(E_Int face, E_Int &np, E_Int *ngon, E_Int *indPG);
E_Int *get_cell(E_Int cell, E_Int &nf, E_Int *nface, E_Int *indPH);
E_Int get_stride(E_Int i, E_Int *indir);
E_Int *get_facets(E_Int i, E_Int *cn, E_Int *indir);
E_Int set_faces_type(AMesh *M);
E_Int set_cells_type(AMesh *M);
E_Int get_reorient(E_Int face, E_Int cell, E_Int normalIn, AMesh *M);
void reorder_cells(AMesh *M);
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
void make_dual_graph(AMesh *M);
void compute_ref_cells_centers(AMesh *M, const std::vector<E_Int> &ref_cells);

// Refine
void get_ref_faces_and_cells(AMesh *, std::vector<E_Int> &ref_faces,
  std::vector<E_Int> &ref_cells);
void resize_data_for_refinement(AMesh *M, size_t nref_cells,
  size_t nref_faces);
void update_external_pe_after_ref(E_Int cell, AMesh *M);
void refine_mesh(AMesh *, const std::vector<E_Int> &ref_faces,
  const std::vector<E_Int> &ref_cells);
void update_patch_faces_after_ref(AMesh *M);
void update_global_cells_after_ref(AMesh *M);
void update_global_faces_after_ref(AMesh *M);

// Topo
E_Int Orient_boundary(AMesh *M);
E_Int Build_own_nei(AMesh *M);
void close_mesh(AMesh *M);
void Right_shift(E_Int *pn, E_Int pos, E_Int size);
E_Int Get_pos(E_Int e, E_Int *pn, E_Int size);
E_Int check_closed_cell(E_Int cell, E_Int *pf, E_Int nf, AMesh *M);
void set_cells_for_2D(AMesh *M);
E_Int check_face_aligned_with_vector(E_Int face, E_Float vec[3], AMesh *M);
E_Int check_canon_cells(AMesh *M);

// Metric
void compute_ref_data(AMesh *M, E_Float **fields, E_Int nfields);
void reconstruct_parent_quad(E_Int face, AMesh *M, E_Int pn[4]);
void smooth_ref_data(AMesh *M);

// Quad
void Q6_get_ordered_data(AMesh *M, E_Int NODE, E_Int reorient,
  E_Int *children, E_Int *local);
void Q6_refine(E_Int quad, AMesh *M);
void Q6_unrefine(E_Int quad, AMesh *M);
void Q9_get_ordered_data(AMesh *M, E_Int NODE, E_Int reorient,
  E_Int *children, E_Int *local);
void Q9_refine(E_Int quad, AMesh *M);
void Q9_unrefine(E_Int quad, AMesh *M);
void reorder_quad(E_Int *local, E_Int *pn, E_Int reorient, E_Int i0);

// Tri
void T6_get_ordered_data(AMesh *M, E_Int NODE, E_Int reorient,
  E_Int *children, E_Int *local);
void T6_refine(E_Int tri, AMesh *);
void T6_unrefine(E_Int tri, AMesh *);
void reorder_tri(E_Int *local, E_Int *pn, E_Int reorient, E_Int i0);

// Hexa
void H18_refine(E_Int hexa, AMesh *M);
void H18_unrefine(E_Int hexa, AMesh *M);
void H27_refine(E_Int hexa, AMesh *M);
void H27_unrefine(E_Int hexa, AMesh *M);
void reorder_hexa(E_Int hexa, E_Int nf, E_Int *pf, AMesh *M);
void make_ref_data_hexa(E_Int hexa, AMesh *M, E_Float *pM, const pDirs &Dirs);
void make_pdirs_hexa(E_Int hexa, AMesh *M, pDirs &Dirs);
E_Int check_canon_hexa(E_Int hexa, AMesh *M);

// Penta
void Pe12_refine(E_Int penta, AMesh *M);
void Pe12_unrefine(E_Int penta, AMesh *M);
void Pe18_refine(E_Int penta, AMesh *M);
void Pe18_unrefine(E_Int penta, AMesh *M);
void reorder_penta(E_Int penta, E_Int nf, E_Int *pf, AMesh *M);
void make_ref_data_penta(E_Int penta, AMesh *M, E_Float *pM, const pDirs &Dirs);
void make_pdirs_penta(E_Int penta, AMesh *M, pDirs &Dirs);
E_Int check_canon_penta(E_Int penta, AMesh *M);

// Tetra
void reorder_tetra(E_Int tetra, E_Int nf, E_Int *pf, AMesh *M);
void unrefine_tetra(E_Int tetra, AMesh *M);
void refine_tetra(E_Int tetra, AMesh *M);
void make_ref_data_tetra(E_Int tetra, AMesh *M, E_Float *pM, const pDirs &Dirs);
void make_pdirs_tetra(E_Int tetra, AMesh *M, pDirs &Dirs);
E_Int check_canon_tetra(E_Int tetra, AMesh *M);

// Pyra
void reorder_pyra(E_Int pyra, E_Int nf, E_Int *pf, AMesh *M);
void refine_pyra(E_Int pyra, AMesh *M);
void unrefine_pyra(E_Int pyra, AMesh *M);
void make_ref_data_pyra(E_Int pyra, AMesh *M, E_Float *pM, const pDirs &Dirs);
void make_pdirs_pyra(E_Int pyra, AMesh *M, pDirs &Dirs);
E_Int check_canon_pyra(E_Int pyra, AMesh *M);

// Unrefine
void get_unref_faces_and_cells(AMesh *, std::vector<E_Int> &unref_faces,
  std::vector<E_Int> &unref_cells);
void unrefine_faces(const std::vector<E_Int> &unref_faces,
  const std::vector<E_Int> &patterns, size_t start, size_t stop, AMesh *M);
void unrefine_cells(const std::vector<E_Int> &unref_cells, size_t start,
  size_t stop, AMesh *M);
void assign_pattern_to_adapt_faces(const std::vector<E_Int> &adapt_faces,
  std::vector<E_Int> &patterns, AMesh *M);
void unrefine_mesh(AMesh *M, std::vector<E_Int> &unref_faces,
  std::vector<E_Int> &unref_cells);
void update_external_pe_after_unref(E_Int cell, AMesh *M);

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

// Comm
void Comm_interface_data_i(AMesh *M, E_Int *data, E_Int stride, E_Int **rbuf);
void Comm_interface_data_f(AMesh *M, E_Float *data, E_Int stride,
  E_Float **rbuf);
void Comm_waitall(AMesh *M);
AMesh *load_balance_mesh(AMesh *M);

// Gradient
E_Float gradient_norm_inf(E_Float G[3]);

// Hessian
E_Float hessian_norm_inf(E_Float H[6]);

// Geom
void compute_principal_vecs(AMesh *M, std::vector<pDirs> &Dirs);

#endif

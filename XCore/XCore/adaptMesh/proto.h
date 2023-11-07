/*    
    Copyright 2013-2023 Onera.

    This file is part of Cassiopee.

    Cassiopee is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Cassiopee is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Cassiopee.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef PROTO_H
#define PROTO_H

#include "struct.h"

/* distribute */
mesh *redistribute_mesh(mesh *);

/* comm */
void comm_waitall(mesh *);
void comm_interface_data_d(mesh *, E_Float *, E_Int, E_Float **);
void comm_interface_data_i(mesh *, E_Int *, E_Int, E_Int **);

/* mesh */
void compute_face_centers(mesh *);
void compute_cell_centers(mesh *);
void compute_face_center(mesh *, E_Int);
void compute_cell_center(mesh *, E_Int);
E_Int get_neighbour(E_Int, E_Int, mesh *);
void mesh_free(mesh *);
E_Float mesh_memsize(mesh *);
void mesh_save_memory(mesh *);
mesh_leaves mesh_get_leaves(mesh *, tree *, tree *, E_Int);

/* metric */
void hessian_to_metric(E_Float *, E_Int);
//void compute_ref_data(mesh *, E_Float *);
void smooth_ref_data(mesh *);
void smooth_ref_data_seq(mesh *);
std::vector<E_Int> compute_canon_info(E_Int, mesh *, E_Int *);
void deduce_nei_ref_data(E_Int, E_Int, E_Int, E_Int, E_Int *);
E_Int is_metric_valid(E_Float *);
E_Int make_ref_data(mesh *, E_Float **, E_Int, E_Float *);
//E_Int apply_freeze_vector(mesh *M, E_Float *);
void compute_principal_vecs(mesh *, std::vector<pDirs> &);
void compute_cells_principal_vecs(mesh *, const std::vector<E_Int> &, std::vector<pDirs> &);

/* topo */
void reorder_hexa(mesh *);
void reorient_skin(mesh *);
void build_own_nei(mesh *);
void topo_init_mesh(mesh *);
E_Int get_reorient(E_Int, E_Int, E_Int, mesh *);
E_Int get_pos(E_Int, E_Int *, E_Int);
void order_quad(E_Int *, E_Int *, E_Int, E_Int);
void right_shift(E_Int *, E_Int, E_Int);

/* math */
void compute_lsq_grad_matrices(mesh *);
E_Float *compute_grad(mesh *, E_Float *);
void compute_lsq_hess_matrices(mesh *);
E_Float *compute_hessian(mesh *, E_Float *);
E_Int BiCGStab(E_Float *, E_Float *, E_Float *, E_Int);
void eigen(E_Float *, E_Float *, E_Float *, E_Float *, E_Float *);
E_Int feq(E_Float, E_Float);
E_Float dot(const E_Float *, const E_Float *, E_Int);
E_Float norm(const E_Float *, E_Int);
void cross(E_Float *, E_Float *, E_Float *);
void symmat_dot_vec(const E_Float *, const E_Float *, E_Float *c);
E_Int compute_hessians_ngon(K_FLD::FldArrayI &, E_Float *, E_Float *,
  E_Float *, E_Int *, E_Int *, E_Float *, const std::vector<E_Float *> &,
  const std::vector<E_Float *> &, std::vector<E_Float *> &);
E_Int compute_gradients_ngon(K_FLD::FldArrayI &, E_Float *, E_Float *,
  E_Float *, E_Int *, E_Int *, E_Float *,
  const std::vector<E_Float *> &, std::vector<E_Float *> &);
void compute_cell_centers_and_vols(K_FLD::FldArrayI &, E_Float *,
  E_Float *, E_Float *, E_Int *, E_Int *, E_Float *, E_Float *,
  E_Float *, E_Float *);
void compute_face_centers_and_areas(K_FLD::FldArrayI &, E_Float *,
  E_Float *, E_Float *, E_Float *, E_Float *);

/* cut */
void cut_face_x(E_Int, mesh *, tree *, E_Int, E_Int);
void cut_face_y(E_Int, mesh *, tree *, E_Int, E_Int);
void cut_face_xy(E_Int, mesh *, tree *);
void cut_cell_x(E_Int, mesh *, tree *, tree *);
void cut_cell_y(E_Int, mesh *, tree *, tree *);
void cut_cell_z(E_Int, mesh *, tree *, tree *);
void cut_cell_xy(E_Int, mesh *, tree *, tree *);
void cut_cell_xz(E_Int, mesh *, tree *, tree *);
void cut_cell_yz(E_Int, mesh *, tree *, tree *);
void cut_cell_xyz(E_Int, mesh *, tree *, tree *);
void resize_data_for_refinement(mesh *, tree *, tree *, E_Int);
std::vector<E_Int> get_ref_cells(mesh *);
E_Int is_cell_to_refine(E_Int *);

/* tree */
tree *tree_new(E_Int, E_Int);
void tree_insert_children(tree *, E_Int, E_Int, E_Int);
void tree_XFREE(tree *);
void tree_resize(tree *, E_Int, E_Int);
E_Int tree_get_nchildren(tree *, E_Int);
E_Int *tree_get_children(tree *, E_Int);
void tree_print(tree *);
void tree_save_memory(tree *, E_Int);
E_Float tree_memsize(tree *);
void tree_get_face_leaves(tree *, E_Int, std::vector<E_Int> &);

#endif

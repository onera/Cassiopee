#ifndef PROTO_H
#define PROTO_H

#ifndef M_PI
#define M_PI       3.14159265358979323846
#endif

#include "struct.h"

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

/* metric */
void hessian_to_metric(E_Float *, mesh *);
std::vector<E_Int> compute_ref_data(mesh *, E_Float *);
void smooth_ref_data(mesh *, std::vector<E_Int> &);
void compute_canon_info(E_Int, mesh *, E_Int *);

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
E_Float dot(E_Float *, E_Float *, E_Int);
void cross(E_Float *, E_Float *, E_Float *);
void symmat_dot_vec(E_Float *, E_Float *, E_Float *c);

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
void resize_data_for_refinement(mesh *, tree *, tree *, E_Int, E_Int);
std::vector<E_Int> get_ref_cells(mesh *, std::vector<E_Int> &, E_Int *, E_Int *);
E_Int is_cell_to_refine(E_Int *);

/* tree */
tree *tree_new(E_Int, E_Int);
void tree_insert_children(tree *, E_Int, E_Int, E_Int);
void tree_free(tree *);
void tree_resize(tree *, E_Int, E_Int);
E_Int tree_get_nchildren(tree *, E_Int);
E_Int *tree_get_children(tree *, E_Int);
void tree_print(tree *);

#endif

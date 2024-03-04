#include "proto.h"
#include <map>
#include <unordered_map>
#include <stack>

static
E_Int compute_bin_point(E_Float x, E_Float y, E_Int ndiv, const BBox &bbox)
{
  E_Int I = x / bbox.xmax * 0.99 * ndiv;
  E_Int J = y / bbox.ymax * 0.99 * ndiv;
  E_Int bin = I + ndiv * J;
  return bin;
}

PyObject *K_XCORE::intersectMesh(PyObject *self, PyObject *args)
{
  PyObject *MASTER, *SLAVE;
  
  if (!PYPARSETUPLE_(args, OO_, &MASTER, &SLAVE)) {
    RAISE("Bad input.");
    return NULL;
  }

  // Check master mesh
  E_Int ret;
  E_Int ni, nj, nk;
  K_FLD::FldArrayF *fm;
  K_FLD::FldArrayI *cnm;
  char *varString, *eltType;
  ret = K_ARRAY::getFromArray3(MASTER, varString, fm, ni, nj, nk, cnm, eltType);

  if (ret <= 0) {
    RAISE("Bad master mesh");
    return NULL;
  }

  if (ret == 1) {
    RAISE("Master mesh is not NGon.");
    RELEASESHAREDS(MASTER, fm);
    return NULL;
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);

  if (posx == -1 || posy == -1 || posz == -1) {
    RELEASESHAREDU(MASTER, fm, cnm);
    RAISE("Bad master point coordinates.");
    return NULL;
  }

  posx++; posy++; posz++;

  E_Float *Xm = fm->begin(posx);
  E_Float *Ym = fm->begin(posy);
  E_Float *Zm = fm->begin(posz);
  E_Int npm = fm->getSize();

  // Check slave mesh
  K_FLD::FldArrayF *fs;
  K_FLD::FldArrayI *cns;
  ret = K_ARRAY::getFromArray3(SLAVE, varString, fs, ni, nj, nk, cns, eltType);

  if (ret <= 0) {
    RAISE("Bad slave mesh");
    return NULL;
  }

  if (ret == 1) {
    RAISE("Slave mesh is not NGon.");
    RELEASESHAREDS(SLAVE, fs);
    RELEASESHAREDU(MASTER, fm, cnm);
    return NULL;
  }

  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);

  if (posx == -1 || posy == -1 || posz == -1) {
    RELEASESHAREDU(MASTER, fm, cnm);
    RELEASESHAREDU(SLAVE, fs, cns);
    RAISE("Bad slave point coordinates.");
    return NULL;
  }

  posx++; posy++; posz++;

  E_Float *Xs = fs->begin(posx);
  E_Float *Ys = fs->begin(posy);
  E_Float *Zs = fs->begin(posz);
  E_Int nps = fs->getSize();

  // Init and orient meshes
  Mesh *M = mesh_init(*cnm, Xm, Ym, Zm, npm);
  Mesh *S = mesh_init(*cns, Xs, Ys, Zs, nps);

  printf("Master cells: %d\n", M->nc);
  printf("Master points: %d\n", M->np);
  printf("Slave cells: %d\n", S->nc);
  printf("Slave points: %d\n", S->np);

  mesh_write(M, "mesh.dat");
  mesh_write(S, "cylinder.dat");

  // Make the bounding boxes
  BBox bboxm = bbox_make(M);
  BBox bboxs = bbox_make(S);

  // Which slave points are inside master bbox?
  std::set<E_Int> bbox_points;

  for (E_Int i = 0; i < S->np; i++) {
    if (bbox_is_point_inside(bboxm, S->x[i], S->y[i], S->z[i])) {
      bbox_points.insert(i);
    }
  }

  printf("Points inside master bbox: %d\n", bbox_points.size());

  point_set_write(bbox_points, S, "bbox_points.dat");

  // Make point cells
  std::vector<std::vector<E_Int>> point_cells = mesh_make_point_cells(S);

  // Isolate bbox cells
  std::set<E_Int> bbox_cells;

  for (auto point : bbox_points) {
    for (auto cell : point_cells[point])
      bbox_cells.insert(cell);
  }

  Mesh *C = mesh_extract_from_cell_set(S, bbox_cells);
  mesh_write(C, "bbox_cells.dat");

  /************/

  // We want to isolate ON_FACES and ON_CELLS
  std::stack<E_Int> stk;
  for (auto cell : bbox_cells) stk.push(cell);

  std::vector<E_Int> ON_FACES(S->nf, 0);
  std::set<E_Int> ON_CELLS;

  while (!stk.empty()) {
    E_Int cell = stk.top();
    stk.pop();

    E_Int nf = -1;
    E_Int *pf = mesh_get_cell(cell, nf, S);
    for (E_Int i = 0; i < nf; i++) {
      E_Int face = pf[i];

      E_Int nei = mesh_get_neighbour(cell, face, S);
      if (nei == -1) continue;

      // Mark nei as ON if it is not a bbox_cell
      if (bbox_cells.find(nei) == bbox_cells.end()) {
        ON_FACES[face] = 1;
        ON_CELLS.insert(nei);
      }
    }
  }

  Mesh *O = mesh_extract_from_cell_set(S, ON_CELLS);
  mesh_write(O, "on_cells.dat");

  std::vector<E_Int> on_face_list;
  for (E_Int i = 0; i < S->nf; i++) {
    if (ON_FACES[i])
      on_face_list.push_back(i);
  }

  Mesh *F = mesh_make_surface_mesh_from_face_list(S, on_face_list);
  mesh_write(F, "on_faces.dat");

  std::set<E_Int> kept_cells;
  for (E_Int i = 0; i < S->nc; i++) {
    if (bbox_cells.find(i) == bbox_cells.end())
      kept_cells.insert(i);
  }

  Mesh *K = mesh_extract_from_cell_set(S, kept_cells);
  mesh_write(K, "kept_cells.dat");

  std::vector<E_Int> base_faces;

  E_Float DIR[3] = {sqrt(2),0,sqrt(2)};

  for (auto face : on_face_list) {
    E_Float n[3];
    geom_compute_face_normal(face, S, n);

    E_Float val = K_MATH::dot(n, DIR, 3);

    //printf("(%f %f %f) -> %f\n", n[0], n[1], n[2], val);
    
    if (val > 0.8)
      base_faces.push_back(face);
  }

  Mesh *B = mesh_make_surface_mesh_from_face_list(S, base_faces);
  mesh_write(B, "base_faces.dat");

  /***********************************************************/

  // Master external faces
  auto efaces = mesh_get_external_faces_indices(M);

  std::vector<E_Int> FACES;
  for (auto face : efaces) {
    E_Int np = -1;
    E_Int *pn = mesh_get_face(face, np, M);
    E_Float Z = 0.0;
    for (E_Int i = 0; i < np; i++)
      Z += M->z[pn[i]];
    if (fabs(Z) <= 1e-6)
      FACES.push_back(face);
  }

  Mesh *MF = mesh_make_surface_mesh_from_face_list(M, FACES);
  mesh_write(MF, "master_patch.dat");

  // All base faces have the same normal: the direction of the projection
  E_Float Proj_dir[3];
  geom_compute_face_normal(base_faces[0], S, Proj_dir);
  printf("Projection direction: (%f, %f, %f)\n", Proj_dir[0],
    Proj_dir[1], Proj_dir[2]);

  // Gather the points to be projected
  std::set<E_Int> proj_points;

  for (auto face : base_faces) {
    E_Int np = -1;
    E_Int *pn = mesh_get_face(face, np, S);
    for (E_Int i = 0; i < np; i++) proj_points.insert(pn[i]);
  }

  printf("Points to be projected: %lu\n", proj_points.size());

  // Project base faces points on master skin polyhedron

  std::map<E_Int, Edge_Hit> point_hit_table;

  for (auto point : proj_points) {
    E_Float px = S->x[point];
    E_Float py = S->y[point];
    E_Float pz = S->z[point];
    E_Float qx = S->x[point] + Proj_dir[0] * bboxm.dmax;
    E_Float qy = S->y[point] + Proj_dir[1] * bboxm.dmax;
    E_Float qz = S->z[point] + Proj_dir[2] * bboxm.dmax;

    Edge_Hit EH;
    EH.T = -1;

    for (auto i = 0; i < FACES.size(); i++) {
      E_Int FACE = FACES[i];
      E_Int np = -1;
      E_Int *pn = mesh_get_face(FACE, np, M);

      assert(np == 4);

      E_Int hit = 0;

      // First triangle
      hit += geom_ray_triangle_intersect(px, py, pz, qx, qy, qz,
                                         M->x[pn[0]], M->y[pn[0]], M->z[pn[0]],
                                         M->x[pn[1]], M->y[pn[1]], M->z[pn[1]],
                                         M->x[pn[2]], M->y[pn[2]], M->z[pn[2]],
                                         EH);
      

      if (hit) {
        EH.T = i;
        point_hit_table[point] = EH;
        break;
      }

      // Second triangle
      hit += geom_ray_triangle_intersect(px, py, pz, qx, qy, qz,
                                         M->x[pn[0]], M->y[pn[0]], M->z[pn[0]],
                                         M->x[pn[2]], M->y[pn[2]], M->z[pn[2]],
                                         M->x[pn[3]], M->y[pn[3]], M->z[pn[3]],
                                         EH);
      

      if (hit) {
        EH.T = i;
        point_hit_table[point] = EH;
        break;
      }
    }

    // Point must hit!
    assert(EH.T != -1);
  }

  // Write edges and intersection points
  point_hits_write(point_hit_table, S, "edges.dat");

  // Write the intersected faces
  std::vector<E_Int> ifaces;
  for (const auto &edge : point_hit_table)
    ifaces.push_back(FACES[edge.second.T]);
  
  Mesh *IF = mesh_make_surface_mesh_from_face_list(M, ifaces);
  mesh_write(IF, "intersected_faces.dat");

  // Make the external faces connectivity
  std::vector<E_Int> fadj, xadj(1, 0);
  for (auto face : FACES) {
    E_Int np = -1;
    E_Int *pn = mesh_get_face(face, np, M);
    xadj.push_back(np);
    for (E_Int i = 0; i < np; i++)
      fadj.push_back(pn[i]);
  }

  for (auto i = 0; i < FACES.size(); i++) xadj[i+1] += xadj[i];

  std::vector<E_Int> fneis;
  K_CONNECT::build_face_neighbourhood(fadj, xadj, fneis);


  /*****************************************************/

  // Build the path of each edge on the master surface
  // Path: all the faces that an edge projection goes through

  std::set<Edge_NO> base_edges;

  std::map<Edge_NO, std::vector<E_Int>> proj_edges_path;

  for (auto face : base_faces) {
    E_Int np = -1;
    E_Int *pn = mesh_get_face(face, np, S);

    for (E_Int i = 0; i < np; i++) {
      E_Int p = pn[i];
      E_Int q = pn[(i+1)%np];

      Edge_NO e(p, q);

      if (base_edges.find(e) != base_edges.end())
        continue;
      
      base_edges.insert(e);

      // First time encountering this edge, do it

      assert(point_hit_table.find(e.q) != point_hit_table.end());

      // End point coordinates
      E_Float proj_qx = point_hit_table[e.q].x;
      E_Float proj_qy = point_hit_table[e.q].y;
      E_Float proj_qz = point_hit_table[e.q].z;

      // Find all the faces that [proj(p), proj(q)] traverses
      E_Int sf = point_hit_table[e.p].T; // start face
      E_Int ef = point_hit_table[e.q].T; // end face

      E_Int found = 0;

      proj_edges_path[e].push_back(FACES[sf]);

      // TODO(Imad): this is dangerous, loop might never terminate due to
      // floating point error...

      while (!found) {
        E_Int np = -1;
        E_Int *pn = mesh_get_face(FACES[sf], np, M);
        // End point lie to the left of all edges to terminate

        E_Int j = 0;
        for (; j < np; j++) {
          E_Int P = pn[j];
          E_Int Q = pn[(j+1)%np];

          if (geom_is_right(M->x[P], M->y[P], M->z[P],
                            M->x[Q], M->y[Q], M->z[Q],
                            proj_qx, proj_qy, proj_qz)) {
            sf = fneis[xadj[sf]+j];
            proj_edges_path[e].push_back(FACES[sf]);
            assert(sf != -1);
            break;
          }
        }

        if (j == np) found = 1;
      }

      assert(sf == ef);
    }
  }

  for (const auto &path : proj_edges_path) {
    for (auto face : path.second)
      printf("%d -> ",  face);
    puts("");


    if (path.second.size() == 6) {
      write_edge_path(M, path.first, path.second, "path.dat");
    }
  }

  return Py_None;
}
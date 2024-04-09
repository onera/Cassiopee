#include "proto.h"
#include <map>
#include <unordered_map>
#include <stack>

PyObject *K_XCORE::removeIntersectingKPlanes(PyObject *self, PyObject *args)
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

  printf("ni: %d\n", ni);
  printf("nj: %d\n", nj);
  printf("nk: %d\n", nk);

  if (ret <= 0) {
    RAISE("Bad slave mesh");
    return NULL;
  }

  if (ret != 1) {
    RAISE("Slave mesh should be structured.");
    RELEASESHAREDB(ret, SLAVE, fs, cns);
    RELEASESHAREDU(MASTER, fm, cnm);
    return NULL;
  }

  posx = K_ARRAY::isCoordinateXPresent(varString);
  posy = K_ARRAY::isCoordinateYPresent(varString);
  posz = K_ARRAY::isCoordinateZPresent(varString);

  if (posx == -1 || posy == -1 || posz == -1) {
    RELEASESHAREDU(MASTER, fm, cnm);
    RELEASESHAREDS(SLAVE, fs);
    RAISE("Slave point coordinates not found.");
    return NULL;
  }

  posx++; posy++; posz++;

  E_Float *Xs = fs->begin(posx);
  E_Float *Ys = fs->begin(posy);
  E_Float *Zs = fs->begin(posz);
  E_Int nps = fs->getSize();
  E_Int nfld = fs->getNfld();
  E_Int nij = ni*nj;

  // Init and orient master mesh
  Mesh *M = mesh_init(*cnm, Xm, Ym, Zm, npm);

  printf("Master cells: %d\n", M->nc);
  printf("Master points: %d\n", M->np);

  BBox bboxm = bbox_make(M);

  // Detect at which k does the slave mesh intersect the master mesh
  // Get projected points coordinates

  printf("Slave points: %d\n", nps);
  printf("Slave nfld: %d\n", nfld);

  E_Int kmax = 0;

  for (E_Int k = nk-1; k >= 0; k--) {
    E_Int inside = 0;
    
    for (E_Int j = 0; j < nj && !inside; j++) {
      for (E_Int i = 0; i < ni; i++) {
        E_Int iq = i + ni*j + nij*k;
        
        // If q in bbox, go to lower level plane
        if (bbox_is_point_inside(bboxm, Xs[iq], Ys[iq], Zs[iq])) {
          inside = 1;
          break;
        }
      }
    }

    if (inside) continue;

    // No points in this plane within the bbox, stop
    assert(inside == 0);
    kmax = k+1;
    break;
  }

  printf("kmax: %d\n", kmax);

  // Build returned array
  PyObject *tpl = K_ARRAY::buildArray3(nfld, varString, ni, nj, kmax);
  E_Float *fnp = K_ARRAY::getFieldPtr(tpl);
  FldArrayF subzone(nij*kmax, nfld, fnp, true);

  for (E_Int dim = 1; dim <= 3; dim++) {
    E_Float* sp = subzone.begin(dim);
    E_Float* fp = fs->begin(dim);

    for (E_Int k = 0; k < kmax; k++) {
      for (E_Int j = 0; j < nj; j++) {
        for (E_Int i = 0; i < ni; i++) {
          E_Int ind = i + j*ni + k*nij;
          sp[ind] = fp[ind];
        }
      }
    }
  }

  // Also return a list of indices of the points to be projected
  npy_intp dims[2];
  dims[0] = (npy_intp)nij;
  dims[1] = 1;

  PyArrayObject *PL = (PyArrayObject *)PyArray_SimpleNew(1, dims, E_NPY_INT);
  E_Int *np = (E_Int *)PyArray_DATA(PL);
  E_Int kshift = (kmax-1)*nij;

  for (E_Int j = 0; j < nj; j++) {
    for (E_Int i = 0; i < ni; i++) {
      E_Int ind = i + j*ni + kshift;
      *np++ = ind;
    }
  }

  np = (E_Int *)PyArray_DATA(PL);
  write_point_list(Xs, Ys, Zs, np, nij, "POINTS");

  /***************************************************/

  // Master external faces
  auto efaces = mesh_get_external_faces_indices(M);

  // TODO(Imad): input master patch instead of geom identification
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
  mesh_write(MF, "master_patch");

  // Make the external faces connectivity
  std::vector<E_Int> fadj, xadj(1, 0);
  for (auto face : FACES) {
    E_Int np = -1;
    E_Int *pn = mesh_get_face(face, np, M);
    xadj.push_back(np);
    assert(np == 4);
    fadj.push_back(pn[0]);
    fadj.push_back(pn[3]);
    fadj.push_back(pn[2]);
    fadj.push_back(pn[1]);
  }

  for (auto i = 0; i < FACES.size(); i++) xadj[i+1] += xadj[i];

  std::vector<E_Int> fneis;
  K_CONNECT::build_face_neighbourhood(fadj, xadj, fneis);

  // Master patch Point-To-Face
  std::unordered_map<E_Int, std::vector<E_Int>> point_faces;

  for (size_t i = 0; i < FACES.size(); i++) {
    for (E_Int j = xadj[i]; j < xadj[i+1]; j++) {
      E_Int point = fadj[j];
      point_faces[point].push_back(i);
    }
  }

  // Make edges
  std::set<Edge_NO> master_edges;

  for (size_t i = 0; i < FACES.size(); i++) {
    for (E_Int j = xadj[i]; j < xadj[i+1]; j++) {
      E_Int p = fadj[j];
      E_Int q = fadj[(j+1)%4];
      Edge_NO e(p, q);
      master_edges.insert(e);
    }
  }

  printf("Master edges: %zu\n", master_edges.size());

  // Project points onto master surface
  std::unordered_map<E_Int, Edge_Hit> point_hit_table;

  for (E_Int i = 0; i < nij; i++) {
    E_Int p = np[i];
    E_Int q = p + nij;

    E_Float px = Xs[p];
    E_Float py = Ys[p];
    E_Float pz = Zs[p];
    E_Float qx = Xs[q];
    E_Float qy = Ys[q];
    E_Float qz = Zs[q];
    E_Float dir[3] = {qx-px, qy-py, qz-pz};
    E_Float NORM = K_MATH::norm(dir, 3);
    E_Float dx = dir[0] / NORM;
    E_Float dy = dir[1] / NORM;
    E_Float dz = dir[2] / NORM;

    // TODO(Imad): improve
    // For now intersect with all master patch triangles

    Edge_Hit EH;
    EH.F = -1;
    EH.T = -1;

    for (size_t fid = 0; fid < FACES.size(); fid++) {
      
      E_Int *pn = &fadj[xadj[fid]];
      E_Int A = pn[0];
      E_Int B = pn[1];
      E_Int C = pn[2];
      E_Int D = pn[3];

      // First triangle ABC

      E_Int hit = geom_ray_triangle_intersect(px, py, pz, dx, dy, dz,
                                         M->x[A], M->y[A], M->z[A],
                                         M->x[B], M->y[B], M->z[B],
                                         M->x[C], M->y[C], M->z[C],
                                         EH);
      
      if (hit) {
        EH.T = 0;
        EH.F = fid;

        E_Float u, v, w;
        u = EH.u;
        v = EH.v;
        w = EH.w;

        // Classify hit points P within a triangle ABC:
        // 1- P is A
        // 2- P is B
        // 3- P is C
        // 4- P in [AB]
        // 5- P in [BC]
        // 6- P in [AC]
        // 7- P inside ABC

        if (FEQ(u, 1) && FEQ(v, 0) && FEQ(w, 0)) {
          EH.is_vertex = 0;
        } else if (FEQ(u, 0) && FEQ(v, 1) && FEQ(w, 0)) {
          EH.is_vertex = 1;
        } else if (FEQ(u, 0) && FEQ(v, 0) && FEQ(w, 1)) {
          EH.is_vertex = 2;
        } else if (!FEQ(u, 0.0) && !FEQ(v, 1.0) && FEQ(w, 0.0)) {
          // P on AB
          EH.on_edge = 0;
        } else if (FEQ(u, 0.0) && !FEQ(v, 1.0) && !FEQ(w, 0.0)) {
          // P on BC
          EH.on_edge = 1;
        }

        point_hit_table[p] = EH;

        break;
      }

      // First triangle CDA

      hit = geom_ray_triangle_intersect(px, py, pz, dx, dy, dz,
                                        M->x[C], M->y[C], M->z[C],
                                        M->x[D], M->y[D], M->z[D],
                                        M->x[A], M->y[A], M->z[A],
                                        EH);
      
      if (hit) {
        EH.T = 1;
        EH.F = fid;

        E_Float u, v, w;
        u = EH.u;
        v = EH.v;
        w = EH.w;

        if (FEQ(u, 1) && FEQ(v, 0) && FEQ(w, 0)) {
          EH.is_vertex = 2;
        } else if (FEQ(u, 0) && FEQ(v, 1) && FEQ(w, 0)) {
          EH.is_vertex = 3;
        } else if (FEQ(u, 0) && FEQ(v, 0) && FEQ(w, 1)) {
          EH.is_vertex = 0;
        } else if (!FEQ(u, 0.0) && !FEQ(v, 1.0) && FEQ(w, 0.0)) {
          // P on CD
          EH.on_edge = 2;
        } else if (FEQ(u, 0.0) && !FEQ(v, 1.0) && !FEQ(w, 0.0)) {
          // P on DA
          EH.on_edge = 3;
        }

        point_hit_table[p] = EH;

        break;
      }
    }

    // Point must hit!
    assert(EH.T != -1);
    assert(EH.F != -1);
  }

  point_hits_write(point_hit_table, "hits");
  puts("ok point hits");

  /***********************************/

  // Master dcel
  dcel Dm(MF);

  // Slave dcel
  std::vector<E_Float> X(nij), Y(nij), Z(nij);
  for (E_Int j = 0; j < nj; j++) {
    for (E_Int i = 0; i < ni; i++) {
      E_Int coord = i + ni*j;
      E_Int p = coord + kshift;
      const auto &EHp = point_hit_table[p];
      X[coord] = EHp.x;
      Y[coord] = EHp.y;
      Z[coord] = EHp.z;
    }
  }

  dcel Ds(ni, nj, X, Y, Z);

  dcel D(Dm, Ds);

  D.resolve();


  PyObject *out = PyList_New(0);
  PyList_Append(out, tpl);
  PyList_Append(out, (PyObject *)PL);
  Py_DECREF(tpl);
  Py_DECREF(PL);
  
  RELEASESHAREDU(MASTER, fm, cnm);
  RELEASESHAREDS(SLAVE, fs);
  // TODO(Imad): free M
  
  return out;
}

PyObject *K_XCORE::intersectMesh(PyObject *self, PyObject *args)
{

  return Py_None;
}

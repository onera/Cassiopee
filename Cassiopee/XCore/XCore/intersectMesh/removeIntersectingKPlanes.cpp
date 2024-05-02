#include "proto.h"
#include <unordered_map>

PyObject *K_XCORE::removeIntersectingKPlanes(PyObject *self, PyObject *args)
{
  PyObject *MASTER, *SLAVE, *PATCH;
  
  if (!PYPARSETUPLE_(args, OOO_, &MASTER, &SLAVE, &PATCH)) {
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
    RELEASESHAREDU(MASTER, fm, cnm);
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

  // Check intersection patch
  E_Int *patch = NULL;
  E_Int patch_size = -1;
  ret = K_NUMPY::getFromNumpyArray(PATCH, patch, patch_size, true);
  if (ret != 1) {
    RELEASESHAREDU(MASTER, fm, cnm);
    RELEASESHAREDS(SLAVE, fs);
    RAISE("Bad master patch.");
    return NULL;
  }

  for (E_Int i = 0; i < patch_size; i++)
    patch[i] -= 1;

  // Init and orient master mesh
  Mesh *M = mesh_init(*cnm, Xm, Ym, Zm, npm);

  BBox bboxm = bbox_make(M);

  /**************************************************************************/

  // Detect at which k does the slave mesh intersect the master mesh
  // Get projected points coordinates

  // Max plane index that doesn't intersection with master bbox (zero-based)
  E_Int kmax = 0;

  for (E_Int k = 0; k < nk; k++) {
    E_Int inside = 0;

    for (E_Int j = 0; j < nj && !inside; j++) {
      for (E_Int i = 0; i < ni; i++) {
        E_Int p = i + ni*j + nij*k;

        // If p in bbox, stop
        if (bbox_is_point_inside(bboxm, Xs[p], Ys[p], Zs[p])) {
          inside = 1;
          kmax = k;
          break;
        }
      }
    }

    if (inside)
      break;
  }

  // Points to be projected nij*(kmax-1) .. nij*kmax
  std::vector<E_Int> np(nij);
  E_Int kshift = (kmax-1)*nij;
  E_Int *ptr = &np[0];
  for (E_Int j = 0; j < nj; j++) {
    for (E_Int i = 0; i < ni; i++) {
      E_Int ind = i + j*ni + kshift;
      *ptr++ = ind;
    }
  }

  write_point_list(Xs, Ys, Zs, &np[0], nij, "POINTS");

  /**************************************************************************/

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

    for (size_t fid = 0; fid < patch_size; fid++) {
      E_Int face = patch[fid];
      E_Int np = -1;
      E_Int *pn = mesh_get_face(face, np, M);
      assert(np == 4);

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

        point_hit_table[p] = EH;

        break;
      }

      // Second triangle CDA

      hit = geom_ray_triangle_intersect(px, py, pz, dx, dy, dz,
                                        M->x[C], M->y[C], M->z[C],
                                        M->x[D], M->y[D], M->z[D],
                                        M->x[A], M->y[A], M->z[A],
                                        EH);
      
      if (hit) {
        EH.T = 1;
        EH.F = fid;

        point_hit_table[p] = EH;

        break;
      }
    }

    // Point must hit!
    assert(EH.T != -1);
    assert(EH.F != -1);
  }

  point_hits_write(point_hit_table, "hits");

  /*************************************************************************/

  // Make out cartesian mesh
  PyObject *tpl;
  nk = kmax + 1; 
  tpl = K_ARRAY::buildArray3(3, "x,y,z", ni, nj, nk, 3);

  K_FLD::FldArrayF *f;
  K_FLD::FldArrayI *c;
  K_ARRAY::getFromArray3(tpl, varString, f, ni, nj, nk, c, eltType);

  E_Float* xt = f->begin(1);
  E_Float* yt = f->begin(2);
  E_Float* zt = f->begin(3);

  // Copy all the points up to kmax
  for (E_Int k = 0; k < kmax; k++) {
    for (E_Int j = 0; j < nj; j++) {
      for (E_Int i = 0; i < ni; i++) {
        E_Int ind = i + j*ni + k*nij;
        xt[ind] = Xs[ind];
        yt[ind] = Ys[ind];
        zt[ind] = Zs[ind];
      }
    }
  }

  // Copy the projected points
  for (E_Int i = 0; i < nij; i++) {
    E_Int p = np[i];
    auto EH = point_hit_table[p];
    E_Float x = EH.x;
    E_Float y = EH.y;
    E_Float z = EH.z;
    xt[p+nij] = x;
    yt[p+nij] = y;
    zt[p+nij] = z;
  }

  // Tag the projected points
  npy_intp dims[2];
  dims[1] = 1;
  dims[0] = (npy_intp)ni*nj*nk;
  PyArrayObject *tag = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
  E_Float *ptag = (E_Float *)PyArray_DATA(tag);
  for (E_Int i = 0; i < nij*kmax; i++) ptag[i] = 0.0;
  for (E_Int i = nij*kmax; i < nij*nk; i++) ptag[i] = 1.0;

  PyObject *out = PyList_New(0);
  PyList_Append(out, tpl);
  PyList_Append(out, (PyObject *)tag);
  RELEASESHAREDS(tpl, f);
  Py_DECREF(tpl);
  Py_DECREF(tag);
  RELEASESHAREDU(MASTER, fm, cnm);
  RELEASESHAREDS(SLAVE, fs);
  // TODO(Imad): free M

  return out;
}
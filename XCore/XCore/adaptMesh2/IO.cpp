#include "Proto.h"

//#define PRINTPATCH

static
PyObject *extract_conformal_mesh(AMesh *M, E_Int closed_mesh)
{
  // Count nface size
  E_Int nface_size = 0;

  M->closed_indPH = (E_Int *)XRESIZE(M->closed_indPH, (M->ncells+1)*sizeof(E_Int));
  M->closed_indPH[0] = 0;

  for (E_Int i = 0; i < M->ncells; i++) {
    E_Int nf = -1;
    E_Int pf[24];
    // TODO(Imad): a function that returns nf only
    get_full_cell(i, M, nf, pf);

    nface_size += nf;
    M->closed_indPH[i+1] = nf;
  }

  for (E_Int i = 0; i < M->ncells; i++)
    M->closed_indPH[i+1] += M->closed_indPH[i];

  assert(M->closed_indPH[M->ncells] == nface_size);

  E_Int ngon_size = -1;

  if (closed_mesh) {
    puts("Exporting closed mesh");
    close_mesh(M);
    ngon_size = M->closed_indPG[M->nfaces];
  } else {
    puts("Exporting non-closed mesh");
    ngon_size = M->indPG[M->nfaces];
  }

  //M->closed_nface = (E_Int *)XRESIZE(M->closed_nface, nface_size * sizeof(E_Int));

  const char *varStringOut = "CoordinateX,CoordinateY,CoordinateZ";

  PyObject *m = K_ARRAY::buildArray3(3, varStringOut, M->npoints, M->ncells,
    M->nfaces, "NGON", ngon_size, nface_size, 3, false, 3);

  FldArrayF *fo;
  FldArrayI *cno;
  K_ARRAY::getFromArray3(m, fo, cno);

  E_Int *ngono  = cno->getNGon();
  E_Int *nfaceo = cno->getNFace();
  E_Int *indPGo = cno->getIndPG();
  E_Int *indPHo = cno->getIndPH();
  
  memcpy(indPHo, M->closed_indPH, (M->ncells+1)*sizeof(E_Int));
  assert(indPHo[M->ncells] == nface_size);

  for (E_Int i = 0; i < M->ncells; i++) {
    E_Int nf = -1;
    E_Int pf[24];
    get_full_cell(i, M, nf, pf);
    E_Int *ptr = &nfaceo[indPHo[i]];
    //E_Int *p = &M->closed_nface[M->closed_indPH[i]];
    for (E_Int j = 0; j < nf; j++) {
      ptr[j] = pf[j]+1;
      //p[j] = pf[j];
    }
  }

  if (closed_mesh) {
    for (E_Int i = 0; i < ngon_size; i++)
      ngono[i] = M->closed_ngon[i] + 1;
    
    memcpy(indPGo, M->closed_indPG, (M->nfaces+1) * sizeof(E_Int));
  } else {
    for (E_Int i = 0; i < ngon_size; i++)
      ngono[i] = M->ngon[i] + 1;
    
    memcpy(indPGo, M->indPG, (M->nfaces+1) * sizeof(E_Int));
  }

  E_Float *px = fo->begin(1);
  E_Float *py = fo->begin(2);
  E_Float *pz = fo->begin(3);
  memcpy(px, M->x, M->npoints * sizeof(E_Float));
  memcpy(py, M->y, M->npoints * sizeof(E_Float));
  memcpy(pz, M->z, M->npoints * sizeof(E_Float));

  assert(K_CONNECT::check_open_cells(*cno, NULL) == 0);
  assert(K_CONNECT::check_overlapping_cells(*cno) == 0);

  // BCs
  PyObject *BCS = PyList_New(0);

  npy_intp dims[2];
  
  for (E_Int i = 0; i < M->nbc; i++) {
    dims[0] = (npy_intp)M->bcsizes[i];
    dims[1] = 1;

    PyArrayObject *PL = (PyArrayObject *)PyArray_SimpleNew(1, dims, E_NPY_INT);
    E_Int *op = M->ptlists[i];
    E_Int *np = (E_Int *)PyArray_DATA(PL);
    for (E_Int j = 0; j < M->bcsizes[i]; j++)
      *np++ = op[j] + 1;

    PyObject *tpl = Py_BuildValue("[Os]", (PyObject *)PL, M->bcnames[i]);
    Py_DECREF(PL);
    PyList_Append(BCS, tpl);
    Py_DECREF(tpl);
  }

#ifdef PRINTPATCH
  // Patch faces as BCs
  for (E_Int i = 0; i < M->npatches; i++) {
    Patch *P = &M->patches[i];

    dims[0] = (npy_intp)P->nf;
    dims[1] = 1;

    PyArrayObject *PL = (PyArrayObject *)PyArray_SimpleNew(1, dims, E_NPY_INT);
    E_Int *op = P->pf;
    E_Int *np = (E_Int *)PyArray_DATA(PL);
    for (E_Int j = 0; j < P->nf; j++)
      *np++ = op[j] + 1;

    char name[64];
    sprintf(name, "proc%d-%d", M->pid, P->nei);

    PyObject *tpl = Py_BuildValue("[Os]", (PyObject *)PL, name);
    Py_DECREF(PL);
    PyList_Append(BCS, tpl);
    Py_DECREF(tpl);
  }
#endif

  // parent elements
  dims[0] = (npy_intp)M->nfaces;
  dims[1] = 1;

  PyArrayObject *OWN = (PyArrayObject *)PyArray_SimpleNew(1, dims, E_NPY_INT);
  PyArrayObject *NEI = (PyArrayObject *)PyArray_SimpleNew(1, dims, E_NPY_INT);
  E_Int *owner = (E_Int *)PyArray_DATA(OWN);
  E_Int *neigh = (E_Int *)PyArray_DATA(NEI);

  memcpy(owner, M->owner, M->nfaces*sizeof(E_Int));
  memcpy(neigh, M->neigh, M->nfaces*sizeof(E_Int));

  PyObject *out = PyList_New(0);
  PyList_Append(out, m);
  PyList_Append(out, BCS);
  PyList_Append(out, (PyObject *)OWN);
  PyList_Append(out, (PyObject *)NEI);
  Py_DECREF(m);
  Py_DECREF(BCS);
  Py_DECREF(OWN);
  Py_DECREF(NEI);

  return out;
}

// Assumes only HEXA for now
static
PyObject *extract_BE_mesh(AMesh *M)
{
  const char *varString = "CoordinateX,CoordinateY,CoordinateZ";

  PyObject *m = K_ARRAY::buildArray3(3, varString, M->npoints, M->ncells,
    "HEXA", false, 3);
  
  FldArrayF *f;
  FldArrayI *cn;

  K_ARRAY::getFromArray3(m, f, cn);

  E_Int *pcn = cn->begin();
  E_Int l = 0;

  for (E_Int i = 0; i < M->ncells; i++) {
    E_Int *pf = get_facets(i, M->nface, M->indPH);

    E_Int pb[4], P[8];

    E_Int bot = pf[0];
    E_Int top = pf[1];
    E_Int lft = pf[2];

    E_Int clvl = M->cellTree->level(i);

    // BOT
    E_Int flvl = M->faceTree->level(bot);
    if (flvl == clvl) {
      E_Int *pn = get_facets(bot, M->ngon, M->indPG);
      for (E_Int j = 0; j < 4; j++) P[j] = pn[j];
      E_Int reorient = get_reorient(bot, i, normalIn_H[0], M);
      if (reorient) std::swap(P[1], P[3]);
    } else {
      reconstruct_parent_quad(bot, M, pb);
      for (E_Int j = 0; j < 4; j++) P[j] = pb[j];
      E_Int reorient = get_reorient(bot, i, normalIn_H[0], M);
      if (reorient) std::swap(P[1], P[3]);
    }

    // LFT
    flvl = M->faceTree->level(lft);
    if (flvl == clvl) {
      E_Int *pn = get_facets(lft, M->ngon, M->indPG);
      E_Int local[4];
      for (E_Int j = 0; j < 4; j++) local[j] = pn[j];
      E_Int reorient = get_reorient(lft, i, normalIn_H[2], M);
      E_Int i0 = Get_pos(P[0], local, 4);
      assert(i0 != -1);
      Right_shift(local, i0, 4);
      if (reorient) std::swap(local[1], local[3]);
      assert(local[0] == P[0]);
      assert(local[1] == P[3]);
      P[7] = local[2];
      P[4] = local[3];
    } else {
      reconstruct_parent_quad(lft, M, pb);
      E_Int reorient = get_reorient(lft, i, normalIn_H[2], M);
      E_Int i0 = Get_pos(P[0], pb, 4);
      assert(i0 != -1);
      Right_shift(pb, i0, 4);
      if (reorient) std::swap(pb[1], pb[3]);
      assert(pb[0] == P[0]);
      assert(pb[1] == P[3]);
      P[7] = pb[2];
      P[4] = pb[3];
    }

    // TOP
    flvl = M->faceTree->level(top);
    if (flvl == clvl) {
      E_Int *pn = get_facets(top, M->ngon, M->indPG);
      E_Int local[4];
      for (E_Int j = 0; j < 4; j++) local[j] = pn[j];
      E_Int reorient = get_reorient(top, i, normalIn_H[1], M);
      E_Int i0 = Get_pos(P[4], local, 4);
      assert(i0 != -1);
      Right_shift(local, i0, 4);
      if (reorient) std::swap(local[1], local[3]);
      assert(local[0] == P[4]);
      assert(local[3] == P[7]);
      P[5] = local[1];
      P[6] = local[2];
    } else {
      reconstruct_parent_quad(top, M, pb);
      E_Int reorient = get_reorient(top, i, normalIn_H[1], M);
      E_Int i0 = Get_pos(P[4], pb, 4);
      assert(i0 != -1);
      Right_shift(pb, i0, 4);
      if (reorient) std::swap(pb[1], pb[3]);
      assert(pb[0] == P[4]);
      assert(pb[3] == P[7]);
      P[5] = pb[1];
      P[6] = pb[2];
    }

    // One-based
    for (E_Int j = 0; j < 8; j++) {
      pcn[l++] = P[j]+1;
    }
  }

  E_Float *px = f->begin(1);
  E_Float *py = f->begin(2);
  E_Float *pz = f->begin(3);

  memcpy(px, M->x, M->npoints*sizeof(E_Float));
  memcpy(py, M->y, M->npoints*sizeof(E_Float));
  memcpy(pz, M->z, M->npoints*sizeof(E_Float));

  return m;
}

#define BCMODE_BE 0
#define BCMODE_NGON 1

static
PyObject *extractBCMesh_NGON(AMesh *M)
{
  // BCs
  PyObject *BCs = PyList_New(0);

  npy_intp dims[2];
  
  for (E_Int i = 0; i < M->nbc; i++) {
    dims[0] = (npy_intp)M->bcsizes[i];
    dims[1] = 1;

    PyArrayObject *PL = (PyArrayObject *)PyArray_SimpleNew(1, dims, E_NPY_INT);
    E_Int *op = M->ptlists[i];
    E_Int *np = (E_Int *)PyArray_DATA(PL);
    for (E_Int j = 0; j < M->bcsizes[i]; j++)
      *np++ = op[j] + 1;

    PyObject *tpl = Py_BuildValue("[Os]", (PyObject *)PL, M->bcnames[i]);
    Py_DECREF(PL);
    PyList_Append(BCs, tpl);
    Py_DECREF(tpl);
  }

  return BCs;
}

static
PyObject *extractBCMesh_BE(AMesh *M)
{
  // Count number of bc faces

  std::vector<E_Int> be_sizes(M->nbc, 0);
  
  for (E_Int i = 0; i < M->nbc; i++) {
    E_Int *ptlist = M->ptlists[i];

    for (E_Int j = 0; j < M->bcsizes[i]; j++) {
      E_Int type = M->faceTree->type(ptlist[i]);
      
      if (type == QUAD) be_sizes[i] += 4;
      else if (type == TRI) be_sizes[i] += 3;
      else assert(0);
    }
  }

  PyObject *BCs = PyList_New(0);

  npy_intp dims[2];

  for (E_Int i = 0; i < M->nbc; i++) {
    dims[0] = (npy_intp)be_sizes[i];
    dims[1] = 1;

    PyArrayObject *PL = (PyArrayObject *)PyArray_SimpleNew(1, dims, E_NPY_INT);
    E_Int *op = M->ptlists[i];
    E_Int *np = (E_Int *)PyArray_DATA(PL);

    for (E_Int j = 0; j < M->bcsizes[i]; j++) {
      E_Int face = op[j];
      E_Int stride = -1;
      E_Int *pn = get_face(face, stride, M->ngon, M->indPG);
      for (E_Int k = 0; k < stride; k++)
        *np++ = pn[k]+1;
    }

    PyObject *tpl = Py_BuildValue("[Os]", (PyObject *)PL, M->bcnames[i]);
    Py_DECREF(PL);
    PyList_Append(BCs, tpl);
    Py_DECREF(tpl);
  }

  return BCs;
}

PyObject *K_XCORE::extractBoundaryMesh(PyObject *self, PyObject *args)
{
  PyObject *AMESH;
  E_Int MODE; // BE or NGON

  if (!PYPARSETUPLE_(args, O_ I_, &AMESH, &MODE)) {
    RAISE("Wrong input.");
    return NULL;
  }

  if (!PyCapsule_IsValid(AMESH, "AMesh")) {
    RAISE("Bad mesh capsule");
    return NULL;
  }

  AMesh *M = (AMesh *)PyCapsule_GetPointer(AMESH, "AMesh");

  assert(MODE == BCMODE_BE || MODE == BCMODE_NGON);

  if (MODE == BCMODE_NGON) return extractBCMesh_NGON(M);

  return extractBCMesh_BE(M);
}

PyObject *K_XCORE::ExtractLeafMesh(PyObject *self, PyObject *args)
{
  PyObject *AMESH;
  E_Int conformize;

  if (!PYPARSETUPLE_(args, O_ I_, &AMESH, &conformize)) {
    RAISE("Wrong input.");
    return NULL;
  }

  if (!PyCapsule_IsValid(AMESH, "AMesh")) {
    RAISE("Bad AMesh capsule.");
    return NULL;
  }

  AMesh *M = (AMesh *)PyCapsule_GetPointer(AMESH, "AMesh");

  if (conformize)
    return extract_conformal_mesh(M, 1);
  else
    return extract_BE_mesh(M);
}

PyObject *K_XCORE::_assignRefDataToAM(PyObject *self, PyObject *args)
{
  PyObject *AMESH, *REF;
  if (!PYPARSETUPLE_(args, OO_, &AMESH, &REF)) {
    RAISE("Bad input.");
    return NULL;
  }

  if (!PyCapsule_IsValid(AMESH, "AMesh")) {
    RAISE("Bad AdaptMesh hook.");
    return NULL;
  }

  AMesh *M = (AMesh *)PyCapsule_GetPointer(AMESH, "AMesh");

  E_Int ret, nfld, size;
  E_Int *ref_data = NULL;
  ret = K_NUMPY::getFromNumpyArray(REF, ref_data, size, nfld, true);
  assert(ret == 1 && nfld == 1 && size == M->ncells);

  M->ref_data = (E_Int *)XRESIZE(M->ref_data, M->ncells*sizeof(E_Int));
  for (E_Int i = 0; i < M->ncells; i++) {
    M->ref_data[i] = (E_Int)ref_data[i]; 
  }

  //smooth_ref_data(M);
  smooth_ref_data(M);

  return Py_None;
}

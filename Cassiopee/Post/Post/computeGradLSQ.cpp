/*    
    Copyright 2013-2026 ONERA.

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
#include "post.h"

static
E_Int parse_pointlists_and_rfields(PyObject *ptlists, PyObject *rfields,
  std::vector<E_Int *> &pfaces, std::vector<E_Int> &npfaces,
  std::vector<std::vector<E_Float *>> &rflds)
{
  //if (PyList_Size(ptlists) == 0) return 0;
  if (ptlists == Py_None) return 0;

  E_Int psize = PyList_Size(ptlists);
  pfaces.resize(psize);
  rflds.resize(psize);
  npfaces.resize(psize);
    
  for (E_Int i = 0; i < psize; i++) 
  {
    PyObject *patch = PyList_GetItem(ptlists, i);
    E_Int ret = K_NUMPY::getFromNumpyArray(patch, pfaces[i], npfaces[i]);
    if (ret != 1)
    {     
      PyErr_SetString(PyExc_TypeError, "computeGradLSQ: incorrect numpy.");
      return 0;
    }
  }

  E_Int rsize = PyList_Size(rfields);
  assert(rsize == psize);

  for (E_Int i = 0; i < rsize; i++) 
  {
    PyObject *flds = PyList_GetItem(rfields, i);
    E_Int nfld = PyList_Size(flds);
    rflds[i].resize(nfld);
    for (E_Int j = 0; j < nfld; j++) 
    {
      PyObject *fld = PyList_GetItem(flds, j);
      E_Int size = -1;
      E_Int ret = K_NUMPY::getFromNumpyArray(fld, rflds[i][j], size);
      if (ret != 1)
      {     
        PyErr_SetString(PyExc_TypeError, "computeGradLSQ: incorrect numpy.");
        return 0;
      }
      assert(size == npfaces[i]);
    }
  }

  return 1;
}

static
void make_A_grad_matrices(K_FLD::FldArrayI &cn, E_Int *owner, E_Int *neigh,
  E_Int *count_neis, E_Float *fcenters, E_Float *cx, E_Float *cy, E_Float *cz,
  E_Float *lsqG)
{
  E_Int ncells = cn.getNElts();
  memset(count_neis, 0, ncells*sizeof(E_Int));

  // Internal faces
  E_Int nfaces = cn.getNFaces();
  E_Int *indPH = cn.getIndPH();

  for (E_Int i = 0; i < nfaces; i++) {
    E_Int nei = neigh[i]-1;
    if (nei == -1) continue;

    E_Int own = owner[i]-1;

    E_Float d[3] = {cx[nei]-cx[own], cy[nei]-cy[own], cz[nei]-cz[own]};

    E_Float *lo = &lsqG[3*(indPH[own] + count_neis[own]++)];
    E_Float *ln = &lsqG[3*(indPH[nei] + count_neis[nei]++)];

    for (E_Int j = 0; j < 3; j++) {
      lo[j] = d[j];
      ln[j] = -d[j];
    }
  }

  // TODO(Imad): boundary contributions
}

static
void correct_A_grad_matrices(K_FLD::FldArrayI &cn, E_Int *owner,
  E_Int *count_neis, E_Float *CX, E_Float *CY, E_Float *CZ,
  const std::vector<E_Int *> &pfaces, const std::vector<E_Int> &npfaces,
  const std::vector<std::vector<E_Float *>> &rfields,
  E_Float *lsqG)
{
  E_Int *indPH = cn.getIndPH();

  for (size_t i = 0; i < pfaces.size(); i++) {
    // Coordinates are the last three vectors of rfields[i]
    E_Int size = rfields[i].size();
    const auto &cx = rfields[i][size-3];
    const auto &cy = rfields[i][size-2];
    const auto &cz = rfields[i][size-1];

    const auto &patch = pfaces[i];
    for (E_Int j = 0; j < npfaces[i]; j++) {
      E_Int face = patch[j]-1;
      assert(face >= 0 && face < cn.getNFaces());
      E_Int own = owner[face]-1;
      assert(own >= 0 && own < cn.getNElts());

      E_Float *lo = &lsqG[3*(indPH[own] + count_neis[own]++)];

      lo[0] = cx[j] - CX[own];
      lo[1] = cy[j] - CY[own];
      lo[2] = cz[j] - CZ[own];
    }
  }
}

static
void deduce_lsq_grad_matrices(K_FLD::FldArrayI &cn, const E_Int *count_neis,
  const E_Float *lsqG, E_Float *lsqGG)
{
  E_Int ncells = cn.getNElts();
  E_Int *indPH = cn.getIndPH();

  for (E_Int i = 0; i < ncells; i++) {
    E_Float *pGG = &lsqGG[9*i];
    const E_Float *pG = &lsqG[3*indPH[i]];
    E_Int idx = 0;
    for (E_Int j = 0; j < 3; j++) {
      for (E_Int k = 0; k < 3; k++) {
        pGG[idx] = 0.0;
        for (E_Int l = 0; l < count_neis[i]; l++)
          pGG[idx] += pG[l*3+j] * pG[l*3+k];
        idx++;
      }
    }
  }
}

static
void make_RHS_vector(K_FLD::FldArrayI &cn, E_Int *count_neis, E_Int *owner,
  E_Int *neigh, const E_Float *Field, E_Float *b)
{
  E_Int ncells = cn.getNElts();
  memset(count_neis, 0, ncells*sizeof(E_Int));

  // Construct b vector
  E_Int nfaces = cn.getNFaces();
  E_Int *indPH = cn.getIndPH();

  for (E_Int i = 0; i < nfaces; i++) {
    E_Int nei = neigh[i]-1;
    if (nei == -1) continue;

    E_Int own = owner[i]-1;
    b[indPH[own] + count_neis[own]++] = Field[nei] - Field[own];
    b[indPH[nei] + count_neis[nei]++] = Field[own] - Field[nei];
  }

  // TODO(Imad): boundary faces contributions
}

static
void correct_RHS_vector(K_FLD::FldArrayI &cn, E_Int *count_neis, E_Int *owner,
  const E_Float *Field, const std::vector<E_Int *> &pfaces,
  const std::vector<E_Int> &npfaces,
  const std::vector<std::vector<E_Float *>> &rfields, E_Int field_idx,
  E_Float *b)
{
  E_Int *indPH = cn.getIndPH();

  for (size_t i = 0; i < pfaces.size(); i++) {
    const auto &patch = pfaces[i];
    const auto &rflds = rfields[i];

    // We're dealing with field of index field_idx
    const auto &rfld = rflds[field_idx];

    for (E_Int j = 0; j < npfaces[i]; j++) {
      E_Int face = patch[j]-1;
      E_Int own = owner[face]-1;
      assert(own >= 0 && own < cn.getNElts());
      b[indPH[own] + count_neis[own]++] = rfld[j] - Field[own];
      assert(count_neis[own] <= 6);
    }
  }
}

static
E_Int make_gradient(K_FLD::FldArrayI &cn, E_Int *count_neis, E_Float *b,
  E_Float *lsqG, E_Float *lsqGG, E_Float *gx, E_Float *gy, E_Float *gz)
{
  E_Int ncells = cn.getNElts();
  E_Int *indPH = cn.getIndPH();
  
  E_Int bad_gradient = 0;
  E_Float B[3];
  E_Float G[3];

  for (E_Int i = 0; i < ncells; i++) {
    // Construct B vector
    E_Float *pG = &lsqG[3*indPH[i]];
    E_Float *pb = &b[indPH[i]];
    for (E_Int j = 0; j < 3; j++) {
      B[j] = 0.0;
      for (E_Int k = 0; k < count_neis[i]; k++)
        B[j] += pb[k] * pG[k*3+j];
    }

    // Solve
    E_Int converged = K_LINEAR::BiCGStab(&lsqGG[9*i], B, 3, G);

    // Copy
    gx[i] = G[0];
    gy[i] = G[1];
    gz[i] = G[2];

    bad_gradient |= !converged;
    if (!converged) {
      printf("Gradient: BiCGStab not converged\n");
      printf("Cell: " SF_D_ "\n", i);
      E_Float *pt = &lsqGG[9*i];
      for (E_Int j = 0; j < 3; j++) {
        for (E_Int k = 0; k < 3; k++) {
          printf("%.3e ", pt[3*j+k]);
        }
        puts("");
      }
      for (E_Int j = 0; j < 3; j++) {
        printf("%.3e ", B[j]);
      }
      puts("");
    }
  }

  return bad_gradient;
}

/* Compute least-squares gradients of input fields on a well-oriented mesh */
PyObject *K_POST::computeGradLSQ(PyObject *self, PyObject *args)
{
  PyObject *arr, *fields, *pe, *cx, *cy, *cz, *fcenters, *ptlists, *rfields;
  if (!PYPARSETUPLE_(args, OOOO_ OOOO_ O_, &arr, &fields, &pe, &cx, &cy, &cz,
    &fcenters, &ptlists, &rfields)) 
  {
    PyErr_SetString(PyExc_ValueError, "computeGradLSQ: wrong input.");
    return NULL;
  }

  E_Int ret;

  // Check fields
  E_Int fsize = PyList_Size(fields);
  if (fsize == 0) return Py_None;

  // Connectivity
  E_Int ni, nj, nk;
  K_FLD::FldArrayF* f; K_FLD::FldArrayI* cn;
  char* varString; char* eltType;
  ret = K_ARRAY::getFromArray3(arr, varString, f, ni, nj, nk, cn, eltType);
  
  if (ret <= 0) 
  {
    PyErr_SetString(PyExc_TypeError, "computeGradLSQ: bad mesh.");
    return NULL;
  }

  if (ret == 1) 
  { 
    PyErr_SetString(PyExc_TypeError, "computeGradLSQ: only for NGons."); 
    RELEASESHAREDS(arr, f);
    return NULL; 
  }

  // Parent elements
  E_Int nfaces = cn->getNFaces();
  E_Int *PE = NULL;
  E_Int size = -1;
  ret = K_NUMPY::getFromNumpyArray(pe, PE, size);
  if (ret != 1 || size != 2*nfaces) 
  {
    RELEASESHAREDU(arr, f, cn);
    PyErr_SetString(PyExc_ValueError, "computeGradLSQ: bad parent elements array.");
    return NULL;
  }
  E_Int *owner = PE;
  E_Int *neigh = PE + nfaces;

  // Cell centers
  E_Int ncells = cn->getNElts();
  E_Float *CX, *CY, *CZ;
  CX = CY = CZ = NULL;
  ret = K_NUMPY::getFromNumpyArray(cx, CX, size);
  ret &= K_NUMPY::getFromNumpyArray(cy, CY, size);
  ret &= K_NUMPY::getFromNumpyArray(cz, CZ, size);
  if (ret != 1 || size != ncells) 
  {
    RELEASESHAREDU(arr, f, cn);
    PyErr_SetString(PyExc_ValueError, "computeGradLSQ: bad cell centers array.");
    return NULL;
  }

  // Face centers
  E_Float *FC = NULL;
  ret = K_NUMPY::getFromNumpyArray(fcenters, FC, size);
  if (ret != 1 || size != 3*nfaces) 
  {
    RELEASESHAREDU(arr, f, cn);
    PyErr_SetString(PyExc_ValueError, "computeGradLSQ: bad face centers array.");
    return NULL;
  }

  // Fields
  std::vector<E_Float *> Fields(fsize);
  for (E_Int i = 0; i < fsize; i++) 
  {
    PyObject *field = PyList_GetItem(fields, i);
    E_Int size = -1;
    ret = K_NUMPY::getFromNumpyArray(field, Fields[i], size);
    assert(size == ncells && ret == 1);
  }

  // Parse pointlists and rfields
  std::vector<E_Int *> pfaces;
  std::vector<std::vector<E_Float *>> rflds;
  std::vector<E_Int> npfaces;
  E_Int parRun = parse_pointlists_and_rfields(ptlists, rfields, pfaces, npfaces, rflds);

  // Make lsq A matrices 
  E_Int sizeNFace = cn->getSizeNFace();
  E_Float *lsqG = (E_Float *)malloc(3*sizeNFace * sizeof(E_Float));
  E_Int *count_neis = (E_Int *)malloc(ncells * sizeof(E_Int));
  make_A_grad_matrices(*cn, owner, neigh, count_neis, FC, CX, CY, CZ, lsqG);

  if (parRun) 
  {
    correct_A_grad_matrices(*cn, owner, count_neis, CX, CY, CZ, pfaces, npfaces,
      rflds, lsqG);
  }

  // Deduce tAA matrices
  E_Float *lsqGG = (E_Float *)malloc(9*ncells * sizeof(E_Float));
  deduce_lsq_grad_matrices(*cn, count_neis, lsqG, lsqGG);

  E_Float *b = (E_Float *)malloc(sizeNFace * sizeof(E_Float));

  PyObject *out = PyList_New(0);
  npy_intp dims[2];
  dims[0] = ncells;
  dims[1] = 1;

  E_Int tsize = parRun ? fsize-3 : fsize;
  for (E_Int i = 0; i < tsize; i++) 
  {
    const auto &Field = Fields[i];

    PyArrayObject *Gx = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    PyArrayObject *Gy = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    PyArrayObject *Gz = (PyArrayObject *)PyArray_SimpleNew(1, dims, NPY_DOUBLE);

    E_Float *gx = (E_Float *)PyArray_DATA(Gx);
    E_Float *gy = (E_Float *)PyArray_DATA(Gy);
    E_Float *gz = (E_Float *)PyArray_DATA(Gz);

    // Make b vectors
    make_RHS_vector(*cn, count_neis, owner, neigh, Field, b);

    if (parRun) 
    {
      correct_RHS_vector(*cn, count_neis, owner, Field, pfaces, npfaces, rflds,
        i, b);
    }

    // Solve for the gradient
    make_gradient(*cn, count_neis, b, lsqG, lsqGG, gx, gy, gz);

    PyObject *Grad = PyList_New(0);
    PyList_Append(Grad, (PyObject *)Gx);
    PyList_Append(Grad, (PyObject *)Gy);
    PyList_Append(Grad, (PyObject *)Gz);
    Py_DECREF(Gx);
    Py_DECREF(Gy);
    Py_DECREF(Gz);

    PyList_Append(out, Grad);
    Py_DECREF(Grad);
  }

  free(lsqG);
  free(lsqGG);
  free(count_neis);
  free(b);

  RELEASESHAREDU(arr, f, cn);
  Py_DECREF(fcenters);
  Py_DECREF(cx);
  Py_DECREF(cy);
  Py_DECREF(cz);
  Py_DECREF(pe);
  
  return out;
}

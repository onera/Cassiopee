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
#include "Array/Array.h"
#include "String/kstring.h"
#include <stdio.h>
#include <string.h>

using namespace K_FLD;

//=============================================================================
/* Build an empty structured array3
   IN: nfld: nbre de champs
   IN: varString: variables string
   IN: ni,nj,nk: number of points in field
   IN: api (1: array, 2: array2, 3: array3)
   OUT: PyObject created. */
//=============================================================================
PyObject* K_ARRAY::buildArray3(E_Int nfld, const char* varString, 
                               E_Int ni, E_Int nj, E_Int nk, E_Int api)
{
  PyObject* tpl;
  IMPORTNUMPY;

  if (api == 1) // Array1
  {
    npy_intp dim[2];
    dim[1] = ni*nj*nk;
    dim[0] = nfld;
    PyArrayObject* a = (PyArrayObject*)PyArray_SimpleNew(2, dim, NPY_DOUBLE);
    tpl = Py_BuildValue("[sOlll]", varString, a, (long)ni, (long)nj, (long)nk);
    Py_DECREF(a);
  }
  else // Array2 ou Array3
  {
    npy_intp dim[3];
    dim[0] = ni; dim[1] = nj; dim[2] = nk;
    PyObject* rake = PyList_New(0);
    for (E_Int n=0; n < nfld; n++)
    {
      PyArrayObject* a = (PyArrayObject*)PyArray_EMPTY(3, dim, NPY_DOUBLE, 1);
      PyList_Append(rake, (PyObject*)a); Py_DECREF(a);
    }
    tpl = Py_BuildValue("[sOlll]", varString, rake, (long)ni, (long)nj, (long)nk);
    Py_DECREF(rake);
  }
  return tpl;
}

//=============================================================================
/* Build an empty NGON array 
   IN: nfld: number of fields
   IN: varString: variable string
   IN: npts: number of vertex
   IN: nelts: number total of elements
   IN: nfaces: number total of faces
   IN: etString: NGON ou NGON*
   IN: center: set to true if field is localised in the centers of
   elements, otherwise let it to false.
   IN: sizeNGon, sizeNFace, nface: connectivity size.
   if sizeNFace == -1, NFACE is not created
   IN: ngonType=1 ou 2 (CGNSv3), ngonType=3 (CGNSv4)
   IN: api=1 (array1, ngonType=1), api=2 (array2, ngonType=2 ou 3), 
   api=3 (array3, ngonType=2 ou 3)
   OUT: PyObject created. */
//=============================================================================
// build empty pour les NGONS
PyObject* K_ARRAY::buildArray3(E_Int nfld, const char* varString,
                               E_Int npts, E_Int nelts, E_Int nfaces,
                               const char* etString,
                               E_Int sizeNGon, E_Int sizeNFace, E_Int ngonType,  
                               E_Bool center, E_Int api)
{
  npy_intp dim[2];
  PyObject* a; PyObject* ac; PyObject* tpl;
  char eltType[12];

  // taille de f
  E_Int fSize = (center) ? nelts : npts;

  IMPORTNUMPY;

  // Corrige eventuellement etString si contradictoire avec center
  if (center) K_ARRAY::starVarString(etString, eltType);
  else K_ARRAY::unstarVarString(etString, eltType);

  // Build array of fields
  if (api == 1) // Array1
  { 
    dim[1] = fSize; dim[0] = nfld;
    a = PyArray_SimpleNew(2, dim, NPY_DOUBLE);
  }
  else // Array2 or Array3
  {
    dim[0] = fSize;
    a = PyList_New(0);
    for (E_Int n=0; n < nfld; n++)
    {
      PyArrayObject* ar = (PyArrayObject*)PyArray_EMPTY(1, dim, NPY_DOUBLE, 1);
      PyList_Append(a, (PyObject*)ar); Py_DECREF(ar);
    }
  } 

  // Build array for connectivity
  if (api == 1) // Array 1
  {
    dim[1] = 4+sizeNGon+sizeNFace; dim[0] = 1;
    ac = PyArray_SimpleNew(2, dim, E_NPY_INT);
    E_Int* data = (E_Int*)PyArray_DATA((PyArrayObject*)ac);
    data[0] = nfaces;
    data[1] = sizeNGon;
    data[sizeNGon+2] = nelts;
    data[sizeNGon+3] = sizeNFace;
  }
  else if (ngonType == 2) // Array2/3 - NGonv3 + indir
  {
    ac = PyList_New(0);
    // ngons - NGON - sizeNGon
    dim[0] = sizeNGon;
    PyObject* ar = PyArray_EMPTY(1, dim, E_NPY_INT, 0);
    PyList_Append(ac, (PyObject*)ar); Py_DECREF(ar);
    // ngons - NFACE - sizeNFace
    dim[0] = sizeNFace;
    ar = PyArray_EMPTY(1, dim, E_NPY_INT, 0);
    PyList_Append(ac, (PyObject*)ar); Py_DECREF(ar);
    // ngons - indPG - nfaces
    dim[0] = nfaces;
    ar = PyArray_EMPTY(1, dim, E_NPY_INT, 0);
    PyList_Append(ac, (PyObject*)ar); Py_DECREF(ar);
    // ngons - indPH - nelts
    dim[0] = nelts;
    ar = PyArray_EMPTY(1, dim, E_NPY_INT, 0);
    PyList_Append(ac, (PyObject*)ar); Py_DECREF(ar);
    // Eventuellement PE - 2*nfaces
    //dim[0] = nfaces; dim[1] = 2;
    //PyObject* ar = PyArray_EMPTY(2, dim, E_NPY_INT, 0);
    //PyList_Append(ac, (PyObject*)ar); Py_DECREF(ar);
  }
  else if (ngonType == 3) // array3 - NGONv4
  {
    ac = PyList_New(0);
    // NGON - sizeNGon
    dim[0] = sizeNGon;
    PyObject* ar = PyArray_EMPTY(1, dim, E_NPY_INT, 0);
    PyList_Append(ac, (PyObject*)ar); Py_DECREF(ar);
    // NFACE - sizeNFace
    dim[0] = sizeNFace;
    ar = PyArray_EMPTY(1, dim, E_NPY_INT, 0);
    PyList_Append(ac, (PyObject*)ar); Py_DECREF(ar);
    // NGON - StartOffset
    dim[0] = nfaces+1;
    ar = PyArray_EMPTY(1, dim, E_NPY_INT, 0);
    E_Int* pt = (E_Int*)PyArray_DATA((PyArrayObject*)ar);
    pt[nfaces] = sizeNGon; 
    PyList_Append(ac, (PyObject*)ar); Py_DECREF(ar);
    // NFACE - startOffset
    dim[0] = nelts+1;
    ar = PyArray_EMPTY(1, dim, E_NPY_INT, 0);
    pt = (E_Int*)PyArray_DATA((PyArrayObject*)ar);
    pt[nelts] = sizeNFace;
    PyList_Append(ac, (PyObject*)ar); Py_DECREF(ar);
    // Eventuellement PE - 2*nfaces
    //dim[0] = nfaces; dim[1] = 2;
    //PyObject* ar = PyArray_EMPTY(2, dim, E_NPY_INT, 0);
    //PyList_Append(ac, (PyObject*)ar); Py_DECREF(ar);
  }
  else
  {
    printf("Warning: buildArray3: invalid api/ngonType. Array not built.\n");
    return NULL;
  }
  tpl = Py_BuildValue("[sOOs]", varString, a, ac, eltType);
  Py_DECREF(a); Py_DECREF(ac);

  return tpl;
}

//=============================================================================
/* Build an NGON array from f unstructured and an empty connectivity
   IN: f: fields
   IN: varString: variable string
   IN: nelts: number total of elements
   IN: nfaces: number total of faces
   IN: etString: NGON ou NGON*
   IN: center: set to true if field is localised in the centers of
   elements, otherwise let it to false.
   IN: sizeNGon, sizeNFace, nface: connectivity size.
   if sizeNFace == -1, NFACE is not created
   IN: ngonType=1 ou 2 (CGNSv3), ngonType=3 (CGNSv4)
   IN: api=1 (array1, ngonType=1), api=2 (array2, ngonType=2 ou 3), 
   api=3 (array3, ngonType=2 ou 3)
   OUT: PyObject created. */
//=============================================================================
PyObject* K_ARRAY::buildArray3(FldArrayF& f, const char* varString,
                               E_Int nelts, E_Int nfaces,
                               const char* etString,
                               E_Int sizeNGon, E_Int sizeNFace, E_Int ngonType,  
                               E_Bool center, E_Int api)
{
  if (K_STRING::cmp(etString, 4, "NGON") != 0)
  {
    printf("Warning: buildArray3: element type is not NGON.");
    return NULL;
  }
  
  if (api == -1) api = f.getApi();
  E_Int nfld = f.getNfld();
  E_Int fSize = f.getSize();

  // Corrige eventuellement etString si contradictoire avec center
  char eltType[256];
  if (center) K_ARRAY::starVarString(etString, eltType);
  else K_ARRAY::unstarVarString(etString, eltType);

  PyObject* tpl = K_ARRAY::buildArray3(
    nfld, varString, fSize, nelts, nfaces, 
    eltType, sizeNGon, sizeNFace, ngonType, center, api
  );
  FldArrayF* f2;
  K_ARRAY::getFromArray3(tpl, f2);

  #pragma omp parallel
  {
    for (E_Int n = 1; n <= nfld; n++)
    {
      E_Float* fp = f.begin(n);
      E_Float* f2p = f2->begin(n);
      #pragma omp for nowait
      for (E_Int i = 0; i < fSize; i++) f2p[i] = fp[i];
    }
  }

  RELEASESHAREDS(tpl, f2);
  return tpl;
}

//=============================================================================
/* Build an empty BE array 
   IN: nfld: number of fields
   IN: varString: variable string
   IN: npts: number of vertex
   IN: nelt: number of elements
   IN: etString: "TRI" ou avec *
   IN: center: set to true if field is localised in the centers of
   elements, otherwise let it to false.
   IN: api=1 (array1), api=2 ou 3 (array2 ou 3)
   OUT: PyObject created. */
//=============================================================================
// build empty pour les single Element (BE)
PyObject* K_ARRAY::buildArray3(E_Int nfld, const char* varString,
                               E_Int npts,
                               E_Int nelts,
                               const char* etString,
                               E_Bool center, E_Int api)
{
  npy_intp dim[2];
  PyObject* a; PyObject* ac; PyObject* tpl;
  char eltType[256];
  // Corrige eventuellement etString si contradictoire avec center
  if (center) K_ARRAY::starVarString(etString, eltType);
  else K_ARRAY::unstarVarString(etString, eltType);

  // taille de f
  E_Int fSize;
  if (center) fSize = nelts;
  else fSize = npts;

  IMPORTNUMPY;

  // Build array of fields
  if (api == 1) // Array1
  { 
    dim[1] = fSize; dim[0] = nfld;
    a = PyArray_SimpleNew(2, dim, NPY_DOUBLE);
  }
  else // Array2 or Array3
  {
    dim[0] = fSize;
    a = PyList_New(0);
    for (E_Int n=0; n < nfld; n++)
    {
      PyArrayObject* ar = (PyArrayObject*)PyArray_EMPTY(1, dim, NPY_DOUBLE, 1);
      PyList_Append(a, (PyObject*)ar); Py_DECREF(ar);
    }
  } 

  // Connectivite
  if (api == 1) // Array1
  {
    E_Int cSize = nelts;
    char st[256]; E_Int dummy; E_Int nvpe;
    eltString2TypeId(eltType, st, nvpe, dummy, dummy);
    dim[1] = cSize; dim[0] = nvpe;
    ac = PyArray_SimpleNew(2, dim, E_NPY_INT);
  }
  else if (api == 2 || api == 3) // Array2 ou 3
  {
    E_Int cSize = nelts;
    char st[256]; E_Int dummy; E_Int nvpe;
    eltString2TypeId(eltType, st, nvpe, dummy, dummy);
    ac = PyList_New(0);
    dim[0] = cSize; dim[1] = nvpe;
    PyObject* ar = PyArray_EMPTY(2, dim, E_NPY_INT, 0);
    PyList_Append(ac, (PyObject*)ar); Py_DECREF(ar);
  }
  else
  {
    PyErr_SetString(PyExc_TypeError, "buildArray: unkown api.");
    return NULL;
  }

  tpl = Py_BuildValue("[sOOs]", varString, a, ac, eltType);
  Py_DECREF(a); Py_DECREF(ac);

  return tpl;
}

//=============================================================================
/* Build an empty ME array 
   IN: nfld: number of fields
   IN: varString: variable string
   IN: npts: number of vertex
   IN: nepc: number of elements for each connect
   IN: etString: "TRI,QUAD" ou avec *
   IN: center: set to true if field is localised in the centers of
   elements, otherwise let it to false.
   IN: api=1 (array1, single connect), api=2 (array2, single connect), 
   api=3 (array3, all connects)

   OUT: PyObject created. */
//=============================================================================
// build empty pour les Multiple Elements (ME)
PyObject* K_ARRAY::buildArray3(E_Int nfld, const char* varString,
                               E_Int npts,
                               std::vector<E_Int>& nepc,
                               const char* etString,
                               E_Bool center, E_Int api)
{
  npy_intp dim[2];
  PyObject* a; PyObject* ac; PyObject* tpl;
  char eltType[256];
  // Corrige eventuellement etString si contradictoire avec center
  if (center) K_ARRAY::starVarString(etString, eltType);
  else K_ARRAY::unstarVarString(etString, eltType);

  // taille de f
  E_Int nelt = 0;
  for (size_t i = 0; i < nepc.size(); i++) nelt += nepc[i];
  E_Int fSize;
  if (center) fSize = nelt;
  else fSize = npts;

  IMPORTNUMPY;

  // Build array of fields
  if (api == 1) // Array1
  { 
    dim[1] = fSize; dim[0] = nfld;
    a = PyArray_SimpleNew(2, dim, NPY_DOUBLE);
  }
  else // Array2 or Array3
  {
    dim[0] = fSize;
    a = PyList_New(0);
    for (E_Int n=0; n < nfld; n++)
    {
      PyArrayObject* ar = (PyArrayObject*)PyArray_EMPTY(1, dim, NPY_DOUBLE, 1);
      PyList_Append(a, (PyObject*)ar); Py_DECREF(ar);
    }
  } 

  // Connectivite
  if (api == 1) // Array1 - force single connect
  {
    E_Int cSize = nelt;
    char st[256]; E_Int dummy; E_Int nvpe;
    eltString2TypeId(eltType, st, nvpe, dummy, dummy);
    dim[1] = cSize; dim[0] = nvpe;
    ac = PyArray_SimpleNew(2, dim, E_NPY_INT);
  }
  else if (api == 2) // Array2 - force single connect
  {
    E_Int cSize = nelt;
    char st[256]; E_Int dummy; E_Int nvpe;
    eltString2TypeId(eltType, st, nvpe, dummy, dummy);
    ac = PyList_New(0);
    dim[0] = cSize; dim[1] = nvpe;
    PyObject* ar = PyArray_EMPTY(2, dim, E_NPY_INT, 0);
    PyList_Append(ac, (PyObject*)ar); Py_DECREF(ar);
  }
  else // Array3
  {
    std::vector<char*> eltTypes;
    K_ARRAY::extractVars(eltType, eltTypes);
    char st[256]; E_Int dummy; E_Int nvpe;
    ac = PyList_New(0);
    for (size_t i = 0; i < eltTypes.size(); i++)
    {
      E_Int cSize = nepc[i];
      eltString2TypeId(eltTypes[i], st, nvpe, dummy, dummy);
      dim[0] = cSize; dim[1] = nvpe;
      PyObject* ar = PyArray_EMPTY(2, dim, E_NPY_INT, 0);
      PyList_Append(ac, (PyObject*)ar); Py_DECREF(ar);
    }
    for (size_t i = 0; i < eltTypes.size(); i++) delete [] eltTypes[i];
  }

  tpl = Py_BuildValue("[sOOs]", varString, a, ac, eltType);
  Py_DECREF(a); Py_DECREF(ac);

  return tpl;
}

// Build an array identical to cn in size (unstructured only)
// but with nfld vars and size vertices. Center and api can be changed.
// if copyConnect is true, performs copy on cn
PyObject* K_ARRAY::buildArray3(E_Int nfld,
                               const char* varString,
                               E_Int npts,
                               FldArrayI& cn,
                               const char* eltType,
                               E_Int center, E_Int api,
                               E_Bool copyConnect)
{
  PyObject* tpl = NULL;
  char eltType2[256];
  // Corrige eventuellement eltType si contradictoire avec center
  if (center > 0) K_ARRAY::starVarString(eltType, eltType2);
  else if (center == 0) K_ARRAY::unstarVarString(eltType, eltType2);
  else // find center from eltType
  { 
    center = 0;
    if (strchr(eltType, '*') != NULL) center = 1;
    strcpy(eltType2, eltType);
  }
  if (K_STRING::cmp(eltType2, "NGON") == 0 || K_STRING::cmp(eltType2, "NGON*") == 0)
  {
    E_Int ngonType = cn.getNGonType();
    E_Int nelts = cn.getNElts();
    E_Int nfaces = cn.getNFaces();
    E_Int sizeNGon = cn.getSizeNGon();
    E_Int sizeNFace = cn.getSizeNFace();
    tpl = K_ARRAY::buildArray3(nfld, varString, npts, nelts, nfaces, 
                               eltType2, sizeNGon, sizeNFace, ngonType,
                               center, api);

    if (copyConnect)
    {
      FldArrayF* f2; FldArrayI* cn2;
      K_ARRAY::getFromArray3(tpl, f2, cn2);
      E_Int *ngon = cn.getNGon(), *ngon2 = cn2->getNGon();
      E_Int *nface = cn.getNFace(), *nface2 = cn2->getNFace();
      E_Int *indPG = NULL, *indPG2 = NULL;
      E_Int *indPH = NULL, *indPH2 = NULL;
      E_Int dim2f = 0, dim2e = 0;
      if (api > 1)
      {
        indPG = cn.getIndPG(); indPG2 = cn2->getIndPG();
        indPH = cn.getIndPH(); indPH2 = cn2->getIndPH();
        if (ngonType == 2) { dim2f = nfaces; dim2e = nelts; }
        else { dim2f = nfaces+1; dim2e = nelts+1; }
      }

      #pragma omp parallel
      {
        #pragma omp for nowait
        for (E_Int i = 0; i < sizeNGon; i++) ngon2[i] = ngon[i];
        #pragma omp for nowait
        for (E_Int i = 0; i < sizeNFace; i++) nface2[i] = nface[i];
        if (api > 1)
        {
          #pragma omp for nowait
          for (E_Int i = 0; i < dim2f; i++) indPG2[i] = indPG[i];
          #pragma omp for
          for (E_Int i = 0; i < dim2e; i++) indPH2[i] = indPH[i];
        }
      }
      RELEASESHAREDU(tpl, f2, cn2);
    }
  }
  else
  {
    E_Int nc = cn.getNConnect();
    std::vector<E_Int> nepc(nc);
    for (E_Int ic = 0; ic < nc; ic++)
    {
      FldArrayI& cm = *(cn.getConnect(ic));
      nepc[ic] = cm.getSize();
    }
    tpl = K_ARRAY::buildArray3(nfld, varString, npts,
                               nepc, eltType2, center, api);

    if (copyConnect)
    {
      FldArrayF* f2; FldArrayI* cn2;  
      K_ARRAY::getFromArray3(tpl, f2, cn2);
      #pragma omp parallel
      {
        for (E_Int ic = 0; ic < nc; ic++)
        { 
          FldArrayI& cm = *(cn.getConnect(ic));
          FldArrayI& cm2 = *(cn2->getConnect(ic));
          #pragma omp for collapse(2) nowait
          for (E_Int i = 0; i < cm.getSize(); i++)
            for (E_Int j = 1; j <= cm.getNfld(); j++)
              cm2(i,j) = cm(i,j);
        }
      }
      RELEASESHAREDU(tpl, f2, cn2);
    }
  }
  return tpl;
}

// Copy from f and cn unstructured
PyObject* K_ARRAY::buildArray3(FldArrayF& f, const char* varString,
                               FldArrayI& cn,  const char* eltType, E_Int api)
{
  if (api == -1) api = f.getApi();
  E_Int nfld = f.getNfld(); E_Int npts = f.getSize();
  if (K_STRING::cmp(eltType, 4, "NGON") == 0)
  {
    E_Int ngonType = cn.getNGonType();
    E_Int nelts = cn.getNElts();
    E_Int nfaces = cn.getNFaces();
    E_Int sizeNGon = cn.getSizeNGon();
    E_Int sizeNFace = cn.getSizeNFace();
    E_Int* ngon = cn.getNGon();
    E_Int* nface = cn.getNFace();
    E_Int* indPG = NULL; E_Int* indPH = NULL;
    if (api > 1) { indPG = cn.getIndPG(); indPH = cn.getIndPH(); }

    E_Bool center = false;
    if (strchr(eltType, '*') != NULL) center = true;
    
    PyObject* tpl = K_ARRAY::buildArray3(
      nfld, varString, npts, nelts, nfaces, 
      eltType, sizeNGon, sizeNFace, ngonType, center, api
    );
    FldArrayF* f2; FldArrayI* cn2;
    K_ARRAY::getFromArray3(tpl, f2, cn2);

    E_Int* ngon2 = cn2->getNGon();
    E_Int* nface2 = cn2->getNFace();
    E_Int* indPG2 = NULL; E_Int* indPH2 = NULL;
    if (api > 1) { indPG2 = cn2->getIndPG(); indPH2 = cn2->getIndPH(); }

    E_Int dim1 = (center) ? nelts : npts;
    E_Int dim2f = (ngonType == 2) ? nfaces : nfaces+1;
    E_Int dim2e = (ngonType == 2) ? nelts : nelts+1;

    #pragma omp parallel
    {
      for (E_Int n = 1; n <= nfld; n++)
      {
        E_Float* fp = f.begin(n);
        E_Float* f2p = f2->begin(n);
        #pragma omp for nowait
        for (E_Int i = 0; i < dim1; i++) f2p[i] = fp[i];
      }
      
      #pragma omp for nowait
      for (E_Int i = 0; i < sizeNGon; i++) ngon2[i] = ngon[i];
      #pragma omp for nowait
      for (E_Int i = 0; i < sizeNFace; i++) nface2[i] = nface[i];

      if (api > 1)
      {
        #pragma omp for nowait
        for (E_Int i = 0; i < dim2f; i++) indPG2[i] = indPG[i];
        #pragma omp for
        for (E_Int i = 0; i < dim2e; i++) indPH2[i] = indPH[i];
      }
    }
    RELEASESHAREDU(tpl, f2, cn2);
    return tpl;
  }
  else // BE/ME
  {
    E_Int nc = cn.getNConnect();
    E_Bool center = false;
    if (strchr(eltType, '*') != NULL) center = true;
    std::vector<E_Int> nepc(nc);
    E_Int ntotElts = 0;
    for (E_Int i = 0; i < nc; i++)
    {
      FldArrayI& cm = *(cn.getConnect(i));
      nepc[i] = cm.getSize(); 
      ntotElts += cm.getSize();
    }
    PyObject* tpl = K_ARRAY::buildArray3(nfld, varString, npts,
                                         nepc, eltType, center, api);
    FldArrayF* f2; FldArrayI* cn2;  
    K_ARRAY::getFromArray3(tpl, f2, cn2);  
    E_Int dim = (center) ? ntotElts : npts;

    #pragma omp parallel
    {
      // copie des champs
      for (E_Int n = 1; n <= nfld; n++)
      {
        E_Float* fp = f.begin(n);
        E_Float* f2p = f2->begin(n);
        #pragma omp for nowait
        for (E_Int i = 0; i < dim; i++) f2p[i] = fp[i];
      }
      
      for (E_Int ic = 0; ic < nc; ic++)
      { 
        FldArrayI& cm = *(cn.getConnect(ic));
        FldArrayI& cm2 = *(cn2->getConnect(ic));
        #pragma omp for collapse(2) nowait
        for (E_Int i = 0; i < cm.getSize(); i++)
        for (E_Int j = 1; j <= cm.getNfld(); j++)
          cm2(i,j) = cm(i,j);
      }
    }
    RELEASESHAREDU(tpl, f2, cn2);
    return tpl;
  }
}

// Build a copy array from f structured
// if api=-1, keep f api
PyObject* K_ARRAY::buildArray3(FldArrayF& f, const char* varString,
                               E_Int ni, E_Int nj, E_Int nk, E_Int api)
{
  if (api == -1) api = f.getApi(); // copie l'api de f
  E_Int nfld = f.getNfld(); E_Int npts = f.getSize();
  PyObject* tpl = K_ARRAY::buildArray3(nfld, varString, ni, nj, nk, api);
  FldArrayF* f2;
  K_ARRAY::getFromArray3(tpl, f2);
  #pragma omp parallel
  {
    for (E_Int n = 1; n <= nfld; n++)
    {
      E_Float* fp = f.begin(n);
      E_Float* f2p = f2->begin(n);
      #pragma omp for nowait
      for (E_Int i = 0; i < npts; i++) f2p[i] = fp[i];
    }
  }
  RELEASESHAREDS(tpl, f2);
  return tpl;
}

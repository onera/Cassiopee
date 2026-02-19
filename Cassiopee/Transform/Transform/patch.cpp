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

// Global operations on several meshes (patch)

# include "transform.h"

using namespace std;
using namespace K_FLD;

// ============================================================================
/* Patch array1 into array2 at position (ii,jj,kk) of array2 */
// ============================================================================
PyObject* K_TRANSFORM::patch(PyObject* self, PyObject* args)
{
  PyObject* array1; PyObject* array2;
  E_Int ii, jj, kk;
  if (!PYPARSETUPLE_(args, OO_ TIII_,
                    &array1, &array2, &ii, &jj, &kk))
  {
    return NULL;
  }
  // Check array
  E_Int im1, jm1, km1, im2, jm2, km2;
  FldArrayF* f1; FldArrayF* f2;
  FldArrayI* cn1; FldArrayI* cn2;
  char *varString1; char *varString2;
  char* eltType1; char* eltType2;
  vector<E_Int> pos1; vector<E_Int> pos2;
  E_Int res1 =
    K_ARRAY::getFromArray3(array1, varString1, f1, im1, jm1, km1,
                           cn1, eltType1);
  E_Int res2 =
    K_ARRAY::getFromArray3(array2, varString2, f2, im2, jm2, km2,
                           cn2, eltType2);
  E_Int ind, ind2;

  if (res1 == 1 && res2 == 1)
  {
    char* varString = new char [strlen(varString1) + strlen(varString2) + 4];
    E_Int res0 =
      K_ARRAY::getPosition(varString1, varString2, pos1, pos2, varString);

    if (res0 == -1)
    {
      RELEASESHAREDS(array1, f1);
      RELEASESHAREDS(array2, f2);
      PyErr_SetString(PyExc_ValueError,
                      "patch: one array is empty.");
      return NULL;
    }
    else if (res0 == 0)
    {
      printf("Warning: patch: arrays have different variables.");
      printf("Only common variables are kept.\n");
    }

    if (im1 > im2 || jm1 > jm2 || km1 > km2)
    {
      RELEASESHAREDS(array1, f1);
      RELEASESHAREDS(array2, f2);
      PyErr_SetString(PyExc_ValueError,
                      "patch: patched element dimensions must be smaller than the element on which it is applied.");
      return NULL;
    }
    if (ii < 1 || ii > im2)
    {
      RELEASESHAREDS(array1, f1);
      RELEASESHAREDS(array2, f2);
      PyErr_SetString(PyExc_ValueError,
                      "patch: position is not in mesh.");
      return NULL;
    }
    if (jj < 1 || jj > jm2)
    {
      RELEASESHAREDS(array1, f1);
      RELEASESHAREDS(array2, f2);
      PyErr_SetString(PyExc_ValueError,
                      "patch: position is not in mesh.");
      return NULL;
    }
    if (kk < 1 || kk > km2)
    {
      RELEASESHAREDS(array1, f1);
      RELEASESHAREDS(array2, f2);
      PyErr_SetString(PyExc_ValueError,
                      "patch: position is not in mesh.");
      return NULL;
    }
    if (im1 != 1 && im2 != 1)
    {
      if (ii + im1 > im2+1)
      {
        RELEASESHAREDS(array1, f1);
        RELEASESHAREDS(array2, f2);
        PyErr_SetString(PyExc_ValueError,
                        "patch: i index is too big.");
        return NULL;
      }
    }
    if (jm1 != 1 && jm2 != 1)
    {
      if (jj + jm1 > jm2+1)
      {
        RELEASESHAREDS(array1, f1);
        RELEASESHAREDS(array2, f2);
        PyErr_SetString(PyExc_ValueError,
                        "patch: j index is too big.");
        return NULL;
      }
    }
    if (km1 != 1 && km2 != 1)
    {
      if (kk + km1 > km2+1)
      {
        RELEASESHAREDS(array1, f1);
        RELEASESHAREDS(array2, f2);
        PyErr_SetString(PyExc_ValueError,
                        "patch: k index is too big.");
        return NULL;
      }
    }
    E_Int api = f2->getApi();
    E_Int nfld = pos1.size();
    vector<E_Float*> fp1(nfld);
    // pointeur sur les champs de array1
    for (E_Int p = 0; p < nfld; p++) { fp1[p] = f1->begin(pos1[p]); }
    // pointeur sur les champs de array2
    vector<E_Float*> fp2(nfld);
    for (E_Int p = 0; p < nfld; p++) { fp2[p] = f2->begin(pos2[p]); }

    for (E_Int k = kk-1; k < kk+km1-1; k++)
      for (E_Int j = jj-1; j < jj+jm1-1; j++)
        for (E_Int i = ii-1; i < ii+im1-1; i++)
        {
          ind = i + j*im2 + k*im2*jm2;
          ind2 = (i-ii+1)+(j-jj+1)*im1+(k-kk+1)*im1*jm1;
          for (E_Int n = 0; n < nfld; n++)
            fp2[n][ind] = fp1[n][ind2];
        }

    // Build array
    PyObject* tpl = K_ARRAY::buildArray3(*f2, varString, im2, jm2, km2, api);
    delete [] varString;
    RELEASESHAREDS(array1, f1);
    RELEASESHAREDS(array2, f2);

    return tpl;
  }
  else if (res1 == 2 || res2 == 2)
  {
    RELEASESHAREDB(res1, array1, f1, cn1);
    RELEASESHAREDB(res2, array2, f2, cn2);
    PyErr_SetString(PyExc_TypeError,
                    "patch: can not be used on an unstructured array.");
    return NULL;
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "patch: invalid array.");
    return NULL;
  }
}

// ============================================================================
/* Patch array1 into array2 at nodes [n1, n2, n3 ...] of array2
   work for structured and unstructured meshes */
// ============================================================================
PyObject* K_TRANSFORM::patch2(PyObject* self, PyObject* args)
{
  PyObject* array1; PyObject* array2;
  PyObject* nodesIndices;
  if (!PYPARSETUPLE_(args, OOO_,
                      &array1, &array2, &nodesIndices))
  {
    return NULL;
  }

  // Check array
  E_Int im1, jm1, km1, im2, jm2, km2;
  FldArrayF* f1; // champ de array1
  FldArrayF* f2; // champ de array2
  FldArrayI* cn1; // connectivite de array1
  FldArrayI* cn2; // connectivite de array2
  char* varString1; // liste des variables contenues dans array1
  char* varString2; // liste des variables contenues dans array2
  char* eltType1; // type d'elements de array1
  char* eltType2; // type d'elements de array2
  E_Int res1 =
    K_ARRAY::getFromArray3(array1, varString1, f1, im1, jm1, km1,
                           cn1, eltType1);
  E_Int res2 =
    K_ARRAY::getFromArray3(array2, varString2, f2, im2, jm2, km2,
                           cn2, eltType2);

  if (res1 != 1 && res1 != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "patch: unknown type of array for array1.");
    return NULL;
  }
  if (res2 != 1 && res2 != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "patch: unknown type of array for array2.");
    return NULL;
  }

  E_Int nbnodes, nf; E_Int* n0;
  K_NUMPY::getFromNumpyArray(nodesIndices, n0, nbnodes, nf);

  if (nbnodes != f1->getSize())
  {
    PyErr_SetString(PyExc_ValueError,
                    "patch: array1 and nodes must have the same size.");
    return NULL;
  }
  E_Int ind2;

  // construction d'un nouveau varString a partir de varString1 et varString2
  vector<E_Int> pos1; vector<E_Int> pos2;
  char* varString = new char [strlen(varString1) + strlen(varString2) + 4];
  E_Int res0 =
    K_ARRAY::getPosition(varString1, varString2, pos1, pos2, varString);

  if (res0 == -1)
  {
    PyErr_SetString(PyExc_ValueError,
                    "patch: one array is empty.");
    return NULL;
  }
  else if (res0 == 0)
  {
    printf("Warning: patch: arrays have different variables.");
    printf("Only common variables are kept.\n");
  }

  E_Int nfld = pos1.size();
  E_Int api = f2->getApi();

  vector<E_Float*> fp1(nfld);
  // pointeur sur les champs de array1
  for (E_Int p = 0; p < nfld; p++) { fp1[p] = f1->begin(pos1[p]); }
  // pointeur sur les champs de array2
  vector<E_Float*> fp2(nfld);
  for (E_Int p = 0; p < nfld; p++) { fp2[p] = f2->begin(pos2[p]); }

  // boucle sur les indices a patcher
  for (E_Int p = 0; p < nfld; p++)
  {
    for (E_Int ind = 0; ind < nbnodes; ind++)
    {
      ind2 = n0[ind]-1;
      fp2[p][ind2] = fp1[p][ind];
    }
  }

  // Build array
  PyObject* tpl = NULL;
  if (res2 == 1)
  {
    tpl = K_ARRAY::buildArray3(*f2, varString, im2, jm2, km2, api);
  }
  else if (res2 == 2)
  {
    tpl = K_ARRAY::buildArray3(*f2, varString, *cn2, eltType2, api);
  }

  Py_DECREF(nodesIndices);
  RELEASESHAREDB(res1, array1, f1, cn1);
  RELEASESHAREDB(res2, array2, f2, cn2);

  delete [] varString;
  return tpl;
}

/*
    Copyright 2013-2025 Onera.

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

# include "transform.h"
# include "Connect/connect.h"
# include "Metis/metis.h"

using namespace K_FLD;
using namespace std;

//=============================================================================
// Split NGon using METIS
// Return a list of element numpys.
//==============================================================================
PyObject* K_TRANSFORM::splitNGon(PyObject* self, PyObject* args)
{
  PyObject* array; E_Int nparts;
  if (!PYPARSETUPLE_(args, O_ I_, &array, &nparts)) return NULL;
  // Check array
  E_Int ni, nj, nk, res;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  res = K_ARRAY::getFromArray3(array, varString,
                               f, ni, nj, nk, cn, eltType);

  if (res != 2)
  {
    if (res == 1) RELEASESHAREDS(array, f);
    PyErr_SetString(PyExc_TypeError,
                    "splitNGon: array is invalid.");
    return NULL;
  }
  if (strcmp(eltType, "TRI")   == 0 || strcmp(eltType, "QUAD") == 0 ||
      strcmp(eltType, "TETRA") == 0 || strcmp(eltType, "HEXA") == 0 ||
      strcmp(eltType, "PENTA") == 0 || strcmp(eltType, "BAR")  == 0 ||
      strcmp(eltType, "PYRA")  == 0 || strcmp(eltType, "NODE") == 0)
  { RELEASESHAREDU(array, f, cn); return array; }
  if (strcmp(eltType, "NGON") != 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "splitNGon: elt type must be NGON.");
    RELEASESHAREDU(array, f, cn); return NULL;
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString); posx++;
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString); posy++;
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString); posz++;

  // Construit le graph
  FldArrayI cFE;
  K_CONNECT::connectNG2FE(*cn, cFE);
  E_Int* cFE1 = cFE.begin(1);
  E_Int* cFE2 = cFE.begin(2);

  E_Int* nface = cn->getNFace();
  E_Int* indPH = cn->getIndPH();
  E_Int nelts = cn->getNElts();
  //printf("nelts=%d\n", nelts);

  E_Int nf;
  E_Int size = 0; // size of adj
  for (E_Int i = 0; i < nelts; i++)
  {
    cn->getElt(i, nf, nface, indPH);
    size += nf;
  }
  //printf("size = %d\n", size);

  E_Int e1, e2, indf;
  idx_t* adj1 = new idx_t [size];
  idx_t* adj = adj1;
  idx_t* xadj = new idx_t [nelts+1];
  size = 0;
  for (E_Int i = 0; i < nelts; i++)
  {
    xadj[i] = size;
    E_Int* face = cn->getElt(i, nf, nface, indPH);

    for (E_Int n = 0; n < nf; n++)
    {
      indf = face[n]-1;
      e1 = cFE1[indf];
      e2 = cFE2[indf];
      //printf("%d - %d %d\n",i+1, e1,e2);
      if (e1 > 0 && e1 != i+1) { adj[size] = e1-1; size++; }
      else if (e2 > 0 && e2 != i+1) { adj[size] = e2-1; size++; }
    }
  }
  xadj[nelts] = size;
  adj = adj1;
  //for (E_Int i = 0; i < nelts+1; i++) printf("%d ",xadj[i]); printf("\n\n");
  //for (E_Int i = 0; i < size; i++) printf("%d ",adj[i]);

  cFE.malloc(0);

  E_Int ncon = 1;
  E_Int objval = 0;
  idx_t* parts = new idx_t [nelts];
  //for (E_Int i = 0; i < nelts; i++) parts[i] = 0; // dbx
  //METIS_PartGraphRecursive(&nelts, &ncon, xadj, adj, NULL, NULL, NULL,
  //                         &nparts, NULL, NULL, NULL, &objval, parts);
  METIS_PartGraphKway(&nelts, &ncon, xadj, adj, NULL, NULL, NULL,
                      &nparts, NULL, NULL, NULL, &objval, parts);
  delete [] xadj; delete [] adj1;

  // output: sortie de la liste des elements pour chaque part
  E_Int p;
  E_Int* partSize = new E_Int [nparts]; // nbre d'elements dans chaque partition
  E_Int** partPtr = new E_Int* [nparts]; // ptr sur les elements

  for (E_Int i = 0; i < nparts; i++) partSize[i] = 0;
  for (E_Int i = 0; i < nelts; i++)
  {
    p = parts[i]; partSize[p] += 1;
  }
  //for (E_Int i = 0; i < nparts; i++) printf("Info: partSize=%d\n", partSize[i]);

  // output numpy of elements
  PyObject* tpl;
  PyObject* l = PyList_New(0);
  for (E_Int i = 0; i < nparts; i++)
  {
    tpl = K_NUMPY::buildNumpyArray(partSize[i], 1, 1, 1);
    partPtr[i] = K_NUMPY::getNumpyPtrI(tpl);
    PyList_Append(l, tpl); Py_DECREF(tpl);
  }

  for (E_Int i = 0; i < nparts; i++) partSize[i] = 0;
  for (E_Int i = 0; i < nelts; i++)
  {
    p = parts[i]; partPtr[p][partSize[p]] = i; partSize[p] += 1;
  }

  delete [] partSize; delete [] partPtr;
  delete [] parts;
  RELEASESHAREDU(array, f, cn);
  return l;
}

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

# include "transform.h"
# include "Connect/connect.h"
# include "Metis/metis.h"

using namespace K_FLD;
using namespace std;

//=============================================================================
// Split Elements using METIS
// Return a list of element numpys.
//==============================================================================
PyObject* K_TRANSFORM::splitElement(PyObject* self, PyObject* args)
{
  PyObject* array; E_Int nparts;
  if (!PYPARSETUPLE_(args, O_ I_, &array, &nparts))
  {
    return NULL;
  }
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
                    "splitElement: array is invalid.");
    return NULL;
  }
  if (strcmp(eltType, "NGON") == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "splitElement: array must be basic elements.");
    RELEASESHAREDU(array, f, cn); return NULL;
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString); posx++;
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString); posy++;
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString); posz++;

  // Construit le graph
  E_Int np = f->getSize();
  E_Int ne = cn->getSize();
  vector< vector<E_Int> > cEEN(ne);
  K_CONNECT::connectEV2EENbrs(eltType, np, *cn, cEEN);

  E_Int size = 0; // size of adj
  for (E_Int i = 0; i < ne; i++)
  {
    vector<E_Int>& voisins = cEEN[i];
    size += voisins.size();
  }
  //printf("size = " SF_D_ "\n", size);
  //printf("size of idx=%d, int=%d\n", sizeof(idx_t), sizeof(E_Int));

  idx_t* adj1 = new idx_t [size];
  idx_t* adj = adj1;
  idx_t* xadj = new idx_t [ne+1];
  size = 0;
  for (E_Int i = 0; i < ne; i++)
  {
    vector<E_Int>& voisins = cEEN[i];
    xadj[i] = size;

    for (size_t n = 0; n < voisins.size(); n++)
    {
      adj[size+n] = voisins[n];
    }
    size += voisins.size();
  }
  xadj[ne] = size;
  adj = adj1;

  E_Int ncon = 1;
  E_Int objval = 0; // retour

  //idx_t options[METIS_NOPTIONS];
  //METIS_SetDefaultOptions(options);
  //options[METIS_OPTION_CONTIG] = 1; // force contiguite
  //options[METIS_OPTION_OBJTYPE] = METIS_OBJTYPE_CUT; // METIS_OBJTYPE_VOL
  //options[METIS_OPTION_MINCONN] = 1; // force min connectivite externe

  idx_t* parts = new idx_t [ne];
  //for (E_Int i = 0; i < ne; i++) parts[i] = 0; // dbx
  //METIS_PartGraphRecursive(&ne, &ncon, xadj, adj, NULL, NULL, NULL,
  //                         &nparts, NULL, NULL, NULL, &objval, parts);
  METIS_PartGraphKway(&ne, &ncon, xadj, adj, NULL, NULL, NULL,
                      &nparts, NULL, NULL, NULL, &objval, parts);
  delete [] xadj; delete [] adj1;

  // output: sortie de la liste des elements pour chaque part
  E_Int p;
  E_Int* partSize = new E_Int [nparts]; // nbre d'elements dans chaque partition
  E_Int** partPtr = new E_Int* [nparts]; // ptr sur les elements

  for (E_Int i = 0; i < nparts; i++) partSize[i] = 0;
  for (E_Int i = 0; i < ne; i++)
  {
    p = parts[i]; partSize[p] += 1;
  }
  for (E_Int i = 0; i < nparts; i++) printf("partSize=" SF_D_ "\n", partSize[i]);

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
  for (E_Int i = 0; i < ne; i++)
  {
    p = parts[i]; partPtr[p][partSize[p]] = i; partSize[p] += 1;
  }

  delete [] partSize; delete [] partPtr;
  delete [] parts;
  RELEASESHAREDU(array, f, cn);
  return l;
}

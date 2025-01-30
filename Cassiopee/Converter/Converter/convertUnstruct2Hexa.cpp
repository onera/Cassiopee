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
#include "converter.h"

using namespace K_FUNC;
using namespace K_FLD;
using namespace std;

// ============================================================================
/* Convert unstructured array to a hexaedrical mesh */
// ============================================================================
PyObject* K_CONVERTER::convertUnstruct2Hexa(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PyArg_ParseTuple(args, "O", &array)) return NULL;

  // Check array
  E_Int nil, njl, nkl;
  FldArrayF* f; FldArrayI* cnl;
  char* varString; char* eltType; string eltType2;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f,
                                     nil, njl, nkl, cnl, eltType);
  if (res != 2) 
  {
    RELEASESHAREDS(array, f);
    PyErr_SetString(PyExc_TypeError, 
                    "convertUnstruct2Hexa: array must be unstructured.");
    return NULL;
  }

  // Acces universel sur BE/ME
  E_Int nc = cnl->getNConnect();
  // Acces universel aux eltTypes
  vector<char*> eltTypes;
  K_ARRAY::extractVars(eltType, eltTypes);

  E_Int npts = f->getSize(), api = f->getApi(), nfld = f->getNfld();
  E_Int nelts, loc = 0;
  vector<E_Int> neltsConn(nc);
  E_Boolean center = false;

  // Boucle sur toutes les connectivites pour determiner les types des
  // nouveaux elements afin de construire la nouvelle connectivite ME
  for (E_Int ic = 0; ic < nc; ic++)
  {
    FldArrayI& cm = *(cnl->getConnect(ic));
    char* eltTypConn = eltTypes[ic];
    if (ic > 0) eltType2.append(","); 
    
    if (strncmp(eltTypConn, "QUAD", 4) == 0 || strncmp(eltTypConn, "HEXA", 4) == 0 ||
        strncmp(eltTypConn, "BAR", 3) == 0)
    {
      eltType2.append(eltTypConn);
    }
    else if (strncmp(eltType, "TRI", 3) == 0) eltType2.append("QUAD");
    else if (strncmp(eltType, "PENTA", 5) == 0 || strncmp(eltType, "TETRA", 5) == 0)
    {
      eltType2.append("HEXA");
    }
    else
    {
      PyErr_SetString(PyExc_TypeError,
                      "convertUnstruct2Hexa: invalid element type.");
      return NULL;
    }

    if (strchr(eltTypConn,'*') != NULL) loc += 1;
    neltsConn[ic] = cm.getSize();
  }

  if (loc != 0 and loc != nc)
  {
    PyErr_SetString(PyExc_TypeError,
                    "convertUnstruct2Hexa: invalid element type - mix of "
                    "cell-centered and nodal fields in ME connectivity");
    return NULL;
  }
  else if (loc == nc) center = true;

  // Build empty ME connectivity
  PyObject* tpl = K_ARRAY::buildArray3(nfld, varString, npts, neltsConn,
                                       eltType2.c_str(), center, api);
  FldArrayF* f2; FldArrayI* cnl2;
  K_ARRAY::getFromArray3(tpl, f2, cnl2);

  // Boucle sur toutes les connectivites pour les remplir
  for (E_Int ic = 0; ic < nc; ic++)
  {
    FldArrayI& cm = *(cnl->getConnect(ic));
    FldArrayI& cm2 = *(cnl2->getConnect(ic));
    char* eltTypConn = eltTypes[ic];
    nelts = cm.getSize();

#pragma omp parallel default(shared) if (nelts > __MIN_SIZE_MEAN__)
    {
      if (strncmp(eltTypConn, "QUAD", 4) == 0)
      {
        // Copy existing connectivity
#pragma omp for
        for (E_Int et = 0; et < nelts; et++)
        {
          cm2(et,1) = cm(et,1);
          cm2(et,2) = cm(et,2);
          cm2(et,3) = cm(et,3);
          cm2(et,4) = cm(et,4);
        }
      }
      else if (strncmp(eltTypConn, "HEXA", 4) == 0)
      {
        // Copy existing connectivity
#pragma omp for
        for (E_Int et = 0; et < nelts; et++)
        {
          cm2(et,1) = cm(et,1);
          cm2(et,2) = cm(et,2);
          cm2(et,3) = cm(et,3);
          cm2(et,4) = cm(et,4);
          cm2(et,5) = cm(et,5);
          cm2(et,6) = cm(et,6);
        }
      }
      else if (strncmp(eltTypConn, "BAR", 3) == 0)
      {
        // Copy existing connectivity
#pragma omp for
        for (E_Int et = 0; et < nelts; et++)
        {
          cm2(et,1) = cm(et,1);
          cm2(et,2) = cm(et,2);
        }
      }
      else if (strncmp(eltType, "TRI", 3) == 0)
      {
#pragma omp for
        for (E_Int et = 0; et < nelts; et++)
        {
          cm2(et,1) = cm(et,1);
          cm2(et,2) = cm(et,2);
          cm2(et,3) = cm(et,3);
          cm2(et,4) = cm(et,3);
        }
      }
      else if (strncmp(eltType, "TETRA", 5) == 0)
      {
#pragma omp for
        for (E_Int et = 0; et < nelts; et++)
        {
          cm2(et,1) = cm(et,1);
          cm2(et,2) = cm(et,2);
          cm2(et,3) = cm(et,3);
          cm2(et,4) = cm(et,3);
          cm2(et,5) = cm(et,4);
          cm2(et,6) = cm(et,4);
          cm2(et,7) = cm(et,4);
          cm2(et,8) = cm(et,4);
        }
      }
      else if (strncmp(eltType, "PENTA", 5) == 0)
      {
#pragma omp for
        for (E_Int et = 0; et < nelts; et++)
        {
          cm2(et,1) = cm(et,1);
          cm2(et,2) = cm(et,2);
          cm2(et,3) = cm(et,3);
          cm2(et,4) = cm(et,3);
          cm2(et,5) = cm(et,4);
          cm2(et,6) = cm(et,5);
          cm2(et,7) = cm(et,6);
          cm2(et,8) = cm(et,6);
        }
      }
    }
  }

  for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];

#pragma omp parallel
  {
    // Copy fields to f2
    for (E_Int n = 1; n <= nfld; n++)
    {
      E_Float* fp = f->begin(n);
      E_Float* f2p = f2->begin(n);
#pragma omp for
      for (E_Int i = 0; i < npts; i++) f2p[i] = fp[i];
    }
  }
  
  RELEASESHAREDU(array, f, cnl);
  RELEASESHAREDU(tpl, f2, cnl2);
  return tpl;
}

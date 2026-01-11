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
#include "converter.h"

using namespace K_FUNC;
using namespace K_FLD;

// ============================================================================
/* Convert an unstructured array to a hexaedrical mesh */
// ============================================================================
PyObject* K_CONVERTER::convertUnstruct2Hexa(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PYPARSETUPLE_(args, O_, &array)) return NULL;

  // Check array
  E_Int ni, nj, nk;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f,
                                     ni, nj, nk, cn, eltType);
  if (res != 2) 
  {
    RELEASESHAREDS(array, f);
    PyErr_SetString(PyExc_TypeError, 
                    "convertUnstruct2Hexa: array must be unstructured.");
    return NULL;
  }
  else if (K_STRING::cmp(eltType, 4, "NGON") == 0)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError, 
                    "convertUnstruct2Hexa: array must be unstructured, not NGON.");
    return NULL;
  }

  // Acces universel sur BE/ME
  E_Int nc = cn->getNConnect();
  std::vector<char*> eltTypes;
  K_ARRAY::extractVars(eltType, eltTypes);

  E_Int npts = f->getSize(), api = f->getApi(), nfld = f->getNfld();
  E_Bool center = strchr(eltType,'*') != NULL;

  // Compute cumulative number of elements per conn. of the input ME (offsets)
  std::vector<E_Int> cumnepc(nc+1); cumnepc[0] = 0;
  for (E_Int ic = 0; ic < nc; ic++)
  {
    K_FLD::FldArrayI& cm = *(cn->getConnect(ic));
    cumnepc[ic+1] = cumnepc[ic] + cm.getSize();
  }

  // Determine the output eltType from the input dimension
  char* eltType2 = new char[K_ARRAY::VARSTRINGLENGTH]; eltType2[0] = '\0';
  E_Int dim = K_CONNECT::getDimME(eltTypes);
  if (dim == 2) strcpy(eltType2, "QUAD");
  else strcpy(eltType2, "HEXA");

  // Build new BE connectivity
  PyObject* tpl = K_ARRAY::buildArray3(nfld, varString, npts, cumnepc[nc],
                                       eltType2, center, api);
  FldArrayF* f2; FldArrayI* cn2;
  K_ARRAY::getFromArray3(tpl, f2, cn2);
  FldArrayI& cm2 = *(cn2->getConnect(0));

  // Boucle sur toutes les connectivites pour les remplir
  #pragma omp parallel
  {
    E_Int ind;
    for (E_Int ic = 0; ic < nc; ic++)
    {
      FldArrayI& cm = *(cn->getConnect(ic));
      char* eltTypConn = eltTypes[ic];
      E_Int nelts = cm.getSize();
      E_Int offset = cumnepc[ic];

      if (K_STRING::cmp(eltTypConn, 4, "QUAD") == 0)
      {
        // Copy existing connectivity
        #pragma omp for nowait
        for (E_Int i = 0; i < nelts; i++)
        {
          ind = i + offset;
          cm2(ind,1) = cm(i,1);
          cm2(ind,2) = cm(i,2);
          cm2(ind,3) = cm(i,3);
          cm2(ind,4) = cm(i,4);
        }
      }
      else if (K_STRING::cmp(eltTypConn, 4, "HEXA") == 0)
      {
        // Copy existing connectivity
        #pragma omp for nowait
        for (E_Int i = 0; i < nelts; i++)
        {
          ind = i + offset;
          cm2(ind,1) = cm(i,1);
          cm2(ind,2) = cm(i,2);
          cm2(ind,3) = cm(i,3);
          cm2(ind,4) = cm(i,4);
          cm2(ind,5) = cm(i,5);
          cm2(ind,6) = cm(i,6);
        }
      }
      else if (K_STRING::cmp(eltTypConn, 3, "BAR") == 0)
      {
        // Copy existing connectivity
        #pragma omp for nowait
        for (E_Int i = 0; i < nelts; i++)
        {
          ind = i + offset;
          cm2(ind,1) = cm(i,1);
          cm2(ind,2) = cm(i,2);
        }
      }
      else if (K_STRING::cmp(eltType, 3, "TRI") == 0)
      {
        #pragma omp for nowait
        for (E_Int i = 0; i < nelts; i++)
        {
          ind = i + offset;
          cm2(ind,1) = cm(i,1);
          cm2(ind,2) = cm(i,2);
          cm2(ind,3) = cm(i,3);
          cm2(ind,4) = cm(i,3);
        }
      }
      else if (K_STRING::cmp(eltType, 5, "TETRA") == 0)
      {
        #pragma omp for nowait
        for (E_Int i = 0; i < nelts; i++)
        {
          ind = i + offset;
          cm2(ind,1) = cm(i,1);
          cm2(ind,2) = cm(i,2);
          cm2(ind,3) = cm(i,3);
          cm2(ind,4) = cm(i,3);
          cm2(ind,5) = cm(i,4);
          cm2(ind,6) = cm(i,4);
          cm2(ind,7) = cm(i,4);
          cm2(ind,8) = cm(i,4);
        }
      }
      else if (K_STRING::cmp(eltType, 5, "PENTA") == 0)
      {
        #pragma omp for nowait
        for (E_Int i = 0; i < nelts; i++)
        {
          ind = i + offset;
          cm2(ind,1) = cm(i,1);
          cm2(ind,2) = cm(i,2);
          cm2(ind,3) = cm(i,3);
          cm2(ind,4) = cm(i,3);
          cm2(ind,5) = cm(i,4);
          cm2(ind,6) = cm(i,5);
          cm2(ind,7) = cm(i,6);
          cm2(ind,8) = cm(i,6);
        }
      }
    }

    // Copy fields to f2
    for (E_Int n = 1; n <= nfld; n++)
    {
      E_Float* fp = f->begin(n);
      E_Float* f2p = f2->begin(n);
      #pragma omp for nowait
      for (E_Int i = 0; i < npts; i++) f2p[i] = fp[i];
    }
  }
  
  delete[] eltType2;
  for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];
  RELEASESHAREDU(array, f, cn);
  RELEASESHAREDU(tpl, f2, cn2);
  return tpl;
}

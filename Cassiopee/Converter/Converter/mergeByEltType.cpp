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
#include <map>

using namespace K_FLD;

// ============================================================================
// Merge an unstructured array by element type such that each eltType is listed
// at most once
// ============================================================================
PyObject* K_CONVERTER::mergeByEltType(PyObject* self, PyObject* args)
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
                    "mergeByEltType: array must be unstructured.");
    return NULL;
  }
  else if (K_STRING::cmp(eltType, 4, "NGON") == 0)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError, 
                    "mergeByEltType: array must be unstructured, not NGON.");
    return NULL;
  }

  E_Int nc = cn->getNConnect();
  if (nc == 1)
  {
    RELEASESHAREDU(array, f, cn);
    return array;
  }

  // Check if at least one element type is repeated twice
  E_Bool isUnique = true;
  std::vector<char*> eltTypes;
  K_ARRAY::extractVars(eltType, eltTypes);

  // An element type is mapped to its position in the output conn.
  E_Int nc2 = 0;
  std::map<std::string, E_Int> eltMap;
  std::map<E_Int, std::pair<E_Int, E_Int> > ioConnMap;
  char* eltType2 = new char[K_ARRAY::VARSTRINGLENGTH]; eltType2[0] = '\0';
  std::vector<E_Int> nepc2;
  std::vector<std::vector<E_Int> > cumnepc2;

  for (E_Int ic = 0; ic < nc; ic++)
  {
    FldArrayI& cm = *(cn->getConnect(ic));
    E_Int nelts = cm.getSize();

    auto result = eltMap.insert(std::make_pair(std::string(eltTypes[ic]), nc2));
    if (result.second)  // newly inserted element
    {
      if (nc2 > 0) strcat(eltType2, ",");
      strcat(eltType2, eltTypes[ic]);
      nepc2.push_back(nelts);
      cumnepc2.push_back({0, nelts});
      ioConnMap[ic] = std::make_pair(nc2, 0);
      nc2++;
    }
    else  // already visited
    {
      isUnique = false;
      E_Int jc2 = result.first->second;
      nepc2[jc2] += nelts;  // accummulate
      cumnepc2[jc2].push_back(nepc2[jc2]);
      ioConnMap[ic] = std::make_pair(jc2, (E_Int)cumnepc2[jc2].size()-2);
    }
  }
    
  if (isUnique)
  {
    delete [] eltType2;
    for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];
    RELEASESHAREDU(array, f, cn);
    return array;
  }

  // Build new ME connectivity
  E_Int npts = f->getSize(), api = f->getApi(), nfld = f->getNfld();
  E_Bool center = strchr(eltType, '*') != NULL;

  PyObject* tpl = K_ARRAY::buildArray3(nfld, varString, npts, nepc2,
                                       eltType2, center, api);
  FldArrayF* f2; FldArrayI* cn2;
  K_ARRAY::getFromArray3(tpl, f2, cn2);

  #pragma omp parallel
  {
    E_Int ind, nelts, nvpe, jc2, offset;
    std::string eltTypeIc;

    for (E_Int ic = 0; ic < nc; ic++)
    {
      FldArrayI& cm = *(cn->getConnect(ic));
      nelts = cm.getSize();
      nvpe = cm.getNfld();

      eltTypeIc = std::string(eltTypes[ic]);
      jc2 = ioConnMap[ic].first;
      FldArrayI& cm2 = *(cn2->getConnect(jc2));
      offset = cumnepc2[jc2][ioConnMap[ic].second];

      // Copy existing connectivity
      #pragma omp for collapse(2) nowait
      for (E_Int i = 0; i < nelts; i++)
      for (E_Int j = 1; j <= nvpe; j++)
      {
        ind = i + offset;
        cm2(ind,j) = cm(i,j);
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
  RELEASESHAREDU(tpl, f2, cn2);
  RELEASESHAREDU(array, f, cn);
  return tpl;
}

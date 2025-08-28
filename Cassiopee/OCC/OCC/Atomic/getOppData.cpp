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

#include "occ.h"
#include <unordered_set>
#include <unordered_map>

// ============================================================================
/* Return the opposite data for an edge on a face */
// ============================================================================
PyObject* K_OCC::getOppData(PyObject* self, PyObject* args)
{
  // array is the face TRI mesh of the opposite face
  // indArray are the indices (starting at 0) in receiver mesh of the edge
  // indArrayD are the indices in donor/opposite mesh of the edge
  // nvR is the number of vertices of the receiver face
  // neD is the number of elements of the donor/opposite face without ghost cells
  PyObject* array;
  PyObject* indArray; PyObject* indArrayD;
  E_Int nvR, neD;
  if (!PYPARSETUPLE_(args, OOO_ II_ ,
                    &array, &indArray, &indArrayD, &nvR, &neD))
    return NULL;

  // Check array
  E_Int nil, njl, nkl;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, nil, njl, nkl, 
                                     cn, eltType);
  
  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "getOppData: invalid array.");
    return NULL;
  }

  if (res == 1)
  {
    RELEASESHAREDS(array, f);
    PyErr_SetString(PyExc_TypeError,
                    "getOppData: only for TRI or QUAD array.");
    return NULL;
  }

  E_Int nvpe, fac; 
  if (strcmp(eltType, "TRI") == 0) { nvpe = 3; fac = 3; }
  else if (strcmp(eltType, "QUAD") == 0) { nvpe = 4; fac = 2; }
  else
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "getOppData: only for TRI or QUAD array.");
    return NULL;
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "getOppData: coordinates missing in TRI array.");
    return NULL;
  }

  // Check inds
  E_Int* inds; E_Int* indsD;
  E_Int size; E_Int nfld;
  E_Int res2 = K_NUMPY::getFromNumpyArray(indArray, inds, size, nfld);
  res2 &= K_NUMPY::getFromNumpyArray(indArrayD, indsD, size, nfld);
  if (res2 == 0)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "getOppData: invalid index array.");
    return NULL;
  }

  // Recherche des vertex voisins et des triangles adjacents
  E_Int nvD = f->getSize();
  // printf("nvR=%d \n", nvR);
  // printf("nvD=%d \n", nvD);
  // printf("neD=%d \n", neD);
  
  // List of elements of faceOpp connected to the edge vertices
  std::vector< std::vector<E_Int> > cVE(nvD);
  K_CONNECT::connectEV2VE(*cn, cVE);

  E_Int ind, indd, inde, indn;
  E_Int vidx = 0; // vertex counter for newly inserted points
  FldArrayI& cm = *(cn->getConnect(0));
  
  std::unordered_map<E_Int, E_Int> vertexMap;
  std::vector<E_Int> vertexInsOrder; // maintains insertion order
  std::vector<E_Int> vertexIdxR; // keeps track of vertex index of receiver face
  TopologyOpt E;
  std::vector<E_Int> elt(nvpe);
  std::unordered_set<TopologyOpt, JenkinsHash<TopologyOpt> > eltSet;
  
  // Insert edge vertices first
  vertexInsOrder.reserve(fac*size); // ballpark
  vertexIdxR.reserve(fac*size);     // ballpark
  for (E_Int i = 0; i < size; i++)
  {
    ind = inds[size-1-i]+1; indd = indsD[i]+1;
    vertexMap.insert(std::make_pair(indd, i+1));
    vertexInsOrder.push_back(indd);
    vertexIdxR.push_back(ind);
  }

  // Insert vertices of ghost cells and ghost cells themselves
  // Loop over all but corner vertices
  for (E_Int i = 1; i < size-1; i++)
  {
    indd = indsD[i];
    const std::vector<E_Int>& elts = cVE[indd];
    for (size_t e = 0; e < elts.size(); e++)
    {
      inde = elts[e];
      if (neD <= inde) continue; // skip ghost cells of opposite face
      for (E_Int j = 0; j < nvpe; j++)
      {
        indn = cm(inde, j+1);
        // Inserted vertex, indn, maps to -1 if it isn't in the map already
        auto resV = vertexMap.insert(std::make_pair(indn, -1));
        if (resV.first->second == -1)
        {
          vidx++; resV.first->second = size+vidx;
          vertexInsOrder.push_back(indn);
          vertexIdxR.push_back(nvR+vidx);
        }
        elt[j] = resV.first->second;
      }
      E.set(elt, nvpe);
      eltSet.insert(E);
    }
  }
  
  vertexInsOrder.shrink_to_fit();
  E_Int vsize = vertexMap.size();
  E_Int esize = eltSet.size();

  // Export receiver and donor numpy point lists
  PyObject* rcvIds = K_NUMPY::buildNumpyArray(vsize, 1, 1, 1);
  E_Int* ridsp = K_NUMPY::getNumpyPtrI(rcvIds);
  PyObject* donorIds = K_NUMPY::buildNumpyArray(vsize-size, 1, 1, 1);
  E_Int* didsp = K_NUMPY::getNumpyPtrI(donorIds);

  // Ajout des triangles supplementaires (ghost cells)
  nfld = f->getNfld(); E_Int api = f->getApi();
  PyObject* tpl2 = K_ARRAY::buildArray3(nfld, varString, vsize, esize,
                                        eltType, false, api);
  FldArrayF* f2; FldArrayI* cn2;
  K_ARRAY::getFromArray3(tpl2, f2, cn2);
  
  #pragma omp parallel if (vsize > __MIN_SIZE_MEAN__)
  {
    E_Int ind1, ind2;
    // Fields
    for (E_Int n = 1; n <= nfld; n++)
    {
      E_Float* fp = f->begin(n);
      E_Float* f2p = f2->begin(n);
      #pragma omp for
      for (E_Int i = 0; i < vsize; i++)
      {
        ind1 = vertexInsOrder[i];
        ind2 = vertexMap[ind1];
        f2p[ind2-1] = fp[ind1-1];
      }
    }

    // Point lists
    #pragma omp for
    for (E_Int i = 0; i < size; i++) ridsp[i] = vertexIdxR[i];

    #pragma omp for
    for (E_Int i = size; i < vsize; i++)
    {
      ridsp[i] = vertexIdxR[i];
      didsp[i-size] = vertexInsOrder[i];
    }
  }

  // Ghost cell connectivity (TRI or QUAD)
  ind = 0;
  FldArrayI& cm2 = *(cn2->getConnect(0));
  for (auto it = eltSet.begin(); it != eltSet.end(); it++)
  {
    const auto vert = it->p_;
    for (E_Int j = 0; j < nvpe; j++)
    {
      cm2(ind,j+1) = vert[j];
    }
    ind++;
  }

  RELEASESHAREDU(tpl2, f2, cn2);
  RELEASESHAREDU(array, f, cn);
  Py_DECREF(indArray); Py_DECREF(indArrayD);

  return Py_BuildValue("OOO", tpl2, rcvIds, donorIds);
}

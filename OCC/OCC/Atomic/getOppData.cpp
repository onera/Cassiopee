/*    
    Copyright 2013-2024 Onera.

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
  // indArray are the indices (1st vertex: 0) in TRI mesh of the edge
  PyObject* array;
  PyObject* indArray;
  if (!PYPARSETUPLE_(args, OO_ , &array, &indArray)) return NULL;

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
                    "getOppData: only for TRI array.");
    return NULL;
  }

  if (strcmp(eltType, "TRI") != 0)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "getOppData: only for TRI array.");
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
  E_Int* inds; E_Int size; E_Int nfld;
  E_Int res2 = K_NUMPY::getFromNumpyArray(indArray, inds, size, nfld, true);
  if (res2 == 0)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "getOppData: invalid index array.");
    return NULL;
  }

  // Recherche des vertex voisins et des triangles adjacents
  E_Int nv = f->getSize();
  //printf("nv=%d \n", nv);
  
  // List of elements of faceOpp connected to the edge vertices
  std::vector< std::vector<E_Int> > cVE(nv);
  K_CONNECT::connectEV2VE(*cn, cVE);

  E_Int ind, inde, indn;
  E_Int vidx = size; // vertex counter starting from max range index of edge
  FldArrayI& cm = *(cn->getConnect(0));
  
  std::unordered_map<E_Int, E_Int> vertexMap;
  std::vector<E_Int> vertexInsOrder; // maintains insertion order
  TopologyOpt E;
  E_Int nvpe = 3; std::vector<E_Int> elt(nvpe);
  std::unordered_set<TopologyOpt, JenkinsHash<TopologyOpt> > eltSet;
  
  // Insert edge vertices in reverse order
  vertexInsOrder.reserve(2*size); // ballpark
  for (E_Int i = 0; i < size; i++)
  {
    ind = inds[i]+1;
    vertexMap.insert(std::make_pair(ind, size-i));
    vertexInsOrder.push_back(ind);
  }

  // Insert vertices of ghost cells and ghost cells themselves
  // Loop over all but corner vertices
  for (E_Int i = 1; i < size-1; i++)
  {
    ind = inds[i];
    const std::vector<E_Int>& elts = cVE[ind];
    for (size_t e = 0; e < elts.size(); e++)
    {
      inde = elts[e];
      for (E_Int j = 0; j < nvpe; j++)
      {
        indn = cm(inde, j+1);
        // Inserted vertex, indn, maps to -1 if it isn't in the map already
        auto resV = vertexMap.insert(std::make_pair(indn, -1));
        if (resV.first->second == -1)
        {
          vidx++; resV.first->second = vidx;
          vertexInsOrder.push_back(indn);
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

  // Ajout des triangles supplementaires (ghost cells)
  nfld = f->getNfld(); E_Int api = f->getApi();
  PyObject* tpl2 = K_ARRAY::buildArray3(nfld, varString, vsize, esize,
                                        "TRI", false, api);
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
  }

  // Connectivity TRI
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
  Py_DECREF(indArray);

  return tpl2;
}
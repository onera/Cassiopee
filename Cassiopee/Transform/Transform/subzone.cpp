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

// subzone
#include <map>
#include <unordered_set>
#include <unordered_map>
#include "transform.h"

using namespace std;
using namespace K_FLD;

// ============================================================================
/* Subzone a structured mesh from imin,imax,... indices */
// ============================================================================
PyObject* K_TRANSFORM::subzoneStruct(PyObject* self, PyObject* args)
{
  E_Int i1, j1, k1; E_Int i2, j2, k2;
  E_Int ind, ind2; PyObject* array;
  if (!PYPARSETUPLE_(args, O_ TIII_ TIII_,
                     &array, &i1, &j1, &k1, &i2, &j2, &k2))
  {
      return NULL;
  }
  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = 
    K_ARRAY::getFromArray3(array, varString, f, im, jm, km, cn, eltType);
  E_Int imjm = im*jm;

  if (res == 1)
  { 
    // Negative -> max indices
    if (i1 < 0) i1 = im+i1+1;
    if (j1 < 0) j1 = jm+j1+1;
    if (k1 < 0) k1 = km+k1+1;
    
    if (i2 < 0) i2 = im+i2+1;
    if (j2 < 0) j2 = jm+j2+1;
    if (k2 < 0) k2 = km+k2+1;
    
    E_Int in = i2-i1+1;
    E_Int jn = j2-j1+1;
    E_Int kn = k2-k1+1;
    E_Int injn = in*jn;
    E_Int nfld = f->getNfld();
    
    // Check
    if (i2 > im || i1 > im || i2 < i1 || i1 < 1 ||
        j2 > jm || j1 > jm || j2 < j1 || j1 < 1 ||
        k2 > km || k1 > km || k2 < k1 || k1 < 1)
    {
      RELEASESHAREDS(array, f);
      printf("Warning: subzone: mesh dimensions are: " SF_D_ " x " SF_D_ " x " SF_D_ "\n",
             im, jm, km);
      PyErr_SetString(PyExc_TypeError,
                      "subzone: wrong choice of index.");
      return NULL;
    }
    
    // Construit l'array resultat
    //E_Int api = f->getApi();
    //if (api == 2) api=3;
    PyObject* tpl= K_ARRAY::buildArray3(nfld, varString, in, jn, kn);
    E_Float* fnp = K_ARRAY::getFieldPtr(tpl);
    FldArrayF subzone0(injn*kn, nfld, fnp, true);

    if (k2-k1 > j2-j1)
    {
#pragma omp parallel shared (im,subzone0,f) private (ind,ind2) if (k2 > 30)
      {
        for (E_Int n = 1; n <= nfld; n++)
        {
          E_Float* sp = subzone0.begin(n);
          E_Float* fp = f->begin(n);
#pragma omp for nowait
          for (E_Int k = k1; k <= k2; k++)
            for (E_Int j = j1; j <= j2; j++)
              for (E_Int i = i1; i <= i2; i++)
              {
                ind  = i-1  + (j-1)*im + (k-1)*imjm;
                ind2 = (i-i1)+(j-j1)*in+(k-k1)*injn;
                sp[ind2] = fp[ind];
              }
        }
      }
    }
    else
    {
#pragma omp parallel shared (im,subzone0,f) private (ind,ind2) if (j2 > 30)
      {
        for (E_Int n = 1; n <= nfld; n++)
        {
          E_Float* sp = subzone0.begin(n);
          E_Float* fp = f->begin(n);
#pragma omp for nowait
          for (E_Int j = j1; j <= j2; j++)
            for (E_Int k = k1; k <= k2; k++)
              for (E_Int i = i1; i <= i2; i++)
              {
                ind  = i-1 + (j-1)*im + (k-1)*imjm;
                ind2 = (i-i1)+(j-j1)*in+(k-k1)*injn;
                sp[ind2] = fp[ind];
              }
        }
      }
    }

    RELEASESHAREDS(array, f);
    return tpl;
  }
  else if (res == 2)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "subzone: cannot be used on an unstructured array.");
    return NULL;
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "subzone: unknown type of array.");
    return NULL;
  }
}     

// ============================================================================
/* Subzone an unstructured mesh given a list of vertex indices */
// ============================================================================
PyObject* K_TRANSFORM::subzoneUnstruct(PyObject* self, PyObject* args)
{
  PyObject* array;
  PyObject* listOfNodes;
  if (!PYPARSETUPLE_(args, OO_, &array, &listOfNodes)) return NULL;

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, im, jm, km, cn, eltType); 
  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "subzone: unknown type of array.");
    return NULL;
  }

  if (res == 1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "subzone: cannot be used on a structured array.");
    RELEASESHAREDS(array, f); return NULL;
  }
  if (strcmp(eltType,"NGON") == 0) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "subzone: type='nodes' not implemented for a NGON array.");
    RELEASESHAREDU(array, f, cn); return NULL;
  }
  FldArrayI indices;
  E_Int ok = K_ARRAY::getFromList(listOfNodes, indices);
  if (ok == 0) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "subzone: 2nd argument must be an integer list or a numpy.");
    RELEASESHAREDU(array, f, cn); return NULL;
  }
  
  E_Int n = indices.getSize();
  E_Int npts = f->getSize(), nfld = f->getNfld(), api = f->getApi();
  E_Int* indicesp = indices.begin();
  FldArrayI tmap(npts); tmap.setAllValuesAt(-1); E_Int* tmapP = tmap.begin();

  if (strcmp(eltType, "NODE") == 0) 
  {
    PyObject* tpl = K_ARRAY::buildArray3(nfld, varString, n, 0, eltType,
                                         false, api);
    FldArrayF* f2; K_ARRAY::getFromArray3(tpl, f2);

    #pragma omp parallel
    {
      E_Int indv;
      
      #pragma omp for
      for (E_Int i = 0; i < n; i++) tmapP[indicesp[i]-1] = i;
      
      // Mapping f -> f2
      for (E_Int eq = 1; eq <= nfld; eq++)
      {
        E_Float* fp = f->begin(eq);
        E_Float* f2p = f2->begin(eq);
        #pragma omp for
        for (E_Int i = 0; i < n; i++)
        {
          indv = indicesp[i]-1;
          f2p[i] = fp[indv];
        }
      }
    }

    RELEASESHAREDS(tpl, f2); RELEASESHAREDU(array, f, cn);
    return tpl;
  }

  // Selectionne les elements subzones - BE/ME
  #pragma omp parallel for
  for (E_Int i = 0; i < n; i++) tmapP[indicesp[i]-1] = i;

  E_Int nc = cn->getNConnect(), nc2 = 0;
  vector<vector<E_Int> > eltList(nc);
  vector<E_Int> nelts2(nc);

  // Tag les points deja inseres dans un element
  // Permet de savoir si un point n'est pas utilise
  FldArrayI tag(n); tag.setAllValuesAtNull();
  E_Int* tagp = tag.begin();

  for (E_Int ic = 0; ic < nc; ic++)
  {
    E_Int compt, indv;
    FldArrayI& cm = *(cn->getConnect(ic));
    E_Int nelts = cn->getSize();
    E_Int nvpe = cn->getNfld();

    for (E_Int i = 0; i < nelts; i++)
    {
      compt = 0;
      for (E_Int j = 1; j <= nvpe; j++)
      {
        indv = cm(i,j)-1;
        if (tmapP[indv] != -1) compt++; 
      }
      if (compt == nvpe)
      {
        // Add element to list and mark vertex as seen
        eltList[ic].push_back(i);
        for (E_Int j = 1; j <= nvpe; j++)
        {
          indv = cm(i,j)-1;
          tagp[tmapP[indv]]++;
        }
      }
    }
    if (eltList[ic].size())
    {
      // Count number of elements in this connectivity
      nelts2[nc2] = eltList[ic].size(); nc2++;
    }
  }
  nelts2.resize(nc2);

  // Build connectivity
  E_Int ierr = 0;
  PyObject* tpl = K_ARRAY::buildArray3(nfld, varString, n, nelts2, eltType, false, api);
  FldArrayF* f2; FldArrayI* cn2;
  K_ARRAY::getFromArray3(tpl, f2, cn2);

  #pragma omp parallel
  {
    E_Int indv, noe, nvpe, ic2 = 0;

    // Check if any non-used vertices
    #pragma omp for
    for (E_Int i = 0; i < n; i++)
    {
      if (tagp[i] == 0) ierr = 1;
    }

    for (E_Int ic = 0; ic < nc; ic++)
    {
      if (!eltList[ic].size()) continue;
      FldArrayI& cm = *(cn->getConnect(ic));
      FldArrayI& cm2 = *(cn2->getConnect(ic2));
      nvpe = cm.getNfld();

      #pragma omp for
      for (E_Int i = 0; i < nelts2[ic2]; i++)
      {
        noe = eltList[ic][i];
        for (E_Int j = 1; j <= nvpe; j++)
        {  
          indv = cm(noe,j)-1;
          cm2(i,j) = tmapP[indv]+1;
        }
      }
      ic2++;
    }

    for (E_Int eq = 1; eq <= nfld; eq++) 
    {
      E_Float* fp = f->begin(eq);
      E_Float* f2p = f2->begin(eq);
      #pragma omp for
      for (E_Int i = 0; i < n; i++)
      {
        indv = indicesp[i]-1;
        f2p[i] = fp[indv];
      }
    }
  }

  if (ierr == 1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "subzoneUnstruct: indices for unstructured subzone must be contiguous.");
    RELEASESHAREDU(array, f, cn);
    return NULL;
  }

  RELEASESHAREDU(tpl, f2, cn2); RELEASESHAREDU(array, f, cn);
  return tpl;
}

// ============================================================================
/* Subzone a unstructured zone: input array is located here at centers       */
// ============================================================================
PyObject* K_TRANSFORM::subzoneUnstructBoth(PyObject* self, PyObject* args)
{
  PyObject *arrayNodes, *arrayCenters;
  PyObject* listOfNodes;
  if (!PYPARSETUPLE_(args, OOO_, &arrayNodes, &arrayCenters, &listOfNodes))
  {
    return NULL;
  }

  FldArrayI indices;
  E_Int ok = K_ARRAY::getFromList(listOfNodes, indices);
  if (ok == 0) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "subzone: 2nd argument must be an integer list or a numpy.");
    return NULL;
  }

  // Check array of nodes
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(arrayNodes, varString, f, im, jm, km, cn, eltType); 
  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "subzone: unknown type of array.");
    return NULL;
  }
  if (res == 1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "subzone: cannot be used on a structured array.");
    RELEASESHAREDS(arrayNodes, f); return NULL;
  }
  if (strcmp(eltType,"NGON") == 0) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "subzone: type='nodes' not implemented for a NGON array.");
    RELEASESHAREDU(arrayNodes, f, cn); return NULL;
  }
  // Check array of centers
  E_Int imc, jmc, kmc;
  FldArrayF* fc; FldArrayI* cnc;
  char* varStringc; char* eltTypec;
  E_Int resc = K_ARRAY::getFromArray3(arrayCenters, varStringc, fc, imc, jmc, kmc, cnc, eltTypec); 
  if (resc != 1 && resc != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "subzone: unknown type of array.");
    RELEASESHAREDU(arrayNodes, f, cn); return NULL;
  }
  if (resc == 1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "subzone: cannot be used on a structured array.");
    RELEASESHAREDU(arrayNodes, f, cn); RELEASESHAREDS(arrayCenters,fc); return NULL;
  }
  if (K_STRING::cmp(eltTypec,"NGON*") == 0) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "subzone: type='nodes' not implemented for a NGON array.");
    RELEASESHAREDU(arrayNodes,f, cn); RELEASESHAREDU(arrayCenters, fc, cnc); 
    return NULL;
  }

  if (cnc->getSize() != cn->getSize())
  {
    PyErr_SetString(PyExc_TypeError,
                    "subzone: sizes of connectivities at nodes and centers are not equal.");
    RELEASESHAREDU(arrayNodes,f, cn); RELEASESHAREDU(arrayCenters, fc, cnc); 
    return NULL;
  }

  PyObject* l = PyList_New(0);
  E_Int n = indices.getSize();
  E_Int npts = f->getSize(), nelts = fc->getSize(), api = f->getApi();
  E_Int nfld = f->getNfld(), nfldc = fc->getNfld();
  E_Int* indicesp = indices.begin();
  FldArrayI tmap(npts); tmap.setAllValuesAt(-1); E_Int* tmapP = tmap.begin();
  
  if (strcmp(eltType, "NODE") == 0) 
  {
    PyObject* tpln = K_ARRAY::buildArray3(nfld, varString, n, 0, eltType,
                                          false, api);
    FldArrayF* f2; K_ARRAY::getFromArray3(tpln, f2);
    PyObject* tplc = K_ARRAY::buildArray3(nfldc, varStringc, n, 0, eltType,
                                          true, api);
    FldArrayF* fc2; K_ARRAY::getFromArray3(tplc, fc2);

    #pragma omp parallel
    {
      E_Int indv;

      #pragma omp for
      for (E_Int i = 0; i < n; i++) tmapP[indicesp[i]-1] = i;
      
      for (E_Int eq = 1; eq <= nfld; eq++) 
      {
        E_Float* fp = f->begin(eq);
        E_Float* f2p = f2->begin(eq);
        #pragma omp for
        for (E_Int i = 0; i < n; i++)
        {
          indv = indicesp[i]-1;
          f2p[i] = fp[indv];
        }
      }

      for (E_Int eq = 1; eq <= nfldc; eq++) 
      {
        E_Float* fcp = fc->begin(eq);
        E_Float* fc2p = fc2->begin(eq);
        #pragma omp for
        for (E_Int i = 0; i < n; i++)
        {
          indv = indicesp[i]-1;
          fc2p[i] = fcp[indv];
        }
      }
    }
    
    RELEASESHAREDS(tpln, f2); PyList_Append(l, tpln); Py_DECREF(tpln);
    RELEASESHAREDS(tplc, fc2); PyList_Append(l, tplc); Py_DECREF(tplc);
    RELEASESHAREDU(arrayNodes, f, cn); RELEASESHAREDU(arrayCenters, fc, cnc);  
    return l;
  }

  E_Int neltstot = 0;
  E_Int nc = cn->getNConnect();
  for (E_Int ic = 0; ic < nc; ic++)
  {
    FldArrayI& cm = *(cn->getConnect(ic));
    neltstot += cm.getSize();
  }
  
  if (nelts != neltstot)
  {
    PyErr_SetString(PyExc_TypeError,
                    "subzone: arrays located at nodes and centers are not consistent.");
    RELEASESHAREDU(arrayNodes, f, cn); RELEASESHAREDU(arrayCenters, fc, cnc); 
    return NULL;
  }

  // Selectionne les elements subzones - BE/ME
  const E_Int numThreads = __NUMTHREADS__;
  vector<vector<vector<E_Int> > > threadEltList(numThreads);

  // Tag les points deja inseres dans un element
  // Permet de savoir si un point n'est pas utilise
  FldArrayI tag(n); tag.setAllValuesAtNull();
  E_Int* tagp = tag.begin();
  
  #pragma omp parallel num_threads(numThreads)
  {
    #pragma omp for
    for (E_Int i = 0; i < n; i++) tmapP[indicesp[i]-1] = i;

    E_Int threadId = __CURRENT_THREAD__;
    vector<vector<E_Int> >& threadIdEltList = threadEltList[threadId];
    threadIdEltList.resize(nc);
  
    for (E_Int ic = 0; ic < nc; ic++)
    {
      E_Int compt, indv;
      FldArrayI& cm = *(cn->getConnect(ic));
      E_Int nelts = cn->getSize();
      E_Int nvpe = cn->getNfld();

      #pragma omp for
      for (E_Int i = 0; i < nelts; i++)
      {
        compt = 0;
        for (E_Int j = 1; j <= nvpe; j++)
        {
          indv = cm(i,j)-1;
          if (tmapP[indv] != -1) compt++; 
        }
        if (compt == nvpe)
        {
          // Add element to list and mark vertex as seen
          threadIdEltList[ic].push_back(i);
          for (E_Int j = 1; j <= nvpe; j++)
          {
            indv = cm(i,j)-1;
            tagp[tmapP[indv]]++;
          }
        }
      }
    }
  }

  // Combine the lists of elements from all threads
  E_Int nc2 = 0;
  vector<E_Int> nelts2(nc);
  vector<vector<E_Int> > eltList(nc);
  for (E_Int ic = 0; ic < nc; ic++)
  {
    for (E_Int t = 0; t < numThreads; t++)
    {
      eltList[ic].insert(eltList[ic].end(), threadEltList[t][ic].begin(), threadEltList[t][ic].end());
    }
    if (eltList[ic].size()) { nelts2[nc2] = eltList[ic].size(); nc2++; }
  }
  nelts2.resize(nc2);

  // Build connectivity
  E_Int ierr = 0;
  PyObject* tpln = K_ARRAY::buildArray3(nfld, varString, n, nelts2, eltType, false, api);
  FldArrayF* f2; FldArrayI* cn2;
  K_ARRAY::getFromArray3(tpln, f2, cn2);

  #pragma omp parallel
  {
    E_Int indv, noe, nvpe, ic2 = 0;

    // Check if any non-used vertices
    #pragma omp for
    for (E_Int i = 0; i < n; i++)
    {
      if (tagp[i] == 0) ierr = 1;
    }

    for (E_Int ic = 0; ic < nc; ic++)
    {
      if (!eltList[ic].size()) continue;
      FldArrayI& cm = *(cn->getConnect(ic));
      FldArrayI& cm2 = *(cn2->getConnect(ic2));
      nvpe = cm.getNfld();

      #pragma omp for
      for (E_Int i = 0; i < nelts2[ic2]; i++)
      {
        noe = eltList[ic][i];
        for (E_Int j = 1; j <= nvpe; j++)
        {  
          indv = cm(noe,j)-1;
          cm2(i,j) = tmapP[indv]+1;
        }
      }
      ic2++;
    }
  }

  if (ierr == 1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "subzoneUnstruct: indices for unstructured subzone must be contiguous.");
    RELEASESHAREDU(arrayNodes, f, cn); RELEASESHAREDU(arrayCenters, fc, cnc);
    return NULL;
  }

  PyObject* tplc = K_ARRAY::buildArray3(nfldc, varStringc, n, *cn2, eltType, 1, api, true);
  FldArrayF* fc2; K_ARRAY::getFromArray3(tplc, fc2);

  #pragma omp parallel
  {
    E_Int ind;

    // Champs aux noeuds
    for (E_Int eq = 1; eq <= nfld; eq++) 
    {
      E_Float* fp = f->begin(eq);
      E_Float* f2p = f2->begin(eq);
      #pragma omp for
      for (E_Int i = 0; i < n; i++)
      {
        ind = indicesp[i]-1;
        f2p[i] = fp[ind];
      }
    }

    // Champs aux centres
    for (E_Int eq = 1; eq <= nfldc; eq++) 
    {
      E_Float* fcp = fc->begin(eq);
      E_Float* fc2p = fc2->begin(eq);
      for (E_Int ic = 0; ic < nc; ic++)
      {
        #pragma omp for
        for (E_Int i = 0; i < nelts2[ic]; i++)
        {
          ind = eltList[ic][i];
          fc2p[i] = fcp[ind];
        }
      }
    }
  }

  RELEASESHAREDU(tpln, f2, cn2); PyList_Append(l, tpln); Py_DECREF(tpln);
  RELEASESHAREDS(tplc, fc2); PyList_Append(l, tplc); Py_DECREF(tplc);
  RELEASESHAREDU(arrayNodes, f, cn); RELEASESHAREDU(arrayCenters, fc, cnc);
  return l;
}

// ============================================================================
/* Subzone an unstructured mesh by element indices */
// ============================================================================
PyObject* K_TRANSFORM::subzoneElements(PyObject* self, PyObject* args)
{
  PyObject* array;
  PyObject* listOfElts;
  if (!PYPARSETUPLE_(args, OO_, &array, &listOfElts)) return NULL;

  // Build element list
  FldArrayI eltList;
  E_Int ret = K_ARRAY::getFromList(listOfElts, eltList);
  if (ret == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "subzone: argument must be a list of element indices (starting from 0).");
    return NULL;
  }

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, im, jm, km, cn, eltType); 
  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "subzone: unknown type of array.");
    return NULL;
  }
  if (res == 1)
  {    
    PyErr_SetString(PyExc_TypeError,
                    "subzone: can not be used on a structured array.");
    RELEASESHAREDS(array, f); return NULL;
  }
  
  PyObject* tpl = NULL;
  E_Int n = eltList.getSize();
  E_Int npts = f->getSize(), nfld = f->getNfld(), api = f->getApi();
  if (api == 2) api = 3;
  
  if (K_STRING::cmp(eltType, "NGON") == 0) // NGON
  {
    E_Int shift = 1; if (api == 3) shift = 0;
    E_Int ngonType = 1; // CGNSv3 compact array1
    if (api == 2) ngonType = 2; // CGNSv3, array2
    else if (api == 3) ngonType = 3; // force CGNSv4, array3

    E_Int *ngon = cn->getNGon(), *indPG = cn->getIndPG();
    E_Int *nface = cn->getNFace(), *indPH = cn->getIndPH();
    E_Int sizeFN = cn->getSizeNGon(), sizeEF = cn->getSizeNFace();
    E_Int nfacesTot = cn->getNFaces();

    E_Int sizeEF2 = 0, sizeFN2 = 0;
    FldArrayI cEFTemp(sizeEF); E_Int* ptrEFTemp = cEFTemp.begin();
    FldArrayI cPHTemp(n);
    FldArrayI indirFaces(nfacesTot); indirFaces.setAllValuesAt(-1);
    E_Int* indirFacesp = indirFaces.begin();
    E_Int indface, indFaceOut, nfaces, pose, posf;
    E_Int nbFacesOut = 0;
    vector<E_Int> origIndicesOfFaces;

    for (E_Int noe = 0; noe < n; noe++)
    {
      // construction de la connectivite elt/faces
      pose = eltList[noe];
      E_Int* elt = cn->getElt(pose, nfaces, nface, indPH);
      ptrEFTemp[0] = nfaces;
      for (E_Int nof = 0; nof < nfaces; nof++)
      {
        indface = elt[nof]-1;
        if (indirFacesp[indface] == -1) 
        {
          indFaceOut = nbFacesOut; 
          indirFacesp[indface] = indFaceOut; 
          nbFacesOut++;
          origIndicesOfFaces.push_back(indface);
        }
        else indFaceOut = indirFacesp[indface];
        ptrEFTemp[nof+shift] = indFaceOut+1;
      }      
      ptrEFTemp += nfaces+shift; sizeEF2 += nfaces+shift; 
      cPHTemp[noe] = nfaces;
    }
    indirFaces.malloc(0); cEFTemp.resize(sizeEF2);

    // construction de la connectivite Faces/Noeuds
    FldArrayI cFNTemp(sizeFN); E_Int* ptrFNTemp = cFNTemp.begin();
    FldArrayI cPGTemp(nbFacesOut);
    FldArrayI indirNodes(npts); indirNodes.setAllValuesAt(-1); E_Int* indirNp = indirNodes.begin();
    E_Int indnode, nbnodes;
    E_Int nUniqueNodes = 0;
    for (E_Int nfe = 0; nfe < nbFacesOut; nfe++)
    {
      posf = origIndicesOfFaces[nfe]; //demarre a 0
      E_Int* face = cn->getFace(posf, nbnodes, ngon, indPG);
      ptrFNTemp[0] = nbnodes;
      for (E_Int p = 0; p < nbnodes; p++)
      {
        indnode = face[p]-1;
        if (indirNp[indnode] == -1) //creation
        {
          indirNp[indnode] = nUniqueNodes+1;
          ptrFNTemp[p+shift] = nUniqueNodes+1;
          nUniqueNodes++;
        }
        else 
        {
          ptrFNTemp[p+shift] = indirNp[indnode];
        }
      }
      ptrFNTemp += nbnodes+shift; sizeFN2 += nbnodes+shift;
      cPGTemp[nfe] = nbnodes;
    }
    origIndicesOfFaces.clear();      

    // construit l'array de sortie
    tpl = K_ARRAY::buildArray3(nfld, varString, nUniqueNodes, n, nbFacesOut,
                               "NGON", sizeFN2, sizeEF2, ngonType, false, api);
    FldArrayF* f2; FldArrayI* cn2;
    K_ARRAY::getFromArray3(tpl, f2, cn2);
    E_Int *ngon2 = cn2->getNGon(), *nface2 = cn2->getNFace();
    E_Int *indPG2 = NULL, *indPH2 = NULL;
    if (api == 2 || api == 3)
    {
      indPG2 = cn2->getIndPG(); indPH2 = cn2->getIndPH();
    }
  
    #pragma omp parallel default(shared)
    {
      E_Int indf;
      for (E_Int eq = 1; eq <= nfld; eq++)
      {
        E_Float* fp = f->begin(eq);
        E_Float* f2p = f2->begin(eq);
        #pragma omp for
        for (E_Int ind = 0; ind < npts; ind++)
        {
          indf = indirNp[ind]-1;
          if (indf > -1) f2p[indf] = fp[ind];
        }
      }

      // reconstruction de la connectivite finale
      #pragma omp for nowait
      for (E_Int i = 0; i < sizeFN2; i++) ngon2[i] = cFNTemp[i];
      #pragma omp for nowait
      for (E_Int i = 0; i < sizeEF2; i++) nface2[i] = cEFTemp[i];

      if (api == 2 || api == 3) // array2 or array3
      {
        #pragma omp for nowait
        for (E_Int i = 0; i < nbFacesOut; i++) indPG2[i] = cPGTemp[i];
        #pragma omp for nowait
        for (E_Int i = 0; i < n; i++) indPH2[i] = cPHTemp[i];
      }
    }

    RELEASESHAREDU(tpl, f2, cn2);
  }
  else // maillage par elements BE/ME
  {
    E_Int nc = cn->getNConnect();
    E_Int binIndex, nelts;
    E_Int nc2 = 0, elOffset = 0, nUniqueNodes = 0;
    vector<E_Int> neltspc(nc+1), neltspc2(nc), nvpe(nc);

    char* eltType2 = new char[50]; strcpy(eltType2, "");
    vector<char*> eltTypes;
    K_ARRAY::extractVars(eltType, eltTypes);
    
    // Compute the cumulative number of elements per connectivity,
    // these are the bounds of the bins that are used to inform on
    // the element type of an element index
    neltspc[0] = elOffset;
    for (E_Int ic = 0; ic < nc; ic++)
    {
      FldArrayI& cm = *(cn->getConnect(ic));
      nvpe[ic] = cm.getNfld();
      nelts = cm.getSize();
      elOffset += nelts;
      neltspc[ic+1] = elOffset;
    }

    // Bin input elements
    vector<vector<E_Int> > binnedEltList(nc);
    for (E_Int i = 0; i < n; i++)
    {
      E_Int noe = eltList[i];
      binIndex = std::upper_bound(neltspc.begin(), neltspc.end(), noe) - neltspc.begin() - 1;
      binnedEltList[binIndex].push_back(noe);
    }

    // Selectionne les elements subzones, calcule le nombre de noeuds uniques
    FldArrayI indirNodes(npts); indirNodes.setAllValuesAt(-1);
    E_Int* indirNodesp = indirNodes.begin();
    vector<E_Int> listOfNodes(npts);

    for (E_Int ic = 0; ic < nc; ic++)
    {
      // Skip empty connectivities
      if (binnedEltList[ic].size())
      {
        strcat(eltType2, eltTypes[ic]); // Build eltType2
        neltspc2[nc2] = binnedEltList[ic].size();
        
        FldArrayI& cm = *(cn->getConnect(ic));
        for (E_Int i = 0; i < neltspc2[nc2]; i++)
        {
          E_Int noe = binnedEltList[ic][i];
          for (E_Int v = 1; v <= nvpe[ic]; v++)
          {
            E_Int indv = cm(noe,v)-1;
            if (indirNodesp[indv] == -1)
            { 
              listOfNodes[nUniqueNodes] = indv; nUniqueNodes++;
              indirNodesp[indv] = nUniqueNodes;
            }
          }
        }

        nc2++;
      }
    }
    neltspc2.resize(nc2); listOfNodes.resize(nUniqueNodes);

    // Create connectivities
    tpl = K_ARRAY::buildArray3(nfld, varString, nUniqueNodes, neltspc2,
                               eltType2, false, api);
    FldArrayF* f2; FldArrayI* cn2;
    K_ARRAY::getFromArray3(tpl, f2, cn2);

    delete[] eltType2;
    for (E_Int ic = 0; ic < nc; ic++) delete [] eltTypes[ic];

    #pragma omp parallel default(shared)
    {
      E_Int indv, noe;
      for (E_Int ic = 0; ic < nc; ic++)
      {
        #pragma omp for
        for (size_t i = 0; i < binnedEltList[ic].size(); i++)
        {
          if (!binnedEltList[ic].size()) continue;
          FldArrayI& cm = *(cn->getConnect(ic));
          FldArrayI& cm2 = *(cn2->getConnect(ic));

          noe = binnedEltList[ic][i];
          for (E_Int v = 1; v <= nvpe[ic]; v++)
          {
            indv = cm(noe,v)-1;
            cm2(i,v) = indirNodesp[indv];
          }
        }
      }
      
      for (E_Int eq = 1; eq <= nfld; eq++)
      {
        E_Float* fp = f->begin(eq);
        E_Float* f2p = f2->begin(eq);
        #pragma omp for
        for (E_Int i = 0; i < nUniqueNodes; i++) f2p[i] = fp[listOfNodes[i]];    
      }
    }

    indirNodes.malloc(0);
    RELEASESHAREDU(tpl, f2, cn2);
  }

  RELEASESHAREDU(array, f, cn);
  return tpl;
}

// ============================================================================
/* Subzone an unstructured mesh by element indices */
// ============================================================================
PyObject* K_TRANSFORM::subzoneElementsBoth(PyObject* self, PyObject* args)
{
  PyObject *arrayNodes, *arrayCenters;
  PyObject* listOfElts;
  if (!PYPARSETUPLE_(args, OOO_, &arrayNodes, &arrayCenters, &listOfElts))
  {
      return NULL;
  }

  // Build element list
  FldArrayI eltList;
  E_Int ret = K_ARRAY::getFromList(listOfElts, eltList);
  if (ret == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "subzone: argument must be a list of element indices (starting from 0).");
    return NULL;
  }

  // Check array of nodes
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(arrayNodes, varString, f, im, jm, km, cn, eltType); 
  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "subzone: unknown type of array.");
    return NULL;
  }
  if (res == 1)
  {    
    PyErr_SetString(PyExc_TypeError,
                    "subzone: can not be used on a structured array.");
    RELEASESHAREDS(arrayNodes, f); return NULL;
  }

  // Check array of centers
  E_Int imc, jmc, kmc;
  FldArrayF* fc; FldArrayI* cnc;
  char* varStringc; char* eltTypec;
  E_Int resc = K_ARRAY::getFromArray3(arrayCenters, varStringc, fc, imc, jmc, kmc, cnc, eltTypec); 
  if (resc != 1 && resc != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "subzone: unknown type of array.");
    RELEASESHAREDU(arrayNodes, f, cn); return NULL;
  }
  if (resc == 1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "subzone: cannot be used on a structured array.");
    RELEASESHAREDU(arrayNodes, f, cn); RELEASESHAREDS(arrayCenters, fc); return NULL;
  }

  PyObject* l = PyList_New(0);
  E_Int n = eltList.getSize();

  E_Int npts = f->getSize(), nfld = f->getNfld(), api = f->getApi();
  E_Int nfldc = fc->getNfld();
  if (api == 2) api = 3;

  if (K_STRING::cmp(eltType, "NGON") == 0) // NGON
  {
    E_Int shift = 1; if (api == 3) shift = 0;
    E_Int ngonType = 1; // CGNSv3 compact array1
    if (api == 2) ngonType = 2; // CGNSv3, array2
    else if (api == 3) ngonType = 3; // force CGNSv4, array3

    E_Int *ngon = cn->getNGon(), *indPG = cn->getIndPG();
    E_Int *nface = cn->getNFace(), *indPH = cn->getIndPH();
    E_Int sizeFN = cn->getSizeNGon(), sizeEF = cn->getSizeNFace();
    E_Int nfacesTot = cn->getNFaces();

    E_Int sizeEF2 = 0, sizeFN2 = 0;
    FldArrayI cEFTemp(sizeEF); E_Int* ptrEFTemp = cEFTemp.begin();
    FldArrayI cPHTemp(n);
    FldArrayI indirFaces(nfacesTot); indirFaces.setAllValuesAt(-1);
    E_Int* indirFacesp = indirFaces.begin();
    E_Int indface, indFaceOut, nfaces, pose, posf;
    E_Int nbFacesOut = 0;
    vector<E_Int> origIndicesOfFaces;

    for (E_Int noe = 0; noe < n; noe++)
    {
      // construction de la connectivite elt/faces
      pose = eltList[noe];
      E_Int* elt = cn->getElt(pose, nfaces, nface, indPH);
      ptrEFTemp[0] = nfaces;
      for (E_Int nof = 0; nof < nfaces; nof++)
      {
        indface = elt[nof]-1;
        if (indirFacesp[indface] == -1) 
        {
          indFaceOut = nbFacesOut; 
          indirFacesp[indface] = indFaceOut; 
          nbFacesOut++;
          origIndicesOfFaces.push_back(indface);
        }
        else indFaceOut = indirFacesp[indface];
        ptrEFTemp[nof+shift] = indFaceOut+1;
      }      
      ptrEFTemp += nfaces+shift; sizeEF2 += nfaces+shift; 
      cPHTemp[noe] = nfaces;
    }
    indirFaces.malloc(0); cEFTemp.resize(sizeEF2);

    // construction de la connectivite Faces/Noeuds
    FldArrayI cFNTemp(sizeFN); E_Int* ptrFNTemp = cFNTemp.begin();
    FldArrayI cPGTemp(nbFacesOut);
    FldArrayI indirNodes(npts); indirNodes.setAllValuesAt(-1); E_Int* indirNp = indirNodes.begin();
    E_Int indnode, nbnodes;
    E_Int nUniqueNodes = 0;
    for (E_Int nfe = 0; nfe < nbFacesOut; nfe++)
    {
      posf = origIndicesOfFaces[nfe]; //demarre a 0
      E_Int* face = cn->getFace(posf, nbnodes, ngon, indPG);
      ptrFNTemp[0] = nbnodes;
      for (E_Int p = 0; p < nbnodes; p++)
      {
        indnode = face[p]-1;
        if (indirNp[indnode] == -1) //creation
        {
          indirNp[indnode] = nUniqueNodes+1;
          ptrFNTemp[p+shift] = nUniqueNodes+1;
          nUniqueNodes++;
        }
        else 
        {
          ptrFNTemp[p+shift] = indirNp[indnode];
        }
      }
      ptrFNTemp += nbnodes+shift; sizeFN2 += nbnodes+shift;
      cPGTemp[nfe] = nbnodes;
    }
    origIndicesOfFaces.clear();

    // construit l'array de sortie
    PyObject* tpln = K_ARRAY::buildArray3(nfld, varString, nUniqueNodes, n, nbFacesOut,
                                          eltType, sizeFN2, sizeEF2, ngonType,
                                          false, api);
    FldArrayF* f2; FldArrayI* cn2;
    K_ARRAY::getFromArray3(tpln, f2, cn2);
    E_Int *ngon2 = cn2->getNGon(), *nface2 = cn2->getNFace();
    E_Int *indPG2 = NULL, *indPH2 = NULL;
    if (api == 2 || api == 3)
    {
      indPG2 = cn2->getIndPG(); indPH2 = cn2->getIndPH();
    }

    // reconstruction de la connectivite finale
    #pragma omp parallel default(shared)
    {
      #pragma omp for nowait
      for (E_Int i = 0; i < sizeFN2; i++) ngon2[i] = cFNTemp[i];
      #pragma omp for nowait
      for (E_Int i = 0; i < sizeEF2; i++) nface2[i] = cEFTemp[i];

      if (api == 2 || api == 3) // array2 or array3
      {
        #pragma omp for nowait
        for (E_Int i = 0; i < nbFacesOut; i++) indPG2[i] = cPGTemp[i];
        #pragma omp for nowait
        for (E_Int i = 0; i < n; i++) indPH2[i] = cPHTemp[i];
      }
    }

    PyObject* tplc = K_ARRAY::buildArray3(nfldc, varStringc, n, *cn2,
                                          eltTypec, 1, api, true);
    FldArrayF* fc2; K_ARRAY::getFromArray3(tplc, fc2);

    #pragma omp parallel default(shared)
    {
      E_Int indf;
      // Champs aux noeuds
      for (E_Int eq = 1; eq <= nfld; eq++)
      {
        E_Float* fp = f->begin(eq);
        E_Float* f2p = f2->begin(eq);
        #pragma omp for
        for (E_Int ind = 0; ind < npts; ind++)
        {
          indf = indirNp[ind]-1;
          if (indf > -1) f2p[indf] = fp[ind];
        }
      }
      // Champs aux centres
      for (E_Int eq = 1; eq <= nfldc; eq++) 
      {
        E_Float* fcp = fc->begin(eq);
        E_Float* fc2p = fc2->begin(eq);
        #pragma omp for
        for (E_Int ind = 0; ind < n; ind++)
        {
          indf = eltList[ind];
          fc2p[ind] = fcp[indf];
        }
      }
    }
    RELEASESHAREDU(tpln, f2, cn2); PyList_Append(l, tpln); Py_DECREF(tpln);
    RELEASESHAREDS(tplc, fc2); PyList_Append(l, tplc); Py_DECREF(tplc);
  }
  else
  {
    E_Int nc = cn->getNConnect();
    E_Int binIndex, nelts;
    E_Int nc2 = 0, elOffset = 0, nUniqueNodes = 0;
    vector<E_Int> neltspc(nc+1), neltspc2(nc), nvpe(nc);

    char* eltType2 = new char[50]; strcpy(eltType2, "");
    vector<char*> eltTypes;
    K_ARRAY::extractVars(eltType, eltTypes);
    
    // Compute the cumulative number of elements per connectivity,
    // these are the bounds of the bins that are used to inform on
    // the element type of an element index
    neltspc[0] = elOffset;
    for (E_Int ic = 0; ic < nc; ic++)
    {
      FldArrayI& cm = *(cn->getConnect(ic));
      nvpe[ic] = cm.getNfld();
      nelts = cm.getSize();
      elOffset += nelts;
      neltspc[ic+1] = elOffset;
    }

    // Bin input elements
    vector<vector<E_Int> > binnedEltList(nc);
    for (E_Int i = 0; i < n; i++)
    {
      E_Int noe = eltList[i];
      binIndex = std::upper_bound(neltspc.begin(), neltspc.end(), noe) - neltspc.begin() - 1;
      binnedEltList[binIndex].push_back(noe);
    }

    // Selectionne les elements subzones, calcule le nombre de noeuds uniques
    FldArrayI indirNodes(npts); indirNodes.setAllValuesAt(-1);
    E_Int* indirNodesp = indirNodes.begin();
    vector<E_Int> listOfNodes(npts);

    for (E_Int ic = 0; ic < nc; ic++)
    {
      // Skip empty connectivities
      if (binnedEltList[ic].size())
      {
        strcat(eltType2, eltTypes[ic]); // Build eltType2
        neltspc2[nc2] = binnedEltList[ic].size();
        
        FldArrayI& cm = *(cn->getConnect(ic));
        for (E_Int i = 0; i < neltspc2[nc2]; i++)
        {
          E_Int noe = binnedEltList[ic][i];
          for (E_Int v = 1; v <= nvpe[ic]; v++)
          {
            E_Int indv = cm(noe,v)-1;
            if (indirNodesp[indv] == -1)
            { 
              listOfNodes[nUniqueNodes] = indv; nUniqueNodes++;
              indirNodesp[indv] = nUniqueNodes;
            }
          }
        }

        nc2++;
      }
    }
    neltspc2.resize(nc2); listOfNodes.resize(nUniqueNodes);
    
    // Create connectivities
    PyObject* tpln = K_ARRAY::buildArray3(nfld, varString, nUniqueNodes, 
                                          neltspc2, eltType2, false, api);
    FldArrayF* f2; FldArrayI* cn2;
    K_ARRAY::getFromArray3(tpln, f2, cn2);
    
    #pragma omp parallel default(shared)
    {
      E_Int indv, noe;
      for (E_Int ic = 0; ic < nc; ic++)
      {
        #pragma omp for
        for (size_t i = 0; i < binnedEltList[ic].size(); i++)
        {
          if (!binnedEltList[ic].size()) continue;
          FldArrayI& cm = *(cn->getConnect(ic));
          FldArrayI& cm2 = *(cn2->getConnect(ic));

          noe = binnedEltList[ic][i];
          for (E_Int v = 1; v <= nvpe[ic]; v++)
          {
            indv = cm(noe,v)-1;
            cm2(i,v) = indirNodesp[indv];
          }
        }
      }
    }

    char* eltType2c = new char[50];
    K_ARRAY::starVarString(eltType2, eltType2c);
    PyObject* tplc = K_ARRAY::buildArray3(nfldc, varStringc, n, *cn2,
                                          eltType2c, 1, api, true);
    FldArrayF* fc2; K_ARRAY::getFromArray3(tplc, fc2);

    delete[] eltType2; delete[] eltType2c;
    for (E_Int ic = 0; ic < nc; ic++) delete [] eltTypes[ic];

    #pragma omp parallel default(shared)
    {  
      // Champs aux noeuds
      for (E_Int eq = 1; eq <= nfld; eq++)
      {
        E_Float* fp = f->begin(eq);
        E_Float* f2p = f2->begin(eq);
        #pragma omp for
        for (E_Int i = 0; i < nUniqueNodes; i++) f2p[i] = fp[listOfNodes[i]];    
      }

      // Champs aux centres
      for (E_Int eq = 1; eq <= nfldc; eq++)
      {
        E_Float* fcp = fc->begin(eq);
        E_Float* fc2p = fc2->begin(eq);
        #pragma omp for
        for (E_Int i = 0; i < n; i++) fc2p[i] = fcp[eltList[i]];    
      }
    }

    indirNodes.malloc(0);
    RELEASESHAREDU(tpln, f2, cn2); PyList_Append(l, tpln); Py_DECREF(tpln);
    RELEASESHAREDS(tplc, fc2); PyList_Append(l, tplc); Py_DECREF(tplc);
  }
  RELEASESHAREDU(arrayNodes, f, cn); RELEASESHAREDU(arrayCenters, fc, cnc);
  return l;
}

// ============================================================================
/* Subzone an unstructured mesh by face indices as global index */
// ============================================================================
PyObject* K_TRANSFORM::subzoneFaces(PyObject* self, PyObject* args)
{
  PyObject* array;
  PyObject* listOfFaces;
  if (!PYPARSETUPLE_(args, OO_, &array, &listOfFaces)) return NULL;

  // Build face list
  FldArrayI faceList;
  E_Int ret = K_ARRAY::getFromList(listOfFaces, faceList);
  if (ret == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "subzone: argument must be a list of face indices (starting from 1).");
    return NULL;
  }

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, im, jm, km, cn, eltType); 
  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "subzone: unknown type of array.");
    return NULL;
  }

  if (res == 1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "subzone: cannot be used on a structured array.");
    RELEASESHAREDS(array, f); return NULL;
  }
  
  E_Int n = faceList.getSize();
  E_Int* faceListp = faceList.begin();
  E_Int nfld = f->getNfld(); E_Int npts = f->getSize();
  E_Int api = f->getApi();

  if (K_STRING::cmp(eltType, "NGON") == 0) // IN: NGON, OUT: NGON
  {
    E_Int nbnodes, fidx;
    E_Int indedgen = 1;
    E_Int npts2 = 0; E_Int sizeEF2 = 0;
    
    // Acces non universel sur les ptrs
    E_Int* ngon = cn->getNGon();
    E_Int* indPG = cn->getIndPG();
    E_Int shift = 1;
    if (api == 3) shift = 0;

    // Calcul du nombre de points et aretes uniques dans la nouvelle
    // connectivite
    vector<E_Int> indirVertices(npts, -1);
    // Les aretes sont hashees pour determiner le nombre unique d'aretes et
    // construire la nouvelle connectivite 2D
    vector<E_Int> edge(2);
    std::pair<E_Int, E_Bool> initEdge(-1, false);
    TopologyOpt E;
    std::unordered_map<TopologyOpt, std::pair<E_Int, E_Bool>, JenkinsHash<TopologyOpt> > edgeMap;

    // Loop over all faces to extract
    for (E_Int i = 0; i < n; i++)
    {
      fidx = faceListp[i]-1;
      E_Int* face = cn->getFace(fidx, nbnodes, ngon, indPG);
      sizeEF2 += nbnodes+shift;

      // Loop over all edges of that face
      for (E_Int j = 0; j < nbnodes; j++)
      {
        edge[0] = face[j]-1; edge[1] = face[(j+1)%nbnodes]-1;
        E.set(edge.data(), 2);
        
        // If indirection table is -1 for the first vertex, then first time
        // this vertex is encountered. The second vertex will be dealt with
        // when considering the next edge
        if (indirVertices[edge[0]] == -1) indirVertices[edge[0]] = ++npts2;

        // If the value associated with E is -1, then first time this edge
        // is encountered, set to current unique edge count
        auto resE = edgeMap.insert(std::make_pair(E, initEdge));
        if (resE.first->second.first == -1)
        {
          resE.first->second.first = indedgen;
          indedgen++;
        }
      }
    }

    // Construction des nouvelles connectivites Elmt/Faces et Face/Noeuds
    E_Int nfaces2 = edgeMap.size();
    E_Int sizeFN2 = (2+shift)*nfaces2;
    E_Int nelts2 = n;
    
    E_Int ngonType = 1; // CGNSv3 compact array1
    if (api == 2) ngonType = 2; // CGNSv3, array2
    else if (api == 3) ngonType = 3; // force CGNSv4, array3
    E_Bool center = false;
    PyObject* tpl = K_ARRAY::buildArray3(nfld, varString, npts2, nelts2,
                                         nfaces2, "NGON", sizeFN2, sizeEF2,
                                         ngonType, center, api);
    FldArrayF* f2; FldArrayI* cn2;
    K_ARRAY::getFromArray3(tpl, f2, cn2);
    E_Int* ngon2 = cn2->getNGon();
    E_Int* nface2 = cn2->getNFace();
    E_Int *indPG2 = NULL, *indPH2 = NULL;
    if (api == 2 || api == 3)
    {
      indPG2 = cn2->getIndPG(); indPH2 = cn2->getIndPH();
    }

    E_Int c1 = 0, c2 = 0; // positions in ngon2 and nface2
    for (E_Int i = 0; i < nelts2; i++)
    {
      fidx = faceListp[i]-1;
      E_Int* face = cn->getFace(fidx, nbnodes, ngon, indPG);
      nface2[c2] = nbnodes;
      if (api == 2 || api == 3) indPH2[i] = nbnodes;

      for (E_Int p = 0; p < nbnodes; p++)
      {
        edge[0] = face[p]-1; edge[1] = face[(p+1)%nbnodes]-1;
        E.set(edge.data(), 2);
        // Find the edge in the edge map
        auto resE = edgeMap.find(E);
        if (not resE->second.second)
        {
          resE->second.second = true;
          ngon2[c1] = 2;
          ngon2[c1+shift] = indirVertices[edge[0]];
          ngon2[c1+1+shift] = indirVertices[edge[1]];
          c1 += 2+shift;
        }
        nface2[c2+p+shift] = resE->second.first;
      }
      c2 += nbnodes+shift;
    }

    #pragma omp parallel
    {
      E_Int indv;
      if (api == 2 || api == 3)
      {
        #pragma omp for
        for(E_Int i = 0; i < nfaces2; i++) indPG2[i] = 2;
      }
    
      for(E_Int eq = 1; eq <= nfld; eq++)
      {
        E_Float* fp = f->begin(eq);
        E_Float* f2p = f2->begin(eq);
        #pragma omp for
        for (E_Int ind = 0; ind < npts; ind++)
        {
          indv = indirVertices[ind]-1;
          if (indv > -1) f2p[indv] = fp[ind];
        }
      }
    }

    RELEASESHAREDU(array, f, cn);
    RELEASESHAREDU(tpl, f2, cn2);
    return tpl;
  }
  else // Basic faces
  {    
    E_Int ierr, vidx, fidx, eidx, fic, nfpe, nvpf, nfaces, nelts;
    
    vector<char*> eltTypes;
    K_ARRAY::extractVars(eltType, eltTypes);
    // facetspc: list of facets per connectivity
    E_Int nc = cn->getNConnect();
    vector<vector<vector<E_Int> > > facetspc(nc);

    // Get dimensionality
    E_Int dim = 3;
    if (strcmp(eltTypes[0], "NODE") == 0) dim = 0;
    else if (strcmp(eltTypes[0], "BAR") == 0) dim = 1;
    else if (strcmp(eltTypes[0], "TRI") == 0 ||
             strcmp(eltTypes[0], "QUAD") == 0) dim = 2;

    // Compute total number of faces and fill facets of ME
    E_Int nfacesTot = 0;
    // nfpc: number of faces per connectivity (cumulative)
    vector<E_Int> nfpc(nc+1); nfpc[0] = 0; 
    
    for (E_Int ic = 0; ic < nc; ic++)
    {
      FldArrayI& cm = *(cn->getConnect(ic));
      nelts = cm.getSize();
      // fill the list of facets for this BE
      ierr = K_CONNECT::getEVFacets(facetspc[ic], eltTypes[ic]);
      if (ierr == 1)
      {
        PyErr_SetString(PyExc_TypeError,
                        "subzone: element type not taken into account.");
        for (E_Int ic = 0; ic < nc; ic++) delete [] eltTypes[ic];
        RELEASESHAREDU(array, f, cn); return NULL;
      }

      // number of faces for this connectivity is the product of the number
      // of elements by the number of facets
      nfaces = nelts*facetspc[ic].size();
      nfpc[ic+1] = nfaces;
      nfacesTot += nfaces;
    }
    //std::cout << "nfacesTot: " << nfacesTot << std::endl;
    for (E_Int ic = 1; ic < nc+1; ic++) nfpc[ic] += nfpc[ic-1]; // cumulated

    // Faces are hashed to build the new 2D connectivity
    TopologyOpt F;
    vector<E_Int> face(4);
    struct FaceAttrs {
      E_Int ic_; E_Int eidx_; E_Int fidx_;

      FaceAttrs(E_Int ic, E_Int eidx, E_Int fidx):
        ic_(ic), eidx_(eidx), fidx_(fidx) {}
    };
    std::map<TopologyOpt, FaceAttrs> faceMap; // std::map for reproducibility

    // Loop over all faces to extract and fill the mapping table
    E_Int npts2 = 0;
    vector<E_Int> indirVertices(npts);

    #pragma omp parallel for if (npts > __MIN_SIZE_MEAN__)
    for (E_Int i = 0; i < npts; i++) indirVertices[i] = -1;

    for (E_Int i = 0; i < n; i++)
    {
      fidx = faceListp[i]-1; // global face indexing
      fic = -1; // connectivity associated to this face index
      for (E_Int j = 0; j < nc; j++)
        if (fidx < nfpc[j+1]) { fic = j; break; }
      if (fic == -1)
      {
        printf("Warning: subzone: face index " SF_D_ " is larger than the "
               "total number of faces: " SF_D_ ". Skipping...\n", fidx+1, nfacesTot);
        continue;
      }
      fidx -= nfpc[fic]; // face indexing relative to this connectivity
      nfpe = facetspc[fic].size();
      eidx = fidx/nfpe; fidx = fidx%nfpe; // relative to the element eidx

      FldArrayI& cm = *(cn->getConnect(fic));

      // Loop over vertices of the facet
      const vector<E_Int>& facet = facetspc[fic][fidx];
      nvpf = facet.size();
      for (E_Int j = 0; j < nvpf; j++)
      {
        vidx = cm(eidx, facet[j])-1;
        if (indirVertices[vidx] == -1) indirVertices[vidx] = ++npts2;
        face[j] = indirVertices[vidx];
      }
      F.set(face.data(), nvpf);
      faceMap.insert({F, FaceAttrs(fic, eidx, fidx)});
    }

    // Build new BE connectivity
    E_Int nelts2 = faceMap.size();
    E_Bool center = false;
    char eltTypeFaces[10];
    if (dim == 1) strcpy(eltTypeFaces, "NODE");
    else if (dim == 2) strcpy(eltTypeFaces, "BAR");
    else if (strcmp(eltTypes[0], "TETRA") == 0) strcpy(eltTypeFaces, "TRI");
    else strcpy(eltTypeFaces, "QUAD");
    PyObject* tpl = K_ARRAY::buildArray3(nfld, varString, npts2, nelts2,
                                         eltTypeFaces, center, api);
    FldArrayF* f2; FldArrayI* cn2;
    K_ARRAY::getFromArray3(tpl, f2, cn2);
    FldArrayI& cm2 = *(cn2->getConnect(0));

    // Copy connectivity to cn2
    eidx = 0;
    E_Int indv;
    for (const auto& face : faceMap)
    {
      nvpf = face.first.n_;
      FaceAttrs fattrs = face.second;
      FldArrayI& cm = *(cn->getConnect(fattrs.ic_));
      const vector<E_Int>& facet = facetspc[fattrs.ic_][fattrs.fidx_];
      
      for (E_Int j = 0; j < nvpf; j++)
      {
        indv = cm(fattrs.eidx_,facet[j])-1;
        cm2(eidx,j+1) = indirVertices[indv];
      }
      eidx++;
    }
    
    #pragma omp parallel if (npts > __MIN_SIZE_MEAN__)
    {
      E_Int indv;
      // Copy fields to f2
      for (E_Int eq = 1; eq <= nfld; eq++)
      {
        E_Float* fp = f->begin(eq);
        E_Float* f2p = f2->begin(eq);
        #pragma omp for
        for (E_Int i = 0; i < npts; i++)
        {
          indv = indirVertices[i]-1;
          if (indv > -1) f2p[indv] = fp[i];
        }
      }
    }

    for (E_Int ic = 0; ic < nc; ic++) delete [] eltTypes[ic];
    RELEASESHAREDU(tpl, f2, cn2);
    RELEASESHAREDU(array, f, cn);
    return tpl;
  }

  RELEASESHAREDU(array, f, cn);
  return NULL;
}

// ============================================================================
/* Subzone a unstructured mesh by face indices - fields located at centers
   are subzoned too */
// ============================================================================
PyObject* K_TRANSFORM::subzoneFacesBoth(PyObject* self, PyObject* args)
{
  PyObject *arrayNodes, *arrayCenters;
  PyObject* listOfFaces;
  if (!PYPARSETUPLE_(args, OOO_, &arrayNodes, &arrayCenters, &listOfFaces))
  {
    return NULL;
  }

  // Build face list
  FldArrayI faceList;
  E_Int ret = K_ARRAY::getFromList(listOfFaces, faceList);
  if (ret == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "subzone: argument must be a list of face indices (starting from 1).");
    return NULL;
  }

  // Check array of nodes
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(arrayNodes, varString, f,
                                     im, jm, km, cn, eltType); 
  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "subzone: unknown type of array.");
    return NULL;
  }
  if (res == 1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "subzone: cannot be used on a structured array.");
    RELEASESHAREDS(arrayNodes, f); return NULL;
  }
  
  // Check array of centers
  E_Int imc, jmc, kmc;
  FldArrayF* fc; FldArrayI* cnc;
  char* varStringc; char* eltTypec;
  E_Int resc = K_ARRAY::getFromArray3(arrayCenters, varStringc, fc,
                                      imc, jmc, kmc, cnc, eltTypec); 
  if (resc != 1 && resc != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "subzone: unknown type of array.");
    RELEASESHAREDU(arrayNodes, f, cn); return NULL;
  }
  if (resc == 1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "subzone: cannot be used on a structured array.");
    RELEASESHAREDU(arrayNodes, f, cn); RELEASESHAREDS(arrayCenters, fc); return NULL;
  }

  E_Int n = faceList.getSize();
  E_Int* faceListp = faceList.begin();
  E_Int nfld = f->getNfld(); E_Int nfldc = fc->getNfld();
  E_Int npts = f->getSize(); E_Int api = f->getApi();
  PyObject* l = PyList_New(0);

  if (K_STRING::cmp(eltType, "NGON") == 0) // IN: NGON, OUT: NGON
  {
    E_Int nbnodes, fidx, nf;
    E_Int indedgen = 1;
    E_Int npts2 = 0; E_Int sizeEF2 = 0;
    
    // Acces non universel sur les ptrs
    E_Int* ngon = cn->getNGon(); E_Int* nface = cn->getNFace();
    E_Int* indPG = cn->getIndPG(); E_Int* indPH = cn->getIndPH();
    E_Int shift = 1;
    E_Int nelts = cn->getNElts(); E_Int nfaces = cn->getNFaces();
    if (api == 3) shift = 0;

    // Calcul du nombre de points et aretes uniques dans la nouvelle
    // connectivite
    vector<E_Int> indirVertices(npts, -1);
    // Les aretes sont hashees pour determiner le nombre unique d'aretes et
    // construire la nouvelle connectivite 2D
    vector<E_Int> edge(2);
    std::pair<E_Int, E_Bool> initEdge(-1, false);
    TopologyOpt E;
    std::unordered_map<TopologyOpt, std::pair<E_Int, E_Bool>, JenkinsHash<TopologyOpt> > edgeMap;

    // Loop over all faces to extract
    for (E_Int i = 0; i < n; i++)
    {
      fidx = faceListp[i]-1;
      E_Int* face = cn->getFace(fidx, nbnodes, ngon, indPG);
      sizeEF2 += nbnodes+shift;

      // Loop over all edges of that face
      for (E_Int j = 0; j < nbnodes; j++)
      {
        edge[0] = face[j]-1; edge[1] = face[(j+1)%nbnodes]-1;
        E.set(edge.data(), 2);
        
        // If indirection table is -1 for the first vertex, then first time
        // this vertex is encountered. The second vertex will be dealt with
        // when considering the next edge
        if (indirVertices[edge[0]] == -1) indirVertices[edge[0]] = ++npts2;

        // If the value associated with E is -1, then first time this edge
        // is encountered, set to current unique edge count
        auto resE = edgeMap.insert(std::make_pair(E, initEdge));
        if (resE.first->second.first == -1)
        {
          resE.first->second.first = indedgen;
          indedgen++;
        }
      }
    }

    // Construction des nouvelles connectivites Elmt/Faces et Face/Noeuds
    E_Int nfaces2 = edgeMap.size();
    E_Int sizeFN2 = (2+shift)*nfaces2;
    E_Int nelts2 = n;
    
    E_Int ngonType = 1; // CGNSv3 compact array1
    if (api == 2) ngonType = 2; // CGNSv3, array2
    else if (api == 3) ngonType = 3; // force CGNSv4, array3
    E_Bool center = false;
    PyObject* tpl = K_ARRAY::buildArray3(nfld, varString, npts2, nelts2,
                                         nfaces2, "NGON", sizeFN2, sizeEF2,
                                         ngonType, center, api);
    FldArrayF* f2; FldArrayI* cn2;
    K_ARRAY::getFromArray3(tpl, f2, cn2);
    E_Int* ngon2 = cn2->getNGon();
    E_Int* nface2 = cn2->getNFace();
    E_Int *indPG2 = NULL, *indPH2 = NULL;
    if (api == 2 || api == 3)
    {
      indPG2 = cn2->getIndPG(); indPH2 = cn2->getIndPH();
    }

    E_Int c1 = 0, c2 = 0; // positions in ngon2 and nface2
    for (E_Int i = 0; i < nelts2; i++)
    {
      fidx = faceListp[i]-1;
      E_Int* face = cn->getFace(fidx, nbnodes, ngon, indPG);
      nface2[c2] = nbnodes;
      if (api == 2 || api == 3) indPH2[i] = nbnodes;

      for (E_Int p = 0; p < nbnodes; p++)
      {
        edge[0] = face[p]-1; edge[1] = face[(p+1)%nbnodes]-1;
        E.set(edge.data(), 2);
        // Find the edge in the edge map
        auto resE = edgeMap.find(E);
        if (not resE->second.second)
        {
          resE->second.second = true;
          ngon2[c1] = 2;
          ngon2[c1+shift] = indirVertices[edge[0]];
          ngon2[c1+1+shift] = indirVertices[edge[1]];
          c1 += 2+shift;
        }
        nface2[c2+p+shift] = resE->second.first;
      }
      c2 += nbnodes+shift;
    }

    // Build parent elements array to compute face-centered data from
    // cell-centered data
    FldArrayI parentElts(nfaces, 2); parentElts.setAllValuesAtNull();
    E_Int* PE1 = parentElts.begin(1); E_Int* PE2 = parentElts.begin(2);
    for (E_Int i = 0; i < nelts; i++)
    {
      E_Int* elem = cn->getElt(i, nf, nface, indPH);
      for (E_Int j = 0; j < nf; j++)
      {
        E_Int ind = elem[j]-1;
        if (PE1[ind] == 0) PE1[ind] = i+1;
        else PE2[ind] = i+1;
      }
    }

    FldArrayF* fc2 = new FldArrayF(nelts2, nfldc);
    #pragma omp parallel
    {
      E_Int indv, indf, etg, etd;
      if (api == 2 || api == 3)
      {
        #pragma omp for
        for(E_Int i = 0; i < nfaces2; i++) indPG2[i] = 2;
      }
    
      for(E_Int eq = 1; eq <= nfld; eq++)
      {
        E_Float* fp = f->begin(eq);
        E_Float* f2p = f2->begin(eq);
        
        #pragma omp for
        for (E_Int ind = 0; ind < npts; ind++)
        {
          indv = indirVertices[ind]-1;
          if (indv > -1) f2p[indv] = fp[ind];
        }
      }

      for(E_Int eq = 1; eq <= nfldc; eq++)
      {
        E_Float* fcp = fc->begin(eq);
        E_Float* fc2p = fc2->begin(eq);

        #pragma omp for
        for (E_Int ind = 0; ind < nelts2; ind++)
        {
          indf = faceListp[ind]-1;
          etg = PE1[indf]; etd = PE2[indf];
          if (etd > 0) fc2p[ind] = 0.5*(fcp[etg-1] + fcp[etd-1]);
          else fc2p[ind] = fcp[etg-1];        
        }
      }
    }

    parentElts.malloc(0);
    RELEASESHAREDU(arrayNodes, f, cn); RELEASESHAREDU(arrayCenters, fc, cnc);

    // Build array
    PyList_Append(l, tpl);
    FldArrayI* cnc2 = new FldArrayI(); *cnc2 = *cn2;
    PyObject* tplc = K_ARRAY::buildArray3(*fc2, varStringc, *cnc2, "NGON", api);
    PyList_Append(l, tplc); Py_DECREF(tplc); delete fc2; delete cnc2;
    RELEASESHAREDU(tpl, f2, cn2);
    return l;
  }
  else // Basic faces
  {
    E_Int ierr, vidx, fidx, eidx, fic, nfpe, nvpf, nfaces, nelts;
    
    vector<char*> eltTypes;
    K_ARRAY::extractVars(eltType, eltTypes);
    // facetspc: list of facets per connectivity
    E_Int nc = cn->getNConnect();
    vector<vector<vector<E_Int> > > facetspc(nc);

    // Get dimensionality
    E_Int dim = 3;
    if (strcmp(eltTypes[0], "NODE") == 0) dim = 0;
    else if (strcmp(eltTypes[0], "BAR") == 0) dim = 1;
    else if (strcmp(eltTypes[0], "TRI") == 0 ||
             strcmp(eltTypes[0], "QUAD") == 0) dim = 2;

    // Compute total number of faces and fill facets of ME
    E_Int nfacesTot = 0;
    // nfpc: number of faces per connectivity (cumulative)
    vector<E_Int> nfpc(nc+1); nfpc[0] = 0; 
    
    for (E_Int ic = 0; ic < nc; ic++)
    {
      FldArrayI& cm = *(cn->getConnect(ic));
      nelts = cm.getSize();
      // fill the list of facets for this BE
      ierr = K_CONNECT::getEVFacets(facetspc[ic], eltTypes[ic]);
      if (ierr == 1)
      {
        PyErr_SetString(PyExc_TypeError,
                        "subzone: element type not taken into account.");
        for (E_Int ic = 0; ic < nc; ic++) delete [] eltTypes[ic];
        RELEASESHAREDU(arrayNodes, f, cn); RELEASESHAREDU(arrayCenters, fc, cnc);
        return NULL;
      }

      // number of faces for this connectivity is the product of the number
      // of elements by the number of facets
      nfaces = nelts*facetspc[ic].size();
      nfpc[ic+1] = nfaces;
      nfacesTot += nfaces;
    }
    for (E_Int ic = 1; ic < nc+1; ic++) nfpc[ic] += nfpc[ic-1]; // cumulated

    // Build parent elements array to compute face-centered data from
    // cell-centered data
    vector<vector<E_Int> > parentElts(2);
    parentElts[0].resize(nfacesTot); parentElts[1].resize(nfacesTot);
    auto pe1 = parentElts[0].begin();
    auto pe2 = parentElts[1].begin();

    // Faces are hashed to build the new 2D connectivity
    TopologyOpt F;
    vector<E_Int> face(4);
    struct FaceAttrs {
      E_Int ic_; E_Int eidx_; E_Int fidx_; E_Int glbFidx_;

      FaceAttrs(E_Int ic, E_Int eidx, E_Int fidx, E_Int glbFidx):
        ic_(ic), eidx_(eidx), fidx_(fidx), glbFidx_(glbFidx) {}
    };
    //std::unordered_map<TopologyOpt, FaceAttrs, JenkinsHash<TopologyOpt> > faceMap;
    std::map<TopologyOpt, FaceAttrs> faceMap; // std::map for reproducibility

    // Loop over all faces to extract and fill the indirection table
    E_Int npts2 = 0;
    vector<E_Int> indirVertices(npts);

    #pragma omp parallel if (npts > __MIN_SIZE_MEAN__)
    {
      #pragma omp for
      for (E_Int i = 0; i < npts; i++) indirVertices[i] = -1;
      #pragma omp for
      for (E_Int i = 0; i < nfacesTot; i++) { pe1[i] = -1; pe2[i] = -1; }
    }

    // Loop over all faces to select
    for (E_Int i = 0; i < n; i++)
    {
      E_Int indf = faceListp[i]-1; // global face indexing
      fic = -1; // connectivity associated to this face index
      for (E_Int j = 0; j < nc; j++)
        if (indf < nfpc[j+1]) { fic = j; break; }
      if (fic == -1)
      {
        printf("Warning: subzone: face index " SF_D_ " is larger than the "
               "total number of faces: " SF_D_ ". Skipping...\n", indf+1, nfacesTot);
        continue;
      }
      fidx = indf - nfpc[fic]; // face indexing relative to this connectivity
      nfpe = facetspc[fic].size();
      eidx = fidx/nfpe; fidx = fidx%nfpe; // relative to the element eidx

      FldArrayI& cm = *(cn->getConnect(fic));

      // Loop over vertices of the facet
      const vector<E_Int>& facet = facetspc[fic][fidx];
      nvpf = facet.size();
      for (E_Int j = 0; j < nvpf; j++)
      {
        vidx = cm(eidx, facet[j])-1;
        if (indirVertices[vidx] == -1) { indirVertices[vidx] = ++npts2; }
        face[j] = indirVertices[vidx];
      }
      F.set(face.data(), nvpf);
      auto resF = faceMap.insert({F, FaceAttrs(fic, eidx, fidx, indf)});
      if (resF.second) pe1[indf] = eidx+1;
      else pe2[resF.first->second.glbFidx_] = eidx+1;
    }

    // Build new BE connectivity
    E_Int nelts2 = faceMap.size();
    E_Bool center = false;
    char eltTypeFaces[10];
    if (dim == 1) strcpy(eltTypeFaces, "NODE");
    else if (dim == 2) strcpy(eltTypeFaces, "BAR");
    else if (strcmp(eltTypes[0], "TETRA") == 0) strcpy(eltTypeFaces, "TRI");
    else strcpy(eltTypeFaces, "QUAD");
    PyObject* tpl = K_ARRAY::buildArray3(nfld, varString, npts2, nelts2,
                                         eltTypeFaces, center, api);
    FldArrayF* f2; FldArrayI* cn2;
    K_ARRAY::getFromArray3(tpl, f2, cn2);
    FldArrayI& cm2 = *(cn2->getConnect(0));

    // Copy connectivity to cn2
    eidx = 0;
    E_Int indv;
    for (const auto& face : faceMap)
    {
      nvpf = face.first.n_;
      FaceAttrs fattrs = face.second;
      FldArrayI& cm = *(cn->getConnect(fattrs.ic_));
      const vector<E_Int>& facet = facetspc[fattrs.ic_][fattrs.fidx_];
      
      for (E_Int j = 0; j < nvpf; j++)
      {
        indv = cm(fattrs.eidx_,facet[j])-1;
        cm2(eidx,j+1) = indirVertices[indv];
      }
      eidx++;
    }

    // To fill fc2, continue populating parentElement by looping over all
    // remaining faces
    // Convert faceList to an unordered_set for fast lookups
    std::unordered_set<E_Int> selectedFacesSet(faceList.begin(), faceList.end());
    for (E_Int i = 0; i < nfacesTot; i++)
    {
      if (selectedFacesSet.find(i+1) == selectedFacesSet.end()) // new face
      {
        fic = -1; // connectivity associated to this face index
        for (E_Int j = 0; j < nc; j++)
          if (i < nfpc[j+1]) { fic = j; break; }
        if (fic == -1)
        {
          printf("Warning: subzone: face index " SF_D_ " is larger than the "
                "total number of faces: " SF_D_ ". Skipping...\n", i+1, nfacesTot);
          continue;
        }
        fidx = i - nfpc[fic]; // face indexing relative to this connectivity
        nfpe = facetspc[fic].size();
        eidx = fidx/nfpe; fidx = fidx%nfpe; // relative to the element eidx

        FldArrayI& cm = *(cn->getConnect(fic));

        // Loop over vertices of the facet
        const vector<E_Int>& facet = facetspc[fic][fidx];
        nvpf = facet.size();
        for (E_Int j = 0; j < nvpf; j++)
        {
          vidx = cm(eidx, facet[j])-1;
          if (indirVertices[vidx] == -1) { indirVertices[vidx] = -(++npts2); }
          face[j] = K_FUNC::E_abs(indirVertices[vidx]);
        }
        F.set(face.data(), nvpf);
        auto resF = faceMap.insert({F, FaceAttrs(fic, eidx, fidx, i)});
        if (resF.second) pe1[i] = eidx+1;
        else pe2[resF.first->second.glbFidx_] = eidx+1;
      }
    }

    FldArrayF* fc2 = new FldArrayF(nelts2, nfldc);
    #pragma omp parallel if (npts > __MIN_SIZE_MEAN__)
    {
      E_Int indv;
      // Copy nodal fields to f2
      for (E_Int eq = 1; eq <= nfld; eq++)
      {
        E_Float* fp = f->begin(eq);
        E_Float* f2p = f2->begin(eq);
        #pragma omp for
        for (E_Int i = 0; i < npts; i++)
        {
          indv = indirVertices[i]-1;
          if (indv > -1) f2p[indv] = fp[i];
        }
      }
    }

    E_Int indf, etg, etd;
    E_Int j = 0;
    for (E_Int eq = 1; eq <= nfldc; eq++)
    {
      E_Float* fcp = fc->begin(eq);
      E_Float* fc2p = fc2->begin(eq);

      for (E_Int i = 0; i < n; i++)
      {
        indf = faceListp[i]-1;
        etg = pe1[indf];
        if (etg > 0) // to skip duplicated faces in the input
        {
          etd = pe2[indf];
          if (etd > 0) fc2p[j] = 0.5*(fcp[etg-1] + fcp[etd-1]);
          else fc2p[j] = fcp[etg-1];
          j++;
        }
      }
    }

    parentElts.clear();
    for (E_Int ic = 0; ic < nc; ic++) delete [] eltTypes[ic];
    RELEASESHAREDU(arrayNodes, f, cn); RELEASESHAREDU(arrayCenters, fc, cnc);

    // Build array
    PyList_Append(l, tpl);
    FldArrayI* cnc2 = new FldArrayI(); *cnc2 = *cn2;
    PyObject* tplc = K_ARRAY::buildArray3(*fc2, varStringc, *cnc2, eltTypeFaces, api);
    PyList_Append(l, tplc); Py_DECREF(tplc); delete fc2; delete cnc2;
    RELEASESHAREDU(tpl, f2, cn2);
    return l;
  }

  RELEASESHAREDU(arrayNodes, f, cn); RELEASESHAREDU(arrayCenters, fc, cnc); 
  return NULL;
}

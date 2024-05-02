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

// subzone 

# include "transform.h"
using namespace std;
using namespace K_FLD;

// ============================================================================
/* Subzone a structured mesh */
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
      printf("Warning: subzone: mesh dimensions are: %d x %d x %d\n", im,jm,km);
      PyErr_SetString(PyExc_TypeError,
                      "subzone: wrong choice of index.");
      return NULL;
    }
    
    // Construit l'array resultat
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
    std::vector<char*> eltTypes;
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
    std::vector<std::vector<E_Int> > binnedEltList(nc);
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
    std::vector<char*> eltTypes;
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
    std::vector<std::vector<E_Int> > binnedEltList(nc);
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
  RELEASESHAREDU(arrayNodes, f, cn); RELEASESHAREDS(arrayCenters, fc);
  return l;
}
// ============================================================================
/* Subzone an unstructured mesh by face indices */
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
    RELEASESHAREDS(array,f); return NULL;
  }
  
  E_Int n = faceList.getSize();
  E_Int* faceListp = faceList.begin();
  E_Int nf, nt2;

  char eltTypeFaces[10];
  int fBAR[] = { 
    1, 
    2 };
  int fTRI[] = {
    1, 2,
    2, 3,
    3, 1 };
  int fQUAD[] = {
    1, 2,
    2, 3,
    3, 4,
    4, 1 };
  int fTETRA[] = {
    1, 3, 2,
    1, 2, 4,
    2, 3, 4,
    3, 1, 4 };
  int fPENTA[] = {
    1, 2, 5, 4,
    2, 3, 6, 5,
    3, 1, 4, 6,
    1, 3, 2, 1,
    4, 5, 6, 4
  };
  int fPYRA[] = {
    1, 4, 3, 2,
    1, 2, 5, 5,
    2, 3, 5, 2,
    3, 4, 5, 3,
    4, 1, 5, 4
  };  
  int fHEXA[] = {
    1, 4, 3, 2,
    1, 2, 6, 5,
    2, 3, 7, 6,
    3, 4, 8, 7,
    1, 5, 8, 4,
    5, 6, 7, 8 };
  
  int* fPtr = NULL; E_Int type = 0;
  nt2 = 0; nf = 0;
  // nf: nbre de face de l'element d'entree
  // eltTypeFaces: type de l'element de sortie (face de l'elt d'entree)
  // nt2: nbre de noeuds de l'element de sortie
  // type=0 (les faces de sortie sont a elements basiques
  // type=1 (entree: elements basiques, sortie: NGON)
  if (K_STRING::cmp(eltType, "BAR") == 0)
  {
    type = 0; nf = 2; strcpy(eltTypeFaces, "NODE"); nt2 = 1; fPtr = fBAR;
  }
  else if (K_STRING::cmp(eltType, "TRI") == 0)
  {
    type = 0; nf = 3; strcpy(eltTypeFaces, "BAR"); nt2 = 2; fPtr = fTRI;
  }
  else if (K_STRING::cmp(eltType, "QUAD") == 0)
  {
    type = 0; nf = 4; strcpy(eltTypeFaces, "BAR"); nt2 = 2; fPtr = fQUAD;
  }
  else if (K_STRING::cmp(eltType, "TETRA") == 0)
  {
    type = 0; nf = 4; strcpy(eltTypeFaces, "TRI"); nt2 = 3; fPtr = fTETRA;
  }
  else if (K_STRING::cmp(eltType, "HEXA") == 0)
  {
    type = 0; nf = 6; strcpy(eltTypeFaces, "QUAD"); nt2 = 4; fPtr = fHEXA;
  }
  else if (K_STRING::cmp(eltType, "PENTA") == 0)
  {
    type = 0; nf = 5; strcpy(eltTypeFaces, "QUAD"); nt2 = 4; fPtr = fPENTA;
  }
  else if (K_STRING::cmp(eltType, "PYRA") == 0)
  {
    type = 0; nf = 5; strcpy(eltTypeFaces, "QUAD"); nt2 = 4; fPtr = fPYRA; 
  }
  else if (K_STRING::cmp(eltType, "NGON") == 0)
  {
    type = 2; strcpy(eltTypeFaces, "NGON");
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "subzone: element type not taken into account.");
    RELEASESHAREDU(array, f, cn); return NULL;
  }
  E_Int nfld = f->getNfld(); E_Int npts = f->getSize();
  if (type == 0) // Basic faces
  {    
    FldArrayI* c2n = new FldArrayI(n,nt2);
    // Selectionne les faces subzonees
    FldArrayI indirNodes(npts); indirNodes.setAllValuesAt(-1);
    E_Int* indirNp = indirNodes.begin();

    E_Int ind, indElt, nof, indnode, indnodef;
    indnodef = 0;
    for (E_Int i = 0; i < n; i++)
    {
      ind = faceListp[i]-1;
      indElt = ind / nf;
      nof = ind - indElt*nf;
      // printf("elt=%d face=%d\n", indElt, nof);
      for (E_Int v = 1; v <= nt2; v++)
      { 
        indnode = (*cn)(indElt, fPtr[nof*nt2+v-1])-1;
        if (indirNp[indnode] == -1) 
        {
          indnodef += 1;
          indirNp[indnode] = indnodef;
          (*c2n)(i,v) = indnodef;
        }
        else 
          (*c2n)(i,v) = indirNp[indnode];
      }
    }
    
    FldArrayF* f2 = new FldArrayF(indnodef, nfld);
    E_Int indf;
    for(E_Int eq = 1; eq <= nfld; eq++)
    {
      E_Float* fp = f->begin(eq);
      E_Float* fnp = f2->begin(eq);
      for (E_Int ind = 0; ind < npts; ind++)
      {
        indf = indirNp[ind]-1;
        if (indf > -1) fnp[indf] = fp[ind];
      }
    }
    // Build array 
    PyObject* tpl = K_ARRAY::buildArray(*f2, varString, 
                                        *c2n, -1, eltTypeFaces);
    RELEASESHAREDU(array,f, cn);
    delete f2; delete c2n;
    return tpl;
  }

  else if (type == 2) // IN: NGON, OUT: NGON
  {
    E_Int* ptr = cn->begin();
    //E_Int sizeFN = ptr[1]; ptr += 2; 
    //E_Int ne = ptr[sizeFN]; // ncnmbre d'elements pour la connectivite cn
    E_Int next1=0; // nbre de faces pour la nouvelle connectivite 
    E_Int next2=0; // nbre d'elements pour la nouvelle connectivite
    E_Int sizeco1=0; // taille de la nouvelle connectivite de faces 
    E_Int sizeco2=0; // taille de la nouvelle connectivite d'elements 
    E_Int nbnodes;

    E_Int pos, ind;
    FldArrayI posFace;
    K_CONNECT::getPosFaces(*cn, posFace);
    E_Int* posFacep = posFace.begin();
    for (E_Int i = 0; i < n; i++)
    {
      ind = faceListp[i]-1;
      pos = posFacep[ind];
      ptr = cn->begin()+pos;
      nbnodes = ptr[0];
      sizeco1 += 3*nbnodes; next1 += nbnodes;
      sizeco2 += nbnodes+1; next2++;
    }

    FldArrayI c2n(sizeco1+sizeco2+4);
    E_Int* coutp = c2n.begin();
    E_Int* ptro1 = coutp+2;
    E_Int* ptro2 = coutp+4+sizeco1;
    FldArrayI indirNodes(npts); indirNodes.setAllValuesAt(-1);
    E_Int* indirNp = indirNodes.begin();
    E_Int indnewface = 1;
    E_Int indvertp1, indvertp2;
    E_Int indvertn = 0;
    vector< vector<E_Int> > indicesFaces(npts); // faces creees associees aux noeuds 

    for (E_Int i = 0; i < n; i++)
    {
      ind = faceListp[i]-1;
      pos = posFacep[ind];
      ptr = cn->begin()+pos;
      nbnodes = ptr[0];
      
      ptro2[0] = nbnodes;
      for (E_Int p = 0; p < nbnodes; p++)
      { 
        ptro1[0] = 2;
        indvertp1 = ptr[p+1]-1;
        indvertp2 = ptr[(p+1)%(nbnodes)+1]-1;
        vector<E_Int>& facesp1 = indicesFaces[indvertp1];
        vector<E_Int>& facesp2 = indicesFaces[indvertp2];

        // creation de la face P1P2 
        E_Int foundFace = -1; 
        if (facesp1.size() > 0 && facesp2.size() > 0) 
        {
          for (size_t fi = 0; fi < facesp1.size(); fi++)
            for (size_t fj = 0; fj < facesp2.size(); fj++)
            {
              if (facesp1[fi] == facesp2[fj] ) 
              {
                foundFace = facesp1[fi]; 
                break;
              }
              if (foundFace != -1) break;
            }
        }
        if (foundFace == -1)
        { 
          if (indirNp[indvertp1] == -1) 
          {
            indvertn += 1;
            indirNp[indvertp1] = indvertn;
            ptro1[1] = indvertn;
          }
          else 
          {
            ptro1[1] = indirNp[indvertp1];
          }
        
          if (indirNp[indvertp2] == -1) 
          {
            indvertn += 1;
            indirNp[indvertp2] = indvertn;
            ptro1[2] = indvertn;
          }
          else 
          {
            ptro1[2] = indirNp[indvertp2];
          }
          ptro1 += 3;

          facesp1.push_back(indnewface); facesp2.push_back(indnewface);
          ptro2[p+1] = indnewface; indnewface++;
        }
        else 
        {
          ptro2[p+1] = foundFace;
        }
      }
      ptro2 += nbnodes+1;
    }
    FldArrayF* f2 = new FldArrayF(indvertn, nfld);

#pragma omp parallel default(shared)
    {
      E_Int indf;
      for (E_Int eq = 1; eq <= nfld; eq++)
      {
        E_Float* fp = f->begin(eq);
        E_Float* fnp = f2->begin(eq);
#pragma omp for
        for (E_Int ind = 0; ind < npts; ind++)
        {
          indf = indirNp[ind]-1;
          if (indf > -1) fnp[indf] = fp[ind];
        }
      }
    }
    // Cree la nouvelle connectivite complete
    E_Int nbFaces = indnewface-1;
    E_Int sizeFN = nbFaces*3;
    E_Int nbElts = next2;
    E_Int sizeEF = sizeco2;
    FldArrayI* cnout = new FldArrayI(sizeFN+sizeEF+4);
    E_Int* cnoutp = cnout->begin();
    ptro1 = c2n.begin()+2;
    ptro2 = c2n.begin()+4+sizeco1;
    cnoutp[0] = nbFaces; // nombre de faces
    cnoutp[1] = sizeFN; // taille du tableau de faces
    cnoutp += 2;
    E_Int* cnoute = cnoutp + sizeFN;
    cnoute[0] = nbElts; 
    cnoute[1] = sizeEF;
    cnoute += 2;

#pragma omp parallel default(shared)
    {
#pragma omp for
      for (E_Int i = 0; i < sizeFN; i++) cnoutp[i] = ptro1[i];
#pragma omp for
      for (E_Int i = 0; i < sizeEF; i++) cnoute[i] = ptro2[i];
    }
    // Build array 
    PyObject* tpl = K_ARRAY::buildArray(*f2, varString, 
                                        *cnout, -1, eltTypeFaces);
    RELEASESHAREDU(array,f, cn);
    delete f2; delete cnout;
    return tpl;
  }
  else if (type == 1)
    PyErr_SetString(PyExc_TypeError,
                    "subzone: not implemented.\n");

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
  if (!PYPARSETUPLE_(args, OOO_,
                        &arrayNodes, &arrayCenters, &listOfFaces))
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
  E_Int res = K_ARRAY::getFromArray(arrayNodes, varString, f, im, jm, km, cn, eltType, true); 
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
    RELEASESHAREDS(arrayNodes,f); return NULL;
  }
  
  // Check array of centers
  E_Int imc, jmc, kmc;
  FldArrayF* fc; FldArrayI* cnc;
  char* varStringc; char* eltTypec;
  E_Int resc = K_ARRAY::getFromArray(arrayCenters, varStringc, fc, imc, jmc, kmc, cnc, eltTypec, true); 
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

  E_Int n = faceList.getSize();
  E_Int* faceListp = faceList.begin();
  E_Int nf, nt2;

  char eltTypeFaces[10];
  int fBAR[] = { 
    1, 
    2 };
  int fTRI[] = {
    1, 2,
    2, 3,
    3, 1 };
  int fQUAD[] = {
    1, 2,
    2, 3,
    3, 4,
    4, 1 };
  int fTETRA[] = {
    1, 3, 2,
    1, 2, 4,
    2, 3, 4,
    3, 1, 4 };

  int fPENTA[] = {
    1, 2, 5, 4,
    2, 3, 6, 5,
    3, 1, 4, 6,
    1, 3, 2, 1,
    4, 5, 6, 4
  };
  int fPYRA[] = {
    1, 4, 3, 2,
    1, 2, 5, 1,
    2, 3, 5, 2,
    3, 4, 5, 3,
    4, 1, 5, 4
  };  

  int fHEXA[] = {
    1, 4, 3, 2,
    1, 2, 6, 5,
    2, 3, 7, 6,
    3, 4, 8, 7,
    1, 5, 8, 4,
    5, 6, 7, 8 };
  
  int* fPtr = NULL; E_Int type = 0;
  nt2 = 0; nf = 0;
  // nf: nbre de faces de l'element d'entree
  // eltTypeFaces: type de l'element de sortie (face de l'elt d'entree)
  // nt2: nbre de noeuds de l'element de sortie
  // type=0 (les faces de sortie est a elements basiques
  // type=1 (entree: elements basiques, sortie: NGONS)
  // type=2 (entree: NGON, sortie: NGON)
  if (K_STRING::cmp(eltType, "BAR") == 0)
  {
    type = 0; nf = 2; strcpy(eltTypeFaces, "NODE"); nt2 = 1; fPtr = fBAR;
  }
  else if (K_STRING::cmp(eltType, "TRI") == 0)
  {
    type = 0; nf = 3; strcpy(eltTypeFaces, "BAR"); nt2 = 2; fPtr = fTRI;
  }
  else if (K_STRING::cmp(eltType, "QUAD") == 0)
  {
    type = 0; nf = 4; strcpy(eltTypeFaces, "BAR"); nt2 = 2; fPtr = fQUAD;
  }
  else if (K_STRING::cmp(eltType, "TETRA") == 0)
  {
    type = 0; nf = 4; strcpy(eltTypeFaces, "TRI"); nt2 = 3; fPtr = fTETRA;
  }
  else if (K_STRING::cmp(eltType, "HEXA") == 0)
  {
    type = 0; nf = 6; strcpy(eltTypeFaces, "QUAD"); nt2 = 4; fPtr = fHEXA;
  }
  else if (K_STRING::cmp(eltType, "PENTA") == 0)
  {
    //type = 1; nf = 5; strcpy(eltTypeFaces, "NGON");
    type = 0; nf = 5; strcpy(eltTypeFaces, "QUAD"); nt2 = 4; fPtr = fPENTA;
  }
  else if (K_STRING::cmp(eltType, "PYRA") == 0)
  {
    //type = 1; nf = 5; strcpy(eltTypeFaces, "NGON");
    type = 0; nf = 5; strcpy(eltTypeFaces, "QUAD"); nt2 = 4; fPtr = fPYRA; 
  }
  else if (K_STRING::cmp(eltType, "NGON") == 0)
  {
    type = 2; strcpy(eltTypeFaces, "NGON");
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "subzone: element type not taken into account.");
    RELEASESHAREDU(arrayNodes,f, cn); RELEASESHAREDU(arrayCenters, fc, cnc); return NULL;
  }
  E_Int nfld = f->getNfld(); E_Int npts = f->getSize();
  E_Int nfldc = fc->getNfld();
  PyObject* l = PyList_New(0);

  if (type == 0) // Basic faces
  { 
    E_Int nvert = cn->getNfld();
    FldArrayF* fc2 = new FldArrayF(n, nfldc);
    FldArrayI* c2n = new FldArrayI(n,nt2);
    // Selectionne les faces subzonees
    FldArrayI indirNodes(npts); indirNodes.setAllValuesAt(-1);
    E_Int* indirNp = indirNodes.begin();
    vector< vector<E_Int> > cVE(npts);
    K_CONNECT::connectEV2VE(*cn, cVE);
    E_Int ind, indElt, nof, indnode, indnodef, sizeVE1, eltV, eltOpp, match, indnode2;
    indnodef = 0;
    for (E_Int i = 0; i < n; i++)
    {
      ind = faceListp[i]-1;
      indElt = ind / nf;
      nof = ind - indElt*nf;
      // printf("elt=%d face=%d\n", indElt, nof);
      
      // Determination du 2e elt partageant cette face
      // on prend le 1er noeud de la face + recup des elts contenant ce noeud
      indnode = (*cn)(indElt, fPtr[nof*nt2])-1;
      vector<E_Int>& cVE1 = cVE[indnode];
      sizeVE1 = cVE1.size();
      eltOpp = -1;
      for (E_Int ve = 0; ve < sizeVE1; ve++)
      {
        eltV = cVE1[ve]; 
        if ( eltV != indElt)
        {
          match = 0;
          for (E_Int v = 1; v <= nt2; v++)
          {
            indnode = (*cn)(indElt,fPtr[nof*nt2+v-1])-1;
            for (E_Int nodei = 1; nodei <= nvert; nodei++)
            {
              indnode2 = (*cn)(eltV, nodei)-1;
              if (indnode == indnode2) { match++; break; }
            }
          }
          if (match == nt2){ eltOpp = eltV; break;}
        }
      }
      if ( eltOpp == -1) //pas d elt voisin trouve
      {
        for (E_Int eq = 1; eq <= nfldc; eq++)
          (*fc2)(i,eq) = (*fc)(indElt,eq); 
      }
      else 
      {
        for (E_Int eq = 1; eq <= nfldc; eq++)
          (*fc2)(i,eq) = 0.5*((*fc)(indElt,eq)+(*fc)(eltOpp,eq)); 
      }
      // Creation de la connectivite pour la face creee
      for (E_Int v = 1; v <= nt2; v++)
      { 
        indnode = (*cn)(indElt, fPtr[nof*nt2+v-1])-1;
        if (indirNp[indnode] == -1) 
        {
          indnodef += 1;
          indirNp[indnode] = indnodef;
          (*c2n)(i,v) = indnodef;
        }
        else 
          (*c2n)(i,v) = indirNp[indnode];
      }
    }
    
    FldArrayF* f2 = new FldArrayF(indnodef, nfld);
    E_Int indf;
    for(E_Int eq = 1; eq <= nfld; eq++)
    {
      E_Float* fp = f->begin(eq);
      E_Float* fnp = f2->begin(eq);
      for (E_Int ind = 0; ind < npts; ind++)
      {
        indf = indirNp[ind]-1;
        if (indf > -1) fnp[indf] = fp[ind];
      }
    }
    // Build array 
    PyObject* tpln = K_ARRAY::buildArray(*f2, varString, 
                                         *c2n, -1, eltTypeFaces);
    PyList_Append(l, tpln); Py_DECREF(tpln); delete f2;
    RELEASESHAREDU(arrayNodes,f, cn);RELEASESHAREDU(arrayCenters, fc, cnc); 
    FldArrayI* cncout = new FldArrayI(); *cncout = *c2n;
    PyObject* tplc = K_ARRAY::buildArray(*fc2, varStringc, 
                                         *cncout, -1, eltTypeFaces, true);
    PyList_Append(l, tplc); Py_DECREF(tplc); delete fc2;  delete cncout; delete c2n;
    return l;
  }
  else if (type == 2) // IN: NGON, OUT: NGON
  {
    E_Int* ptr = cn->begin();
    E_Int nfacesTot = ptr[0];
    E_Int sizeFN = ptr[1]; ptr += 2; 
    E_Int ne = ptr[sizeFN]; // nombre d'elements pour la connectivite cn
    E_Int nfaces2=0; // nbre de faces pour la nouvelle connectivite 
    E_Int nelts2=0; // nbre d'elements pour la nouvelle connectivite
    E_Int sizeFN2=0; // taille de la nouvelle connectivite de faces 
    E_Int sizeEF2=0; // taille de la nouvelle connectivite d'elements 

    FldArrayI posFace;
    K_CONNECT::getPosFaces(*cn, posFace);
    E_Int* posFacep = posFace.begin();
    E_Int ind, pos, nf, nbnodes;
    for (E_Int i = 0; i < n; i++)
    {
      ind = faceListp[i]-1;
      pos = posFacep[ind];
      ptr = cn->begin()+pos;
      nbnodes = ptr[0];
      sizeFN2 += 3*nbnodes; nfaces2 += nbnodes;
      sizeEF2 += nbnodes+1; nelts2++;
    }

    FldArrayI c2n(sizeFN2+sizeEF2+4);
    E_Int* coutp = c2n.begin();
    E_Int* ptro1 = coutp+2;
    E_Int* ptro2 = coutp+4+sizeFN2;
    FldArrayI indirNodes(npts); indirNodes.setAllValuesAt(-1);
    E_Int* indirNp = indirNodes.begin();
    E_Int indnewface = 1;
    E_Int indvertp1, indvertp2;
    E_Int indvertn = 0;
    vector< vector<E_Int> > indicesFaces(npts);// faces creees associees aux noeuds 
    
    E_Int* ptrNF = cn->begin()+sizeFN+4;
    FldArrayI parentElts(nfacesTot,2); parentElts.setAllValuesAtNull();
    E_Int* PE1 = parentElts.begin(1);
    E_Int* PE2 = parentElts.begin(2);
    for (E_Int i = 0; i < ne; i++)
    {
      nf = ptrNF[0];
      for (E_Int j = 1; j <= nf; j++)
      {
        ind = ptrNF[j]-1;
        if (PE1[ind] == 0) PE1[ind] = i+1;
        else PE2[ind] = i+1;
      }
      ptrNF += nf+1;
    }
    FldArrayF* fc2 = new FldArrayF(n, nfldc);

#pragma omp parallel default(shared)
    {
      E_Int indf, etg, etd;
      for (E_Int eq = 1; eq <= nfldc; eq++)
      {
        E_Float* fp2l = fc2->begin(eq);
        E_Float* fpl = fc->begin(eq);
#pragma omp for
        for (E_Int i = 0; i < n; i++)
        {
          indf = faceListp[i]-1;
          etg = PE1[indf]; etd = PE2[indf];
          if (etg > 0 && etd > 0 )
            fp2l[i] = 0.5*(fpl[etg-1]+fpl[etd-1]);
          else if (etg==0 && etd>0)
            fp2l[i] = fpl[etd-1];        
          else if (etg>0 && etd==0) 
            fp2l[i] = fpl[etg-1];        
        }
      }
    }
    parentElts.malloc(0);

    for (E_Int i = 0; i < n; i++)
    {
      ind = faceListp[i]-1;
      pos = posFacep[ind];
      ptr = cn->begin()+pos;
      nbnodes = ptr[0];
      
      ptro2[0] = nbnodes;
      for (E_Int p = 0; p < nbnodes; p++)
      { 
        ptro1[0] = 2;
        indvertp1 = ptr[p+1]-1;
        indvertp2 = ptr[(p+1)%(nbnodes)+1]-1;
        vector<E_Int>& facesp1 = indicesFaces[indvertp1];
        vector<E_Int>& facesp2 = indicesFaces[indvertp2];

        //creation de la face P1P2 
        E_Int foundFace = -1; 
        if (facesp1.size() > 0 && facesp2.size() > 0) 
        {
          for (size_t fi = 0; fi < facesp1.size(); fi++)
            for (size_t fj = 0; fj < facesp2.size(); fj++)
            {
              if (facesp1[fi] == facesp2[fj] ) 
              {
                foundFace = facesp1[fi]; 
                break;
              }
              if (foundFace != -1) break;
            }        
        }
        if (foundFace == -1)
        { 
          if (indirNp[indvertp1] == -1) 
          {
            indvertn += 1;
            indirNp[indvertp1] = indvertn;
            ptro1[1] = indvertn;
          }
          else 
          {
            ptro1[1] = indirNp[indvertp1];
          }
        
          if (indirNp[indvertp2] == -1) 
          {
            indvertn += 1;
            indirNp[indvertp2] = indvertn;
            ptro1[2] = indvertn;
          }
          else 
          {
            ptro1[2] = indirNp[indvertp2];
          }
          ptro1 += 3;

          facesp1.push_back(indnewface); facesp2.push_back(indnewface);
          ptro2[p+1] = indnewface; indnewface++;
        }
        else 
        {
          ptro2[p+1] = foundFace;
        }
      }
      ptro2 += nbnodes+1;
    }
    posFace.malloc(0);
    for (E_Int i = 0; i < npts; i++)
      indicesFaces[i].clear();
    indicesFaces.clear();
    FldArrayF* f2 = new FldArrayF(indvertn, nfld);

#pragma omp parallel default(shared)
    {
      E_Int indf;
      for (E_Int eq = 1; eq <= nfld; eq++)
      {
        E_Float* fp = f->begin(eq);
        E_Float* fnp = f2->begin(eq);
#pragma omp for
        for (E_Int ind = 0; ind < npts; ind++)
        {
          indf = indirNp[ind]-1;
          if (indf > -1) fnp[indf] = fp[ind];
        }
      }
    }
    // Cree la nouvelle connectivite complete
    E_Int nbFaces = indnewface-1;
    sizeFN = nbFaces*3;
    E_Int nbElts = nelts2;
    E_Int sizeEF = sizeEF2;
    FldArrayI* cnout = new FldArrayI(sizeFN+sizeEF+4);
    E_Int* cnoutp = cnout->begin();
    ptro1 = c2n.begin()+2;
    ptro2 = c2n.begin()+4+sizeFN2;
    cnoutp[0] = nbFaces; // nombre de faces
    cnoutp[1] = sizeFN; // taille du tableau de faces
    cnoutp += 2;
    E_Int* cnoute = cnoutp + sizeFN;
    cnoute[0] = nbElts; 
    cnoute[1] = sizeEF;
    cnoute += 2;

#pragma omp parallel default(shared)
    {
#pragma omp for
      for (E_Int i = 0; i < sizeFN; i++) cnoutp[i] = ptro1[i];
#pragma omp for
      for (E_Int i = 0; i < sizeEF; i++) cnoute[i] = ptro2[i];
    }
    // Build array 
    PyObject* tpln = K_ARRAY::buildArray(*f2, varString, 
                                         *cnout, -1, eltTypeFaces);
    PyList_Append(l, tpln); Py_DECREF(tpln); delete f2;
    RELEASESHAREDU(arrayNodes,f, cn); RELEASESHAREDU(arrayCenters, fc, cnc); 
    FldArrayI* cncout = new FldArrayI(); *cncout = *cnout;
    PyObject* tplc = K_ARRAY::buildArray(*fc2, varStringc, 
                                         *cncout, -1, eltTypeFaces, true);
    PyList_Append(l, tplc); Py_DECREF(tplc); delete fc2; 
    delete cnout; delete cncout;
    return l;
  }
  else if (type == 1)
    PyErr_SetString(PyExc_TypeError,
                    "subzone: not implemented.\n");

  RELEASESHAREDU(arrayNodes,f, cn); RELEASESHAREDU(arrayCenters, fc, cnc); 
  return NULL;
}

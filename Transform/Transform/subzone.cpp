/*    
    Copyright 2013-2019 Onera.

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
  if (!PYPARSETUPLEI(args,
                    "O(lll)(lll)", "O(iii)(iii)",
                    &array, &i1, &j1, &k1, &i2, &j2, &k2))
  {
      return NULL;
  }
  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = 
    K_ARRAY::getFromArray(array, varString, f, im, jm, km, cn, eltType, true); 
  E_Int imjm = im*jm;

  if (res == 1)
  { 
    // Negative -> max indices
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
    PyObject* tpl= K_ARRAY::buildArray(nfld, varString, in, jn, kn);
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
  if (!PyArg_ParseTuple(args, "OO",
                        &array, &listOfNodes))
  {
      return NULL;
  }

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray(array, varString, f, im, jm, km, cn, eltType, true); 
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
  if (strcmp(eltType,"NGON") == 0) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "subzone: type='nodes' not implemented for a NGON array.");
    RELEASESHAREDU(array,f, cn); return NULL;
  }
  FldArrayI indices;
  E_Int ok = K_ARRAY::getFromList(listOfNodes, indices);
  if (ok == 0) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "subzone: 2nd argument must be an integer list or a numpy.");
    RELEASESHAREDU(array,f,cn); return NULL;
  }
  E_Int n = indices.getSize();

  // Mapping f -> f2
  E_Int npts = f->getSize();
  E_Int* indicesp = indices.begin();
  FldArrayI tmap(npts); tmap.setAllValuesAt(-1); E_Int* tmapP = tmap.begin();
  for (E_Int i = 0; i < n; i++) tmapP[indicesp[i]-1] = i;

  E_Int nfld = f->getNfld();

  if (strcmp(eltType, "NODE") == 0) 
  {
    PyObject* tpl = K_ARRAY::buildArray(nfld, varString, n, 0, -1, eltType);
    E_Float* fnp = K_ARRAY::getFieldPtr(tpl);

    FldArrayF f2(n, nfld, fnp, true);
    for (E_Int v = 1; v <= nfld; v++)
    {
      E_Float* fpl = f->begin(v);
      E_Float* f2pl = f2.begin(v);
      for (E_Int i = 0; i < n; i++)
      {
        E_Int ind = indicesp[i]-1;
        f2pl[i] = fpl[ind];
      }
    }
    RELEASESHAREDU(array,f,cn);
    return tpl;
  }
  // Selectionne les elements subzones
  E_Int nelts = cn->getSize();
  E_Int nvert = cn->getNfld();
  E_Int compt = 0; E_Int indv;
  vector<E_Int> vectOfElts;

  for (E_Int et = 0; et < nelts; et++)
  {
    compt = 0;
    for (E_Int nov = 1; nov <= nvert; nov++)
    {
      indv = (*cn)(et,nov)-1;
      if (tmapP[indv] != -1) compt++; 
    }
    if (compt == nvert) vectOfElts.push_back(et);
    //if (compt > 0) vectOfElts.push_back(et); //CB: extraction des qu'un node
  }
  
  E_Int nelts2 = vectOfElts.size();
  FldArrayI tag(n); // tag les points deja inseres dans un element
  tag.setAllValuesAtNull();
  E_Int* tagp = tag.begin();

  E_Int csize = nelts2*nvert; 
  PyObject* tpl = K_ARRAY::buildArray(nfld, varString,
                                      n, nelts2,-1, eltType, false, csize);
  E_Float* fnp = K_ARRAY::getFieldPtr(tpl);
  FldArrayF f2(n, nfld, fnp, true);
  for (E_Int v = 1; v <= nfld; v++) 
  {
    E_Float* fpl = f->begin(v);
    E_Float* f2pl = f2.begin(v);
    for (E_Int i = 0; i < n; i++)
    {
      E_Int ind = indicesp[i]-1;
      f2pl[i] = fpl[ind];
    }
  }
  E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
  FldArrayI cn2(nelts2, nvert, cnnp, true);

  E_Int indp, indn;
  for (E_Int nov = 1; nov <= nvert; nov++)
  {
    E_Int* cnl = cn->begin(nov);
    E_Int* cnl2 = cn2.begin(nov);
    for (E_Int noet = 0; noet < nelts2; noet++)
    {  
      indp = cnl[vectOfElts[noet]]-1;
      indn = tmapP[indp]; tagp[indn]++;
      cnl2[noet] = indn+1;
    }
  }
  for (E_Int i = 0; i < n; i++)
  {
    if (tagp[i] == 0)
    {
      PyErr_SetString(PyExc_TypeError,
                      "subzoneUnstruct: indices for unstructured subzone must be contiguous.");
      RELEASESHAREDU(array,f,cn); 
      return NULL;
    }
  }
  RELEASESHAREDU(array,f,cn);
  return tpl;
}
// ============================================================================
/* Subzone a unstructured zone: input array is located here at centers       */
// ============================================================================
PyObject* K_TRANSFORM::subzoneUnstructBoth(PyObject* self, PyObject* args)
{
  PyObject *arrayNodes, *arrayCenters;
  PyObject* listOfNodes;
  if (!PyArg_ParseTuple(args, "OOO",
                        &arrayNodes, &arrayCenters, &listOfNodes))
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
  E_Int n = indices.getSize();

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
  if (strcmp(eltType,"NGON") == 0) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "subzone: type='nodes' not implemented for a NGON array.");
    RELEASESHAREDU(arrayNodes,f, cn); return NULL;
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
  if (K_STRING::cmp(eltTypec,"NGON*") == 0) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "subzone: type='nodes' not implemented for a NGON array.");
    RELEASESHAREDU(arrayNodes,f, cn); RELEASESHAREDU(arrayCenters, fc, cnc); 
    return NULL;
  }

  E_Int npts = f->getSize();
  E_Int* indicesp = indices.begin();
  FldArrayI tmap(npts); tmap.setAllValuesAt(-1); E_Int* tmapP = tmap.begin();
  for (E_Int i = 0; i < n; i++) tmapP[indicesp[i]-1] = i;
  E_Int nfld = f->getNfld();
  E_Int nfldc = fc->getNfld();
  E_Int nelts = fc->getSize();
  if (cnc->getSize() != cn->getSize())
  {
    PyErr_SetString(PyExc_TypeError,
                    "subzone: sizes of connectivities at nodes and centers are not equal.");
    RELEASESHAREDU(arrayNodes,f, cn); RELEASESHAREDU(arrayCenters, fc, cnc); 
    return NULL;
  }
  PyObject* l = PyList_New(0);

  if (strcmp(eltType, "NODE") == 0) 
  {
    FldArrayF* fout = new FldArrayF(n,nfld);
    FldArrayI* cnout = new FldArrayI(); 
    for (E_Int v = 1; v <= nfld; v++) 
    {
      E_Float* fpl = f->begin(v);
      E_Float* f2pl = fout->begin(v);

      for (E_Int i = 0; i < n; i++)
      {
        E_Int ind = indicesp[i]-1;
        f2pl[i] = fpl[ind];
      }
    }
    PyObject* tpl1 = K_ARRAY::buildArray(*fout, varString,*cnout, -1, eltType);
    PyList_Append(l,tpl1); Py_DECREF(tpl1);
    delete fout; delete cnout;

    FldArrayF* fcout = new FldArrayF(n,nfldc);
    FldArrayI* cncout = new FldArrayI();
    for (E_Int v = 1; v <= nfldc; v++) 
    {
      E_Float* fcpl = fc->begin(v);
      E_Float* fc2pl = fcout->begin(v);
      for (E_Int i = 0; i < n; i++)
      {
        E_Int ind = indicesp[i]-1;
        fc2pl[i] = fcpl[ind];
      }
    }
    PyObject* tpl2 = K_ARRAY::buildArray(*fcout, varStringc,*cncout,-1, eltTypec, true);
    PyList_Append(l,tpl2); Py_DECREF(tpl2);
    delete fcout; delete cncout;

    RELEASESHAREDU(arrayNodes,f, cn); RELEASESHAREDU(arrayCenters, fc, cnc);  
    return l;
  }
  if (nelts != cn->getSize())
  {
    PyErr_SetString(PyExc_TypeError,
                    "subzone: arrays located at nodes and centers are not consistent.");
    RELEASESHAREDU(arrayNodes,f, cn); RELEASESHAREDU(arrayCenters, fc, cnc); 
    return NULL;
  }
  // Selectionne les elements subzones
  E_Int nvert = cn->getNfld();
  E_Int compt = 0;
  vector<E_Int> vectOfElts;
  for (E_Int et = 0; et < nelts; et++)
  {
    compt = 0;
    for (E_Int nov = 1; nov <= nvert; nov++)
    {
      E_Int indv = (*cn)(et,nov)-1;
      if (tmapP[indv] != -1) compt++; 
    }
    if (compt == nvert) vectOfElts.push_back(et);
  }
  E_Int nelts2 = vectOfElts.size();
  FldArrayI tag(n); // tag les points deja inseres dans un element
  tag.setAllValuesAtNull();
  E_Int* tagp = tag.begin();

  E_Int csize = nelts2*nvert; 
  PyObject* tpl = K_ARRAY::buildArray(nfld, varString,
                                      n, nelts2,-1, eltType, false, csize);
  E_Float* fnp = K_ARRAY::getFieldPtr(tpl);
  FldArrayF f2(n, nfld, fnp, true);
  for (E_Int v = 1; v <= nfld; v++) 
  {
    E_Float* fpl = f->begin(v);
    E_Float* f2pl = f2.begin(v);
    for (E_Int i = 0; i < n; i++)
    {
      E_Int ind = indicesp[i]-1;
      f2pl[i] = fpl[ind];
    }
  }
  E_Int* cnnp = K_ARRAY::getConnectPtr(tpl);
  FldArrayI cn2(nelts2, nvert, cnnp, true);

  for (E_Int nov = 1; nov <= nvert; nov++)
  {
    E_Int* cnl = cn->begin(nov);
    E_Int* cnl2 = cn2.begin(nov);
    for (E_Int noet = 0; noet < nelts2; noet++)
    {  
      E_Int indp = cnl[vectOfElts[noet]]-1;
      E_Int indn = tmapP[indp]; tagp[indn]++;
      cnl2[noet] = indn+1;
    }
  }
  for (E_Int i = 0; i < n; i++)
  {
    if (tagp[i] == 0)
    {
      PyErr_SetString(PyExc_TypeError,
                      "subzone: indices for unstructured subzone must be contiguous.");
      RELEASESHAREDU(arrayNodes,f,cn); RELEASESHAREDU(arrayCenters,fc,cnc); 
      return NULL;
    }
  }
  // Centers
  PyObject* tplc = K_ARRAY::buildArray(nfldc, varStringc,
                                       nelts2, nelts2, -1, eltTypec);
  E_Int* cnnpc = K_ARRAY::getConnectPtr(tplc);
  K_KCORE::memcpy__(cnnpc, cn2.begin(), nelts2*nvert);
  E_Float* fcnp = K_ARRAY::getFieldPtr(tplc);
  FldArrayF fc2(nelts2, nfldc, fcnp, true);
  for (E_Int v = 1; v <= nfldc; v++) 
  {
    E_Float* fpl = fc->begin(v);
    E_Float* f2pl = fc2.begin(v);
    for (E_Int noet = 0; noet < nelts2; noet++)
    {
      f2pl[noet] = fpl[vectOfElts[noet]];
    }
  }
  PyList_Append(l,tpl); Py_DECREF(tpl);
  PyList_Append(l,tplc); Py_DECREF(tplc);
  RELEASESHAREDU(arrayNodes,f,cn);RELEASESHAREDU(arrayCenters,fc,cnc);
  return l;
}
// ============================================================================
/* Subzone an unstructured mesh by element indices */
// ============================================================================
PyObject* K_TRANSFORM::subzoneElements(PyObject* self, PyObject* args)
{
  PyObject* array;
  PyObject* listOfElts;
  if (!PyArg_ParseTuple(args, "OO",
                        &array, &listOfElts))
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

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray(array, varString, f, im, jm, km, cn, eltType, true); 
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
  E_Int npts = f->getSize(); E_Int nfld = f->getNfld();
  E_Int n = eltList.getSize();
  E_Int* eltListp = eltList.begin();
  E_Int nt = cn->getNfld();
  PyObject* tpl = NULL;
  if (K_STRING::cmp(eltType, "NGON") == 0) // NGONS
  {
    FldArrayI posFace;
    K_CONNECT::getPosFaces(*cn, posFace);
    E_Int* posFacep = posFace.begin();
    FldArrayI posElt;
    K_CONNECT::getPosElts(*cn, posElt);

    E_Int* ptr = cn->begin();
    E_Int nfacesTot = ptr[0];
    E_Int sizeFN = ptr[1]; ptr += sizeFN+2;
    //E_Int neltsTot = ptr[0];
    E_Int sizeEF = ptr[1];
    E_Int* ptrEF = NULL;
    E_Int* ptrFN = NULL;

    FldArrayI cEFTemp(sizeEF); E_Int* ptrEFTemp = cEFTemp.begin();
    E_Int sizeEF2 = 0;
    E_Int sizeFN2 = 0;
    FldArrayI indirFaces(nfacesTot); indirFaces.setAllValuesAt(-1); E_Int* indirFacesp = indirFaces.begin();
    E_Int indface, indFaceOut, nfaces, pose, posf, elt;
    E_Int nbFacesOut = 0;
    vector<E_Int> origIndicesOfFaces;

    for (E_Int noe = 0; noe < n; noe++)
    {
      elt = eltListp[noe];
      pose = posElt[elt];
      
      //construction de la connectivite elt/faces
      ptrEF = cn->begin()+pose;
      nfaces = ptrEF[0];
      ptrEFTemp[0] = nfaces;
      for (E_Int nof = 0; nof < nfaces; nof++)
      {
        indface = ptrEF[nof+1]-1;
        if (indirFacesp[indface] == -1) 
        {
          posf = posFacep[indface];
          ptrFN = cn->begin()+posf;
          sizeFN2 += ptrFN[0]+1;
          
          indFaceOut = nbFacesOut; 
          indirFacesp[indface] = indFaceOut; 
          nbFacesOut++;
          origIndicesOfFaces.push_back(indface);
        }
        else indFaceOut = indirFacesp[indface];
        ptrEFTemp[nof+1] = indFaceOut+1;
      }      
      ptrEFTemp += nfaces+1; sizeEF2 += nfaces+1; 
    }
    indirFaces.malloc(0); posElt.malloc(0); cEFTemp.resize(sizeEF2);

    // construction de la connectivite Faces/Noeuds
    FldArrayI cFNTemp(sizeFN); E_Int* ptrFNTemp = cFNTemp.begin(); 
    FldArrayI indirNodes(npts); indirNodes.setAllValuesAt(-1); E_Int* indirNp = indirNodes.begin();
    E_Int indnode, indfaceo, nbnodes;
    E_Int numNode = 0;
    for (E_Int nfe = 0 ; nfe < nbFacesOut; nfe++)
    {
      indfaceo = origIndicesOfFaces[nfe];//demarre a 0
      posf = posFacep[indfaceo];
      ptrFN = cn->begin()+posf;
      nbnodes = ptrFN[0];
      ptrFNTemp[0] = nbnodes;
      for (E_Int p = 1; p <= nbnodes; p++)
      {
        indnode = ptrFN[p]-1;
        if (indirNp[indnode] == -1) //creation
        {
          indirNp[indnode] = numNode+1;
          ptrFNTemp[p] = numNode+1;
          numNode++;
        }
        else 
        {
          ptrFNTemp[p] = indirNp[indnode];
        }
      }
      ptrFNTemp += nbnodes+1;
    }
    posFace.malloc(0); origIndicesOfFaces.clear();      

    // construit l'array de sortie
    E_Int csize = sizeFN2+sizeEF2+4;
    char eltType[10]; strcpy(eltType,"NGON");
    tpl= K_ARRAY::buildArray(nfld, varString, numNode, csize, 8, eltType, 
                             false, csize);
    
    E_Float* fnp = K_ARRAY::getFieldPtr(tpl);
    FldArrayF f2(numNode,nfld, fnp, true);
  
#pragma omp parallel default(shared)
    {
      E_Int indf;
      for (E_Int eq = 1; eq <= nfld; eq++)
      {
        E_Float* fp = f->begin(eq);
        E_Float* fnp = f2.begin(eq);
#pragma omp for
        for (E_Int ind = 0; ind < npts; ind++)
        {
          indf = indirNp[ind]-1;
          if (indf > -1) fnp[indf] = fp[ind];
        }
      }
    }
    // reconstruction de la connectivite finale
    E_Int* cnp = K_ARRAY::getConnectPtr(tpl);
    cnp[0] = nbFacesOut; cnp[1] = sizeFN2;
    cnp += 2;
    E_Int* cne = cnp + sizeFN2;
    cne[0] = n; 
    cne[1] = sizeEF2;
    cne += 2;
    E_Int* ptro1 = cFNTemp.begin();  
    E_Int* ptro2 = cEFTemp.begin(); 
    
#pragma omp parallel default(shared)
    {
#pragma omp for
      for (E_Int i = 0; i < sizeFN2; i++) cnp[i] = ptro1[i];
#pragma omp for
      for (E_Int i = 0; i < sizeEF2; i++) cne[i] = ptro2[i];
    }
  }
  else // maillage par elements
  {
    // Pour les elts simples, on prend seulement la connectivite des elements en entree
    FldArrayI* c2n = new FldArrayI(n, nt);
    E_Int et, indv;
    FldArrayI indirNodes(npts); indirNodes.setAllValuesAt(-1); E_Int* indirNodesp = indirNodes.begin();
    E_Int nbnodes = 0;
    vector<E_Int> listOfNodes;
    
    // Selectionne les elements subzones
    for (E_Int i = 0; i < n; i++)
    {
      et = eltListp[i];
      for (E_Int v = 1; v <= nt; v++)
      {
        indv = (*cn)(et,v)-1;
        if (indirNodesp[indv] == -1) 
        { 
          listOfNodes.push_back(indv);
          nbnodes++;
          indirNodesp[indv] = nbnodes;
          (*c2n)(i,v) = nbnodes;
        }
        else  
          (*c2n)(i,v) = indirNodesp[indv];
      }
    }
  
    indirNodes.malloc(0);
    FldArrayF* fout = new FldArrayF(nbnodes,nfld);
#pragma omp parallel default(shared)
    {
      for (E_Int eq = 1; eq <= nfld; eq++)
      {
        E_Float* fp = f->begin(eq);
        E_Float* fp2 = fout->begin(eq);
#pragma omp for
        for (E_Int i = 0; i < nbnodes; i++)
          fp2[i] = fp[listOfNodes[i]];    
      }
    }
    listOfNodes.clear();
    
    // Build array 
    tpl = K_ARRAY::buildArray(*fout, varString, *c2n, -1, eltType);
    delete fout; delete c2n;
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
  if (!PyArg_ParseTuple(args, "OOO",
                        &arrayNodes, &arrayCenters, &listOfElts))
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
                    "subzone: can not be used on a structured array.");
    RELEASESHAREDS(arrayNodes, f); return NULL;
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

  E_Int npts = f->getSize(); 
  E_Int nfld = f->getNfld();
  E_Int nfldc = fc->getNfld();
  E_Int n = eltList.getSize();
  E_Int* eltListp = eltList.begin();
  E_Int nt = cn->getNfld();
  PyObject* l = PyList_New(0);

  if (K_STRING::cmp(eltType, "NGON") == 0) // NGONS
  {
    FldArrayI posFace;
    K_CONNECT::getPosFaces(*cn, posFace);
    E_Int* posFacep = posFace.begin();
    FldArrayI posElt;
    K_CONNECT::getPosElts(*cn, posElt);

    E_Int* ptr = cn->begin();
    E_Int nfacesTot = ptr[0];
    E_Int sizeFN = ptr[1]; ptr += sizeFN+2;
    //E_Int neltsTot = ptr[0];
    E_Int sizeEF = ptr[1];
    E_Int* ptrEF = NULL;
    E_Int* ptrFN = NULL;

    FldArrayI cEFTemp(sizeEF); E_Int* ptrEFTemp = cEFTemp.begin(); 
    E_Int sizeEF2 = 0;
    E_Int sizeFN2 = 0;
    FldArrayI indirFaces(nfacesTot); indirFaces.setAllValuesAt(-1); E_Int* indirFacesp = indirFaces.begin();
    E_Int indface, indFaceOut, nfaces, pose, posf, elt;
    E_Int nbFacesOut = 0;
    vector<E_Int> origIndicesOfFaces;
    for (E_Int noe = 0; noe < n; noe++)
    {
      elt = eltListp[noe];
      pose = posElt[elt];
      
      //construction de la connectivite elt/faces
      ptrEF = cn->begin()+pose;
      nfaces = ptrEF[0];
      ptrEFTemp[0] = nfaces;
      for (E_Int nof = 0; nof < nfaces; nof++)
      {
        indface = ptrEF[nof+1]-1;
        if (indirFacesp[indface] == -1) 
        {
          posf = posFacep[indface];
          ptrFN = cn->begin()+posf;
          sizeFN2 += ptrFN[0]+1;
          
          indFaceOut = nbFacesOut; 
          indirFacesp[indface]=indFaceOut; 
          nbFacesOut++;
          origIndicesOfFaces.push_back(indface);
        }
        else indFaceOut = indirFacesp[indface];
        ptrEFTemp[nof+1] = indFaceOut+1;
      }      
      ptrEFTemp+= nfaces+1; sizeEF2 += nfaces+1; 
    }
    indirFaces.malloc(0);posElt.malloc(0); cEFTemp.resize(sizeEF2);

    // construction de la connectivite Faces/Noeuds
    FldArrayI cFNTemp(sizeFN); E_Int* ptrFNTemp = cFNTemp.begin(); 
    FldArrayI indirNodes(npts); indirNodes.setAllValuesAt(-1); E_Int* indirNp = indirNodes.begin();
    E_Int indnode, indfaceo, nbnodes;
    E_Int numNode = 0;
    for (E_Int nfe = 0 ; nfe < nbFacesOut; nfe++)
    {
      indfaceo = origIndicesOfFaces[nfe];//demarre a 0
      posf = posFacep[indfaceo];
      ptrFN = cn->begin()+posf;
      nbnodes = ptrFN[0];
      ptrFNTemp[0] = nbnodes;
      for (E_Int p = 1; p <= nbnodes; p++)
      {
        indnode = ptrFN[p]-1;
        if ( indirNp[indnode] == -1) //creation
        {
          indirNp[indnode] = numNode+1;
          ptrFNTemp[p] = numNode+1;
          numNode++;
        }
        else 
        {
          ptrFNTemp[p] = indirNp[indnode];
        }
      }
      ptrFNTemp+= nbnodes+1;
    }
    posFace.malloc(0); origIndicesOfFaces.clear();      

    // construit l array de sortie
    E_Int csize = sizeFN2+sizeEF2+4;
    char eltType[10]; strcpy(eltType,"NGON");
    FldArrayF* f2 = new FldArrayF(numNode,nfld);
  
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
    // reconstruction de la connectivite finale
    FldArrayI* c2n = new FldArrayI(csize);
    E_Int* cnp = c2n->begin();
    cnp[0] = nbFacesOut; cnp[1] = sizeFN2; cnp += 2;
    E_Int* cne = cnp + sizeFN2;
    cne[0] = n; cne[1] = sizeEF2; cne += 2;
    E_Int* ptro1 = cFNTemp.begin();  
    E_Int* ptro2 = cEFTemp.begin(); 
    
#pragma omp parallel default(shared)
    {
#pragma omp for
      for (E_Int i = 0; i < sizeFN2; i++) cnp[i] = ptro1[i];
#pragma omp for
      for (E_Int i = 0; i < sizeEF2; i++) cne[i] = ptro2[i];
    }
    //champs aux centres des elts
    FldArrayF* fc2 = new FldArrayF(n, nfldc);
    for (E_Int v = 1; v <= nfldc; v++) 
    {
      E_Float* fcpl = fc->begin(v);
      E_Float* fc2pl = fc2->begin(v);
      for (E_Int i = 0; i < n; i++)
      {
        elt = eltListp[i];
        fc2pl[i] = fcpl[elt];
      }
    }
    FldArrayI* cncout = new FldArrayI();  *cncout = *c2n;
    PyObject* tpl1 = K_ARRAY::buildArray(*f2, varString, *c2n, -1, eltType, true); delete f2;
    PyObject* tpl2 = K_ARRAY::buildArray(*fc2, varStringc, *cncout, -1, eltType, true); delete fc2;
    PyList_Append(l, tpl1); Py_DECREF(tpl1); 
    PyList_Append(l, tpl2); Py_DECREF(tpl2); 
    delete cncout; delete c2n; 
  }
  else
  {
    // Pour les elts simples, on prend seulement la connectivite des elements en entree
    FldArrayI* c2n = new FldArrayI(n, nt);
    E_Int et, indv;
    FldArrayI indirNodes(npts); indirNodes.setAllValuesAt(-1); E_Int* indirNodesp = indirNodes.begin();
    E_Int nbnodes = 0;
    vector<E_Int> listOfNodes;
    
    // Selectionne les elements subzones
    for (E_Int i = 0; i < n; i++)
    {
      et = eltListp[i];
      for (E_Int v = 1; v <= nt; v++)
      {
        indv = (*cn)(et,v)-1;
        if ( indirNodesp[indv] == -1) 
        { 
          listOfNodes.push_back(indv);
          nbnodes++;
          indirNodesp[indv] = nbnodes;
          (*c2n)(i,v) = nbnodes;
        }
        else  
          (*c2n)(i,v) = indirNodesp[indv];
      }
    }
  
    indirNodes.malloc(0);
    FldArrayF* fout = new FldArrayF(nbnodes,nfld);
#pragma omp parallel default(shared)
    {
      for (E_Int eq = 1; eq <= nfld; eq++)
      {
        E_Float* fp = f->begin(eq);
        E_Float* fp2 = fout->begin(eq);
#pragma omp for
        for (E_Int i = 0; i < nbnodes; i++)
          fp2[i] = fp[listOfNodes[i]];    
      }
    }
    listOfNodes.clear();
    
    // Build array of nodes
    PyObject* tpl1 = K_ARRAY::buildArray(*fout, varString, *c2n, -1, eltType);
    delete fout; 
    PyList_Append(l, tpl1); Py_DECREF(tpl1);

    //champs aux centres des elts
    PyObject* tpl2 = K_ARRAY::buildArray(nfldc, varStringc, n, c2n->getSize(), -1, eltType, true); 
    E_Float* fcp = K_ARRAY::getFieldPtr(tpl2);
    FldArrayF fc2(n, nfldc, fcp, true);
    E_Int* cnpc = K_ARRAY::getConnectPtr(tpl2);
    K_KCORE::memcpy__(cnpc, c2n->begin(), c2n->getSize()*c2n->getNfld());
    for (E_Int v = 1; v <= nfldc; v++) 
    {
      E_Float* fcpl = fc->begin(v);
      E_Float* fc2pl = fc2.begin(v);
      for (E_Int i = 0; i < n; i++)
      {
        et = eltListp[i];
        fc2pl[i] = fcpl[et];
      }
    }
    PyList_Append(l, tpl2); Py_DECREF(tpl2);
    delete c2n;
  }
  RELEASESHAREDU(arrayNodes, f, cn);
  return l;
}
// ============================================================================
/* Subzone a unstructured mesh by face indices */
// ============================================================================
PyObject* K_TRANSFORM::subzoneFaces(PyObject* self, PyObject* args)
{
  PyObject* array;
  PyObject* listOfFaces;
  if (!PyArg_ParseTuple(args, "OO",
                        &array, &listOfFaces))
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

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray(array, varString, f, im, jm, km, cn, eltType, true); 
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
    type = 1; nf = 5; strcpy(eltTypeFaces, "NGON");
  }
  else if (K_STRING::cmp(eltType, "PYRA") == 0)
  {
    type = 1; nf = 5; strcpy(eltTypeFaces, "NGON");
  }
  else if (K_STRING::cmp(eltType, "NGON") == 0)
  {
    type = 2; strcpy(eltTypeFaces, "NGON");
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "subzone: element type not taken into account.");
    RELEASESHAREDU(array,f, cn); return NULL;
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
    vector< vector<E_Int> > indicesFaces(npts);// faces creees associees aux noeuds 

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

  RELEASESHAREDU(array,f, cn);
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
  if (!PyArg_ParseTuple(args, "OOO",
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
    type = 1; nf = 5; strcpy(eltTypeFaces, "NGON");
  }
  else if (K_STRING::cmp(eltType, "PYRA") == 0)
  {
    type = 1; nf = 5; strcpy(eltTypeFaces, "NGON");
  }
  else if (K_STRING::cmp(eltType, "NGON") == 0)
  {
    type = 2; strcpy(eltTypeFaces, "NGON");
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "subzone: element type not taken into account.");
    RELEASESHAREDU(arrayNodes,f, cn); RELEASESHAREDU(arrayCenters, fc, cnc);  return NULL;
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
    PyObject* tpl1 = K_ARRAY::buildArray(*f2, varString, 
                                         *c2n, -1, eltTypeFaces);
    PyList_Append(l, tpl1); Py_DECREF(tpl1); delete f2;
    RELEASESHAREDU(arrayNodes,f, cn);RELEASESHAREDU(arrayCenters, fc, cnc); 
    FldArrayI* cncout = new FldArrayI(); *cncout = *c2n;
    PyObject* tpl2 = K_ARRAY::buildArray(*fc2, varStringc, 
                                         *cncout, -1, eltTypeFaces, true);
    PyList_Append(l, tpl2); Py_DECREF(tpl2); delete fc2;  delete cncout; delete c2n;
    return l;
  }
  else if (type == 2) // IN: NGON, OUT: NGON
  {
    E_Int* ptr = cn->begin();
    E_Int nfacesTot = ptr[0];
    E_Int sizeFN = ptr[1]; ptr += 2; 
    E_Int ne = ptr[sizeFN]; // nombre d'elements pour la connectivite cn
    E_Int next1=0; // nbre de faces pour la nouvelle connectivite 
    E_Int next2=0; // nbre d'elements pour la nouvelle connectivite
    E_Int sizeco1=0; // taille de la nouvelle connectivite de faces 
    E_Int sizeco2=0; // taille de la nouvelle connectivite d'elements 

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
    PyObject* tpl1 = K_ARRAY::buildArray(*f2, varString, 
                                         *cnout, -1, eltTypeFaces);
    PyList_Append(l, tpl1); Py_DECREF(tpl1); delete f2;
    RELEASESHAREDU(arrayNodes,f, cn); RELEASESHAREDU(arrayCenters, fc, cnc); 
    FldArrayI* cncout = new FldArrayI(); *cncout = *cnout;
    PyObject* tpl2 = K_ARRAY::buildArray(*fc2, varStringc, 
                                         *cncout, -1, eltTypeFaces, true);
    PyList_Append(l, tpl2); Py_DECREF(tpl2); delete fc2; 
    delete cnout; delete cncout;
    return l;
  }
  else if (type == 1)
    PyErr_SetString(PyExc_TypeError,
                    "subzone: not implemented.\n");

  RELEASESHAREDU(arrayNodes,f, cn); RELEASESHAREDU(arrayCenters, fc, cnc); 
  return NULL;
}

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
# include "connector.h"

//=============================================================================
/* Met le cellN a 2 pour les points situes sur une frontiere overlap sur depth
   rangees. Peut etre localise en centres ou en noeuds. array doit deja
   etre en centres ou en noeuds
*/
//=============================================================================
PyObject* K_CONNECTOR::applyBCOverlapStruct(PyObject* self, PyObject* args)
{
  PyObject *array;
  E_Int imin, imax, jmin, jmax, kmin, kmax;
  E_Int depth; E_Int loc; E_Int cellNInterpValue;
  char* cellNName; 
  if (!PYPARSETUPLE_(args, O_ TIII_ TIII_ III_ S_,
                     &array, &imin, &jmin, &kmin, &imax, &jmax, &kmax, 
                     &depth, &loc, &cellNInterpValue, &cellNName))
    return NULL;
  
  E_Int shift = 0;
  if (loc == 0) shift = 1; // loc='nodes'

  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, im, jm, km, 
                                     cn, eltType);
  if (res != 1) 
  {    
    PyErr_SetString(PyExc_TypeError, 
                    "applyBCOverlaps: 1st argument must be structured.");
    RELEASESHAREDB(res, array, f, cn); return NULL;
  }

  E_Float interpolatedValue = cellNInterpValue;

  // verif cellN 
  E_Int posc = K_ARRAY::isNamePresent(cellNName,varString);
  if (posc == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "applyBCOverlaps: 1st arg must contain cellN variable.");
    RELEASESHAREDS(array, f); return NULL;
  }
  posc++;
  
  if (loc == 0) // nodes
  {
    if (imin < 1 || imax > im || jmin < 1 || jmax > jm || kmin < 1 || kmax > km)
    {
      PyErr_SetString(PyExc_TypeError,
                      "applyBCOverlaps: indices of structured window are not valid.");
      RELEASESHAREDS(array, f); return NULL;
    }
  }
  else // centers
  {
    if (imin < 1 || imax > im+1 || jmin < 1 || jmax > jm+1 || kmin < 1 || kmax > km+1)
    {
      PyErr_SetString(PyExc_TypeError,
                      "applyBCOverlaps: indices of structured window are not valid.");
      RELEASESHAREDS(array, f); return NULL;
    }
  }
  
  if (imin == imax)
  {
    if (imin == 1) {imin = 1; imax = K_FUNC::E_min(depth,im);}
    else { imin =  K_FUNC::E_max(imin-depth+shift,1); imax = imax-1+shift;}
  }
  else 
  {imax = imax-1+shift;}
  
  if (jmin == jmax )
  {
    if (jmin == 1) {jmin = 1; jmax = K_FUNC::E_min(depth, jm);}
    else {jmin = K_FUNC::E_max(jmin-depth+shift,1); jmax = jmax-1+shift;}
  }
  else 
  {jmax = jmax-1+shift;} 
    
  if (kmin == kmax)
  {
    if (kmin == 1) { kmin = 1 ; kmax = K_FUNC::E_min(depth, km);}
    else {kmin = K_FUNC::E_max(kmin-depth+shift,1); kmax = kmax-1+shift;}
  }
  else {kmax = kmax-1+shift;}

  E_Int imjm = im*jm;
  E_Float* cellNt = f->begin(posc);
#pragma omp parallel default(shared)
  {
    E_Int ind;
# pragma omp for
    for (E_Int k = kmin; k <= kmax; k++)
      for (E_Int j = jmin; j <= jmax; j++)
        for (E_Int i = imin; i <= imax; i++)
        {
          ind = (i-1) + (j-1)* im + (k-1)*imjm;
          if (cellNt[ind] != 0.) cellNt[ind] = interpolatedValue;
        }
  }
  RELEASESHAREDS(array, f);
  Py_INCREF(Py_None);
  return Py_None;
}

//=============================================================================
/* Met le cellN a 2 pour les centres des elements dans un voisinage de depth 
   elements d'une face definissant une frontiere overlap 
*/
//=============================================================================
PyObject* K_CONNECTOR::applyBCOverlapsNG(PyObject* self, PyObject* args)
{
  PyObject *array, *faceList;
  E_Int depth; E_Int loc; E_Int cellNInterpValue;
  char* cellNName;
  if (!PYPARSETUPLE_(args, OO_ III_ S_,
                    &array, &faceList, &depth, &loc, 
                    &cellNInterpValue, &cellNName))
  {
    return NULL;
  }

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, im, jm, km, 
                                     cn, eltType); 
  if (res != 2) 
  {    
    PyErr_SetString(PyExc_TypeError, 
                    "applyBCOverlaps: 1st argument not valid.");
    RELEASESHAREDB(res, array, f, cn); return NULL;
  }
  
  // verif cellN 
  E_Int posc = K_ARRAY::isNamePresent(cellNName, varString);
  if (posc == -1) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "applyBCOverlaps: 1st arg must contain cellN variable.");
    RELEASESHAREDU(array, f, cn); return NULL;
  }
  posc++;
  
  // Get list of face indices defined as a numpy array
  FldArrayI* indicesF;
  E_Int ret = K_NUMPY::getFromPointList(faceList, indicesF);
  if (ret != 1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "applyBCOverlaps: 2nd arg must be a numpy array of ints.");
    RELEASESHAREDU(array, f, cn); return NULL;
  }

  E_Float interpolatedValue = cellNInterpValue;
  
  PyObject* tpl;
  E_Int npts, nelts, sizeFN, nov1;
  E_Int nvoisins, nvertexV, neltsV, etv;
  E_Int* ptr;
  FldArrayI cFE;
  E_Int nfacesBC = indicesF->getSize()*indicesF->getNfld();
  std::vector<E_Int> voisins; std::vector<E_Int> voisinsL;  
 
  if (strcmp(eltType, "NGON") != 0 && strcmp(eltType, "NGON*") != 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "applyBCOverlaps: unstructured array must be a NGON or NGON*.");
    RELEASESHAREDN(faceList, indicesF);
    RELEASESHAREDU(array, f, cn); return NULL;
  }

  E_Int* indicesFp = indicesF->begin();
  if (loc == 0) // loc='nodes'
  {
    npts = f->getSize();
      
    ptr = cn->begin();
    sizeFN = ptr[1]; ptr += sizeFN+2;
    nelts = ptr[0];
    ptr += 2;

    FldArrayI tag(nelts); tag.setAllValuesAtNull();
    E_Int api = f->getApi();
    tpl = K_ARRAY::buildArray3(*f, varString,
                               *cn, "NGON", api);
    FldArrayF* out;
    K_ARRAY::getFromArray3(tpl, out);
    E_Float* cellNp = out->begin(posc);

    FldArrayI posFaces;
    K_CONNECT::getPosFaces(*cn, posFaces);

    ptr = cn->begin();
    // identify vertices in tagged face 
    for (E_Int nof = 0; nof < nfacesBC; nof++)
    {
      E_Int numFace = indicesFp[nof]-1;
      E_Int* ptrfn = ptr+posFaces[numFace];
      for (E_Int nov = 0; nov < ptrfn[0]; nov++)
      {
        nov1 = ptrfn[nov+1];
        cellNp[nov1-1] = interpolatedValue;
        voisins.push_back(nov1-1);
      }
    }
    //std::unique(voisins.begin(), voisins.end()); 

    // depth voisinage
    if (depth > 1)
    {        
      std::vector< std::vector<E_Int> > cVN(npts);
      K_CONNECT::connectNG2VNbrs(*cn, cVN);
      for (E_Int d = 1; d < depth; d++)
      {
        nvoisins = voisins.size();
        for (E_Int no = 0; no < nvoisins; no++)
        {
          E_Int vertex = voisins[no];

          //parcours de ses voisins
          std::vector<E_Int> vertexN = cVN[vertex];
          nvertexV = vertexN.size();
          for (E_Int novv = 0; novv < nvertexV; novv++)
          {
            nov1 = vertexN[novv]-1;
            cellNp[nov1] = interpolatedValue;
            voisinsL.push_back(nov1);
          }
        }
        //std::unique(voisinsL.begin(), voisinsL.end());
        voisins.clear(); voisins = voisinsL; voisinsL.clear();
      }
    }
    RELEASESHAREDS(tpl, out);
  } // cellN at nodes for NGON 
  else 
  {
    K_CONNECT::connectNG2FE(*cn, cFE);
    nelts = f->getSize(); // nombre total d'elements
        
    E_Int api = f->getApi();
    tpl = K_ARRAY::buildArray3(*f, varString,
                               *cn, "NGON", api);
    
    FldArrayF* out;
    K_ARRAY::getFromArray3(tpl, out);
    E_Float* cellNp = out->begin(posc);

    E_Int* cFE1 = cFE.begin(1);
    E_Int* cFE2 = cFE.begin(2);
      
    // identify elements with tagged face 
    for (E_Int nof = 0; nof < nfacesBC; nof++)
    {
      E_Int numFace = indicesFp[nof]-1;
      E_Int e1 = cFE1[numFace];
      E_Int e2 = cFE2[numFace];
      if (e1 > 0) 
      {
        cellNp[e1-1] = 2.;
        voisins.push_back(e1-1);
      }
      if (e2 > 0)
      {
        cellNp[e2-1] = 2.;
        voisins.push_back(e2-1);
      }
    }
    //std::unique(voisins.begin(), voisins.end());
    
    // depth elements
    if (depth > 1)
    {
      std::vector< std::vector<E_Int> > cEEN(nelts);
      K_CONNECT::connectFE2EENbrs(cFE,cEEN);
      for (E_Int d = 1; d < depth; d++)
      {
        nvoisins = voisins.size();
        for (E_Int noe = 0; noe < nvoisins; noe++)
        {
          E_Int et = voisins[noe];
          std::vector<E_Int>& eltsV = cEEN[et]; // demarrent a 0
          neltsV = eltsV.size();
          for (E_Int vv = 0; vv < neltsV; vv++)
          {
            etv = eltsV[vv];
            cellNp[etv] = interpolatedValue;
            voisinsL.push_back(etv);
          }
        }
        //std::unique(voisinsL.begin(), voisinsL.end());
        voisins.clear(); voisins = voisinsL; voisinsL.clear();
      }
    }
    RELEASESHAREDS(tpl, out);  
  }// End cellN at elements
  
  RELEASESHAREDN(faceList, indicesF);
  RELEASESHAREDU(array, f, cn);
  return tpl;
}

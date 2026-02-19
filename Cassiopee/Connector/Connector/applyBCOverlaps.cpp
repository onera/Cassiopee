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
# pragma omp for collapse(3)
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
  E_Int npts, nelts, nov1;
  E_Int nvoisins, nvertexV, neltsV, etv, nv;
  FldArrayI cFE;
  E_Int nfacesBC = indicesF->getSize()*indicesF->getNfld();
  std::vector<E_Int> voisins; std::vector<E_Int> voisinsL;
  voisins.reserve(nfacesBC*4);  // ballpark
  voisinsL.reserve(nfacesBC*4);
 
  if (K_STRING::cmp(eltType, 4, "NGON") != 0)
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
    E_Int api = f->getApi();
    nelts = cn->getNElts();
    E_Int* ngon = cn->getNGon(); E_Int* indPG = cn->getIndPG();

    FldArrayI tag(nelts); tag.setAllValuesAtNull();
    
    tpl = K_ARRAY::buildArray3(*f, varString, *cn, "NGON", api);
    FldArrayF* f2;
    K_ARRAY::getFromArray3(tpl, f2);
    E_Float* f2p = f2->begin(posc);

    // identify vertices in tagged face 
    for (E_Int i = 0; i < nfacesBC; i++)
    {
      E_Int numFace = indicesFp[i]-1;
      E_Int* face = cn->getFace(numFace, nv, ngon, indPG);
      for (E_Int j = 0; j < nv; j++)
      {
        nov1 = face[j];
        f2p[nov1-1] = interpolatedValue;
        voisins.push_back(nov1-1);
      }
    }
    std::sort(voisins.begin(), voisins.end());
    voisins.erase(
        std::unique(voisins.begin(), voisins.end()),
        voisins.end()
    );

    // depth voisinage
    if (depth > 1)
    {
      std::vector<std::vector<E_Int> > cVN(npts);
      K_CONNECT::connectNG2VNbrs(*cn, cVN);
      for (E_Int d = 1; d < depth; d++)
      {
        nvoisins = voisins.size();
        for (E_Int i = 0; i < nvoisins; i++)
        {
          E_Int vertex = voisins[i];

          //parcours de ses voisins
          const std::vector<E_Int>& vertexN = cVN[vertex];
          nvertexV = vertexN.size();
          for (E_Int j = 0; j < nvertexV; j++)
          {
            nov1 = vertexN[j]-1;
            f2p[nov1] = interpolatedValue;
            voisinsL.push_back(nov1);
          }
        }

        std::sort(voisinsL.begin(), voisinsL.end());
        voisinsL.erase(
            std::unique(voisinsL.begin(), voisinsL.end()),
            voisinsL.end()
        );
        voisins = std::move(voisinsL);
      }
    }
    RELEASESHAREDS(tpl, f2);
  }
  else 
  {
    K_CONNECT::connectNG2FE(*cn, cFE);
    nelts = f->getSize(); // nombre total d'elements

    E_Int api = f->getApi();
    tpl = K_ARRAY::buildArray3(*f, varString, *cn, "NGON", api);
    FldArrayF* f2;
    K_ARRAY::getFromArray3(tpl, f2);
    E_Float* f2p = f2->begin(posc);

    E_Int* cFE1 = cFE.begin(1);
    E_Int* cFE2 = cFE.begin(2);

    // identify elements with tagged face 
    for (E_Int i = 0; i < nfacesBC; i++)
    {
      E_Int numFace = indicesFp[i]-1;
      E_Int e1 = cFE1[numFace];
      E_Int e2 = cFE2[numFace];
      if (e1 > 0) 
      {
        f2p[e1-1] = 2.;
        voisins.push_back(e1-1);
      }
      if (e2 > 0)
      {
        f2p[e2-1] = 2.;
        voisins.push_back(e2-1);
      }
    }

    std::sort(voisins.begin(), voisins.end());
    voisins.erase(
        std::unique(voisins.begin(), voisins.end()),
        voisins.end()
    );
    
    // depth elements
    if (depth > 1)
    {
      std::vector<std::vector<E_Int> > cEEN(nelts);
      K_CONNECT::connectFE2EENbrs(cFE, cEEN);
      for (E_Int d = 1; d < depth; d++)
      {
        nvoisins = voisins.size();
        for (E_Int i = 0; i < nvoisins; i++)
        {
          E_Int et = voisins[i];
          std::vector<E_Int>& eltsV = cEEN[et]; // demarrent a 0
          neltsV = eltsV.size();
          for (E_Int j = 0; j < neltsV; j++)
          {
            etv = eltsV[j];
            f2p[etv] = interpolatedValue;
            voisinsL.push_back(etv);
          }
        }

        std::sort(voisinsL.begin(), voisinsL.end());
        voisinsL.erase(
            std::unique(voisinsL.begin(), voisinsL.end()),
            voisinsL.end()
        );
        voisins = std::move(voisinsL);
      }
    }
    RELEASESHAREDS(tpl, f2);  
  }
  
  RELEASESHAREDN(faceList, indicesF);
  RELEASESHAREDU(array, f, cn);
  return tpl;
}

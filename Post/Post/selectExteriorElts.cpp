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

// selectExteriorElts

# include <stdio.h>
# include <string.h>
# include "post.h"
# include "Connect/merge.h"

using namespace std;
using namespace K_FLD;

//=============================================================================
/* Selectionne les elts exterieurs (frontieres) d'un array */
// ============================================================================
PyObject* K_POST::selectExteriorElts(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PyArg_ParseTuple(args, "O", &array)) return NULL;
  
  char* varString; char* eltType;
  FldArrayF* f; FldArrayI* cn;
  E_Int ni, nj, nk; 
  E_Int res = K_ARRAY::getFromArray(array, varString, f, ni, nj, nk, 
                                    cn, eltType, true);
  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "selectExteriorElts: array is invalid.");
    return NULL;
  }
  if (res == 1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "exteriorElts: only for unstructured arrays.");
    RELEASESHAREDB(res, array, f, cn); return NULL;
  }
  //E_Int posx = K_ARRAY::isCoordinateXPresent(varString)+1;
  //E_Int posy = K_ARRAY::isCoordinateYPresent(varString)+1;
  //E_Int posz = K_ARRAY::isCoordinateZPresent(varString)+1;

  PyObject* t;
  if (strcmp(eltType, "NGON") == 0)
  {
    t = selectExteriorEltsNGon(*f, *cn, varString);
  }
  else 
  {
    /*
    if (posx == 0 || posy == 0 || posz == 0 || 
        strcmp(eltType, "BAR") == 0 || strcmp(eltType, "PYRA") == 0 ||
        strcmp(eltType, "PENTA") == 0)
      t = selectExteriorEltsBasic(*f, *cn, eltType, varString);
    else
      t = selectExteriorEltsBasic2(*f, *cn, eltType, varString, 
                                   posx, posy, posz);
    */
    t = selectExteriorEltsBasic(*f, *cn, eltType, varString);
  }
  RELEASESHAREDU(array, f, cn); 
  return t;
}

//=============================================================================
// Recherche topologique des elements exterieurs
// La partie tres chere est la construction EV2EENbrs
//==============================================================================
PyObject* K_POST::selectExteriorEltsBasic(FldArrayF& f, FldArrayI& cn, 
                                          char* eltType, char* varString)
{
  E_Int nelts = cn.getSize();
  vector< vector<E_Int> > cEEN(nelts);
  K_CONNECT::connectEV2EENbrs(eltType, f.getSize(), cn, cEEN);

  // Nombre de voisins pour un elt interieur
  unsigned int nvoisins = 0;
  if (strcmp(eltType, "BAR") == 0) nvoisins = 2;
  else if (strcmp(eltType, "QUAD") == 0) nvoisins = 4;
  else if (strcmp(eltType, "TRI") == 0) nvoisins = 3;
  else if (strcmp(eltType, "HEXA") == 0) nvoisins = 6;
  else if (strcmp(eltType, "TETRA") == 0) nvoisins = 4;
  else if (strcmp(eltType, "PYRA") == 0) nvoisins = 5;
  else if (strcmp(eltType, "PENTA") == 0) nvoisins = 5;
  
  E_Int nv = cn.getNfld();
  E_Int nthreads = __NUMTHREADS__;
  E_Int net = nelts/nthreads+1;
  E_Int** indirs = new E_Int* [nthreads];
  E_Int* nes = new E_Int [nthreads];
  for (E_Int i = 0; i < nthreads; i++) indirs[i] = new E_Int [net];
  E_Int* prev = new E_Int [nthreads];

#pragma omp parallel default(shared)
  {
    E_Int  ithread = __CURRENT_THREAD__;
    E_Int* indir = indirs[ithread];
    E_Int ne = 0;

#pragma omp for
    for (E_Int e = 0; e < nelts; e++)
    {
      if (cEEN[e].size() != nvoisins)
      { indir[ne] = e; ne++; }
    }
    nes[ithread] = ne;
  }
 
  // Nombre d'elements totaux
  E_Int netot = 0;
  for (E_Int i = 0; i < nthreads; i++) { prev[i] = netot; netot += nes[i]; }

  // Connectivite
  FldArrayI cnn(netot, nv);

#pragma omp parallel default(shared)
  {
    E_Int  ithread = __CURRENT_THREAD__;
    E_Int* indir = indirs[ithread];
    E_Int ne = nes[ithread];
    E_Int p = prev[ithread];
    for (E_Int nfld = 1; nfld <= nv; nfld++)
    {
      E_Int* cnnp = cnn.begin(nfld);
      E_Int* cnp = cn.begin(nfld);
      for (E_Int e = 0; e < ne; e++) cnnp[e+p] = cnp[indir[e]];
    }
  }

  // delete
  for (E_Int i = 0; i < nthreads; i++) delete [] indirs[i];
  delete [] indirs; delete [] nes; delete [] prev;

  PyObject* t = K_ARRAY::buildArray(f, varString, cnn, -1, eltType);
  return t;
}

//=============================================================================
// Recherche geometrique des elements exterieurs
// Moins rapide que la recherche topologique
//==============================================================================
PyObject* K_POST::selectExteriorEltsBasic2(FldArrayF& f, FldArrayI& cn, 
                                           char* eltType, char* varString,
                                           E_Int posx, E_Int posy, E_Int posz)
{
  E_Int nfaces = 0; E_Int nvertex = 0;
  if (strcmp(eltType, "TRI") == 0)
  { nfaces = 3; nvertex = 2; }
  else if (strcmp(eltType, "QUAD") == 0) 
  { nfaces = 4; nvertex = 2; }
  else if (strcmp(eltType, "TETRA") == 0) 
  { nfaces = 4;  nvertex = 3; }
  else if (strcmp(eltType, "BAR") == 0) 
  { nfaces = 2; nvertex = 3; }
  else if (strcmp(eltType, "HEXA") == 0) 
  { nfaces = 6; nvertex = 4; }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "exteriorElts: not implemented for this type of elements.");
    return NULL;
  }
 
  E_Int ne = cn.getSize();
  E_Int nv = cn.getNfld();
  E_Float* fx = f.begin(posx);
  E_Float* fy = f.begin(posy);
  E_Float* fz = f.begin(posz);

  E_Int nfacesmax = ne*nfaces;
  FldArrayI face(nfaces, nvertex);
  K_POST::buildFaceInfo2(cn, face);

  FldArrayF interfaceCenters(nfacesmax, 3);
  E_Float* interfacex = interfaceCenters.begin(1);
  E_Float* interfacey = interfaceCenters.begin(2);
  E_Float* interfacez = interfaceCenters.begin(3);
  E_Float fvertex = 1./nvertex;

#pragma omp parallel default(shared)
  {
    E_Float xm, ym, zm;
    E_Int ind, ff;
#pragma omp for
    for (E_Int et = 0; et < ne; et++)
    {
      for (E_Int j = 0; j < nfaces; j++)
      {
        xm = 0.; ym = 0.; zm = 0.;
        for (E_Int k = 0; k < nvertex; k++) 
        {
          ind = cn(et, face[j+k*nfaces])-1;
          xm += fx[ind]; ym += fy[ind]; zm += fz[ind];
        }
        xm = xm * fvertex; ym = ym * fvertex; zm = zm * fvertex;
        ff = j+nfaces*et;
        interfacex[ff] = xm;
        interfacey[ff] = ym;
        interfacez[ff] = zm;
      }
    }
  }

  ArrayAccessor<FldArrayF> coordAcc(interfaceCenters, 1, 2, 3);
  K_SEARCH::KdTree<FldArrayF> globalKdt(coordAcc);
  
  short* tag = new short [nfacesmax];
  for (E_Int i = 0; i < nfacesmax; i++) tag[i] = 0;
  E_Int nfaceExt = 0;

#pragma omp parallel default(shared)
  {
    E_Float pt[3]; E_Int ff; E_Float dx, dy, dz;
    E_Int ind, indp;
    for (E_Int et = 0; et < ne; et++)
    {
      for (E_Int j = 0; j < nfaces; j++)
      {
        ff = et*nfaces+j;
        pt[0] = interfacex[ff];
        pt[1] = interfacey[ff];
        pt[2] = interfacez[ff]; 
        ind = globalKdt.getClosest(pt);
        indp = globalKdt.getClosest(ind);
        dx = interfacex[indp]-interfacex[ind];
        dy = interfacey[indp]-interfacey[ind];
        dz = interfacez[indp]-interfacez[ind];
        if (dx*dx + dy*dy + dz*dz > 1.e-24) // exterieur
        {
          tag[ff] = 1; //nfaceExt++;
        }
      }
    }
  }
  
  for (E_Int i = 0; i < nfacesmax; i++) nfaceExt += tag[i];

  interfaceCenters.malloc(0);
  FldArrayI cnn(nfaceExt, nv);
  E_Int ee = 0;
  for (E_Int et = 0; et < ne; et++)
  {
    for (E_Int j = 0; j < nfaces; j++)
    {
      if (tag[et*nfaces+j] == 1) // exterieur
      {
        for (E_Int k = 1; k <= nv; k++) cnn(ee, k) = cn(et, k);
        ee++;
      }
    }
  }
  PyObject* t = K_ARRAY::buildArray(f, varString, cnn, -1, eltType);
  return t;
}

//=============================================================================
// Recherche geometrique des elements exterieurs
// Utilisation de merge, moins rapide que la recherche topologique
//==============================================================================
PyObject* K_POST::selectExteriorEltsBasic3(FldArrayF& f, FldArrayI& cn, 
                                           char* eltType, char* varString,
                                           E_Int posx, E_Int posy, E_Int posz)
{
  E_Int nfaces = 0; E_Int nvertex = 0;
  if (strcmp(eltType, "TRI") == 0)
  { nfaces = 3; nvertex = 2; }
  else if (strcmp(eltType, "QUAD") == 0) 
  { nfaces = 4; nvertex = 2; }
  else if (strcmp(eltType, "TETRA") == 0) 
  { nfaces = 4;  nvertex = 3; }
  else if (strcmp(eltType, "BAR") == 0) 
  { nfaces = 2; nvertex = 3; }
  else if (strcmp(eltType, "HEXA") == 0) 
  { nfaces = 6; nvertex = 4; }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "exteriorElts: not implemented for this type of elements.");
    return NULL;
  }
 
  E_Int ne = cn.getSize();
  E_Int nv = cn.getNfld();
  E_Float* fx = f.begin(posx);
  E_Float* fy = f.begin(posy);
  E_Float* fz = f.begin(posz);

  E_Int nfacesmax = ne*nfaces;
  FldArrayI face(nfaces, nvertex);
  K_POST::buildFaceInfo2(cn, face);

  FldArrayF interfaceCenters(nfacesmax, 3);
  E_Float* interfacex = interfaceCenters.begin(1);
  E_Float* interfacey = interfaceCenters.begin(2);
  E_Float* interfacez = interfaceCenters.begin(3);
  E_Float fvertex = 1./nvertex;

#pragma omp parallel default(shared)
  {
    E_Float xm, ym, zm;
    E_Int ind, ff;
#pragma omp for
    for (E_Int et = 0; et < ne; et++)
    {
      for (E_Int j = 0; j < nfaces; j++)
      {
        xm = 0.; ym = 0.; zm = 0.;
        for (E_Int k = 0; k < nvertex; k++) 
        {
          ind = cn(et, face[j+k*nfaces])-1;
          xm += fx[ind]; ym += fy[ind]; zm += fz[ind];
        }
        xm = xm * fvertex; ym = ym * fvertex; zm = zm * fvertex;
        ff = j+nfaces*et;
        interfacex[ff] = xm;
        interfacey[ff] = ym;
        interfacez[ff] = zm;
      }
    }
  }

  ArrayAccessor<FldArrayF> coordAcc(interfaceCenters, 1, 2, 3);
  vector<E_Int> ids;
  ::merge(coordAcc, 1.e-12, ids);

  interfaceCenters.setAllValuesAtNull(); // store int/ext tag
  E_Int nfaceExt = 0;
  E_Int idsSize = ids.size();
  for (E_Int i = 0; i < idsSize; i++)
  {
    if (ids[i] != i) // merged=interior
    {
      interfacex[i] = 1.; nfaceExt++;
    }
  }

  FldArrayI cnn(nfaceExt, nv);
  E_Int ee = 0; E_Int ff = 0; 
  for (E_Int et = 0; et < ne; et++)
  {
    for (E_Int j = 0; j < nfaces; j++)
    {
      if (interfacex[ff] > 0.5) // exterieur
      {
        for (E_Int k = 1; k <= nv; k++) cnn(ee, k) = cn(et, k);
        ee++;
      }
      ff++;
    }
  }
  PyObject* t = K_ARRAY::buildArray(f, varString, cnn, -1, eltType);
  return t;
}

//=============================================================================
// Recherche topologique
//=============================================================================
PyObject* K_POST::selectExteriorEltsNGon(FldArrayF& f, FldArrayI& cn, 
                                         char* varString)
{
  FldArrayI cFE; K_CONNECT::connectNG2FE(cn, cFE);
  E_Int* cFE1 = cFE.begin(1);
  E_Int* cFE2 = cFE.begin(2);
  E_Int* ptr = cn.begin();
  E_Int nfacesTot = ptr[0];
  E_Int sizeFN = ptr[1]; ptr += sizeFN+2;
  E_Int neltsTot = ptr[0];
  E_Int sizeEF = ptr[1];
  
  FldArrayI posFace;
  K_CONNECT::getPosFaces(cn, posFace);
  E_Int* posFacep = posFace.begin();
  
  E_Int nbnodes, nbfaces, extElt, indface;
  E_Int posf, e1, e2;
  E_Int* ptrFN = cn.begin()+2;
  E_Int* ptrEF = cn.begin()+4+sizeFN;

  FldArrayI connectEFTemp(sizeEF); E_Int* cEFTemp = connectEFTemp.begin(); 
  FldArrayI indirFaces(nfacesTot); indirFaces.setAllValuesAt(-1); E_Int* indirFacesp = indirFaces.begin();
  E_Int indFaceExt;
  vector<E_Int> origIndicesOfExternalFaces;
  E_Int sizeFN2 = 0; E_Int sizeEF2 = 0;
  E_Int nbFacesExt = 0; E_Int nbEltsExt = 0;
  for (E_Int et = 0; et < neltsTot; et++)
  {
    nbfaces = ptrEF[0]; // nbre de faces de l'element
    extElt = -1;

    for (E_Int nf = 0; nf < nbfaces; nf++)
    {
      indface = ptrEF[nf+1]-1;
      e1 = cFE1[indface];  // element voisin 1
      e2 = cFE2[indface];  // element voisin 2
      if (e2 == 0 && e1 > 0) { extElt = e1-1; break; }
      else if (e1 == 0 && e2 > 0) { extElt = e2-1; break; }
    }
    if (nbfaces == 1) { extElt = cFE1[ptrEF[1]-1]-1; } // Ghost cells
    if (extElt != -1)
    {
      // construction de la connectivite Elts/Faces des elts externes
      cEFTemp[0] = nbfaces;
      for (E_Int nof = 0; nof < nbfaces; nof++)
      {
        indface = ptrEF[nof+1]-1;        
        if (indirFacesp[indface] == -1) 
        {
          posf = posFacep[indface];
          ptrFN = cn.begin()+posf;
          sizeFN2 += ptrFN[0]+1;

          indFaceExt = nbFacesExt; 
          indirFacesp[indface] = indFaceExt; 
          nbFacesExt++;
          origIndicesOfExternalFaces.push_back(indface);
        }
        else indFaceExt = indirFacesp[indface];

        cEFTemp[nof+1] = indFaceExt+1;
      }
      cEFTemp += nbfaces+1; sizeEF2 += nbfaces+1; nbEltsExt += 1;
    }
    ptrEF += nbfaces+1;
  }
  cFE.malloc(0); indirFaces.malloc(0);
  // on connait sizeEF2, nbEltsExt, nbFacesExt (sans doublons)
  connectEFTemp.resize(sizeEF2);
  // construction de la connectivite Faces/Noeuds
  FldArrayI connectFNTemp(sizeFN); E_Int* cFNTemp = connectFNTemp.begin(); 
  E_Int npts = f.getSize(); E_Int nfld = f.getNfld();

  FldArrayI indirNodes(npts); indirNodes.setAllValuesAt(-1); E_Int* indirNp = indirNodes.begin();
  E_Int indnode, indfaceo;
  E_Int numNode = 0;
  for (E_Int nfe = 0 ; nfe < nbFacesExt; nfe++)
  {
    indfaceo = origIndicesOfExternalFaces[nfe]; //demarre a 0
    posf = posFacep[indfaceo];
    ptrFN = cn.begin()+posf;
    nbnodes = ptrFN[0];
    cFNTemp[0] = nbnodes;
    for (E_Int p = 1; p <= nbnodes; p++)
    {
      indnode = ptrFN[p]-1;
      if (indirNp[indnode] == -1) //creation
      {
        indirNp[indnode] = numNode+1;
        cFNTemp[p] = numNode+1;
        numNode++;
      }
      else
      {
        cFNTemp[p] = indirNp[indnode];
      }
    }
    cFNTemp += nbnodes+1;
  }
  posFace.malloc(0); origIndicesOfExternalFaces.clear();

  // construit l'array de sortie
  E_Int csize = sizeFN2+sizeEF2+4;
  PyObject* tpl= K_ARRAY::buildArray(nfld, varString, numNode, csize, 8, 
                                     "NGON", false, csize); 

  E_Float* fnp = K_ARRAY::getFieldPtr(tpl);
  FldArrayF f2(numNode, nfld, fnp, true);
  
#pragma omp parallel default(shared)
  {
    E_Int indf;
    for (E_Int eq = 1; eq <= nfld; eq++)
    {
      E_Float* fp = f.begin(eq);
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
  cnp[0] = nbFacesExt; cnp[1] = sizeFN2;
  cnp += 2;
  E_Int* cne = cnp + sizeFN2;
  cne[0] = nbEltsExt; 
  cne[1] = sizeEF2;
  cne += 2;
  E_Int* ptro1 = connectFNTemp.begin();  
  E_Int* ptro2 = connectEFTemp.begin(); 

#pragma omp parallel default(shared)
  {
#pragma omp for
    for (E_Int i = 0; i < sizeFN2; i++) cnp[i] = ptro1[i];
#pragma omp for
    for (E_Int i = 0; i < sizeEF2; i++) cne[i] = ptro2[i];
  }
  return tpl;
}

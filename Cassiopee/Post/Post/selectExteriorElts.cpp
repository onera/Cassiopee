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
# include "post.h"
# include "Nuga/include/merge.h"

using namespace K_FLD;

//=============================================================================
/* Selectionne les elts exterieurs (frontieres) d'un array */
// ============================================================================
PyObject* K_POST::selectExteriorElts(PyObject* self, PyObject* args)
{
  PyObject* array;
  if (!PYPARSETUPLE_(args, O_, &array)) return NULL;
  
  char* varString; char* eltType;
  FldArrayF* f; FldArrayI* cn;
  E_Int ni, nj, nk; 
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, ni, nj, nk, 
                                     cn, eltType);
  if (res == 1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "exteriorElts: only for unstructured arrays.");
    RELEASESHAREDB(res, array, f, cn); return NULL;
  }
  else if (res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "selectExteriorElts: array is invalid.");
    return NULL;
  }

  PyObject* tpl;
  if (strcmp(eltType, "NGON") == 0)
  {
    tpl = selectExteriorEltsNGon(*f, *cn, varString);
  }
  else 
  {
    tpl = selectExteriorEltsBasic(*f, *cn, eltType, varString);
  }
  RELEASESHAREDU(array, f, cn); 
  return tpl;
}

//=============================================================================
// Recherche topologique des elements exterieurs
// La partie tres chere est la construction EV2EENbrs
//==============================================================================
PyObject* K_POST::selectExteriorEltsBasic(FldArrayF& f, FldArrayI& cn, 
                                          char* eltType, char* varString)
{
  E_Int nelts = cn.getSize();
  std::vector<std::vector<E_Int> > cEEN(nelts);
  K_CONNECT::connectEV2EENbrs(eltType, f.getSize(), cn, cEEN);

  // Nombre de voisins pour un elt interieur
  size_t nvoisins = 0;
  if (strcmp(eltType, "BAR") == 0) nvoisins = 2;
  else if (strcmp(eltType, "QUAD") == 0) nvoisins = 4;
  else if (strcmp(eltType, "TRI") == 0) nvoisins = 3;
  else if (strcmp(eltType, "HEXA") == 0) nvoisins = 6;
  else if (strcmp(eltType, "TETRA") == 0) nvoisins = 4;
  else if (strcmp(eltType, "PYRA") == 0) nvoisins = 5;
  else if (strcmp(eltType, "PENTA") == 0) nvoisins = 5;
  
  E_Int api = f.getApi();
  E_Int nv = cn.getNfld();
  E_Int nthreads = __NUMTHREADS__;
  E_Int net = nelts/nthreads+1;
  E_Int** indirs = new E_Int* [nthreads];
  E_Int* nes = new E_Int [nthreads];
  for (E_Int i = 0; i < nthreads; i++) indirs[i] = new E_Int [net];
  E_Int* prev = new E_Int [nthreads];

#pragma omp parallel
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

  #pragma omp parallel
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

  PyObject* tpl = K_ARRAY::buildArray3(f, varString, cnn, eltType, api);
  return tpl;
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
 
  E_Int api = f.getApi();
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
  PyObject* tpl = K_ARRAY::buildArray3(f, varString, cnn, eltType, api);
  return tpl;
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
 
  E_Int api = f.getApi();
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

  #pragma omp parallel
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
  std::vector<E_Int> ids;
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
  PyObject* tpl = K_ARRAY::buildArray3(f, varString, cnn, eltType, api);
  return tpl;
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

  E_Int npts = f.getSize(); E_Int nfld = f.getNfld(); E_Int api = f.getApi();
  E_Int* ngon = cn.getNGon(); E_Int* indPG = cn.getIndPG();
  E_Int* nface = cn.getNFace(); E_Int* indPH = cn.getIndPH();
  E_Int nfaces = cn.getNFaces(); E_Int nelts = cn.getNElts();
  E_Int sizeFN = cn.getSizeNGon(); E_Int sizeEF = cn.getSizeNFace();
  E_Int ngonType = cn.getNGonType();
  E_Int shift = 1; if (ngonType == 3) shift = 0;
  E_Bool hasCnOffsets = (ngonType == 2 || ngonType == 3);
  
  E_Int nf, nv, extElt, indface;
  E_Int pos, e1, e2;

  FldArrayI indirFaces(nfaces); indirFaces.setAllValuesAt(-1);
  FldArrayI nface2Temp(sizeEF), indPH2Temp(0);
  if (hasCnOffsets) { indPH2Temp.resize(nelts+1); indPH2Temp[0] = 0; }

  E_Int indFaceExt;
  std::vector<E_Int> origIndicesOfExtFaces;
  origIndicesOfExtFaces.reserve(E_Int(nfaces/10));  // random guess
  E_Int sizeFN2 = 0; E_Int sizeEF2 = 0;
  E_Int nExtFaces = 0; E_Int nExtElts = 0;

  for (E_Int i = 0; i < nelts; i++)
  {
    extElt = -1;
    E_Int* elt = cn.getElt(i, nf, nface, indPH);
    
    if (nf == 1) { extElt = cFE1[elt[0]-1] - 1; }  // Ghost cell
    else
    {
      for (E_Int f = 0; f < nf; f++)
      {
        indface = elt[f] - 1;
        e1 = cFE1[indface];  // element voisin 1
        e2 = cFE2[indface];  // element voisin 2
        if (e2 == 0 && e1 > 0) { extElt = e1 - 1; break; }
        else if (e1 == 0 && e2 > 0) { extElt = e2 - 1; break; }
      }
    }
    
    if (extElt != -1)  // exterior element found
    {
      // construction de la connectivite Elts/Faces des elts externes
      nface2Temp[sizeEF2] = nf;
      if (hasCnOffsets) indPH2Temp[nExtElts+1] = indPH2Temp[nExtElts] + nf;
      for (E_Int f = 0; f < nf; f++)
      {
        indface = elt[f] - 1;
        if (indirFaces[indface] == -1) 
        {
          cn.getFace(indface, nv, ngon, indPG);
          sizeFN2 += nv + shift;
          indFaceExt = nExtFaces; nExtFaces++;
          indirFaces[indface] = indFaceExt;
          origIndicesOfExtFaces.push_back(indface);
        }
        else indFaceExt = indirFaces[indface];

        nface2Temp[sizeEF2+f+shift] = indFaceExt + 1;
      }
      sizeEF2 += nf + shift; nExtElts++;
    }
  }
  cFE.malloc(0); indirFaces.malloc(0);

  // on connait sizeEF2, nExtElts, nExtFaces (sans doublons)
  nface2Temp.resize(sizeEF2);

  // construction de la connectivite Faces/Noeuds
  E_Int npts2 = 0;
  E_Int indnode, indface2;
  FldArrayI indirNodes(npts); indirNodes.setAllValuesAt(-1);
  FldArrayI ngon2Temp(sizeFN), indPG2Temp(0);
  if (hasCnOffsets) { indPG2Temp.resize(nExtFaces+1); indPG2Temp[0] = 0; }

  pos = 0;
  for (E_Int f = 0; f < nExtFaces; f++)
  {
    indface2 = origIndicesOfExtFaces[f];  // starts at 0
    E_Int* face = cn.getFace(indface2, nv, ngon, indPG);
    ngon2Temp[pos] = nv;
    if (hasCnOffsets) indPG2Temp[f+1] = indPG2Temp[f] + nv;
    for (E_Int p = 0; p < nv; p++)
    {
      indnode = face[p] - 1;
      if (indirNodes[indnode] == -1)  // creation
      {
        indirNodes[indnode] = npts2 + 1;
        ngon2Temp[pos+p+shift] = npts2 + 1;
        npts2++;
      }
      else ngon2Temp[pos+p+shift] = indirNodes[indnode];
    }
    pos += nv + shift;
  }
  origIndicesOfExtFaces.clear();

  // Reconstruction de la connectivite finale
  PyObject* tpl = K_ARRAY::buildArray3(
    nfld, varString, npts2, nExtElts, nExtFaces,
    "NGON", sizeFN2, sizeEF2, ngonType, false, api
  ); 
  FldArrayF* f2; FldArrayI* cn2;
  K_ARRAY::getFromArray3(tpl, f2, cn2);

  E_Int *ngon2 = cn2->getNGon(), *nface2 = cn2->getNFace();
  E_Int *indPG2 = NULL, *indPH2 = NULL;
  if (hasCnOffsets) { indPG2 = cn2->getIndPG(); indPH2 = cn2->getIndPH(); }

  #pragma omp parallel
  {
    E_Int indf;
    // Copy compressed fields to f2
    for (E_Int eq = 1; eq <= nfld; eq++)
    {
      E_Float* fp = f.begin(eq);
      E_Float* f2p = f2->begin(eq);
      #pragma omp for nowait
      for (E_Int ind = 0; ind < npts; ind++)
      {
        indf = indirNodes[ind] - 1;
        if (indf > -1) f2p[indf] = fp[ind];
      }
    }

    // Copy compressed connectivity to cn2
    #pragma omp for nowait
    for (E_Int i = 0; i < sizeFN2; i++) ngon2[i] = ngon2Temp[i];
    #pragma omp for nowait
    for (E_Int i = 0; i < sizeEF2; i++) nface2[i] = nface2Temp[i];

    if (hasCnOffsets)
    {
      #pragma omp for nowait
      for (E_Int i = 0; i < nExtFaces; i++) indPG2[i] = indPG2Temp[i];
      #pragma omp for
      for (E_Int i = 0; i < nExtElts; i++) indPH2[i] = indPH2Temp[i];
    }
  }

  RELEASESHAREDU(tpl, f2, cn2);
  return tpl;
}

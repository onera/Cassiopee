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
# include "post.h"
# include "String/kstring.h"
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
    tpl = selectExteriorEltsME(*f, *cn, eltType, varString);
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

//=============================================================================
// Recherche topologique des elements exterieurs utilisant la connectivite
// EV2NNbrs
//==============================================================================
PyObject* K_POST::selectExteriorEltsME(FldArrayF& f, FldArrayI& cn, 
                                       char* eltType, char* varString)
{
  E_Int nc = cn.getNConnect();
  E_Int nfld = f.getNfld();
  E_Int api = f.getApi();
  E_Int npts = f.getSize();
  std::vector<char*> eltTypes;
  K_ARRAY::extractVars(eltType, eltTypes);

  // Compute total number of elements across all connectivities, ntotElts
  std::vector<E_Int> nepc(nc);
  std::vector<E_Int> cumnepc(nc+1); cumnepc[0] = 0;  // cumulative number of elts per conn.
  for (E_Int ic = 0; ic < nc; ic++)
  {
    K_FLD::FldArrayI& cm = *(cn.getConnect(ic));
    E_Int nelts = cm.getSize();
    nepc[ic] = nelts;
    cumnepc[ic+1] = cumnepc[ic] + nelts;
  }
  E_Int ntotElts = cumnepc[nc];

  // Compute number of neighbour elements of internal elements, that is the
  // number of faces per element, nfpe
  std::vector<E_Int> nfpe;
  E_Int ierr = K_CONNECT::getNFPE(nfpe, eltType, true);
  if (ierr != 0) return NULL;

  // Build the element -> number of neighbour elements connectivity
  std::vector<E_Int> cENN(ntotElts);
  K_CONNECT::connectEV2NNbrs(eltType, npts, cn, cENN);

  // Uniform chunks (schedule: static) with at most 'net' elements per thread
  E_Int nthreads = __NUMTHREADS__;
  E_Int net = ntotElts/nthreads + nc;
  // Thread-related arrays are prefixed with 't'. For each thread:
  //  - tindir: maps element indices from new to old ME
  E_Int** tindir = new E_Int* [nthreads];
  //  - tnextepc: number of exterior elements found in each connectivity
  E_Int** tnextepc = new E_Int* [nthreads];
  //  - toffset: cumulative number of exterior elements found in each connectivity
  E_Int** toffset = new E_Int* [nthreads];
  for (E_Int i = 0; i < nthreads; i++)
  {
    tindir[i] = new E_Int [net];
    tnextepc[i] = new E_Int [nc];
    toffset[i] = new E_Int [nc];
  }

  // Number of elements per connectivity of the output ME
  // ('tmp_' is uncompressed: same number of connectivities as the input ME)
  std::vector<E_Int> tmp_nepc2(nc, 0);

  // In a first pass, tag vertex indices that belong to exterior elements
  std::vector<E_Int> vindir(npts, 0);

  #pragma omp parallel
  {
    E_Int indv;
    E_Int eidx;  // global element index
    E_Int nneis;  // number of neighbours of element eidx
    E_Int nextElts = 0;  // number of exterior elements found in all conn. of that thread
    E_Int nextEltsIc;  // number of exterior elements found in a given conn. of that thread
    // Local thread-related arrays are prefixed with 'loc_t'
    E_Int tid = __CURRENT_THREAD__;
    E_Int* loc_tindir = tindir[tid];
    E_Int* loc_tnextepc = tnextepc[tid];
    std::vector<E_Int> loc_ttmp_nepc2(nc, 0);

    for (E_Int ic = 0; ic < nc; ic++)
    {
      FldArrayI& cm = *(cn.getConnect(ic));
      E_Int nvpe = cm.getNfld();
      nextEltsIc = 0;

      #pragma omp for schedule(static)
      for (E_Int i = 0; i < nepc[ic]; i++)
      {
        eidx = cumnepc[ic] + i;
        nneis = cENN[eidx];
        if (nneis != nfpe[ic])  // exterior element found
        {
          loc_tindir[nextElts] = i; nextElts++; nextEltsIc++;
          loc_ttmp_nepc2[ic]++;

          // Tag vertices as exterior vertices
          for (E_Int j = 1; j <= nvpe; j++)
          {
            indv = cm(i, j) - 1;
            vindir[indv] = 1;
          }
        }
      }

      loc_tnextepc[ic] = nextEltsIc;
    }

    #pragma omp critical
    {
      for (E_Int ic = 0; ic < nc; ic++) tmp_nepc2[ic] += loc_ttmp_nepc2[ic];
    }
  }

  // Transform the exterior vertex mask of zeros and ones into a vertex map
  // from old to new connectivities, and get the number of unique exterior
  // vertices, npts2
  E_Int npts2 = K_CONNECT::prefixSum(vindir);

  // Compute thread element offsets in the output ME for each connectivity
  // toffset is a cumulative tnextepc over all conns
  // Used to build cm2 using multiple threads
  {
    E_Int* loc_toffset = toffset[0];
    for (E_Int ic = 0; ic < nc; ic++) loc_toffset[ic] = 0;
  }
  
  for (E_Int i = 1; i < nthreads; i++)
  {
    E_Int* tnextepcm1 = tnextepc[i-1];
    E_Int* loc_toffset = toffset[i];
    E_Int* toffsetm1 = toffset[i-1];
    for (E_Int ic = 0; ic < nc; ic++)
      loc_toffset[ic] = toffsetm1[ic] + tnextepcm1[ic];
  }

  // Free memory
  cENN.clear(); cENN.shrink_to_fit();

  // Build new eltType from connectivities that have at least one element
  E_Int nc2 = 0;
  char* eltType2 = new char[K_ARRAY::VARSTRINGLENGTH];
  eltType2[0] = '\0';
  for (E_Int ic = 0; ic < nc; ic++)
  {
    if (tmp_nepc2[ic] > 0)
    {
      nc2++;
      if (eltType2[0] == '\0') strcpy(eltType2, eltTypes[ic]);
      else
      {
        strcat(eltType2, ",");
        strcat(eltType2, eltTypes[ic]);
      }
    }
  }

  // Compress the number of elements per connectivity of the output ME, ie,
  // drop connectivities containing no exterior elements
  std::vector<E_Int> nepc2(nc2);
  nc2 = 0;
  for (E_Int ic = 0; ic < nc; ic++)
  {
    if (tmp_nepc2[ic] > 0) { nepc2[nc2] = tmp_nepc2[ic]; nc2++; }
  }
 
  // Build new connectivity
  PyObject* tpl = K_ARRAY::buildArray3(nfld, varString, npts2,
                                       nepc2, eltType2, false, api);
  FldArrayF* f2; FldArrayI* cn2;
  K_ARRAY::getFromArray3(tpl, f2, cn2);

  #pragma omp parallel
  {
    E_Int ic2, indv, inde, nelts, nvpe;
    E_Int offR;  // cumulative element offset of a given conn. to Read from cm
    E_Int offW;  // cumulative element offset of a given conn. to Write into cm2
    E_Int tid = __CURRENT_THREAD__;
    E_Int* loc_tindir = tindir[tid];
    E_Int* loc_tnextepc = tnextepc[tid];
    E_Int* loc_toffset = toffset[tid];

    // Fields
    for (E_Int n = 1; n <= nfld; n++)
    {
      E_Float* fp = f.begin(n);
      E_Float* f2p = f2->begin(n);
      #pragma omp for
      for (E_Int i = 0; i < npts; i++)
      {
        indv = vindir[i];
        if (indv > 0) f2p[indv-1] = fp[i];
      }
    }

    // Connectivity
    ic2 = 0;
    offR = 0;
    for (E_Int ic = 0; ic < nc; ic++)
    {
      if (tmp_nepc2[ic] == 0) continue;  // no exterior elements in this conn, skip
      FldArrayI& cm = *(cn.getConnect(ic));
      FldArrayI& cm2 = *(cn2->getConnect(ic2));
      nelts = loc_tnextepc[ic];
      offW = loc_toffset[ic];
      nvpe = cm.getNfld();
      for (E_Int i = 0; i < nelts; i++)
      {
        inde = loc_tindir[i+offR];
        for (E_Int j = 1; j <= nvpe; j++)
        {
          indv = cm(inde, j) - 1;
          cm2(offW+i, j) = vindir[indv];
        }
      }
      ic2++;
      offR += nelts;
    }
  }

  // Free memory
  for (E_Int i = 0; i < nthreads; i++)
  {
    delete [] tindir[i];
    delete [] tnextepc[i];
    delete [] toffset[i];
  }
  delete [] tindir; delete [] tnextepc; delete [] toffset;

  RELEASESHAREDU(tpl, f2, cn2);
  delete [] eltType2;
  for (size_t ic = 0; ic < eltTypes.size(); ic++) delete [] eltTypes[ic];
  return tpl;
}

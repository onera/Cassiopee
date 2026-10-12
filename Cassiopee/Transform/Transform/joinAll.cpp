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
# include "transform.h"
# include <stdio.h>

using namespace std;
using namespace K_FLD;
using namespace K_FUNC;
using namespace K_CONST;

// ============================================================================
/* Join all arrays and their arrays at centers */
// ============================================================================
PyObject* K_TRANSFORM::joinAll(PyObject* self, PyObject* args)
{
  PyObject *arrays, *arraysc=NULL; E_Float tol;
  // Check different signatures
  if (!PYPARSETUPLE_(args, O_ R_, &arrays, &tol))
  {
    PyErr_Clear();
    if (!PYPARSETUPLE_(args, OO_ R_, &arrays, &arraysc, &tol)) return NULL;
  }

  // Check arrays for fields located at nodes
  vector<E_Int> res;
  vector<char*> structVarString; vector<char*> unstructVarString;
  vector<FldArrayF*> structF; vector<FldArrayF*> unstructF;
  vector<E_Int> ni; vector<E_Int> nj; vector<E_Int> nk;
  vector<FldArrayI*> cn; vector<char*> eltType;
  vector<PyObject*> objs, obju;
  K_ARRAY::getFromArrays(arrays, res, structVarString,
                         unstructVarString, structF,
                         unstructF, ni, nj, nk,
                         cn, eltType, objs, obju,
                         true, true, true, false, true);
  E_Int nu = unstructF.size();

  // Check arrays for fields located at centers
  E_Int nuc = 0, nfldc = 0;
  vector<E_Int> resc;
  vector<char*> structVarStringc; vector<char*> unstructVarStringc;
  vector<FldArrayF*> structFc; vector<FldArrayF*> unstructFc;
  vector<E_Int> nic, njc, nkc;
  vector<FldArrayI*> cnc; vector<char*> eltTypec;
  vector<PyObject*> objsc, objuc;
  if (arraysc != NULL)
  {
    K_ARRAY::getFromArrays(arraysc, resc, structVarStringc,
			                     unstructVarStringc, structFc,
			                     unstructFc, nic, njc, nkc,
			                     cnc, eltTypec, objsc, objuc,
			                     false, false, false, false, true);
    nuc = unstructFc.size();
  }

  // Fusion des zones non-structures
  PyObject* tpln = NULL;
  if ((nu == 0 && nuc == 0) || (nuc > 0 && nu != nuc))
  {
    for (E_Int k = 0; k < nu; k++) RELEASESHAREDU(obju[k], unstructF[k], cn[k]);
    for (E_Int k = 0; k < nuc; k++) RELEASESHAREDU(objuc[k], unstructFc[k], cnc[k]);
    PyErr_SetString(PyExc_ValueError,
                    "joinAll: number of arrays at nodes and centers must "
                    "be equal.");
    return NULL;
  }

  const char* eltRef = NULL;
  E_Int missed = 0;
  E_Int nuref = 0;
  E_Int nc = 0, dimRef = -1, dim = -1;
  char newEltType[K_ARRAY::VARSTRINGLENGTH]; newEltType[0] = '\0';
  // Counters for all arrays
  E_Int npts = 0;
  E_Int nfaces = 0, neltsNGON = 0, sizeFN = 0, sizeEF = 0;
  E_Int shift = 1;
  // Table d'indirection des connectivites pour les ME
  // Masque pour le NGON
  vector<vector<E_Int> > indir(nu);

  for (E_Int k = 0; k < nu; k++)
  {
    E_Int nptsk = unstructF[k]->getSize();
    if (nptsk == 0) { indir[k].push_back(-2); missed++; continue; }  // skip empty NODE conn.
    if (eltRef == NULL) { nuref = k; eltRef = eltType[k]; }

    if (K_STRING::cmp(eltType[k], "NGON") == 0)
    {
      // La connectivite fusionee ne doit avoir que des NGONs
      if (K_STRING::cmp(eltRef, "NGON") != 0)
      {
        indir[k].push_back(-2); missed++; continue;
      }
      npts += nptsk;
      neltsNGON += cn[k]->getNElts();
      nfaces += cn[k]->getNFaces();
      sizeFN += cn[k]->getSizeNGon();
      sizeEF += cn[k]->getSizeNFace();
      indir[k].push_back(1);
    }
    else if (K_STRING::cmp(eltRef, "NGON") != 0)
    {
      // Calcul du nombre d'elt types dans la connectivite ME fusionee
      // et de leur identite
      vector<char*> eltTypesk;
      K_ARRAY::extractVars(eltType[k], eltTypesk);

      E_Int nck = cn[k]->getNConnect();
      for (E_Int ic = 0; ic < nck; ic++)
      {
        char* eltTypConn = eltTypesk[ic];
        dim = K_CONNECT::getDimME(eltTypConn);
        // Check dimensionality: merge if identical
        if (dimRef == -1) dimRef = dim;
        else if (dim != dimRef) { indir[k].push_back(-2); missed++; continue; }
        npts += nptsk;
        // Add default value in mapping table
        indir[k].push_back(-1);
        // Concatenate elttypes, discard duplicates
        if (strstr(newEltType, eltTypConn) == NULL)
        {
          strcat(newEltType, eltTypConn); strcat(newEltType, ",");
          nc += 1;
        }
      }

      for (size_t ic = 0; ic < eltTypesk.size(); ic++) delete [] eltTypesk[ic];
    }
    else { indir[k].push_back(-2); missed++; }
  }

  E_Int nfld = unstructF[nuref]->getNfld();
  E_Int api = (nc > 1) ? 3 : unstructF[nuref]->getApi();

  if (npts == 0)  // all input arrays are empty, create an empty conn.
  {
    printf("Warning: joinAll: all arrays are empty.\n");

    for (E_Int k = 0; k < nu; k++) RELEASESHAREDU(obju[k], unstructF[k], cn[k]);
    for (E_Int k = 0; k < nuc; k++) RELEASESHAREDU(objuc[k], unstructFc[k], cnc[k]);

    // The eltType of the output empty conn. is conserved if all input eltTypes
    // are the same, otherwise NODE is chosen.
    const char* eltType2 = NULL;
    eltRef = eltType[0];
    if (K_STRING::cmp(eltRef, "NGON") == 0) eltType2 = "NODE";
    else
    {
      for (E_Int k = 1; k < nu; k++)
      {
        if (K_STRING::cmp(eltType[k], eltRef) != 0) { eltType2 = "NODE"; break; }
      }
      if (eltType2 == NULL) eltType2 = eltRef;
    }

    tpln = K_ARRAY::buildArray3(nfld, unstructVarString[nuref],
                                0, 0, eltType2, false, api);
    if (arraysc == NULL) return tpln;
    else
    {
      PyObject* l = PyList_New(0);
      PyObject* tplc = K_ARRAY::buildArray3(nfldc, unstructVarStringc[nuref],
                                            0, 0, eltType2, false, api);
      PyList_Append(l, tpln); Py_DECREF(tpln);
      PyList_Append(l, tplc); Py_DECREF(tplc);
      return l;
    }
  }

  if (missed > 0)
    printf("Warning: joinAll: some arrays cannot be joined: different mesh "
           "types or empty array(s) found.\n");

  // Build unstructured connectivity
  E_Int ngonType = -1;
  E_Int ntotElts = 0;
  E_Int res2 = 2;
  FldArrayF* f2; FldArrayI* cn2;

  // Remplissage table d'indirection et nombre d'elements par eltType
  vector<E_Int> nepc(nc);                   // number of elements per output conn. (contiguous)
  vector<E_Int> cumnepc(nc, 0);             // cumulative number of elements per output conn. (contiguous)
  vector<vector<E_Int> > cumnepcOrig(nu);   // cumulative number of elements per input list of conn. (may not be contiguous)

  if (K_STRING::cmp(eltRef, "NGON") == 0)
  {
    strcpy(newEltType, "NGON");
    ntotElts = neltsNGON;
    ngonType = cn[nuref]->getNGonType();
    if (ngonType == 3) shift = 0;
    tpln = K_ARRAY::buildArray3(nfld, unstructVarString[nuref], npts, neltsNGON,
                                nfaces, newEltType, sizeFN, sizeEF,
                                ngonType, false, api);
    K_ARRAY::getFromArray3(tpln, f2, cn2);
  }
  else if (dimRef == 0)
  {
    strcpy(newEltType, "NODE");
    tpln = K_ARRAY::buildArray3(nfld, unstructVarString[nuref], npts, 0,
                                newEltType, false, api);
    cn2 = new FldArrayI();
    res2 = 1; K_ARRAY::getFromArray3(tpln, f2);
  }
  else
  {
    // Remove trailing comma in newEltType
    E_Int len = strlen(newEltType);
    newEltType[len-1] = '\0';

    if (nc == 2 && dimRef == 3) // TODO
    {
      // HEXA & TETRA cannot be joined in a conformal mesh, skipping the last
      // one of the two
      if (strstr(newEltType, "HEXA") != NULL && strstr(newEltType, "TETRA") != NULL)
      {
        len = strchr(newEltType, ',') - newEltType + 1;
        newEltType[len-1] = '\0';
        nc = 1;
        printf("Warning: joinAll: joining HEXA and TETRA would result in a "
               "non-conformal mesh. Keeping %s only.\n", newEltType);
      }
    }

    vector<char*> newEltTypes;
    K_ARRAY::extractVars(newEltType, newEltTypes);

    // Remplissage table d'indirection et nombre d'elements par eltType
    for (E_Int k = 0; k < nu; k++)
    {
      vector<char*> eltTypesk;
      K_ARRAY::extractVars(eltType[k], eltTypesk);

      E_Int nck = cn[k]->getNConnect();
      cumnepcOrig[k].resize(nck+1); cumnepcOrig[k][0] = 0;

      for (E_Int ic = 0; ic < nck; ic++)
      {
        if (indir[k][ic] == -2) continue; // skip
        for (E_Int icglb = 0; icglb < nc; icglb++)
        {
          if (K_STRING::cmp(newEltTypes[icglb], eltTypesk[ic]) == 0)
          {indir[k][ic] = icglb; break;}
        }
        if (indir[k][ic] < 0) continue; // skip
        FldArrayI& cmkic = *(cn[k]->getConnect(ic));
        E_Int neltskic = cmkic.getSize();
        nepc[indir[k][ic]] += neltskic;
        cumnepcOrig[k][ic+1] = cumnepcOrig[k][ic] + neltskic;
        ntotElts += neltskic;
      }
      for (size_t ic = 0; ic < eltTypesk.size(); ic++) delete [] eltTypesk[ic];
    }
    for (size_t ic = 0; ic < newEltTypes.size(); ic++) delete [] newEltTypes[ic];

    for (E_Int ic = 1; ic < nc; ic++) cumnepc[ic] = cumnepc[ic-1] + nepc[ic-1];

    tpln = K_ARRAY::buildArray3(nfld, unstructVarString[nuref], npts, nepc,
                                newEltType, false, api);
    K_ARRAY::getFromArray3(tpln, f2, cn2);
  }

  // Nouveaux champs aux centres (la connectivite sera identique a cn2)
  FldArrayF* fc = NULL;
  if (arraysc != NULL)
  {
    E_Bool compact = (api == 1) ? true : false;
    nfldc = unstructFc[nuref]->getNfld();
    fc = new FldArrayF(ntotElts, nfldc, compact);
  }

  // Acces non universel sur les ptrs NGON
  E_Int *ngon = NULL, *nface = NULL, *indPG = NULL, *indPH = NULL;
  if (K_STRING::cmp(eltRef, "NGON") == 0)
  {
    ngon = cn2->getNGon(); nface = cn2->getNFace();
    if (ngonType == 2 || ngonType == 3)
    {
      indPG = cn2->getIndPG(); indPH = cn2->getIndPH();
    }
  }

  #pragma omp parallel
  {
    E_Int offsetSizeFN = 0, offsetSizeEF = 0;
    E_Int offV = 0, offsetFaces = 0, offsetElts = 0;
    std::vector<E_Int> offE(nc, 0);  // element offset per output conn.
    E_Int nptsk, nck, neltskic, neltsk, nfacesk, sizeFNk, sizeEFk, ic2;
    E_Int offOrig, offDst;

    for (E_Int k = 0; k < nu; k++)
    {
      // Skip if the ref elt type is NGON and if current elt type is not NGON
      // NB: Dissimilar BE elt types can be combined to form ME
      if (K_STRING::cmp(eltRef, "NGON") == 0 and indir[k][0] < 0) continue;

      nptsk = unstructF[k]->getSize();
      // Copie des champs aux noeuds
      for (E_Int n = 1; n <= nfld; n++)
      {
        E_Float* fkn = unstructF[k]->begin(n);
        E_Float* fn = f2->begin(n);
        #pragma omp for nowait
        for (E_Int i = 0; i < nptsk; i++) fn[i+offV] = fkn[i];
      }

      if (K_STRING::cmp(eltRef, "NGON") == 0)
      {
        // Copie des champs aux centres
        for (E_Int n = 1; n <= nfldc; n++)
        {
          E_Float* fckn = unstructFc[k]->begin(n);
          E_Float* fcn = fc->begin(n);
          #pragma omp for nowait
          for (E_Int i = 0; i < unstructFc[k]->getSize(); i++)
            fcn[i+offsetElts] = fckn[i];
        }

        neltsk = cn[k]->getNElts();
        nfacesk = cn[k]->getNFaces();
        sizeFNk = cn[k]->getSizeNGon();
        sizeEFk = cn[k]->getSizeNFace();

        // Ajout de la connectivite NGON k
        E_Int *ngonk = cn[k]->getNGon(), *nfacek = cn[k]->getNFace();
        E_Int *indPGk = NULL, *indPHk = NULL;

        #pragma omp for nowait
        for (E_Int i = 0; i < sizeFNk; i++)
          ngon[i+offsetSizeFN] = ngonk[i] + offV;
        #pragma omp for nowait
        for (E_Int i = 0; i < sizeEFk; i++)
          nface[i+offsetSizeEF] = nfacek[i] + offsetFaces;

        if (ngonType == 2 || ngonType == 3)
        {
          indPGk = cn[k]->getIndPG(); indPHk = cn[k]->getIndPH();
          #pragma omp for nowait
          for (E_Int i = 0; i < nfacesk; i++)
            indPG[i+offsetFaces] = indPGk[i] + offsetSizeFN;
          #pragma omp for
          for (E_Int i = 0; i < neltsk; i++)
            indPH[i+offsetElts] = indPHk[i] + offsetSizeEF;
        }

        // Increment NGON offsets
        offsetFaces += nfacesk;
        offsetElts += neltsk;
        offsetSizeFN += sizeFNk;
        offsetSizeEF += sizeEFk;
      }
      else if (dim != 0)  // skip NODE
      {
        // Ajout de la connectivite BE/ME k
        nck = cn[k]->getNConnect();
        for (E_Int ic = 0; ic < nck; ic++)
        {
          ic2 = indir[k][ic];  // output conn. position
          if (ic2 < 0) continue;  // skip
          FldArrayI& cmkic = *(cn[k]->getConnect(ic));
          FldArrayI& cm = *(cn2->getConnect(ic2));
          neltskic = cmkic.getSize();

          #pragma omp for
          for (E_Int i = 0; i < neltskic; i++)
            for (E_Int j = 1; j <= cmkic.getNfld(); j++)
              cm(i+offE[ic2],j) = cmkic(i,j) + offV;

          // Copie des champs aux centres pour cette connectivite d'input
          offOrig = cumnepcOrig[k][ic];
          offDst = cumnepc[ic2] + offE[ic2];

          for (E_Int n = 1; n <= nfldc; n++)
          {
            E_Float* fckn = unstructFc[k]->begin(n);
            E_Float* fcn = fc->begin(n);
            #pragma omp for nowait
            for (E_Int i = 0; i < neltskic; i++)
              fcn[i+offDst] = fckn[i+offOrig];
          }

          // Increment offsets
          offE[ic2] += neltskic;
        }
      }
      offV += nptsk;
    }
  }

  // NGON: Correction for number of vertices per face and number of faces per
  // element for all but the first valid array nuref
  if (K_STRING::cmp(eltRef, "NGON") == 0 and shift == 1)
  {
    E_Int offsetSizeFN = cn[nuref]->getSizeNGon();
    E_Int offsetSizeEF = cn[nuref]->getSizeNFace();
    for (E_Int k = 0; k < nu; k++)
    {
      if (indir[k][0] < 0 || k == nuref) continue;
      E_Int ind = 0;
      E_Int *ngonk = cn[k]->getNGon(), *nfacek = cn[k]->getNFace();
      for (E_Int i = 0; i < cn[k]->getNFaces(); i++)
      {
        ngon[offsetSizeFN+ind] = ngonk[ind];
        ind += ngonk[ind]+1;
      }

      ind = 0;
      for (E_Int i = 0; i < cn[k]->getNElts(); i++)
      {
        nface[offsetSizeEF+ind] = nfacek[ind];
        ind += nfacek[ind]+1;
      }

      offsetSizeFN += cn[k]->getSizeNGon();
      offsetSizeEF += cn[k]->getSizeNFace();
    }
  }

  for (E_Int k = 0; k < nu; k++) RELEASESHAREDU(obju[k], unstructF[k], cn[k]);
  for (E_Int k = 0; k < nuc; k++) RELEASESHAREDU(objuc[k], unstructFc[k], cnc[k]);

  E_Int posx = K_ARRAY::isCoordinateXPresent(unstructVarString[nuref])+1;
  E_Int posy = K_ARRAY::isCoordinateYPresent(unstructVarString[nuref])+1;
  E_Int posz = K_ARRAY::isCoordinateZPresent(unstructVarString[nuref])+1;

  PyObject* tpln2 = NULL;
  if (posx > 0 && posy > 0 && posz > 0)
  {
    // Do not remove degenerated nor duplicated elements - AMR
    tpln2 = K_CONNECT::V_cleanConnectivity(
      unstructVarString[nuref], *f2, *cn2, newEltType, tol,
      true, true, true, false, true, false
    );
  }

  if (arraysc == NULL)
  {
    RELEASESHAREDB(res2, tpln, f2, cn2);
    if (res2 == 1) delete cn2;
    if (tpln2 == NULL) return tpln;
    else { Py_DECREF(tpln); return tpln2; };
  }
  else
  {
    PyObject* l = PyList_New(0);
    PyObject* tplc = NULL;
    char newEltTypec[K_ARRAY::VARSTRINGLENGTH];
    K_ARRAY::starVarString(newEltType, newEltTypec);

    if (tpln2 == NULL)
    {
      tplc = K_ARRAY::buildArray3(*fc, unstructVarStringc[nuref],
                                  *cn2, newEltTypec, api);
      PyList_Append(l, tpln); PyList_Append(l, tplc);
    }
    else
    {
      PyList_Append(l, tpln2);
      FldArrayF* fout; FldArrayI* cnout;
      K_ARRAY::getFromArray3(tpln2, fout, cnout);
      tplc = K_ARRAY::buildArray3(*fc, unstructVarStringc[nuref],
                                  *cnout, newEltTypec, api);
      PyList_Append(l, tplc);
      RELEASESHAREDU(tpln2, fout, cnout); 
      Py_DECREF(tpln2);
    }
    
    RELEASESHAREDB(res2, tpln, f2, cn2);
    if (res2 == 1) delete cn2;
    Py_DECREF(tpln); Py_DECREF(tplc); delete fc;
    return l;
  }
}

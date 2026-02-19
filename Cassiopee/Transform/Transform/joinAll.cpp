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
/* Join all arrays */
// ============================================================================
PyObject* K_TRANSFORM::joinAll(PyObject* self, PyObject* args)
{
  PyObject* arrays; E_Float tol;
  if (!PYPARSETUPLE_(args, O_ R_, &arrays, &tol)) return NULL;

  // Check arrays
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

  // Fusion des zones non-structures
  PyObject* tpl = NULL;
  E_Int nu = unstructF.size();
  if (nu == 0) return NULL;
  char* eltRef = eltType[0];
  E_Int missed = 0;

  E_Int nc = 0, dimRef = -1, dim;
  char newEltType[K_ARRAY::VARSTRINGLENGTH]; newEltType[0] = '\0';
  // Counters for all arrays
  E_Int npts = 0;
  E_Int nfaces = 0, neltsNGON = 0, sizeFN = 0, sizeEF = 0;
  // Table d'indirection des connectivites ME
  vector<vector<E_Int> > indir(nu);

  for (E_Int k = 0; k < nu; k++)
  {
    npts += unstructF[k]->getSize();

    if (strcmp(eltType[k], "NGON") == 0)
    {
      // La connectivite fusionee ne doit avoir que des NGONs
      if (strcmp(eltRef, "NGON") == 0)
      {
        neltsNGON += cn[k]->getNElts();
        nfaces += cn[k]->getNFaces();
        sizeFN += cn[k]->getSizeNGon();
        sizeEF += cn[k]->getSizeNFace();
      }
      else missed++;
    }
    else if (strcmp(eltRef, "NGON") != 0)
    {
      // Calcul du nombre d'elt types dans la connectivite ME fusionee
      // et de leur identite
      vector<char*> eltTypesk;
      K_ARRAY::extractVars(eltType[k], eltTypesk);

      E_Int nck = cn[k]->getNConnect();
      for (E_Int ic = 0; ic < nck; ic++)
      {
        char* eltTypConn = eltTypesk[ic];
        // Check dimensionality: allow merge if identical
        if (dimRef == -1)
        {
          if (strcmp(eltTypesk[ic], "NODE") == 0) dimRef = 0;
          else if (strcmp(eltTypesk[ic], "BAR") == 0) dimRef = 1;
          else if (strcmp(eltTypesk[ic], "TRI") == 0 or
                   strcmp(eltTypesk[ic], "QUAD") == 0) dimRef = 2;
          else dimRef = 3;
        }
        else
        {
          if (strcmp(eltTypesk[ic], "NODE") == 0) dim = 0;
          else if (strcmp(eltTypesk[ic], "BAR") == 0) dim = 1;
          else if (strcmp(eltTypesk[ic], "TRI") == 0 or
                   strcmp(eltTypesk[ic], "QUAD") == 0) dim = 2;
          else dim = 3;
          if (dim != dimRef) { indir[k].push_back(-2); missed++; continue; }
        }
        // Add default value in mapping table
        indir[k].push_back(-1);
        // Concatenate elttypes, discard duplicates
        if (strstr(newEltType, eltTypConn) == NULL)
        {
          strcat(newEltType, eltTypConn); strcat(newEltType, ",");
          nc += 1;
        }
      }

      for (size_t ic = 0; ic < eltTypesk.size(); ic++)
        delete [] eltTypesk[ic];
    }
    else missed++;
  }

  if (missed > 0)
    printf("Warning: joinAll: some arrays cannot be joined: different element "
           "types.\n");

  // Build unstructured connectivity
  E_Int nfld = unstructF[0]->getNfld();
  E_Int api = unstructF[0]->getApi();
  vector<E_Int> neltsME;

  if (strcmp(eltRef, "NGON") == 0)
  {
    strcat(newEltType, "NGON");
    E_Int ngonType = cn[0]->getNGonType();
    tpl = K_ARRAY::buildArray3(nfld, unstructVarString[0], npts, neltsNGON,
                               nfaces, newEltType, sizeFN, sizeEF,
                               ngonType, false, api);
  }
  else
  {
    // Remove trailing comma in newEltType
    E_Int len = strlen(newEltType);
    newEltType[len-1] = '\0';

    if (nc == 2 && dimRef == 3)
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
    neltsME.resize(nc); neltsME.assign(nc, 0);
    for (E_Int k = 0; k < nu; k++)
    {
      vector<char*> eltTypesk;
      K_ARRAY::extractVars(eltType[k], eltTypesk);

      E_Int nck = cn[k]->getNConnect();
      for (E_Int ic = 0; ic < nck; ic++)
      {
        if (indir[k][ic] == -2) continue; // skip
        for (E_Int icglb = 0; icglb < nc; icglb++)
        {
          if (strcmp(newEltTypes[icglb], eltTypesk[ic]) == 0)
          {indir[k][ic] = icglb; break;}
        }
        if (indir[k][ic] < 0) continue; // skip
        FldArrayI& cmkic = *(cn[k]->getConnect(ic));
        neltsME[indir[k][ic]] += cmkic.getSize();
      }

      for (size_t ic = 0; ic < eltTypesk.size(); ic++)
        delete [] eltTypesk[ic];
    }
    for (size_t ic = 0; ic < newEltTypes.size(); ic++)
        delete [] newEltTypes[ic];
    tpl = K_ARRAY::buildArray3(nfld, unstructVarString[0], npts, neltsME,
                               newEltType, false, api);
  }
  FldArrayF* f; FldArrayI* cno;
  K_ARRAY::getFromArray3(tpl, f, cno);

  // Acces non universel sur les ptrs NGON
  E_Int *ngon = NULL, *nface = NULL, *indPG = NULL, *indPH = NULL;
  E_Int ngonType = cno->getNGonType();
  if (strcmp(eltRef, "NGON") == 0)
  {
    ngon = cno->getNGon(); nface = cno->getNFace();
    if (ngonType == 2 || ngonType == 3)
    {
      indPG = cno->getIndPG(); indPH = cno->getIndPH();
    }
  }

  #pragma omp parallel
  {
    E_Int offsetSizeFN = 0, offsetSizeEF = 0;
    E_Int offsetPts = 0, offsetFaces = 0, offsetElts = 0;
    for (E_Int k = 0; k < nu; k++)
    {
      // Skip if the ref elt type is NGON and if current elt type is not NGON
      // NB: Dissimilar BE elt types can be combined to form ME
      if (strcmp(eltRef, "NGON") == 0 and strcmp(eltRef, eltType[k]) != 0)
        continue;

      E_Int nptsk = unstructF[k]->getSize();
      // Copie des champs aux noeuds
      for (E_Int n = 1; n <= nfld; n++)
      {
        E_Float* fkn = unstructF[k]->begin(n);
        E_Float* fn = f->begin(n);
        #pragma omp for
        for (E_Int i = 0; i < nptsk; i++) fn[i+offsetPts] = fkn[i];
      }

      if (strcmp(eltRef, "NGON") == 0)
      {
        E_Int neltsk = cn[k]->getNElts();
        E_Int nfacesk = cn[k]->getNFaces();
        E_Int sizeFNk = cn[k]->getSizeNGon();
        E_Int sizeEFk = cn[k]->getSizeNFace();

        // Ajout de la connectivite NGON k
        E_Int *ngonk = cn[k]->getNGon(), *nfacek = cn[k]->getNFace();
        E_Int *indPGk = NULL, *indPHk = NULL;

        #pragma omp for
        for (E_Int i = 0; i < sizeFNk; i++)
          ngon[i+offsetSizeFN] = ngonk[i] + offsetPts;
        #pragma omp for
        for (E_Int i = 0; i < sizeEFk; i++)
          nface[i+offsetSizeEF] = nfacek[i] + offsetFaces;

        if (ngonType == 2 || ngonType == 3)
        {
          indPGk = cn[k]->getIndPG(); indPHk = cn[k]->getIndPH();
          #pragma omp for
          for (E_Int i = 0; i < nfacesk; i++)
            indPG[i+offsetFaces] = indPGk[i] + offsetFaces;
          #pragma omp for
          for (E_Int i = 0; i < neltsk; i++)
            indPH[i+offsetElts] = indPHk[i] + offsetElts;
        }

        // Increment NGON offsets
        offsetFaces += nfacesk;
        offsetElts += neltsk;
        offsetSizeFN += sizeFNk;
        offsetSizeEF += sizeEFk;
      }
      else
      {
        // Ajout de la connectivite BE/ME k
        E_Int nck = cn[k]->getNConnect();
        for (E_Int ic = 0; ic < nck; ic++)
        {
          if (indir[k][ic] < 0) continue; // skip
          FldArrayI& cmkic = *(cn[k]->getConnect(ic));
          FldArrayI& cm = *(cno->getConnect(indir[k][ic]));
          E_Int neltskic = cmkic.getSize();

          #pragma omp for
          for (E_Int i = 0; i < neltskic; i++)
            for (E_Int j = 1; j <= cmkic.getNfld(); j++)
              // Add offsets
              cm(i+offsetElts,j) = cmkic(i,j) + offsetPts;

          // Increment ME offsets
          offsetElts += neltskic;
        }
      }
      offsetPts += nptsk;
    }
  }

  // NGON: Correction for number of vertices per face and number of faces per
  // element for all but the first array
  if (strcmp(eltRef, "NGON") == 0 and api != 3)
  {
    E_Int offsetSizeFN = cn[0]->getSizeNGon();
    E_Int offsetSizeEF = cn[0]->getSizeNFace();
    for (E_Int k = 1; k < nu; k++)
    {
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

  for (E_Int k = 0; k < nu; k++)
    RELEASESHAREDU(obju[k], unstructF[k], cn[k]);

  E_Int posx = K_ARRAY::isCoordinateXPresent(unstructVarString[0])+1;
  E_Int posy = K_ARRAY::isCoordinateYPresent(unstructVarString[0])+1;
  E_Int posz = K_ARRAY::isCoordinateZPresent(unstructVarString[0])+1;
  if (posx > 0 && posy > 0 && posz > 0)
  {
    K_CONNECT::cleanConnectivity(posx, posy, posz, tol, newEltType,
                                 *f, *cno);
    PyObject* tpl2 = K_ARRAY::buildArray3(*f, unstructVarString[0], *cno, newEltType, api);
    // PyObject* tpl2 = K_CONNECT::V_cleanConnectivity(
    //   unstructVarString[0], *f, *cno, newEltType, tol);
    RELEASESHAREDU(tpl, f, cno); Py_DECREF(tpl);
    return tpl2;
  }
  else
  {
    RELEASESHAREDU(tpl, f, cno);
    return tpl;
  }
}

// ============================================================================
/* Join all arrays and their arrays at centers */
// ============================================================================
PyObject* K_TRANSFORM::joinAllBoth(PyObject* self, PyObject* args)
{
  PyObject *arrays, *arraysc; E_Float tol;
  if (!PYPARSETUPLE_(args, OO_ R_, &arrays, &arraysc, &tol)) return NULL;

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

  // Check arrays for fields located at centers
  vector<E_Int> resc;
  vector<char*> structVarStringc; vector<char*> unstructVarStringc;
  vector<FldArrayF*> structFc; vector<FldArrayF*> unstructFc;
  vector<E_Int> nic, njc, nkc;
  vector<FldArrayI*> cnc; vector<char*> eltTypec;
  vector<PyObject*> objsc, objuc;
  K_ARRAY::getFromArrays(arraysc, resc, structVarStringc,
			                   unstructVarStringc, structFc,
			                   unstructFc, nic, njc, nkc,
			                   cnc, eltTypec, objsc, objuc,
			                   false, false, false, false, true);

  // Fusion des zones non-structures
  PyObject* tpln = NULL;
  E_Int nu = unstructF.size(); E_Int nuc = unstructFc.size();
  if (nu != nuc or (nu == 0 && nuc == 0))
  {
    for (E_Int k = 0; k < nu; k++) RELEASESHAREDU(obju[k], unstructF[k], cn[k]);
    for (E_Int k = 0; k < nuc; k++) RELEASESHAREDU(objuc[k], unstructFc[k], cnc[k]);
    PyErr_SetString(PyExc_ValueError,
                    "joinAllBoth: number of arrays at nodes and centers must "
                    "be equal.");
    return NULL;
  }

  char* eltRef = NULL;
  eltRef = eltType[0];
  E_Int missed = 0;

  E_Int nc = 0, dimRef = -1, dim;
  char newEltType[K_ARRAY::VARSTRINGLENGTH]; newEltType[0] = '\0';
  // Counters for all arrays
  E_Int npts = 0;
  E_Int nfaces = 0, neltsNGON = 0, sizeFN = 0, sizeEF = 0;
  // Table d'indirection des connectivites ME
  vector<vector<E_Int> > indir(nu);

  for (E_Int k = 0; k < nu; k++)
  {
    npts += unstructF[k]->getSize();

    if (strcmp(eltType[k], "NGON") == 0)
    {
      // La connectivite fusionee ne doit avoir que des NGONs
      if (strcmp(eltRef, "NGON") == 0)
      {
        neltsNGON += cn[k]->getNElts();
        nfaces += cn[k]->getNFaces();
        sizeFN += cn[k]->getSizeNGon();
        sizeEF += cn[k]->getSizeNFace();
      }
      else missed++;
    }
    else if (strcmp(eltRef, "NGON") != 0)
    {
      // Calcul du nombre d'elt types dans la connectivite ME fusionee
      // et de leur identite
      vector<char*> eltTypesk;
      K_ARRAY::extractVars(eltType[k], eltTypesk);

      E_Int nck = cn[k]->getNConnect();
      for (E_Int ic = 0; ic < nck; ic++)
      {
        char* eltTypConn = eltTypesk[ic];
        // Check dimensionality: allow merge if identical
        if (dimRef == -1)
        {
          if (strcmp(eltTypesk[ic], "NODE") == 0) dimRef = 0;
          else if (strcmp(eltTypesk[ic], "BAR") == 0) dimRef = 1;
          else if (strcmp(eltTypesk[ic], "TRI") == 0 or
                   strcmp(eltTypesk[ic], "QUAD") == 0) dimRef = 2;
          else dimRef = 3;
        }
        else
        {
          if (strcmp(eltTypesk[ic], "NODE") == 0) dim = 0;
          else if (strcmp(eltTypesk[ic], "BAR") == 0) dim = 1;
          else if (strcmp(eltTypesk[ic], "TRI") == 0 or
                   strcmp(eltTypesk[ic], "QUAD") == 0) dim = 2;
          else dim = 3;
          if (dim != dimRef) { indir[k].push_back(-2); missed++; continue; }
        }
        // Add default value in mapping table
        indir[k].push_back(-1);
        // Concatenate elttypes, discard duplicates
        if (strstr(newEltType, eltTypConn) == NULL)
        {
          strcat(newEltType, eltTypConn); strcat(newEltType, ",");
          nc += 1;
        }
      }

      for (size_t ic = 0; ic < eltTypesk.size(); ic++)
        delete [] eltTypesk[ic];
    }
    else missed++;
  }

  if (missed > 0)
    printf("Warning: joinAll: some arrays cannot be joined: different element "
           "types.\n");

  // Build unstructured connectivity
  E_Int nfld = unstructF[0]->getNfld();
  E_Int nfldc = unstructFc[0]->getNfld();
  E_Int api = unstructF[0]->getApi();
  vector<E_Int> neltsME; E_Int neltstot = 0;

  if (strcmp(eltRef, "NGON") == 0)
  {
    strcpy(newEltType, "NGON");
    neltstot = neltsNGON;
    E_Int ngonType = cn[0]->getNGonType();
    tpln = K_ARRAY::buildArray3(nfld, unstructVarString[0], npts, neltsNGON,
                                nfaces, newEltType, sizeFN, sizeEF,
                                ngonType, false, api);
  }
  else
  {
    // Remove trailing comma in newEltType
    E_Int len = strlen(newEltType);
    newEltType[len-1] = '\0';

    vector<char*> newEltTypes;
    K_ARRAY::extractVars(newEltType, newEltTypes);

    // Remplissage table d'indirection et nombre d'elements par eltType
    neltsME.resize(nc); neltsME.assign(nc, 0);
    for (E_Int k = 0; k < nu; k++)
    {
      vector<char*> eltTypesk;
      K_ARRAY::extractVars(eltType[k], eltTypesk);

      E_Int nck = cn[k]->getNConnect();
      for (E_Int ic = 0; ic < nck; ic++)
      {
        if (indir[k][ic] == -2) continue; // skip
        for (E_Int icglb = 0; icglb < nc; icglb++)
        {
          if (strcmp(newEltTypes[icglb], eltTypesk[ic]) == 0)
          {indir[k][ic] = icglb; break;}
        }
        if (indir[k][ic] < 0) continue; // skip
        FldArrayI& cmkic = *(cn[k]->getConnect(ic));
        neltsME[indir[k][ic]] += cmkic.getSize();
        neltstot += cmkic.getSize();
      }

      for (size_t ic = 0; ic < eltTypesk.size(); ic++)
        delete [] eltTypesk[ic];
    }
    for (size_t ic = 0; ic < newEltTypes.size(); ic++)
        delete [] newEltTypes[ic];
    tpln = K_ARRAY::buildArray3(nfld, unstructVarString[0], npts, neltsME,
                                newEltType, false, api);
  }
  FldArrayF* f; FldArrayI* cno;
  K_ARRAY::getFromArray3(tpln, f, cno);

  // Nouveaux champs aux centres (la connectivite sera identique a cno)
  E_Bool compact = false;
  if (api == 1) compact = true;
  FldArrayF* fc = new FldArrayF(neltstot, nfldc, compact);

  // Acces non universel sur les ptrs NGON
  E_Int *ngon = NULL, *nface = NULL, *indPG = NULL, *indPH = NULL;
  E_Int ngonType = cno->getNGonType();
  if (strcmp(eltRef, "NGON") == 0)
  {
    ngon = cno->getNGon(); nface = cno->getNFace();
    if (ngonType == 2 || ngonType == 3)
    {
      indPG = cno->getIndPG(); indPH = cno->getIndPH();
    }
  }

  #pragma omp parallel
  {
    E_Int offsetSizeFN = 0, offsetSizeEF = 0;
    E_Int offsetPts = 0, offsetFaces = 0, offsetElts = 0;
    for (E_Int k = 0; k < nu; k++)
    {
      // Skip if the ref elt type is NGON and if current elt type is not NGON
      // NB: Dissimilar BE elt types can be combined to form ME
      if (strcmp(eltRef, "NGON") == 0 and strcmp(eltRef, eltType[k]) != 0)
        continue;

      E_Int nptsk = unstructF[k]->getSize();
      // Copie des champs aux noeuds
      for (E_Int n = 1; n <= nfld; n++)
      {
        E_Float* fkn = unstructF[k]->begin(n);
        E_Float* fn = f->begin(n);
        #pragma omp for nowait
        for (E_Int i = 0; i < nptsk; i++) fn[i+offsetPts] = fkn[i];
      }

      // Copie des champs aux centres
      for (E_Int n = 1; n <= nfldc; n++)
      {
        E_Float* fckn = unstructFc[k]->begin(n);
        E_Float* fcn = fc->begin(n);
        #pragma omp for nowait
        for (E_Int i = 0; i < unstructFc[k]->getSize(); i++)
          fcn[i+offsetElts] = fckn[i];
      }

      if (strcmp(eltRef, "NGON") == 0)
      {
        E_Int neltsk = cn[k]->getNElts();
        E_Int nfacesk = cn[k]->getNFaces();
        E_Int sizeFNk = cn[k]->getSizeNGon();
        E_Int sizeEFk = cn[k]->getSizeNFace();

        // Ajout de la connectivite NGON k
        E_Int *ngonk = cn[k]->getNGon(), *nfacek = cn[k]->getNFace();
        E_Int *indPGk = NULL, *indPHk = NULL;

        #pragma omp for nowait
        for (E_Int i = 0; i < sizeFNk; i++)
          ngon[i+offsetSizeFN] = ngonk[i] + offsetPts;
        #pragma omp for nowait
        for (E_Int i = 0; i < sizeEFk; i++)
          nface[i+offsetSizeEF] = nfacek[i] + offsetFaces;

        if (ngonType == 2 || ngonType == 3)
        {
          indPGk = cn[k]->getIndPG(); indPHk = cn[k]->getIndPH();
          #pragma omp for nowait
          for (E_Int i = 0; i < nfacesk; i++)
            indPG[i+offsetFaces] = indPGk[i] + offsetFaces;
          #pragma omp for
          for (E_Int i = 0; i < neltsk; i++)
            indPH[i+offsetElts] = indPHk[i] + offsetElts;
        }

        // Increment NGON offsets
        offsetFaces += nfacesk;
        offsetElts += neltsk;
        offsetSizeFN += sizeFNk;
        offsetSizeEF += sizeEFk;
      }
      else
      {
        // Ajout de la connectivite BE/ME k
        E_Int nck = cn[k]->getNConnect();
        for (E_Int ic = 0; ic < nck; ic++)
        {
          if (indir[k][ic] < 0) continue; // skip
          FldArrayI& cmkic = *(cn[k]->getConnect(ic));
          FldArrayI& cm = *(cno->getConnect(indir[k][ic]));
          E_Int neltskic = cmkic.getSize();

          #pragma omp for
          for (E_Int i = 0; i < neltskic; i++)
            for (E_Int j = 1; j <= cmkic.getNfld(); j++)
              // Add offsets
              cm(i+offsetElts,j) = cmkic(i,j) + offsetPts;

          // Increment ME offsets
          offsetElts += neltskic;
        }
      }
      offsetPts += nptsk;
    }
  }

  // NGON: Correction for number of vertices per face and number of faces per
  // element for all but the first array
  if (strcmp(eltRef, "NGON") == 0 and api != 3)
  {
    E_Int offsetSizeFN = cn[0]->getSizeNGon();
    E_Int offsetSizeEF = cn[0]->getSizeNFace();
    for (E_Int k = 1; k < nu; k++)
    {
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

  for (E_Int k = 0; k < nu; k++)
  {
    RELEASESHAREDU(obju[k], unstructF[k], cn[k]);
    RELEASESHAREDU(objuc[k], unstructFc[k], cnc[k]);
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(unstructVarString[0])+1;
  E_Int posy = K_ARRAY::isCoordinateYPresent(unstructVarString[0])+1;
  E_Int posz = K_ARRAY::isCoordinateZPresent(unstructVarString[0])+1;

  PyObject* l = PyList_New(0);
  if (posx > 0 && posy > 0 && posz > 0)
  {
    K_CONNECT::cleanConnectivity(posx, posy, posz, tol, newEltType, *f, *cno);
    // PyObject* tpln2 = K_CONNECT::V_cleanConnectivity(
    //  unstructVarString[0], *f, *cno, newEltType, tol);
    // PyList_Append(l, tpln2); Py_DECREF(tpln2);
  }
  // else
  // {
  //  PyList_Append(l, tpln);
  // }
  PyObject* tpln2 = K_ARRAY::buildArray3(*f, unstructVarString[0], *cno, newEltType, api);
  PyList_Append(l, tpln2); Py_DECREF(tpln2);

  char newEltTypec[K_ARRAY::VARSTRINGLENGTH];
  K_ARRAY::starVarString(newEltType, newEltTypec);
  PyObject* tplc = K_ARRAY::buildArray3(*fc, unstructVarStringc[0],
                                        *cno, newEltTypec, api);
  PyList_Append(l, tplc); Py_DECREF(tplc); delete fc;
  RELEASESHAREDU(tpln, f, cno);
  return l;
}

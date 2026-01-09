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

using namespace K_FLD;
using namespace std;

//=============================================================================
/*
   splitSharpEdges:
   Decoupe un array TRI, QUAD, BAR en partie lisses (angle entre elements
   inferieur a alphaRef)
   La connectivite doit etre propre.
*/
//=============================================================================
PyObject* K_TRANSFORM::splitSharpEdges(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Float alphaRef;
  E_Float dirVect[3];
  dirVect[0] = 0.; dirVect[1] = 0.; dirVect[2] = 1.;

  if (!PYPARSETUPLE_(args, O_ R_, &array, &alphaRef))
  {
    return NULL;
  }

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res =
    K_ARRAY::getFromArray3(array, varString, f, im, jm, km, cn, eltType);

  if (res == 1)
  {
    RELEASESHAREDS(array, f);
    PyErr_SetString(PyExc_TypeError,
                    "splitSharpEdges: cannot be used on a structured array.");
    return NULL;
  }
  if (res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "splitSharpEdges: unknown type of array.");
    return NULL;
  }
  if (K_STRING::cmp(eltType, "NODE") == 0)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "splitSharpEdges: cannot be used on a NODE-array.");
    return NULL;
  }
  if (K_STRING::cmp(eltType, "PYRA") == 0)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "splitSharpEdges: cannot be used on a PYRA-array.");
    return NULL;
  }
  if (K_STRING::cmp(eltType, "PENTA") == 0)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "splitSharpEdges: cannot be used on a PENTA-array.");
    return NULL;
  }
  if (K_STRING::cmp(eltType, "HEXA") == 0)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "splitSharpEdges: cannot be used on a HEXA-array.");
    return NULL;
  }
  if (K_STRING::cmp(eltType, "MIX") == 0)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "splitSharpEdges: cannot be used on a MIX-array.");
    return NULL;
  }
  E_Int type = 0;
  if (K_STRING::cmp(eltType, "BAR") == 0) type = 0;
  else if (K_STRING::cmp(eltType, "TRI") == 0) type = 1;
  else if (K_STRING::cmp(eltType, "QUAD") == 0) type = 2;
  else if (K_STRING::cmp(eltType, "NGON") == 0) type = 3;
  else
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "splitSharpEdges: can not be used on this type of array.");
    return NULL;
  }

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);

  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "splitSharpEdges: can't find coordinates in array.");
    return NULL;
  }
  posx++; posy++; posz++;

  PyObject* out;
  if (type == 3)
    out = splitSharpEdgesNGon(f, cn, varString,
                              posx, posy, posz, alphaRef, dirVect);
  else out = splitSharpEdgesBasics(f, cn, eltType, varString,
                                   posx, posy, posz, type, alphaRef, dirVect);

  RELEASESHAREDU(array, f, cn);
  return out;
}

//=============================================================================
// Split sharp edges pour un array a elements basiques (BAR,TRI,QUAD)
//=============================================================================
PyObject* K_TRANSFORM::splitSharpEdgesBasics(
  FldArrayF* f, FldArrayI* cn,
  char* eltType, char* varString, E_Int posx, E_Int posy, E_Int posz,
  E_Int type, E_Float alphaRef, E_Float* dirVect)
{
  E_Float* x = f->begin(posx);
  E_Float* y = f->begin(posy);
  E_Float* z = f->begin(posz);

  vector< vector<E_Int> > cEEN(cn->getSize());
  K_CONNECT::connectEV2EENbrs(eltType, f->getSize(), *cn, cEEN);

  E_Int api = f->getApi();
  E_Int nt = cn->getNfld();
  E_Int ne = cn->getSize(); // nbre d'elements
  E_Int nev = 0; // nbre d'elements deja visites
  char* isVisited = (char*)calloc(ne, sizeof(char)); // elt deja visite?
  E_Int* mustBeVisited = (E_Int*)malloc(ne * sizeof(E_Int));
  E_Int mbv, p, i, ie, elt, curr;
  unsigned int iv;
  vector<FldArrayI*> components;

  E_Float alpha;
  E_Int ind;
  E_Float pts1[4][3]; E_Float pts2[4][3];

  mbv = 0;
  while (nev < ne)
  {
    // Recherche le premier elt pas encore visite
    for (p = 0; (isVisited[p] != 0); p++);

    // C'est un nouveau composant
    FldArrayI* c = new FldArrayI(ne, nt);

    mustBeVisited[mbv] = p;
    mbv++; nev++;
    isVisited[p] = 1;
    curr = 0;

    while (mbv > 0)
    {
      mbv--;
      elt = mustBeVisited[mbv];
      for (i = 0; i < nt; i++)
      {
        ind = (*cn)(elt,i+1)-1;
        pts1[i][0] = x[ind]; pts1[i][1] = y[ind]; pts1[i][2] = z[ind];
        (*c)(curr,i+1) = ind+1;
      }
      curr++;

      for (iv = 0; iv < cEEN[elt].size(); iv++)
      {
        ie = cEEN[elt][iv];
        if (isVisited[ie] == 0)
        {
          // Calcul de alpha
          for (i = 0; i < nt; i++)
          {
            ind = (*cn)(ie,i+1)-1;
            pts2[i][0] = x[ind]; pts2[i][1] = y[ind]; pts2[i][2] = z[ind];
          }
          switch (type)
          {
            case 0:
              alpha = K_COMPGEOM::getAlphaAngleBetweenBars(
              pts1[0], pts1[1],
              pts2[0], pts2[1], dirVect);
              break;

            case 1:
              alpha = K_COMPGEOM::getAlphaAngleBetweenTriangles(
                pts1[0], pts1[1], pts1[2],
                pts2[0], pts2[1], pts2[2]); break;

            case 2:
              alpha = K_COMPGEOM::getAlphaAngleBetweenQuads(
                pts1[0], pts1[1], pts1[2], pts1[3],
                pts2[0], pts2[1], pts2[2], pts2[3]); break;

            default: alpha = 180.; break;
          }

          if (alpha == -1000. || K_FUNC::E_abs(alpha-180.) <= alphaRef)
          {
            mustBeVisited[mbv] = ie;
            mbv++; nev++;
            isVisited[ie] = 1;
          }
        }
      }
    }
    c->reAllocMat(curr, nt);
    components.push_back(c);
  }
  free(isVisited);
  free(mustBeVisited);

  // Formation des arrays de sortie + cleanConnectivity
  PyObject* tpl;
  PyObject* l = PyList_New(0);
  E_Int size = components.size();

  for (i = 0; i < size; i++)
  {
    FldArrayF* f0 = new FldArrayF(*f);
    FldArrayF& fp = *f0;
    FldArrayI& cnp = *components[i];
    K_CONNECT::cleanConnectivity(posx, posy, posz, 1.e-10, eltType,
                                 fp, cnp);
    tpl = K_ARRAY::buildArray3(fp, varString, cnp, eltType, api);
    delete &fp; delete &cnp;
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
    // tpl = K_CONNECT::V_cleanConnectivity(varString, *f, *components[i], eltType, 1.e-10);
    // PyList_Append(l, tpl); Py_DECREF(tpl);
    // delete components[i];
  }
  return l;
}

//=============================================================================
// Split sharp edges pour un array NGON de dim=2 ou dim=1
//==============================================================================
PyObject* K_TRANSFORM::splitSharpEdgesNGon(
  FldArrayF* f, FldArrayI* cn, char* varString,
  E_Int posx, E_Int posy, E_Int posz, E_Float alphaRef, E_Float* dirVect)
{
  E_Float* x = f->begin(posx);
  E_Float* y = f->begin(posy);
  E_Float* z = f->begin(posz);

  E_Int api = f->getApi();
  E_Int* ptr = cn->begin();
  E_Int sf = ptr[1];
  E_Int ne = ptr[2+sf]; // nbre d'elements

  FldArrayI cFE;
  K_CONNECT::connectNG2FE(*cn, cFE);

  E_Int* ptrElts = ptr + (sf+4);
  E_Int se = ptr[3+sf]; // taille connectivite elements/faces
  E_Int nev = 0; // nbre d'elements deja visites
  char* isVisited = (char*)calloc(ne, sizeof(char)); // elt deja visite?
  E_Int* mustBeVisited = (E_Int*)malloc(ne * sizeof(E_Int));
  E_Int mbv, p, i, ie, elt, curr, necurr, lt, iv;
  vector<FldArrayI*> components; // connectivite EF locale

  FldArrayI pos; K_CONNECT::getPosElts(*cn, pos);
  FldArrayI posFaces; K_CONNECT::getPosFaces(*cn, posFaces);

  // Recupere la dim en se basant sur la premiere face
  E_Int dim = 3;
  dim = ptr[2];
  dim = min(dim, E_Int(3));

  // Commence par calculer alpha
  E_Int nfaces = ptr[0];
  FldArrayF alphat(nfaces);
  E_Float* alphap = alphat.begin();
  E_Int* cFE1 = cFE.begin(1);
  E_Int* cFE2 = cFE.begin(2);
#pragma omp parallel
  {
  E_Int e1, e2;
  vector<E_Int> indices;
  vector<E_Float*> pts1; vector<E_Float*> pts2;
  E_Float alpha;
  E_Int ind, nvert;

#pragma omp for
  for (E_Int i = 0; i < nfaces; i++)
  {
    e1 = cFE1[i]-1; e2 = cFE2[i]-1;
    if (e1 == -1) alpha = -1000.;
    else if (e2 == -1) alpha = -1000.;
    else
    {
      K_CONNECT::getVertexIndices(cn->begin(), posFaces.begin(), pos[e1], indices);
      nvert = indices.size();
      pts1.reserve(nvert);
      for (E_Int k = 0; k < nvert; k++)
      {
        E_Float* pt = new E_Float[3];
        ind = indices[k]-1;
        pt[0] = x[ind]; pt[1] = y[ind]; pt[2] = z[ind];
        pts1.push_back(pt);
      }
      K_CONNECT::getVertexIndices(cn->begin(), posFaces.begin(), pos[e2], indices);
      nvert = indices.size();
      pts2.reserve(nvert);
      for (E_Int k = 0; k < nvert; k++)
      {
        E_Float* pt = new E_Float[3];
        ind = indices[k]-1;
        pt[0] = x[ind]; pt[1] = y[ind]; pt[2] = z[ind];
        pts2.push_back(pt);
      }

      if (dim == 2)
      {
        alpha = K_COMPGEOM::getAlphaAngleBetweenPolygons(pts1, pts2);
        //printf("mbv=%d, alpha=%f %f\n", mbv, alpha, alphaRef);
      }
      else if (dim == 1)
      {
        if (pts1.size() == 2 && pts2.size() == 2)
        {
          alpha = K_COMPGEOM::getAlphaAngleBetweenBars(
                    pts1[0], pts1[1],
                    pts2[0], pts2[1], dirVect);
        }
        else alpha = 180.;
      }
      else alpha = 180.;

      nvert = pts1.size();
      for (E_Int k = 0; k < nvert; k++) delete [] pts1[k];
      pts1.clear();
      nvert = pts2.size();
      for (E_Int k = 0; k < nvert; k++) delete [] pts2[k];
      pts2.clear();
    }
    alphap[i] = alpha;
  }
  }

  E_Int e1, e2, nf;
  E_Float alpha;

  // split
  mbv = 0;
  while (nev < ne)
  {
    // Recherche le premier elt pas encore visite
    for (p = 0; (isVisited[p] != 0); p++);

    // C'est un nouveau composant
    FldArrayI* c = new FldArrayI(se+1);
    E_Int* pc = c->begin(); // current pointer

    mustBeVisited[mbv] = p;
    mbv++; nev++;
    isVisited[p] = 1;
    curr = 0; necurr = 0;

    while (mbv > 0)
    {
      mbv--;
      elt = mustBeVisited[mbv];
      // copie elt
      ptrElts = &ptr[pos[elt]];
      lt = ptrElts[0];
      pc[0] = lt;
      for (i = 1; i <= pc[0]; i++) pc[i] = ptrElts[i];
      pc += lt+1;
      curr += lt+1; necurr++;

      for (iv = 0; iv < lt; iv++)
      {
        nf = ptrElts[iv+1]-1;
        e1 = cFE1[nf]-1; e2 = cFE2[nf]-1;
        if (e1 == elt) ie = e2;
        else ie = e1;

        if (ie != -1 && isVisited[ie] == 0)
        {
          alpha = alphap[nf];
          if (alpha == -1000. || K_FUNC::E_abs(alpha-180.) <= alphaRef)
          {
            mustBeVisited[mbv] = ie;
            mbv++; nev++;
            isVisited[ie] = 1;
          }
        }
      }
    }

    pc[0] = necurr; // sentinelle
    c->reAlloc(curr+1);
    components.push_back(c);
  }

  free(isVisited);
  free(mustBeVisited);

  // Formation des arrays de sortie + cleanConnectivity
  PyObject* tpl;
  PyObject* l = PyList_New(0);
  E_Int size = components.size();

  for (i = 0; i < size; i++)
  {
    FldArrayF* f0 = new FldArrayF(*f);
    FldArrayF& fp = *f0;
    // nouvelle connectivite
    E_Int si = components[i]->getSize()-1;
    E_Int* comp = components[i]->begin();
    E_Int size = sf+4+si;
    FldArrayI cnp(size);
    E_Int* cnpp = cnp.begin();
    for (E_Int j = 0; j < sf+2; j++) cnpp[j] = ptr[j];
    cnpp += sf+2;
    cnpp[0] = comp[si]; cnpp[1] = si; cnpp += 2;
    for (E_Int j = 0; j < si; j++) cnpp[j] = comp[j];

    K_CONNECT::cleanConnectivityNGon(posx, posy, posz, 1.e-10,
                                     fp, cnp);
    if (dim == 1) // il semble que dans ce cas, il faut l'appeler 2 fois
      K_CONNECT::cleanConnectivityNGon(posx, posy, posz, 1.e-10,
                                       fp, cnp);
    cnp.setNGonType(1);
    tpl = K_ARRAY::buildArray3(fp, varString, cnp, "NGON", api);
    delete &fp; delete components[i];
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }
  return l;
}

//=============================================================================
// Split sharp edges pour un array NGON de dim=2 ou dim=1
// Input : liste d'index et une sortie liste d'index splittee
// Chaque index est associe a un element
//==============================================================================
PyObject* K_TRANSFORM::splitSharpEdgesList(PyObject* self, PyObject* args)
{
  PyObject* array; PyObject* arrayI;
  E_Float alphaRef;
  if (!PYPARSETUPLE_(args, OO_ R_,
                    &array, &arrayI, &alphaRef))
  {
      return NULL;
  }

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res =
    K_ARRAY::getFromArray3(array, varString, f, im, jm, km, cn, eltType);

  if (res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "splitSharpEdgesList: unknown type of array.");
    return NULL;
  }
  if (res == 1)
  {
    RELEASESHAREDS(array, f);
    PyErr_SetString(PyExc_TypeError,
                    "splitSharpEdgesList: can not be used on a structured array.");
    return NULL;
  }
  if (res == 2 && strcmp(eltType, "NGON") != 0)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "splitSharpEdgesList: only for NGON array.");
    return NULL;
  }

  FldArrayI* indexI;
  res = K_NUMPY::getFromNumpyArray(arrayI, indexI);
  if (res == 0)
  {
    PyErr_SetString(PyExc_TypeError,
                    "splitSharpEdgesList: index numpy is invalid.");
    return NULL;
  }

  E_Float dirVect[3];
  dirVect[0] = 0.; dirVect[1] = 0.; dirVect[2] = 1.;
  E_Int* index = indexI->begin();

  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);

  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "splitSharpEdgesList: can't find coordinates in array.");
    return NULL;
  }
  posx++; posy++; posz++;

  E_Float* x = f->begin(posx);
  E_Float* y = f->begin(posy);
  E_Float* z = f->begin(posz);

  E_Int* ptr = cn->begin();
  E_Int sf = ptr[1];
  E_Int ne = ptr[2+sf]; // nbre d'elements

  FldArrayI cFE;
  K_CONNECT::connectNG2FE(*cn, cFE);

  E_Int* ptrElts = ptr + (sf+4);
  E_Int se = ptr[2+sf]; // nbre d'elements
  E_Int nev = 0; // nbre d'elements deja visites
  char* isVisited = (char*)calloc(ne, sizeof(char)); // elt deja visite?
  E_Int* mustBeVisited = (E_Int*)malloc(ne * sizeof(E_Int));
  E_Int mbv, p, i, ie, elt/*, curr*/, necurr, lt, iv;
  vector<FldArrayI*> components; // liste des index par bloc

  FldArrayI pos; K_CONNECT::getPosElts(*cn, pos);
  FldArrayI posFaces; K_CONNECT::getPosFaces(*cn, posFaces);

  // Recupere la dim en se basant sur la premiere face
  E_Int dim = 3;
  dim = ptr[2];
  dim = min(dim, E_Int(3));

  // Commence par calculer alpha
  E_Int nfaces = ptr[0];
  FldArrayF alphat(nfaces);
  E_Float* alphap = alphat.begin();
  E_Int* cFE1 = cFE.begin(1);
  E_Int* cFE2 = cFE.begin(2);
#pragma omp parallel
  {
  E_Int e1, e2;
  vector<E_Int> indices;
  vector<E_Float*> pts1; vector<E_Float*> pts2;
  E_Float alpha;
  E_Int ind, nvert;

#pragma omp for
  for (E_Int i = 0; i < nfaces; i++)
  {
    e1 = cFE1[i]-1; e2 = cFE2[i]-1;
    if (e1 == -1) alpha = -1000.;
    else if (e2 == -1) alpha = -1000.;
    else
    {
      K_CONNECT::getVertexIndices(cn->begin(), posFaces.begin(), pos[e1], indices);
      nvert = indices.size();
      pts1.reserve(nvert);
      for (E_Int k = 0; k < nvert; k++)
      {
        E_Float* pt = new E_Float[3];
        ind = indices[k]-1;
        pt[0] = x[ind]; pt[1] = y[ind]; pt[2] = z[ind];
        pts1.push_back(pt);
      }
      K_CONNECT::getVertexIndices(cn->begin(), posFaces.begin(), pos[e2], indices);
      nvert = indices.size();
      pts2.reserve(nvert);
      for (E_Int k = 0; k < nvert; k++)
      {
        E_Float* pt = new E_Float[3];
        ind = indices[k]-1;
        pt[0] = x[ind]; pt[1] = y[ind]; pt[2] = z[ind];
        pts2.push_back(pt);
      }

      if (dim == 2)
      {
        alpha = K_COMPGEOM::getAlphaAngleBetweenPolygons(pts1, pts2);
        //printf("mbv=%d, alpha=%f %f\n", mbv, alpha, alphaRef);
      }
      else if (dim == 1)
      {
        if (pts1.size() == 2 && pts2.size() == 2)
        {
          alpha = K_COMPGEOM::getAlphaAngleBetweenBars(
                    pts1[0], pts1[1],
                    pts2[0], pts2[1], dirVect);
        }
        else alpha = 180.;
      }
      else alpha = 180.;

      nvert = pts1.size();
      for (E_Int k = 0; k < nvert; k++) delete [] pts1[k];
      pts1.clear();
      nvert = pts2.size();
      for (E_Int k = 0; k < nvert; k++) delete [] pts2[k];
      pts2.clear();
    }
    alphap[i] = alpha;
  }
  }

  E_Int e1, e2, nf;
  E_Float alpha;

  // split
  mbv = 0;
  while (nev < ne)
  {
    // Recherche le premier elt pas encore visite
    for (p = 0; (isVisited[p] != 0); p++);

    // C'est un nouveau composant (morceau de liste)
    FldArrayI* c = new FldArrayI(se+1);
    E_Int* pc = c->begin(); // current pointer

    mustBeVisited[mbv] = p;
    mbv++; nev++;
    isVisited[p] = 1;
    necurr = 0;

    while (mbv > 0)
    {
      mbv--;
      elt = mustBeVisited[mbv];

      // copie index de l'element
      ptrElts = &ptr[pos[elt]];
      lt = ptrElts[0];
      (*pc) = index[elt]; pc += 1;
      necurr++;

      for (iv = 0; iv < lt; iv++)
      {
        nf = ptrElts[iv+1]-1;
        e1 = cFE1[nf]-1; e2 = cFE2[nf]-1;
        if (e1 == elt) ie = e2;
        else ie = e1;

        if (ie != -1 && isVisited[ie] == 0)
        {
          alpha = alphap[nf];
          if (alpha == -1000. || K_FUNC::E_abs(alpha-180.) <= alphaRef)
          {
            mustBeVisited[mbv] = ie;
            mbv++; nev++;
            isVisited[ie] = 1;
          }
        }
      }
    }

    c->reAlloc(necurr);
    components.push_back(c);
  }

  free(isVisited);
  free(mustBeVisited);

  // Formation des arrays de sortie + cleanConnectivity
  PyObject* tpl;
  PyObject* l = PyList_New(0);
  E_Int size = components.size();

  for (i = 0; i < size; i++)
  {
    FldArrayI* c = components[i];
    E_Int nd = c->getSize();
    E_Int* pc = c->begin();
    tpl = K_NUMPY::buildNumpyArray(nd, 1, 1, 1);
    E_Int* pos = K_NUMPY::getNumpyPtrI(tpl);
    for (E_Int k = 0; k < nd; k++) pos[k] = pc[k];
    delete components[i];
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }
  RELEASESHAREDU(array, f, cn);
  RELEASESHAREDN(arrayI, indexI);
  return l;
}

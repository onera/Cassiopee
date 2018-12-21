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

# include "generator.h"
# include "Fld/ArrayAccessor.h"
# include "Connect/merge.h"
# include "Search/KdTree.h"

using namespace K_FUNC;
using namespace K_FLD;
using namespace K_CONST;
using namespace std;

// ============================================================================
/* Close a list of meshes at epsilon */
// ============================================================================
PyObject* K_GENERATOR::closeAllMeshes(PyObject* self, PyObject* args)
{
  PyObject* arrays;
  E_Float eps;
  
  if (!PYPARSETUPLEF(args, "Od", "Of", &arrays, &eps))   
    return NULL;

  // Extract infos from arrays
  vector<E_Int> res;
  vector<char*> structVarString; vector<char*> unstrVarString;
  vector<FldArrayF*> structF; vector<FldArrayF*> unstrF;
  vector<E_Int> nit; vector<E_Int> njt; vector<E_Int> nkt;
  vector<FldArrayI*> cnt; vector<char*> eltType;
  vector<PyObject*> objs, obju;
  E_Boolean skipNoCoord = true;
  E_Boolean skipStructured = false;
  E_Boolean skipUnstructured = false;
  E_Boolean skipDiffVars = false;

  E_Int isOk = K_ARRAY::getFromArrays(
    arrays, res, structVarString, unstrVarString,
    structF, unstrF, nit, njt, nkt, cnt, eltType, objs, obju, 
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured);

  E_Int size;
  if (isOk == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "close: invalid list of arrays.");
    size = structF.size();
    for (E_Int v = 0; v < size; v++) delete structF[v];
    size = unstrF.size();
    for (E_Int v = 0; v < size; v++) delete unstrF[v];
    size = cnt.size();
    for (E_Int v = 0; v < size; v++) delete cnt[v];
    return NULL;
  }
  
  // Position des coordonnees
  vector<E_Int> posxs; vector<E_Int> posys; vector<E_Int> poszs;
  E_Int posx, posy, posz;
  size = structF.size();
  for (E_Int v = 0; v < size; v++)
  {
    posx = K_ARRAY::isCoordinateXPresent(structVarString[v]);
    posy = K_ARRAY::isCoordinateYPresent(structVarString[v]);
    posz = K_ARRAY::isCoordinateZPresent(structVarString[v]);
    posxs.push_back(posx+1);
    posys.push_back(posy+1);
    poszs.push_back(posz+1);
  }
  vector<E_Int> posxu; vector<E_Int> posyu; vector<E_Int> poszu;
  size = unstrF.size();
  for (E_Int v = 0; v < size; v++)
  {
    posx = K_ARRAY::isCoordinateXPresent(unstrVarString[v]);
    posy = K_ARRAY::isCoordinateYPresent(unstrVarString[v]);
    posz = K_ARRAY::isCoordinateZPresent(unstrVarString[v]);
    posxu.push_back(posx+1);
    posyu.push_back(posy+1);
    poszu.push_back(posz+1);
  }

  // Structure
  closeAllStructuredMeshes(structF, nit, njt, nkt, posxs, posys, poszs, eps);
  
  // Non structure
  size = unstrF.size();
  for (E_Int v = 0; v < size; v++)
  {
    FldArrayF& f0 = *unstrF[v]; 
    closeUnstructuredMesh(posxu[v], posyu[v], poszu[v], 
                          eps, eltType[v], f0, *cnt[v]);
  }

  PyObject* l = PyList_New(0);
  PyObject* tpl;
  size = structF.size();
  for (E_Int v = 0; v < size; v++)
  {
    FldArrayF& f0 = *structF[v];
    tpl = K_ARRAY::buildArray(f0, structVarString[v], nit[v], njt[v], nkt[v]);
    delete &f0;
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }
  size = unstrF.size();
  for (E_Int v = 0; v < size; v++)
  {
    FldArrayF& f0 = *unstrF[v];
    tpl = K_ARRAY::buildArray(f0, unstrVarString[v], *cnt[v], -1, eltType[v],
                              false);
    delete &f0; delete cnt[v];
    PyList_Append(l, tpl);
    Py_DECREF(tpl);
  }
  return l;
}

// ============================================================================
/* Close a mesh at epsilon */
// ============================================================================
PyObject* K_GENERATOR::closeMesh(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Float eps;
  
  if (!PYPARSETUPLEF(args, "Od", "Of", &array, &eps)) return NULL;

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = 
    K_ARRAY::getFromArray(array, varString, f, im, jm, km, cn, eltType);
  if (res == 1)
  {
    E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
    E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
    E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
    if (posx == -1 || posy == -1 || posz == -1)
    {
      delete f;
      PyErr_SetString(PyExc_TypeError,
                      "close: can't find coordinates in array.");
      return NULL;
    }
    posx++; posy++; posz++;
    closeStructuredMesh(f->begin(posx), f->begin(posy), f->begin(posz), im, jm, km, eps);

    PyObject* tpl = K_ARRAY::buildArray(*f, varString, im, jm, km); 
    delete f;
    return tpl;
  }
  else if (res == 2)
  { 
     if (K_STRING::cmp(eltType, "NODE*") == 0 ||
         K_STRING::cmp(eltType, "BAR*") == 0 ||
         K_STRING::cmp(eltType, "TRI*") == 0 ||
         K_STRING::cmp(eltType, "QUAD*") == 0 ||
         K_STRING::cmp(eltType, "TETRA*") ==0 ||
         K_STRING::cmp(eltType, "PYRA*") == 0 ||
         K_STRING::cmp(eltType, "PENTA*") == 0 ||
         K_STRING::cmp(eltType, "HEXA*") == 0 ||
         K_STRING::cmp(eltType, "MIX*") == 0)
     {
       delete f; delete cn;
       PyErr_SetString(PyExc_TypeError,
                       "close: array must be defined at vertices.");
       return NULL;
     }
     E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
     E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
     E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
     if (posx == -1 || posy == -1 || posz == -1)
     {
       delete f;
       PyErr_SetString(PyExc_TypeError,
                       "close: can't find coordinates in array.");
       return NULL;
     }
     posx++; posy++; posz++;
     
     closeUnstructuredMesh(posx, posy, posz, eps, eltType, *f, *cn);
     PyObject* tpl = 
       K_ARRAY::buildArray(*f, varString, *cn, -1, eltType, false);
     delete f; delete cn;
     return tpl;
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "close: unrecognised type of array.");
    return NULL;
  }
}

//=============================================================================
/* Close an unstructured mesh. */
//=============================================================================
void K_GENERATOR::closeUnstructuredMesh(E_Int posx, E_Int posy, E_Int posz,
                                        E_Float eps, char* eltType,
                                        FldArrayF& f, FldArrayI& cn)
{
  //E_Int nv = cn.getNfld();
  //if (nv == 2)
  //  closeBARMesh(posx, posy, posz, f, cn);
  
  K_CONNECT::cleanConnectivity(posx, posy, posz, eps, eltType, f, cn);
}

//=============================================================================
// Traitement special pour les BAR-array uniquement dans le plan (x,y):
// On supprime en plus les definitions redondantes des BARs (quand deux BARs
// se superposent)
//=============================================================================
void K_GENERATOR::closeBARMesh(E_Int posx, E_Int posy, E_Int posz,
                               FldArrayF& f, FldArrayI& cn)
{
  // Verifie que la BAR est dans le plan x,y
  E_Float xmin, ymin, zmin, xmax, ymax, zmax;
  K_COMPGEOM::boundingBox(posx, posy, posz,
                          f,
                          xmin, ymin, zmin,
                          xmax, ymax, zmax);
  if (fEqualZero(zmax-zmin) == false) return;

  E_Int i = 0;
  while (i < cn.getSize())
  {
    E_Int ret = closeBARMeshElt(posx, posy, posz,
                                f, cn, i);
    if (ret == 0) i++;
  }
}
//=============================================================================
E_Int K_GENERATOR::closeBARMeshElt(E_Int posx, E_Int posy, E_Int posz,
                                   FldArrayF& f, FldArrayI& cn, E_Int i)
{
  // Segmentation si necessaire de l'elt i
  E_Float* xt = f.begin(posx);
  E_Float* yt = f.begin(posy);
  E_Float* zt = f.begin(posz);
  E_Float A1[3];
  E_Float A2[3];
  E_Float B1[3];
  E_Float B2[3];
  E_Float pi0[3];
  E_Float pi1[3];
  E_Float d1, d2, d3, d4, d5;
  E_Int indA1, indA2, indB1, indB2;
  E_Int res = 0;

  E_Int ne = cn.getSize();
  E_Int nv = f.getSize();
  E_Int nfield = f.getNfld();
  E_Int nemax = 2*ne;
  E_Int nvmax = 2*nv;
  FldArrayF f2 = f;
  f2.reAllocMat(nvmax, nfield);
  FldArrayI cn2(nemax, 2);
  for (E_Int j = 0; j < i; j++)
  {
    cn2(j,1) = cn(j,1); cn2(j,2) = cn(j,2);
  }

  E_Int elt = i; // element courant dans cn2 

  // Loop
  indA1 = cn(i,1)-1;
  indA2 = cn(i,2)-1;
  A1[0] = xt[indA1]; A1[1] = yt[indA1]; A1[2] = zt[indA1];
  A2[0] = xt[indA2]; A2[1] = yt[indA2]; A2[2] = zt[indA2];

  cn2(elt, 1) = indA1+1;

  for (E_Int j = 0; j < ne; j++)
  {
    if (j != i)
    {
      indB1 = cn(j,1)-1;
      indB2 = cn(j,2)-1;
      B1[0] = xt[indB1]; B1[1] = yt[indB1]; B1[2] = zt[indB1];
      B2[0] = xt[indB2]; B2[1] = yt[indB2]; B2[2] = zt[indB2];
      
      res = K_COMPGEOM::intersect2Segments(A1, A2,
                                           B1, B2,
                                           pi0, pi1);
      
      if (res == 1)
      {
        d1 = (pi0[0]-A1[0])*(pi0[0]-A1[0])
          +(pi0[1]-A1[1])*(pi0[1]-A1[1])
          +(pi0[2]-A1[2])*(pi0[2]-A1[2]);
        d2 = (pi0[0]-A2[0])*(pi0[0]-A2[0])
          +(pi0[1]-A2[1])*(pi0[1]-A2[1])
          +(pi0[2]-A2[2])*(pi0[2]-A2[2]);
        
        if (fEqualZero(d1) == false && fEqualZero(d2) == false)
        {
          E_Float a = sqrt(d1) / (sqrt(d1) + sqrt(d2));
          for (E_Int nf = 1; nf <= nfield; nf++)
            f2(nv,nf) = a*f(indA1,nf) + (1.-a)*f(indA2,nf);
          f2(nv,posx) = pi0[0]; 
          f2(nv,posy) = pi0[1]; 
          f2(nv,posz) = pi0[2]; nv++;
          cn2(elt, 2) = nv; elt++;
          cn2(elt, 1) = nv;
          goto end;
        }
        else res = 0;
      }
      else if (res == 2)
      {
        d1 = (A1[0]-A2[0])*(A1[0]-A2[0])
          +(A1[1]-A2[1])*(A1[1]-A2[1])
          +(A1[2]-A2[2])*(A1[2]-A2[2]);
        d2 = (A1[0]-pi0[0])*(A1[0]-pi0[0])
          +(A1[1]-pi0[1])*(A1[1]-pi0[1])
          +(A1[2]-pi0[2])*(A1[2]-pi0[2]);
        d3 = (A1[0]-pi1[0])*(A1[0]-pi1[0])
          +(A1[1]-pi1[1])*(A1[1]-pi1[1])
          +(A1[2]-pi1[2])*(A1[2]-pi1[2]);
        d4 = (A2[0]-pi0[0])*(A2[0]-pi0[0])
          +(A2[1]-pi0[1])*(A2[1]-pi0[1])
          +(A2[2]-pi0[2])*(A2[2]-pi0[2]);
        d5 = (A2[0]-pi1[0])*(A2[0]-pi1[0])
          +(A2[1]-pi1[1])*(A2[1]-pi1[1])
          +(A2[2]-pi1[2])*(A2[2]-pi1[2]);
        
        if (fEqualZero(d2) == false && 
            fEqualZero(d3) == false &&
            fEqualZero(d4) == false &&
            fEqualZero(d5) == false &&
            d2 < d1 - 1.e-12 && d3 < d1 - 1.e-12) // insere les deux
        {
          if (d2 < d3)
          {
            E_Float a = sqrt(d2) / sqrt(d1);
            for (E_Int nf = 1; nf <= nfield; nf++)
              f2(nv,nf) = a*f(indA1,nf) + (1.-a)*f(indA2,nf);
            f2(nv,posx) = pi0[0]; f2(nv,posy) = pi0[1]; f2(nv,posz) = pi0[2]; nv++;
            cn2(elt, 2) = nv; elt++;
            a = sqrt(d3) / sqrt(d1);
            for (E_Int nf = 1; nf <= nfield; nf++)
              f2(nv,nf) = a*f(indA1,nf) + (1.-a)*f(indA2,nf);
            f2(nv,posx) = pi1[0]; f2(nv,posy) = pi1[1]; f2(nv,posz) = pi1[2]; nv++;
            cn2(elt, 1) = nv;
            goto end;
          }
          else
          {
            E_Float a = sqrt(d3) / sqrt(d1);
            for (E_Int nf = 1; nf <= nfield; nf++)
              f2(nv,nf) = a*f(indA1,nf) + (1.-a)*f(indA2,nf);
            f2(nv,posx) = pi1[0]; f2(nv,posy) = pi1[1]; f2(nv,posz) = pi1[2]; nv++;
            cn2(elt, 2) = nv; elt++;
            a = sqrt(d2) / sqrt(d1);
            for (E_Int nf = 1; nf <= nfield; nf++)
              f2(nv,nf) = a*f(indA1,nf) + (1.-a)*f(indA2,nf);
            f2(nv,posx) = pi0[0]; f2(nv,posy) = pi0[1]; f2(nv,posz) = pi0[2]; nv++;
            cn2(elt, 1) = nv;
            goto end;
          }
        }
        else if (fEqualZero(d2) == false && 
                 fEqualZero(d4) == false && 
                 d2 < d1 - 1.e-12) // insere pi0
        {
          E_Float a = sqrt(d2) / sqrt(d1);
          for (E_Int nf = 1; nf <= nfield; nf++)
            f2(nv,nf) = a*f(indA1,nf) + (1.-a)*f(indA2,nf);
          f2(nv,posx) = pi0[0]; f2(nv,posy) = pi0[1]; f2(nv,posz) = pi0[2]; nv++;
          cn2(elt, 2) = nv; elt++;
          cn2(elt, 1) = nv;
          goto end;
        }
        else if (fEqualZero(d3) == false && 
                 fEqualZero(d5) == false && 
                 d3 < d1 - 1.e-12) // insere pi1
        {
          E_Float a = sqrt(d3) / sqrt(d1);
          for (E_Int nf = 1; nf <= nfield; nf++)
            f2(nv,nf) = a*f(indA1,nf) + (1.-a)*f(indA2,nf);
          f2(nv,posx) = pi1[0]; f2(nv,posy) = pi1[1]; f2(nv,posz) = pi1[2]; nv++;
          cn2(elt, 2) = nv; elt++;
          cn2(elt, 1) = nv;
          goto end;
        }
        else
          res = 0;
      }
    }
  }
   
  end:
  cn2(elt, 2) = indA2+1; elt++;
  
  for (E_Int j = i+1; j < ne; j++)
  {
    cn2(elt,1) = cn(j,1); cn2(elt,2) = cn(j,2); elt++;
  }
  f2.reAllocMat(nv, f.getNfld());
  cn2.reAllocMat(elt, 2);
  
  f = f2;
  cn = cn2;

  return res;
}
//=============================================================================
/* closeAllStructuredMeshes */
//=============================================================================
void K_GENERATOR::closeAllStructuredMeshes(
  vector<FldArrayF*>& structF,
  vector<E_Int>& nit, vector<E_Int>& njt, vector<E_Int>& nkt,
  vector<E_Int>& posxt, vector<E_Int>& posyt, vector<E_Int>& poszt,
  E_Float eps)
{
  // Create a kdtree with all the external windows of all zones
  E_Int nzones = structF.size();

  FldArrayF dmin(nzones);
  FldArrayI indmin(nzones,2); // stockage indtmin, indmin par zone

  E_Int sizemax = 0;
  E_Int sizemaxloc = 0;
  for (E_Int v1 = 0; v1 < nzones; v1++)
  {    
    if (nkt[v1] > 1) sizemaxloc=2*njt[v1]*nkt[v1]+2*nit[v1]*nkt[v1]+2*nit[v1]*njt[v1];
    else sizemaxloc = 2*njt[v1]+2*nit[v1];
    sizemax += sizemaxloc;
  }
  // Creation du kdtree et des tableaux d'indirection
  FldArrayF ftemp(sizemax,3);
  E_Float* xp = ftemp.begin(1);
  E_Float* yp = ftemp.begin(2);
  E_Float* zp = ftemp.begin(3);
  FldArrayI indirB(sizemax); FldArrayI indirI(sizemax);
  
  E_Int ind, shifti, shiftj, shiftk;
  E_Int indt = 0; // index in kdtree
  for (E_Int noz = 0; noz < nzones; noz++)
  {    
    E_Int ni = nit[noz]; E_Int nj = njt[noz]; E_Int nk = nkt[noz];
    E_Int ninj = ni*nj;
    E_Float* xt = structF[noz]->begin(posxt[noz]);
    E_Float* yt = structF[noz]->begin(posyt[noz]);
    E_Float* zt = structF[noz]->begin(poszt[noz]);

    // i = 1
    shifti = 0;
    for (E_Int k = 0; k < nk; k++)
      for (E_Int j = 0; j < nj; j++)
      {
        ind = shifti + j*ni + k*ninj;
        xp[indt] = xt[ind]; yp[indt] = yt[ind]; zp[indt] = zt[ind];
        indirB[indt] = noz; indirI[indt] = ind;
        indt++;
      }
    // i = ni
    shifti = ni-1;
    for (E_Int k = 0; k < nk; k++)
      for (E_Int j = 0; j < nj; j++)
      {
        ind = shifti + j*ni + k*ninj;
        xp[indt] = xt[ind]; yp[indt] = yt[ind]; zp[indt] = zt[ind];
        indirB[indt] = noz; indirI[indt] = ind;
        indt++;
      }    
    // j = 1
    shiftj = 0;
    for (E_Int k = 0; k < nk; k++)
      for (E_Int i = 0; i < ni; i++)
      {
        ind = i + shiftj*ni + k*ninj;
        xp[indt] = xt[ind]; yp[indt] = yt[ind]; zp[indt] = zt[ind];
        indirB[indt] = noz; indirI[indt] = ind;
        indt++;
      }
    // j = nj
    shiftj = nj-1;
    for (E_Int k = 0; k < nk; k++)
      for (E_Int i = 0; i < ni; i++)
      {
        ind = i + shiftj*ni + k*ninj;
        xp[indt] = xt[ind]; yp[indt] = yt[ind]; zp[indt] = zt[ind];
        indirB[indt] = noz; indirI[indt] = ind;
        indt++;
      }
    if ( nk > 1)
    {
      // k = 1
      shiftk = 0;
      for (E_Int j = 0; j < nj; j++)
        for (E_Int i = 0; i < ni; i++)
        {
          ind = i + j*ni + shiftk*ninj;
          xp[indt] = xt[ind]; yp[indt] = yt[ind]; zp[indt] = zt[ind];
          indirB[indt] = noz; indirI[indt] = ind;
          indt++;
        }
      // k = nk
      shiftk = nk-1;
      for (E_Int j = 0; j < nj; j++)
        for (E_Int i = 0; i < ni; i++)
        {
          ind = i + j*ni + shiftk*ninj;
          xp[indt] = xt[ind]; yp[indt] = yt[ind]; zp[indt] = zt[ind];
          indirB[indt] = noz; indirI[indt] = ind;
          indt++;
        }    
    }
  }
  ArrayAccessor<FldArrayF> coordAcc(ftemp, 1, 2, 3);
  K_SEARCH::KdTree<FldArrayF> globalKdt(coordAcc);
  FldArrayI tag(sizemax); tag.setAllValuesAtNull();
  E_Float eps10 = 10*eps;
  E_Float bbmin[3], bbmax[3];
  vector<E_Int> candidates;
  E_Float dx, dy, dz;
  for (E_Int indt = 0; indt < sizemax; indt++)
  {
    E_Int noz = indirB[indt]; E_Int ind = indirI[indt];
    
    E_Float* xt = structF[noz]->begin(posxt[noz]);
    E_Float* yt = structF[noz]->begin(posyt[noz]);
    E_Float* zt = structF[noz]->begin(poszt[noz]);
    bbmin[0] = xt[ind]-eps10; bbmin[1] = yt[ind]-eps10; bbmin[2] = zt[ind]-eps10;
    bbmax[0] = xt[ind]+eps10; bbmax[1] = yt[ind]+eps10; bbmax[2] = zt[ind]+eps10;
    candidates.clear();
    globalKdt.getInBox(bbmin, bbmax, candidates);
    //dmin = K_CONST::E_MAX_FLOAT;
    //indmin = -1; indtmin = -1;
    dmin.setAllValuesAt(K_CONST::E_MAX_FLOAT);
    indmin.setAllValuesAt(-1);
    for (unsigned int cr = 0; cr < candidates.size(); cr++)
    {
      E_Int indkr = candidates[cr];
      E_Int nozr = indirB[indkr]; E_Int indr = indirI[indkr];
      if (nozr != noz)
      {
        E_Float* xr = structF[nozr]->begin(posxt[nozr]);
        E_Float* yr = structF[nozr]->begin(posyt[nozr]);
        E_Float* zr = structF[nozr]->begin(poszt[nozr]);
        dx = xr[indr]-xt[ind];
        dy = yr[indr]-yt[ind];
        dz = zr[indr]-zt[ind];
        if (K_FUNC::fEqualZero(dx,eps) == true &&
            K_FUNC::fEqualZero(dy,eps) == true &&
            K_FUNC::fEqualZero(dz,eps) == true && 
            tag[indt] == 0 && tag[indkr] == 0)
        {
          if (dx*dx+dy*dy+dz*dz < dmin[nozr])
          {
            dmin[nozr] = dx*dx+dy*dy+dz*dz;
            indmin(nozr,1) = indkr; indmin(nozr,2) = indr;
          }
        }
      }
    }
    
    E_Int accu = 1;
    E_Float xmean = xt[ind];
    E_Float ymean = yt[ind];
    E_Float zmean = zt[ind];
    for (E_Int nz = 0; nz < nzones; nz++)
    {
      E_Int im = indmin(nz,2);
      if (im != -1)
      {
        accu++;
        E_Float* xm = structF[nz]->begin(posxt[nz]);
        E_Float* ym = structF[nz]->begin(posyt[nz]);
        E_Float* zm = structF[nz]->begin(poszt[nz]);
        xmean += xm[im];
        ymean += ym[im];
        zmean += zm[im];
      }
    }

    if (accu > 1)
    {
      E_Float w = 1./accu;
      xmean = xmean*w; ymean = ymean*w; zmean = zmean*w;
      xt[ind] = xmean; yt[ind] = ymean; zt[ind] = zmean; tag[indt] = 1;
      for (E_Int nz = 0; nz < nzones; nz++)
      {
        E_Int im = indmin(nz,2);
        if (im != -1)
        {
          E_Float* xm = structF[nz]->begin(posxt[nz]);
          E_Float* ym = structF[nz]->begin(posyt[nz]);
          E_Float* zm = structF[nz]->begin(poszt[nz]);
          xm[im] = xmean; ym[im] = ymean; zm[im] = zmean;
          tag[indmin(nz,1)] = 1;
        }
      }
    }

  } // tous les pts
}
//=============================================================================
/* closeAllStructuredMeshes */
//=============================================================================
void K_GENERATOR::closeStructuredMesh(E_Float* xt, E_Float* yt, E_Float* zt,
                                      E_Int ni, E_Int nj, E_Int nk, E_Float eps)
{
  E_Int sizemax = 0;
  if (nk > 1) sizemax = 2*nj*nk+2*ni*nk+2*ni*nj;
  else sizemax = 2*ni+2*nj;
  // Creation du kdtree et des tableaux d indirection
  FldArrayI indirI(sizemax);
  FldArrayF ftemp(sizemax,3);
  E_Float* xp = ftemp.begin(1);
  E_Float* yp = ftemp.begin(2);
  E_Float* zp = ftemp.begin(3);
  E_Int ind;
  E_Int indt = 0; // index in kdtree
  
  E_Int ninj = ni*nj;
  // i = 1
  E_Int shifti = 0;
  for (E_Int k = 0; k < nk; k++)
    for (E_Int j = 0; j < nj; j++)
    {
      ind = shifti + j*ni + k*ninj;
      xp[indt] = xt[ind]; yp[indt] = yt[ind]; zp[indt] = zt[ind];
      indirI[indt] = ind;
      indt++;
    }
    // i = ni
  shifti = ni-1;
  for (E_Int k = 0; k < nk; k++)
    for (E_Int j = 0; j < nj; j++)
    {
      ind = shifti + j*ni + k*ninj;
      xp[indt] = xt[ind]; yp[indt] = yt[ind]; zp[indt] = zt[ind];
      indirI[indt] = ind;
      indt++;
    }    
  // j = 1
  E_Int shiftj = 0;
  for (E_Int k = 0; k < nk; k++)
    for (E_Int i = 0; i < ni; i++)
    {
      ind = i + shiftj*ni + k*ninj;
      xp[indt] = xt[ind]; yp[indt] = yt[ind]; zp[indt] = zt[ind];
      indirI[indt] = ind;
      indt++;
    }
  // j = nj
  shiftj = nj-1;
  for (E_Int k = 0; k < nk; k++)
    for (E_Int i = 0; i < ni; i++)
    {
      ind = i + shiftj*ni + k*ninj;
      xp[indt] = xt[ind]; yp[indt] = yt[ind]; zp[indt] = zt[ind];
      indirI[indt] = ind;
      indt++;
    }
  
  if ( nk > 1 )
  {
    // k = 1
    E_Int shiftk = 0;
    for (E_Int j = 0; j < nj; j++)
      for (E_Int i = 0; i < ni; i++)
      {
        ind = i + j*ni + shiftk*ninj;
        xp[indt] = xt[ind]; yp[indt] = yt[ind]; zp[indt] = zt[ind];
        indirI[indt] = ind;
        indt++;
      }
    // k = nk
    shiftk = nk-1;
    for (E_Int j = 0; j < nj; j++)
      for (E_Int i = 0; i < ni; i++)
      {
        ind = i + j*ni + shiftk*ninj;
        xp[indt] = xt[ind]; yp[indt] = yt[ind]; zp[indt] = zt[ind];
        indirI[indt] = ind;
        indt++;
      }    
  }
  
  ArrayAccessor<FldArrayF> coordAcc(ftemp, 1, 2, 3);
  vector<E_Int> newID(sizemax);
  E_Int nbMerges = merge(coordAcc, eps, newID);

  if ( nbMerges > 0 )
  {
    for (E_Int indt = 0; indt < sizemax; indt++)
    {
      E_Int indtopp = newID[indt];
      if (indtopp != indt )
      {
        E_Int ind = indirI[indt]; E_Int indopp = indirI[indtopp];          
        xt[ind] = xt[indopp];
        yt[ind] = yt[indopp];
        zt[ind] = zt[indopp];
      }  
    }
  }
}

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

# include "generator.h"
# include "Nuga/include/ArrayAccessor.h"
# include "Nuga/include/merge.h"
# include "Nuga/include/KdTree.h"

using namespace K_FUNC;
using namespace K_FLD;
using namespace K_CONST;
using namespace std;

// ============================================================================
/* Close a mesh at epsilon */
// ============================================================================
PyObject* K_GENERATOR::closeMesh(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Float eps;
  E_Bool rmOverlappingPts = true;
  E_Bool rmOrphanPts = true;
  E_Bool rmDuplicatedFaces = true;
  E_Bool rmDuplicatedElts = true;
  E_Bool rmDegeneratedFaces = true;
  E_Bool rmDegeneratedElts = true;
  E_Bool exportIndirPts = false;
  
  if (!PYPARSETUPLE_(args, O_ R_ BBB_ BBBB_,
                    &array, &eps, &rmOverlappingPts, &rmOrphanPts,
                    &rmDuplicatedFaces, &rmDuplicatedElts,
                    &rmDegeneratedFaces, &rmDegeneratedElts,
                    &exportIndirPts)) return NULL;

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;

  E_Int res = K_ARRAY::getFromArray3(array, varString, f, im, jm, km, cn, eltType);
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDB(res, array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "close: can't find coordinates in array.");
    return NULL;
  }
  posx++; posy++; posz++;

  if (res == 1)
  {
    closeStructuredMesh(f->begin(posx), f->begin(posy), f->begin(posz), im, jm, km, eps);
    PyObject* tpl = K_ARRAY::buildArray3(*f, varString, im, jm, km); 
    RELEASESHAREDS(array, f);
    return tpl;
  }
  else if (res == 2)
  { 
    if (strchr(eltType, '*') != NULL)
    {
      RELEASESHAREDU(array, f, cn);
      PyErr_SetString(PyExc_TypeError,
                      "close: array must be defined at vertices.");
      return NULL;
    }

    PyObject* tpl = K_CONNECT::V_cleanConnectivity(
        varString, *f, *cn, eltType, eps,
        rmOverlappingPts, rmOrphanPts, rmDuplicatedFaces, rmDuplicatedElts,
        rmDegeneratedFaces, rmDegeneratedElts, exportIndirPts);

    RELEASESHAREDU(array, f, cn);
    if (tpl == NULL) return array;
    else return tpl;
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "close: unrecognised type of array.");
    return NULL;
  }
}

// ============================================================================
PyObject* K_GENERATOR::closeMeshLegacy(PyObject* self, PyObject* args)
{
  PyObject* array;
  E_Float eps;
  E_Int removeDegen = 0;
  
  if (!PYPARSETUPLE_(args, O_ R_ I_, &array, &eps, &removeDegen)) return NULL;

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;

  E_Int res = K_ARRAY::getFromArray3(array, varString, f, im, jm, km, cn, eltType);
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1)
  {
    RELEASESHAREDS(array, f);
    PyErr_SetString(PyExc_TypeError,
                    "close: can't find coordinates in array.");
    return NULL;
  }
  posx++; posy++; posz++;

  if (res == 1)
  {
    closeStructuredMesh(f->begin(posx), f->begin(posy), f->begin(posz), im, jm, km, eps);
    PyObject* tpl = K_ARRAY::buildArray3(*f, varString, im, jm, km); 
    RELEASESHAREDS(array, f);
    return tpl;
  }
  else if (res == 2)
  { 
    if (strchr(eltType, '*') != NULL)
    {
      RELEASESHAREDU(array, f, cn);
      PyErr_SetString(PyExc_TypeError,
                      "close: array must be defined at vertices.");
      return NULL;
    }

    closeUnstructuredMesh(posx, posy, posz, eps, eltType, *f, *cn, removeDegen);
    PyObject* tpl = K_ARRAY::buildArray3(*f, varString, *cn, eltType);
    RELEASESHAREDU(array, f, cn);
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
/* Close a structured mesh. */
//=============================================================================
void K_GENERATOR::closeStructuredMesh(E_Float* xt, E_Float* yt, E_Float* zt,
                                      E_Int ni, E_Int nj, E_Int nk, E_Float eps)
{
  E_Int ninj = ni*nj;
  E_Int nink = ni*nk;
  E_Int njnk = nj*nk;
  E_Int sizemax = 0;
  if (nk > 1) sizemax = 2*(njnk + nink + ninj);
  else sizemax = 2*(ni + nj);

  // Creation du kdtree et des tableaux d indirection
  FldArrayI indirI(sizemax);
  FldArrayF ftemp(sizemax,3);
  E_Float* xp = ftemp.begin(1);
  E_Float* yp = ftemp.begin(2);
  E_Float* zp = ftemp.begin(3);
  
  #pragma omp parallel
  {
    E_Int ind, indt; // indt: index in kdtree
    // i = 1
    E_Int shifti = 0;
    #pragma omp for
    for (E_Int k = 0; k < nk; k++)
      for (E_Int j = 0; j < nj; j++)
      {
        ind = shifti + j*ni + k*ninj;
        indt = k*nj + j;
        xp[indt] = xt[ind]; yp[indt] = yt[ind]; zp[indt] = zt[ind];
        indirI[indt] = ind;
      }
    // i = ni
    shifti = ni-1;
    #pragma omp for
    for (E_Int k = 0; k < nk; k++)
      for (E_Int j = 0; j < nj; j++)
      {
        ind = shifti + j*ni + k*ninj;
        indt = njnk + k*nj + j;
        xp[indt] = xt[ind]; yp[indt] = yt[ind]; zp[indt] = zt[ind];
        indirI[indt] = ind;
      }    
    // j = 1
    E_Int shiftj = 0;
    #pragma omp for
    for (E_Int k = 0; k < nk; k++)
      for (E_Int i = 0; i < ni; i++)
      {
        ind = i + shiftj*ni + k*ninj;
        indt = 2*njnk + k*ni + i;
        xp[indt] = xt[ind]; yp[indt] = yt[ind]; zp[indt] = zt[ind];
        indirI[indt] = ind;
      }
    // j = nj
    shiftj = nj-1;
    #pragma omp for
    for (E_Int k = 0; k < nk; k++)
      for (E_Int i = 0; i < ni; i++)
      {
        ind = i + shiftj*ni + k*ninj;
        indt = 2*njnk + nink + k*ni + i;
        xp[indt] = xt[ind]; yp[indt] = yt[ind]; zp[indt] = zt[ind];
        indirI[indt] = ind;
      }
    
    if ( nk > 1 )
    {
      // k = 1
      E_Int shiftk = 0;
      #pragma omp for
      for (E_Int j = 0; j < nj; j++)
        for (E_Int i = 0; i < ni; i++)
        {
          ind = i + j*ni + shiftk*ninj;
          indt = 2*(njnk + nink) + j*ni + i;
          xp[indt] = xt[ind]; yp[indt] = yt[ind]; zp[indt] = zt[ind];
          indirI[indt] = ind;
        }
      // k = nk
      shiftk = nk-1;
      #pragma omp for
      for (E_Int j = 0; j < nj; j++)
        for (E_Int i = 0; i < ni; i++)
        {
          ind = i + j*ni + shiftk*ninj;
          indt = 2*(njnk + nink) + ninj + j*ni + i;
          xp[indt] = xt[ind]; yp[indt] = yt[ind]; zp[indt] = zt[ind];
          indirI[indt] = ind;
        }
    }
  }
  
  ArrayAccessor<FldArrayF> coordAcc(ftemp, 1, 2, 3);
  vector<E_Int> newID(sizemax);
  E_Int nbMerges = merge(coordAcc, eps, newID);

  E_Int ind, indtopp, indopp;
  if (nbMerges > 0)
  {
    for (E_Int i = 0; i < sizemax; i++)
    {
      indtopp = newID[i];
      if (indtopp != i)
      {
        ind = indirI[i]; indopp = indirI[indtopp];          
        xt[ind] = xt[indopp];
        yt[ind] = yt[indopp];
        zt[ind] = zt[indopp];
      }  
    }
  }
}

//=============================================================================
/* Close an unstructured mesh. */
//=============================================================================
void K_GENERATOR::closeUnstructuredMesh(E_Int posx, E_Int posy, E_Int posz,
                                        E_Float eps, char* eltType,
                                        FldArrayF& f, FldArrayI& cn, E_Int removeDegen)
{
  //E_Int nv = cn.getNfld();
  //if (nv == 2)
  //  closeBARMesh(posx, posy, posz, f, cn);
  
  K_CONNECT::cleanConnectivity(posx, posy, posz, eps, eltType, f, cn, removeDegen);
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

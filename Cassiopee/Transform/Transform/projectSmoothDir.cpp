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

# include "transform.h"
# include "Nuga/include/BbTree.h"

using namespace K_FLD;
using namespace std;

// ============================================================================
/* Project along direction dir with smoothing */
// ============================================================================
PyObject* K_TRANSFORM::projectSmoothDir(PyObject* self, PyObject* args)
{
  PyObject* arrays;
  PyObject* array2;
  E_Float nx, ny, nz;
  E_Int oriented;

  if (!PYPARSETUPLE_(args, OO_ TRRR_ I_,
                    &arrays, &array2, &nx, &ny, &nz, &oriented))
  {
      return NULL;
  }
  // Extract infos from arrays
  vector<E_Int> resl;
  vector<char*> structVarString; vector<char*> unstrVarString;
  vector<FldArrayF*> structF; vector<FldArrayF*> unstrF;
  vector<E_Int> nit; vector<E_Int> njt; vector<E_Int> nkt;
  vector<FldArrayI*> cnt; vector<char*> eltType;
  vector<PyObject*> objst, objut;
  E_Boolean skipNoCoord = true;
  E_Boolean skipStructured = false;
  E_Boolean skipUnstructured = true;
  E_Boolean skipDiffVars = true;
  E_Int isOk = K_ARRAY::getFromArrays(
    arrays, resl, structVarString, unstrVarString,
    structF, unstrF, nit, njt, nkt, cnt, eltType, objst, objut, 
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);
  E_Int nu = objut.size(); E_Int ns = objst.size();
  if (isOk == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "projectSmoothDir: invalid list of arrays.");
    for (E_Int nos = 0; nos < ns; nos++)
      RELEASESHAREDS(objst[nos], structF[nos]);
    for (E_Int nos = 0; nos < nu; nos++)
      RELEASESHAREDU(objut[nos], unstrF[nos], cnt[nos]);
    return NULL;
  }
  E_Int posx1, posy1, posz1;
  vector<E_Int> posxs; vector<E_Int> posys; vector<E_Int> poszs;
  vector<E_Int> posxu; vector<E_Int> posyu; vector<E_Int> poszu;
  for (E_Int nos = 0; nos < ns; nos++)
  {
    posx1 = K_ARRAY::isCoordinateXPresent(structVarString[nos]); posx1++;
    posy1 = K_ARRAY::isCoordinateYPresent(structVarString[nos]); posy1++;
    posz1 = K_ARRAY::isCoordinateZPresent(structVarString[nos]); posz1++;
    posxs.push_back(posx1); posys.push_back(posy1); poszs.push_back(posz1); 
  }
  for (E_Int nou = 0; nou < nu; nou++)
  {
    posx1 = K_ARRAY::isCoordinateXPresent(unstrVarString[nou]); posx1++;
    posy1 = K_ARRAY::isCoordinateYPresent(unstrVarString[nou]); posy1++;
    posz1 = K_ARRAY::isCoordinateZPresent(unstrVarString[nou]); posz1++;
    posxu.push_back(posx1); posyu.push_back(posy1); poszu.push_back(posz1); 
  }

  // projection surfaces
  E_Int im2, jm2, km2;
  FldArrayF* f2; FldArrayI* cn2;
  char* varString2; char* eltType2;
  E_Int res2 = K_ARRAY::getFromArray(array2, varString2, 
                                     f2, im2, jm2, km2, cn2, eltType2, true); 
  if (res2 != 2)
  {
    for (E_Int nos = 0; nos < ns; nos++)
      RELEASESHAREDS(objst[nos], structF[nos]);
    for (E_Int nos = 0; nos < nu; nos++)
      RELEASESHAREDU(objut[nos], unstrF[nos], cnt[nos]);
    RELEASESHAREDB(res2, array2, f2, cn2);
    PyErr_SetString(PyExc_TypeError,
                    "projectSmoothDir: array2 must be unstructured.");
    return NULL;
  }
  if (strcmp(eltType2, "TRI") != 0)
  {
    for (E_Int nos = 0; nos < ns; nos++)
      RELEASESHAREDS(objst[nos], structF[nos]);
    for (E_Int nos = 0; nos < nu; nos++)
      RELEASESHAREDU(objut[nos], unstrF[nos], cnt[nos]);
    RELEASESHAREDB(res2, array2, f2, cn2);

    PyErr_SetString(PyExc_TypeError,
                    "projectSmoothDir: array2 must be a TRI array.");
    return NULL;
  }

  E_Int posx2 = K_ARRAY::isCoordinateXPresent(varString2);
  E_Int posy2 = K_ARRAY::isCoordinateYPresent(varString2);
  E_Int posz2 = K_ARRAY::isCoordinateZPresent(varString2);
   
  if (posx2 == -1 || posy2 == -1 || posz2 == -1)
  {
    for (E_Int nos = 0; nos < ns; nos++)
      RELEASESHAREDS(objst[nos], structF[nos]);
    for (E_Int nos = 0; nos < nu; nos++)
      RELEASESHAREDU(objut[nos], unstrF[nos], cnt[nos]);
    RELEASESHAREDB(res2, array2, f2, cn2);

    PyErr_SetString(PyExc_TypeError,
                    "projectSmoothDir: can't find coordinates in array2.");
    return NULL;
  }
  posx2++; posy2++; posz2++;

  // Projete
  vector<FldArrayF*> structFields;
  E_Int precond = 0;
  E_Int im1, jm1, km1;
  for (E_Int nos = 0; nos < ns; nos++)
  {
    posx1 = posxs[nos]; posy1 = posys[nos]; posz1 = poszs[nos];
    FldArrayF* f = new FldArrayF(*structF[nos]);
    im1 = nit[nos]; jm1 = njt[nos]; km1 = nkt[nos];
    // Reorder eventuel de array1 pour avoir un tableau ni x nj
    if (km1 != 1 && im1 == 1)
      K_CONNECT::reorderStructField(im1, jm1, km1, *f, 3, 2, -1);
    
    else if (km1 != 1 && jm1 == 1)
      K_CONNECT::reorderStructField(im1, jm1, km1, *f, 1, 3, -2);

    if (precond == 0)
      projectSmoothDirWithoutPrecond(
        nx, ny, nz, im1, jm1, km1, cn2->getSize(), *cn2,
        f2->begin(posx2), f2->begin(posy2), f2->begin(posz2), 
        f->begin(posx1), f->begin(posy1), f->begin(posz1), oriented);
    else 
      projectSmoothDirWithPrecond(
        nx, ny, nz, im1, jm1, km1, cn2->getSize(), *cn2,
        f2->begin(posx2), f2->begin(posy2), f2->begin(posz2), 
        f->begin(posx1), f->begin(posy1), f->begin(posz1), oriented);
    structFields.push_back(f);
  }

  RELEASESHAREDU(array2, f2, cn2);
  for (E_Int nos = 0; nos < ns; nos++)
    RELEASESHAREDS(objst[nos], structF[nos]);
  for (E_Int nou = 0; nou < nu; nou++)
    RELEASESHAREDU(objut[nou], unstrF[nou], cnt[nou]);

  // Build arrays
  PyObject* l = PyList_New(0);
  PyObject* tpl;    
  for (E_Int nos = 0; nos < ns; nos++)
  {
    tpl = K_ARRAY::buildArray(*structFields[nos], structVarString[nos],
                              nit[nos], njt[nos], nkt[nos]);
    delete structFields[nos];
    PyList_Append(l, tpl); Py_DECREF(tpl);
  }

  return l;
}

//=============================================================================
/* projection avec preconditionnement par bboxtree */
//=============================================================================
void K_TRANSFORM::projectSmoothDirWithPrecond(
  E_Float nx, E_Float ny, E_Float nz,
  E_Int im1, E_Int jm1, E_Int km1,
  E_Int nelts2, FldArrayI& cn2,
  E_Float* fx2, E_Float* fy2, E_Float* fz2,
  E_Float* fx, E_Float* fy, E_Float* fz, E_Int oriented)
{
  E_Int npts = im1*jm1*km1;
  FldArrayIS tag = FldArrayIS(npts); tag.setAllValuesAtNull();
  short* tagp = tag.begin();
  E_Float tol = 1.e-6;
  typedef K_SEARCH::BoundingBox<3>  BBox3DType;

  E_Int* cn2p1 = cn2.begin(1); E_Int* cn2p2 = cn2.begin(2); E_Int* cn2p3 = cn2.begin(3);

  // Creation de la bboxtree
  vector<BBox3DType*> boxes(nelts2);// liste des bbox de ts les elements de a2
  K_FLD::FldArrayF bbox(nelts2,6);// xmin, ymin, zmin, xmax, ymax, zmax
  K_COMPGEOM::boundingBoxOfUnstrCells(cn2, fx2, fy2, fz2, bbox);
  E_Float minB[3];  E_Float maxB[3];
  E_Float* xminp = bbox.begin(1); E_Float* xmaxp = bbox.begin(4);
  E_Float* yminp = bbox.begin(2); E_Float* ymaxp = bbox.begin(5);
  E_Float* zminp = bbox.begin(3); E_Float* zmaxp = bbox.begin(6);
  for (E_Int et = 0; et < nelts2; et++)
  {
    minB[0] = xminp[et]; minB[1] = yminp[et]; minB[2] = zminp[et];
    maxB[0] = xmaxp[et]; maxB[1] = ymaxp[et]; maxB[2] = zmaxp[et]; 
    boxes[et] = new BBox3DType(minB, maxB);
  }
  // Build the box tree.
  K_SEARCH::BbTree3D bbtree(boxes);
    
  // Algorithme de projection
  E_Float p[3]; E_Float pr1[3]; E_Float pr2[3]; E_Float pi[3];
  E_Float p0[3]; E_Float p1[3]; E_Float p2[3];
  E_Float dist; E_Float distc; 
  E_Int ret; E_Int ind1, ind2, ind3;
  E_Float dx, dy, dz, normp, ps;
  vector<E_Int> indicesBB;
  E_Int nbboxes, et;

  for (E_Int ind = 0; ind < npts; ind++)
  {
    p[0] = fx[ind]; p[1] = fy[ind]; p[2] = fz[ind];
    pr1[0] = p[0]; pr1[1] = p[1]; pr1[2] = p[2];
    pr2[0] = p[0] + nx; pr2[1] = p[1] + ny; pr2[2] = p[2] + nz;
    bbtree.getIntersectingBoxes(pr1, pr2, indicesBB, tol);
    
    distc = 1e6;
    nbboxes = indicesBB.size();
    for (E_Int noe = 0; noe < nbboxes; noe++)
    {
      et = indicesBB[noe];
      ind1 = cn2p1[et]-1; ind2 = cn2p2[et]-1; ind3 = cn2p3[et]-1;
      p0[0] = fx2[ind1]; p0[1] = fy2[ind1]; p0[2] = fz2[ind1];
      p1[0] = fx2[ind2]; p1[1] = fy2[ind2]; p1[2] = fz2[ind2];
      p2[0] = fx2[ind3]; p2[1] = fy2[ind3]; p2[2] = fz2[ind3];        
      ret = K_COMPGEOM::intersectRayTriangle(p0, p1, p2,
                                             pr1, pr2,
                                             pi);

      if (ret == 1)
      {
        dx = pi[0]-p[0]; dy = pi[1]-p[1]; dz = pi[2]-p[2];
        dist = dx*dx + dy*dy + dz*dz;
        if ( oriented != 0 )
        {
          //distp = sqrt(dist);//PP'
          normp = sqrt(nx*nx+ny*ny+nz*nz);//normale
          ps = nx*dx+ny*dy+nz*dz/(dist*normp);
          if ( ps > 0. ) 
          {
            if (dist < distc) 
            {fx[ind] = pi[0]; fy[ind] = pi[1]; fz[ind] = pi[2]; distc = dist; tagp[ind] = 1;}
          }
        }
        else 
        {
          if (dist < distc) 
          {fx[ind] = pi[0]; fy[ind] = pi[1]; fz[ind] = pi[2]; distc = dist; tagp[ind] = 1;}
        }
      }
    }
    indicesBB.clear();
  }

  for (E_Int et = 0; et < nelts2; et++) delete boxes[et];
  boxes.clear();

  smoothUnprojectedPts(im1, jm1, km1, npts, fx, fy, fz, tagp);

}
//=============================================================================
/* projection sans preconditionnement par bboxtree*/
//=============================================================================
void K_TRANSFORM::projectSmoothDirWithoutPrecond(
  E_Float nx, E_Float ny, E_Float nz,
  E_Int im1, E_Int jm1, E_Int km1,
  E_Int nelts2, FldArrayI& cn2,
  E_Float* fx2, E_Float* fy2, E_Float* fz2,
  E_Float* fx, E_Float* fy, E_Float* fz, E_Int oriented)
{
  E_Int npts = im1*jm1*km1;
  FldArrayIS tag = FldArrayIS(npts); tag.setAllValuesAtNull();
  short* tagp = tag.begin();

  E_Int* cn2p1 = cn2.begin(1); E_Int* cn2p2 = cn2.begin(2); E_Int* cn2p3 = cn2.begin(3);
  // Algorithme de projection
  E_Float p[3]; E_Float pr1[3]; E_Float pr2[3]; E_Float pi[3];
  E_Float p0[3]; E_Float p1[3]; E_Float p2[3];
  E_Float dist; E_Float distc; 
  E_Int ret; E_Int ind1, ind2, ind3;
  E_Float dx, dy, dz, normp, ps;
  for (E_Int ind = 0; ind < npts; ind++)
  {
    p[0] = fx[ind]; p[1] = fy[ind]; p[2] = fz[ind];
    pr1[0] = p[0]; pr1[1] = p[1]; pr1[2] = p[2];
    pr2[0] = p[0] + nx; pr2[1] = p[1] + ny; pr2[2] = p[2] + nz;
    distc = 1e6;
    for (E_Int e = 0; e < nelts2; e++)
    {
      ind1 = cn2p1[e]-1;
      ind2 = cn2p2[e]-1;
      ind3 = cn2p3[e]-1;
      p0[0] = fx2[ind1]; p0[1] = fy2[ind1]; p0[2] = fz2[ind1];
      p1[0] = fx2[ind2]; p1[1] = fy2[ind2]; p1[2] = fz2[ind2];
      p2[0] = fx2[ind3]; p2[1] = fy2[ind3]; p2[2] = fz2[ind3];        
     
      ret = K_COMPGEOM::intersectRayTriangle(p0, p1, p2,
                                             pr1, pr2,
                                             pi);
      if (ret == 1)
      {
        dx = pi[0]-p[0]; dy = pi[1]-p[1]; dz = pi[2]-p[2];
        dist = dx*dx + dy*dy + dz*dz;
        if ( oriented != 0 )
        {
          //distp = sqrt(dist);//PP'
          normp = sqrt(nx*nx+ny*ny+nz*nz);//normale
          ps = nx*dx+ny*dy+nz*dz/(dist*normp);
          if ( ps > 0. ) 
          {
            if (dist < distc) 
              {fx[ind] = pi[0]; fy[ind] = pi[1]; fz[ind] = pi[2]; distc = dist;tagp[ind] = 1;}
          }
        }
        else 
        {
          if (dist < distc) 
          {fx[ind] = pi[0]; fy[ind] = pi[1]; fz[ind] = pi[2]; distc = dist;tagp[ind] = 1;}
        } 
      }
    }
  }
  smoothUnprojectedPts(im1, jm1, km1, npts, fx, fy, fz, tagp);
}
//=============================================================================
// lisse les pts non projetes
//===========================================================================
void K_TRANSFORM::smoothUnprojectedPts(
  E_Int im1, E_Int jm1, E_Int km1, E_Int npts,
  E_Float* fx, E_Float* fy, E_Float* fz, short* tagp)
{
  // Lisse les points non projetes
  FldArrayF coord(npts, 3);
  E_Float* cx = coord.begin(1);
  E_Float* cy = coord.begin(2);
  E_Float* cz = coord.begin(3);
  E_Float eps = 0.5;
  E_Int niter = 10000;
  E_Int ni = im1; E_Int nj = jm1;
  E_Int nij = ni * nj;
  E_Int i, j, k;
  E_Float deltax, deltay, deltaz;
  E_Int indNP, indNM, indPN, indMN;
  E_Float diff = 1.e6;
  E_Int nit = 0;
  E_Float eps4 = 0.25 * eps; E_Float eps1 = 1.-eps;
  E_Int ni1 = ni-1; E_Int nj1 = nj-1;
  E_Float diffx, diffy, diffz;

  while (nit < niter && diff > 1.e-10)
  {
    diff = 0;
    for (E_Int ind = 0; ind < npts; ind++)
    {
      if (tagp[ind] == 0) // point non projete
      {
        k = ind / nij;
        j = (ind - k*nij) / ni;
        i = ind - j*ni - k*nij; 
        indPN = K_FUNC::E_min(i+1, ni1) + j*ni + k*nij;
        indMN = K_FUNC::E_max(i-1, 0) + j*ni + k*nij;
        indNP = i + K_FUNC::E_min(j+1, nj1)*ni + k*nij;
        indNM = i + K_FUNC::E_max(j-1, 0)*ni + k*nij;
        deltax = fx[indPN] + fx[indMN] + fx[indNP] + fx[indNM];
        deltay = fy[indPN] + fy[indMN] + fy[indNP] + fy[indNM];
        deltaz = fz[indPN] + fz[indMN] + fz[indNP] + fz[indNM];
        cx[ind] = fx[ind]*eps1 + eps4 * deltax;
        cy[ind] = fy[ind]*eps1 + eps4 * deltay;
        cz[ind] = fz[ind]*eps1 + eps4 * deltaz;
        diffx = fx[ind] - cx[ind]; 
        diffy = fy[ind] - cy[ind]; 
        diffz = fz[ind] - cz[ind]; 

        diff = K_FUNC::E_max(diff, diffx*diffx + diffy*diffy + diffz*diffz);
      }
    }
    for (E_Int ind = 0; ind < npts; ind++)
    {
      if (tagp[ind] == 0) // point non projete
      {
        fx[ind] = cx[ind]; fy[ind] = cy[ind]; fz[ind] = cz[ind];
      }
    }
    nit++;
  }
  //printf("niter=%d, diff=%f\n", niter, diff);
}

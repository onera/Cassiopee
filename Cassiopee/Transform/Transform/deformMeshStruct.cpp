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

# include <stdio.h>
# include "transform.h"
using namespace std;
using namespace K_FLD;
using namespace K_SEARCH;

#define AMORTI \
  x0 = xw[indw]; y0 = yw[indw]; z0 = zw[indw];                          \
  dist = (xt[ind]-x0)*(xt[ind]-x0)+(yt[ind]-y0)*(yt[ind]-y0)+(zt[ind]-z0)*(zt[ind]-z0); \
  beta0 = dist/(eps+delta);                                             \
  theta = exp(-beta*beta0);                                             \
  dxt[ind] = theta*dx0; dyt[ind] = theta*dy0; dzt[ind] = theta*dz0;

#define BLOCK \
  pt[0] = xt[ind]; pt[1] = yt[ind]; pt[2] = zt[ind];                    \
  indw = kdt.getClosest(pt);                                            \
  dx0 = xw[indw]-pt[0]; dy0 = yw[indw]-pt[1]; dz0 = zw[indw]-pt[2];     \
  dist2 = dx0*dx0+dy0*dy0+dz0*dz0;                                      \
  distmaxloc = K_FUNC::E_max(distmaxloc, dist2);                        \
  ptrIndices[noi]=indw; noi++;                                          \

// ============================================================================
/* Computes the deformation vector for the 4/6 borders of a mesh             */
// ============================================================================
PyObject* K_TRANSFORM::deformMeshStruct(PyObject* self, 
                                        PyObject* args)
{
  PyObject* array; PyObject* deltas;
  E_Float beta;
  if (!PYPARSETUPLE_(args, OO_ R_, &array, &deltas, &beta))
    return NULL;

  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res = K_ARRAY::getFromArray3(array, varString, f, im, jm, km, cn, eltType);
  
  if (res == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "deformMeshStruct: invalid array.");
    return NULL;
  }
  if ( res == 2) 
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "deformMeshStruct: array must be structured.");
    return NULL;
  }

  E_Int posx1 = K_ARRAY::isCoordinateXPresent(varString); posx1++;
  E_Int posy1 = K_ARRAY::isCoordinateYPresent(varString); posy1++;
  E_Int posz1 = K_ARRAY::isCoordinateZPresent(varString); posz1++;

  // Extract infos from surface arrays for which delta is known
  vector<E_Int> resls;
  vector<char*> structVarStrings; vector<char*> unstrVarStrings;
  vector<FldArrayF*> structFs; vector<FldArrayF*> unstrFs;
  vector<E_Int> nits; vector<E_Int> njts; vector<E_Int> nkts;
  vector<FldArrayI*> cnts;
  vector<char*> eltTypes;
  vector<PyObject*> objsts, objuts;
  E_Bool skipNoCoord = true;
  E_Bool skipStructured = true;
  E_Bool skipUnstructured = false; 
  E_Bool skipDiffVars = true;
  E_Int isOk = K_ARRAY::getFromArrays(
    deltas, resls, structVarStrings, unstrVarStrings,
    structFs, unstrFs, nits, njts, nkts, cnts, eltTypes, objsts, objuts, 
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);

  E_Int nwalls = objuts.size();
  if (isOk == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "deformMeshStruct: invalid list of arrays.");
    RELEASESHAREDS(array, f);
    for (E_Int nos = 0; nos < nwalls; nos++)
      RELEASESHAREDU(objuts[nos], unstrFs[nos], cnts[nos]);
    return NULL;
  }
  if (nwalls == 0) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "deformMeshStruct: no valid surface provided.");
    RELEASESHAREDS(array, f);
    for (E_Int nos = 0; nos < nwalls; nos++)
      RELEASESHAREDU(objuts[nos], unstrFs[nos], cnts[nos]);
    return NULL;
  }
  vector<E_Int> posxw; vector<E_Int> posyw; vector<E_Int> poszw;
  vector<E_Int> posdxw; vector<E_Int> posdyw; vector<E_Int> posdzw;
  E_Int posdx, posdy, posdz;
  E_Int found = 1;
  for (E_Int nou = 0; nou < nwalls; nou++)
  {
    posx1 = K_ARRAY::isCoordinateXPresent(unstrVarStrings[nou]); posx1++;
    posy1 = K_ARRAY::isCoordinateYPresent(unstrVarStrings[nou]); posy1++;
    posz1 = K_ARRAY::isCoordinateZPresent(unstrVarStrings[nou]); posz1++;
    posdx = K_ARRAY::isNamePresent("dx",unstrVarStrings[nou]); posdx++;
    if (posdx < 1) { found = 0; break; }
    posdy = K_ARRAY::isNamePresent("dy",unstrVarStrings[nou]); posdy++;
    if (posdx < 1) { found = 0; break; }
    posdz = K_ARRAY::isNamePresent("dz",unstrVarStrings[nou]); posdz++;
    if (posdx < 1) { found = 0; break; }
    posxw.push_back(posx1); posyw.push_back(posy1); poszw.push_back(posz1); 
    posdxw.push_back(posdx); posdyw.push_back(posdy); posdzw.push_back(posdz); 
  }
  if (found == 0) 
  {
    PyErr_SetString(PyExc_TypeError,
                    "deformMeshStruct: surface arrays must all contain dx, dy, dz variables.");
    RELEASESHAREDS(array, f);
    for (E_Int nos = 0; nos < nwalls; nos++)
      RELEASESHAREDU(objuts[nos], unstrFs[nos], cnts[nos]);
    return NULL;
  }

  // Creation du kdtree
  E_Int nptsmax = 0;
  for (E_Int v = 0; v < nwalls; v++) nptsmax += unstrFs[v]->getSize();
  FldArrayF* surfaces = new FldArrayF(nptsmax, 6);
  E_Float* xw =  surfaces->begin(1);
  E_Float* yw =  surfaces->begin(2);
  E_Float* zw =  surfaces->begin(3);
  E_Float* dxw = surfaces->begin(4);
  E_Float* dyw = surfaces->begin(5);
  E_Float* dzw = surfaces->begin(6);
  E_Int c = 0;
  for (E_Int v = 0; v < nwalls; v++)
  {
    FldArrayF* fieldv = unstrFs[v];
    E_Int posxv = posxw[v]; E_Int posyv = posyw[v]; E_Int poszv = poszw[v];
    E_Int posdxv = posdxw[v]; E_Int posdyv = posdyw[v]; E_Int posdzv = posdzw[v];
    E_Float* xw0 = fieldv->begin(posxv); E_Float* dxw0 = fieldv->begin(posdxv);
    E_Float* yw0 = fieldv->begin(posyv); E_Float* dyw0 = fieldv->begin(posdyv);
    E_Float* zw0 = fieldv->begin(poszv); E_Float* dzw0 = fieldv->begin(posdzv);
    E_Int nptsw = fieldv->getSize();
    for (E_Int i = 0; i < nptsw; i++)
    {
      xw[c] = xw0[i]; yw[c] = yw0[i]; zw[c] = zw0[i]; 
      dxw[c] = dxw0[i]; dyw[c] = dyw0[i]; dzw[c] = dzw0[i];
      c++;
    }
  } //fin kdtree
  
  ArrayAccessor<FldArrayF> coordAcc(*surfaces, 1,2,3);
  KdTree<FldArrayF> kdt(coordAcc);

  // Build array
  E_Int nfld = 3;
  E_Int npts = f->getSize();
  PyObject* tpl = K_ARRAY::buildArray(nfld, "dx,dy,dz", im, jm, km);
  E_Float* fp = K_ARRAY::getFieldPtr(tpl);
  FldArrayF fout(npts, nfld, fp, true); fout.setAllValuesAtNull();
  E_Float* dxt = fout.begin(1);
  E_Float* dyt = fout.begin(2);
  E_Float* dzt = fout.begin(3);

  E_Float pt[3]; 

  E_Float* xt = f->begin(posx1);
  E_Float* yt = f->begin(posy1);
  E_Float* zt = f->begin(posz1);
  E_Float beta0, dist, delta, theta, x0, y0, z0, dx0, dy0, dz0, dist2, distmaxloc;
  E_Float distmax = -K_CONST::E_MAX_FLOAT;
  E_Int shift0, shift1;
  E_Int ind, indw, dir, ok, noi;
  E_Int* ptrIndices;
  // On utilise une exponentielle en puissance de 2
  // On attend 70% d'amortissement a une distance de beta*deformation locale 
  //E_Float beta = 7;
  //E_Int factor = 8;
  E_Float eps = 1.e-10;
  beta = 1./(beta*beta); //beta = 1./pow(beta, factor);

  if (km == 1)//2D
  {
    vector<FldArrayI> vectOfIndices(4);
    FldArrayI indIMIN(jm); vectOfIndices[0] = indIMIN; 
    FldArrayI indIMAX(jm); vectOfIndices[1] = indIMAX; 
    FldArrayI indJMIN(im); vectOfIndices[2] = indJMIN; 
    FldArrayI indJMAX(im); vectOfIndices[3] = indJMAX; 
    //j=1
    distmaxloc = -K_CONST::E_MAX_FLOAT;
    ptrIndices=vectOfIndices[2].begin();
    noi=0;
    for (E_Int i = 0; i < im; i++)
    { 
      ind = i;
      BLOCK;
    }
    distmax = distmaxloc; dir = 3;
    
    //i=1
    distmaxloc = -K_CONST::E_MAX_FLOAT;
    ptrIndices=vectOfIndices[0].begin();
    ok = 1; noi=0;
    for (E_Int j = 0; j < jm; j++)
    {      
      ind = j*im;
      BLOCK;
      if (distmaxloc > distmax)  {ok = -1; break;}    
    }
    if (distmaxloc < distmax && ok == 1) {distmax = distmaxloc; dir=1;}
    
    //i=imax
    distmaxloc = -K_CONST::E_MAX_FLOAT;
    ok = 1; noi=0;
    ptrIndices=vectOfIndices[1].begin();
    for (E_Int j = 0; j < jm; j++)
    {      
      ind = im-1+j*im;
      BLOCK;
      if (distmaxloc > distmax)  {ok = -1; break;}    
    }
    if (distmaxloc < distmax  && ok == 1) {distmax = distmaxloc; dir=2;}
    
    //j=jmax
    distmaxloc = -K_CONST::E_MAX_FLOAT;
    ok = 1; noi=0;
    ptrIndices=vectOfIndices[3].begin();
    for (E_Int i = 0; i < im; i++)
    { 
      ind = i+(jm-1)*im;
      BLOCK;
      if (distmaxloc > distmax)  {ok = -1; break;}    
    }
    if (distmaxloc < distmax && ok == 1) {distmax = distmaxloc; dir=4;}
    
    //cleaning
    for (E_Int no = 1; no < 5; no++)
    {
      if (no != dir)  vectOfIndices[no-1].malloc(0);
    }
    ptrIndices = vectOfIndices[dir-1].begin();
    
    // deformation vector computation for each node 
    // detection of the nearest surface point
    if (dir == 1 || dir == 2)
    {
      noi = 0;
      if (dir == 1) { shift0 = 0; shift1 = im-1;}
      else { shift0 = im-1; shift1 = 0;}

      for(E_Int j=0; j<jm; j++)
      {
        ind = shift0+j*im;
        indw = ptrIndices[noi];
        dx0 = dxw[indw]; dy0 = dyw[indw]; dz0 = dzw[indw];//delta vector of nearest surface pt
        delta = dx0*dx0+dy0*dy0+dz0*dz0;//norm of the delta vector

        //Pt on border i=const
        AMORTI;
        
        //Pt on opposite border
        ind=shift1+j*im;
        AMORTI;
        
        // Lateral borders - pas tres joli...
        if ( j == 0 || j == jm-1)
        {
          for (E_Int i = 1; i < im-1; i++)
          {
            ind = i + j*im;
            AMORTI;
          }
        }        
        noi++;
      }
    }//dir=1 or 2
    else if (dir == 3 || dir == 4)
    {
      noi = 0;
      if (dir == 3) {shift0 = 0; shift1 = (jm-1)*im;}
      else { shift0 = (jm-1)*im; shift1 = 0;}

      for(E_Int i=0; i < im; i++)
      {
        ind = i+shift0;
        indw = ptrIndices[noi];
        dx0 = dxw[indw]; dy0 = dyw[indw]; dz0 = dzw[indw];//delta vector of nearest surface pt
        delta = dx0*dx0+dy0*dy0+dz0*dz0;//norm of the delta vector

        //Pt on border j=const
        AMORTI;
        
        //Pt on opposite border
        ind=i+shift1;
        AMORTI;

        // Lateral borders - pas tres joli...
        if (i == 0 || i == im-1)
        {
          for (E_Int j = 1; j < jm-1; j++)
          {
            ind = i + j*im;
            AMORTI;
          }
        }        
        noi++;
      }
    }
  }
  else
  {
    vector<FldArrayI> vectOfIndices(6);
    FldArrayI indIMIN(jm*km); vectOfIndices[0] = indIMIN; 
    FldArrayI indIMAX(jm*km); vectOfIndices[1] = indIMAX; 
    FldArrayI indJMIN(im*km); vectOfIndices[2] = indJMIN; 
    FldArrayI indJMAX(im*km); vectOfIndices[3] = indJMAX; 
    FldArrayI indKMIN(im*jm); vectOfIndices[4] = indKMIN; 
    FldArrayI indKMAX(im*jm); vectOfIndices[5] = indKMAX; 
    printf(" coucou 1\n");

    //k=1
    distmaxloc = -K_CONST::E_MAX_FLOAT;
    ptrIndices=vectOfIndices[4].begin();
    noi=0;
    for(E_Int j=0; j < jm; j++)
    { 
      for (E_Int i = 0; i < im; i++)
      { 
        ind = i+j*im;
        BLOCK;
      }
    } 
    distmax = distmaxloc; dir = 5;

    //k=km
    distmaxloc = -K_CONST::E_MAX_FLOAT;
    ptrIndices=vectOfIndices[5].begin();
    ok = 1; noi=0;
    E_Int imjm = im*jm;
    for (E_Int j = 0; j < jm; j++)
    { 
      for (E_Int i = 0; i < im; i++)
      {      
        ind = i+j*im+(km-1)*imjm;
        BLOCK;
        if (distmaxloc > distmax) {ok = -1; break;}
      }
    } 
    if (distmaxloc < distmax && ok == 1) {distmax = distmaxloc; dir=6;}

    //i=1
    distmaxloc = -K_CONST::E_MAX_FLOAT;
    ptrIndices=vectOfIndices[0].begin();
    ok = 1; noi=0;
    for (E_Int k = 0; k < km; k++)
    { 
      for (E_Int j = 0; j < jm; j++)
      {      
        ind = j*im+k*imjm;
        BLOCK;
        if (distmaxloc > distmax) {ok = -1; break;}
      }
    } 
    if (distmaxloc < distmax && ok == 1) {distmax = distmaxloc; dir=1;}

    //i=im
    distmaxloc = -K_CONST::E_MAX_FLOAT;
    ptrIndices=vectOfIndices[1].begin();
    ok = 1; noi=0;
    for (E_Int k = 0; k < km; k++)
    { 
      for (E_Int j = 0; j < jm; j++)
      {      
        ind = im-1+j*im+k*imjm;
        BLOCK;
        if (distmaxloc > distmax) {ok = -1; break;}
      }
    } 
    if (distmaxloc < distmax && ok == 1) {distmax = distmaxloc; dir=2;}

    //j=1
    distmaxloc = -K_CONST::E_MAX_FLOAT;
    ok = 1; noi=0;
    ptrIndices=vectOfIndices[2].begin();
    for (E_Int k = 0; k < km; k++)
      for (E_Int i = 0; i < im; i++)
      { 
        ind = i+ k*imjm;
        BLOCK;
        if (distmaxloc > distmax)  {ok = -1; break;}
      }
    if (distmaxloc < distmax && ok == 1) {distmax = distmaxloc; dir=3;}

    //j=jmax
    distmaxloc = -K_CONST::E_MAX_FLOAT;
    ok = 1; noi=0;
    ptrIndices=vectOfIndices[3].begin();
    for (E_Int k = 0; k < km; k++)
      for (E_Int i = 0; i < im; i++)
      { 
        ind = i+(jm-1)*im +k*imjm;
        BLOCK;
        if (distmaxloc > distmax)  {ok = -1; break;}
      }
    if (distmaxloc < distmax && ok == 1) {distmax = distmaxloc; dir=4;}
    //cleaning
    for (E_Int no = 1; no < 7; no++)
    {
      if (no != dir)  vectOfIndices[no-1].malloc(0);
    }
    ptrIndices = vectOfIndices[dir-1].begin();
    // deformation vector computation for each node 
    // detection of the nearest surface point
    if (dir == 1 || dir == 2)
    {
      noi = 0;
      if (dir == 1) { shift0 = 0; shift1 = im-1;}
      else { shift0 = im-1; shift1 = 0;}

      for (E_Int k=0; k<km; k++)
        for (E_Int j=0; j<jm; j++)
        {
          ind = shift0+j*im+k*imjm;
          indw = ptrIndices[noi];
          dx0 = dxw[indw]; dy0 = dyw[indw]; dz0 = dzw[indw];//delta vector of nearest surface pt
          delta = dx0*dx0+dy0*dy0+dz0*dz0;//norm of the delta vector
          
          //Pt on border i=const
          AMORTI;
          
          //Pt on opposite border
          ind=shift1+j*im+k*imjm;
          AMORTI;
        
          // Lateral borders - pas tres joli...
          if ( j == 0 || j == jm-1)
          {            
            for (E_Int i = 1; i < im-1; i++)
            {
              ind = i + j*im+k*imjm;
              AMORTI;
            }
          }        
          if (k == 0 || k == km-1)
          {            
            for (E_Int i = 1; i < im-1; i++)
            {
              ind = i + j*im +k*imjm;
              AMORTI;
            }
          } 
          noi++;
        }
    }//dir=1 or 2
    else if (dir == 3 || dir == 4)
    {
      noi = 0;
      if ( dir == 3) {shift0 = 0; shift1 = (jm-1)*im;}
      else { shift0 = (jm-1)*im; shift1 = 0;}
      for(E_Int k=0; k<km; k++)
        for(E_Int i=0; i<im; i++)
        {
          ind = i+shift0+k*imjm;
          indw = ptrIndices[noi];
          dx0 = dxw[indw]; dy0 = dyw[indw]; dz0 = dzw[indw];//delta vector of nearest surface pt
          delta = dx0*dx0+dy0*dy0+dz0*dz0;//norm of the delta vector

          //Pt on border j=const
          AMORTI;
        
          //Pt on opposite border
          ind=i+shift1+k*imjm;
          AMORTI;
          
          // Lateral borders - pas tres joli...
          if ( i == 0 || i == im-1)
          {
            for (E_Int j = 1; j < jm-1; j++)
            {
              ind = i + j*im + k*imjm;
              AMORTI;
            }
          }    
          if (k == 0 || k == km-1)
          {
            for (E_Int j = 1; j < jm-1; j++)
            {
              ind = i + j*im + k*imjm;
             AMORTI;
            }
          }            
          noi++;
      }
    }//dir=3 or 4
    else if (dir == 5 || dir == 6)
    {
      noi = 0;
      if (dir == 5) {shift0 = 0; shift1 = (km-1)*imjm;}
      else { shift0 = (km-1)*imjm; shift1 = 0;}
      for(E_Int j=0; j<jm; j++)
        for(E_Int i=0; i<im; i++)
        {
          ind = i+j*im+shift0;
          indw = ptrIndices[noi];
          dx0 = dxw[indw]; dy0 = dyw[indw]; dz0 = dzw[indw];//delta vector of nearest surface pt
          delta = dx0*dx0+dy0*dy0+dz0*dz0;//norm of the delta vector

          //Pt on border k=const
          AMORTI;
        
          //Pt on opposite border
          ind=i+shift1+j*im;
          AMORTI;
          
          // Lateral borders - pas tres joli...
          if (i == 0 || i == im-1)
          {
            for (E_Int k = 1; k < km-1; k++)
            {
              ind = i + j*im + k*imjm;
              AMORTI;
            }
          }    
          if (j == 0 || j == jm-1)
          {
            for (E_Int k = 1; k < km-1; k++)
            {
              ind = i + j*im + k*imjm;
              AMORTI;
            }
          }            
        noi++;
      }
    }//dir=3 or 4

  } //3D
  
  delete surfaces;
  // Release memory
  RELEASESHAREDS(array, f);
  for (E_Int nou = 0; nou < nwalls; nou++)
    RELEASESHAREDU(objuts[nou], unstrFs[nou], cnts[nou]);
  return tpl;
}

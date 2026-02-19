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
//  nullifyVectorAtBCFace: set Q.n = 0 on the BCFace (Q may be u)

# include "converter.h"
using namespace K_FLD;
using namespace std;

# define BLOCKN \
  d13x = xt[ind3]-xt[ind1]; d24x = xt[ind4]-xt[ind2];\
  d13y = yt[ind3]-yt[ind1]; d24y = yt[ind4]-yt[ind2];\
  d13z = zt[ind3]-zt[ind1]; d24z = zt[ind4]-zt[ind2];\
  nx = 0.5 * (d13y*d24z - d13z*d24y);\
  ny = 0.5 * (d13z*d24x - d13x*d24z);\
  nz = 0.5 * (d13x*d24y - d13y*d24x);\
  normn = 1./K_FUNC::E_max(1.e-10,sqrt(nx*nx+ny*ny+nz*nz));\
  nx = nx*normn;\
  ny = ny*normn;\
  nz = nz*normn;\
  vn = ptrVx[indv]*nx+ptrVy[indv]*ny+ptrVz[indv]*nz;\
  ptrBCFX[noindint] = ptrVx[indv] - vn*nx;\
  ptrBCFY[noindint] = ptrVy[indv] - vn*ny;\
  ptrBCFZ[noindint] = ptrVz[indv] - vn*nz;\
  noindint++;\

//=============================================================================
/* nullify a vector Q at BC faces */
//=============================================================================
PyObject* K_CONVERTER::nullifyVectorAtBCFaceStruct(PyObject* self, PyObject* args)
{
  E_Int imin, imax, jmin, jmax, kmin, kmax;
  PyObject *zone, *dataBCX, *dataBCY, *dataBCZ;
  E_Int Loc;
  char *vectx, *vecty, *vectz; 
  char *GridCoordinates, *FlowSolutionNodes, *FlowSolutionCenters;
  if (!PYPARSETUPLE_(args, OOOO_ IIII_ III_ SSSS_ SS_, 
                    &zone, &dataBCX, &dataBCY, &dataBCZ, 
                    &imin, &imax, &jmin, &jmax, &kmin, &kmax, &Loc,
                    &vectx, &vecty, &vectz, &GridCoordinates, &FlowSolutionNodes, &FlowSolutionCenters)) return NULL;

  vector<PyArrayObject*> hook;
  E_Int im, jm, km, cnSize, cnNfld;
  char* varString; char* eltType;
  vector<E_Float*> fields; vector<E_Int> locs;
  vector<E_Int*> cn;
  K_PYTREE::getFromZone(zone, 1, Loc, varString, fields, locs, im, jm, km, 
                        cn, cnSize, cnNfld, eltType, hook, GridCoordinates, 
                        FlowSolutionNodes, FlowSolutionCenters);
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);
  if (posx == -1 || posy == -1 || posz == -1) 
  {
    PyErr_SetString(PyExc_TypeError, 
                    "nullifyVectorAtBCFaceStruct: coordinates not found in zone.");
    RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL);
    return NULL; 
  }
  E_Float* xt = fields[posx];
  E_Float* yt = fields[posy];
  E_Float* zt = fields[posz];

  E_Int posvx = K_ARRAY::isNamePresent(vectx,varString);
  E_Int posvy = K_ARRAY::isNamePresent(vecty,varString);
  E_Int posvz = K_ARRAY::isNamePresent(vectz,varString);
  if (posvx == -1 || posvy == -1 || posvz == -1) 
  {
    PyErr_SetString(PyExc_TypeError, 
                    "nullifyVectorAtBCFaceStruct: vector not found in zone.");
    RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL);
    return NULL; 
  }
  E_Float* ptrVx = fields[posvx];
  E_Float* ptrVy = fields[posvy];
  E_Float* ptrVz = fields[posvz];

  // flow field at BC (Vertex or Faces)
  FldArrayF *fXInt, *fYInt, *fZInt;
  K_NUMPY::getFromNumpyArray(dataBCX, fXInt);
  K_NUMPY::getFromNumpyArray(dataBCY, fYInt);
  K_NUMPY::getFromNumpyArray(dataBCZ, fZInt);
  E_Float* ptrBCFX = fXInt->begin();
  E_Float* ptrBCFY = fYInt->begin();
  E_Float* ptrBCFZ = fZInt->begin();

  E_Int imjm = im*jm;
  E_Int noindint = 0;
  E_Int indv, shift;
  E_Float d13x, d13y, d13z, d24x, d24y, d24z, vn, normn, nx, ny, nz;
  E_Int ind1, ind2, ind3, ind4;
  if ( Loc == 0)// vertex
  {
    if ( imin == imax )
    {
      for (E_Int k = kmin-1; k < kmax; k++)
        for (E_Int j = jmin-1; j < jmax; j++)
        {
          indv = imin-1+j*im+k*imjm;
          ind1 = indv; ind2 = ind1+im; ind3 = ind2+imjm; ind4 = ind1+imjm;
          BLOCKN;
        }
    }
    else if ( jmin == jmax)
    {
      for (E_Int k = kmin-1; k < kmax; k++)
        for (E_Int i = imin-1; i < imax; i++)
        {
          indv = i+(jmin-1)*im+k*imjm;
          ind1 = indv; ind2 = ind1+imjm; ind3 = ind2+1; ind4 = ind1+1;
          BLOCKN;
        }
    }
    else
    {
      for (E_Int j = jmin-1; j < jmax; j++)
        for (E_Int i = imin-1; i < imax; i++)
        {
          indv = i+j*im+(kmin-1)*imjm;
          ind1 = indv; ind2 = ind1+1; ind3 = ind2+im; ind4 = ind1+im;
          BLOCKN;
        }      
    }
  }
  else //Face centers
  {
    nx = 1.; ny = 0.; nz= 0.;
    E_Int im1 = K_FUNC::E_max(1, im-1); 
    E_Int jm1 = K_FUNC::E_max(1, jm-1); 
    E_Int km1 = K_FUNC::E_max(1, km-1); 
    E_Int im1jm1 = im1*jm1;
    if (imin == imax)
    {
      if (imin == 1) shift = 0;
      else shift= im1-1;  
      for (E_Int k = kmin-1; k < kmax-1; k++)
        for (E_Int j = jmin-1; j < jmax-1; j++)
        {
          indv=shift+j*im1+k*im1jm1;
          ind1 = imin-1+j*im+k*imjm; ind2 = ind1+im; ind3 = ind2+imjm; ind4 = ind1+imjm;
          BLOCKN;
        }
    }    
    else if ( jmin == jmax)
    {
      if (jmin == 1) shift = 0;
      else shift= (jm1-1)*im1;  
      for (E_Int k = kmin-1; k < kmax-1; k++)
        for (E_Int i = imin-1; i < imax-1; i++)
        {
          indv=i+shift+k*im1jm1;
          ind1 = i+(jmin-1)*im+k*imjm; ind2 = ind1+imjm; ind3 = ind2+1; ind4 = ind1+1;
          BLOCKN;
        }
    }      
    
    else if ( kmin == kmax)
    {
      if (kmin == 1) shift = 0;
      else shift= (km1-1)*im1jm1;
      for (E_Int j = jmin-1; j < jmax-1; j++)
        for (E_Int i = imin-1; i < imax-1; i++)
        {
          indv=i+j*im1+shift;
          ind1 = i+j*im+(kmin-1)*imjm; ind2 = ind1+1; ind3 = ind2+im; ind4 = ind1+im;
          BLOCKN;
        }      
    }
  }
  RELEASESHAREDN(dataBCX, fXInt);
  RELEASESHAREDN(dataBCY, fYInt);
  RELEASESHAREDN(dataBCZ, fZInt);
  RELEASESHAREDZ(hook, (char*)NULL, (char*)NULL);
  Py_INCREF(Py_None);
  return Py_None;  
}

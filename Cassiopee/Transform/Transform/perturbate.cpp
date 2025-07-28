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

using namespace std;
using namespace K_FUNC;
using namespace K_FLD;

// ============================================================================
/* Perturbate array coordinates describing a mesh */
/* Perturbate a mesh randomly with given radius */
// ============================================================================
PyObject* K_TRANSFORM::perturbate(PyObject* self, PyObject* args)
{
  E_Float radius; PyObject* array;
  E_Int dim;
  if (!PYPARSETUPLE_(args, O_ R_ I_,
                    &array, &radius, &dim))
  {
    return NULL;
  }
  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  E_Int res;
  res = K_ARRAY::getFromArray3(array, varString, f, im, jm, km, cn, eltType); 

  E_Int posx, posy, posz;
  E_Float rx, ry, rz;
  E_Float rxi, ryi, rzi, rxj, ryj, rzj, rxk, ryk,rzk;
  E_Int ind, indi, indj, indk;
  E_Float l; 

  if (res == 1)
  {
    posx = K_ARRAY::isCoordinateXPresent(varString);
    posy = K_ARRAY::isCoordinateYPresent(varString);
    posz = K_ARRAY::isCoordinateZPresent(varString);
    if (posx == -1 || posy == -1 || posz == -1)
    {
      RELEASESHAREDS(array, f);
      PyErr_SetString(PyExc_TypeError,
                      "perturbate: can't find coordinates in array.");
      return NULL;
    }
    posx++; posy++; posz++;
    
    E_Float* fx = f->begin(posx);
    E_Float* fy = f->begin(posy);
    E_Float* fz = f->begin(posz);

    // Random perturbation
    E_LONG idum = -1;
    // Boundaries are not modified
    if (km == 2 || km == 1)
    {
      for (E_Int j = 1; j < jm-1; j++)
        for (E_Int i = 1; i < im-1; i++)
        {
          rx = K_NOISE::stdRand(&idum);
          ry = 0.; rz = 0.;
          if (dim > 1) ry = K_NOISE::stdRand(&idum);
          if (dim > 2) rz = K_NOISE::stdRand(&idum);

          for (E_Int k = 0; k < km; k++)
          {
            ind  = i + j*im + k*im*jm;
            indi = ind + 1;
            indj = ind + im;
            rxi = fx[indi]-fx[ind];
            ryi = fy[indi]-fy[ind];
            rxj = fx[indj]-fx[ind];
            ryj = fy[indj]-fy[ind];
            fx[ind] += radius*(rx*rxi+ry*rxj);
            fy[ind] += radius*(rx*ryi+ry*ryj);
          }
        }
    }
    else
    {
      for (E_Int k = 1; k < km-1; k++)
        for (E_Int j = 1; j < jm-1; j++)
          for (E_Int i = 1; i < im-1; i++)
          {
            rx = K_NOISE::stdRand(&idum)-0.5;
            ry = 0.; rz = 0.;
            if (dim > 1) ry = K_NOISE::stdRand(&idum)-0.5;
            if (dim > 2) rz = K_NOISE::stdRand(&idum)-0.5;  
            
            ind = i + j*im + k*im*jm;
            indi = i+1 + j*im + k*im*jm;
            indj = i + (j+1)*im + k*im*jm;
            indk = i + j*im + (k+1)*im*jm;
            rxi = fx[indi]-fx[ind];
            ryi = fy[indi]-fy[ind];
            rzi = fz[indi]-fz[ind];
            
            rxj = fx[indj]-fx[ind];
            ryj = fy[indj]-fy[ind];
            rzj = fz[indj]-fz[ind];
           
            rxk = fx[indk]-fx[ind];
            ryk = fy[indk]-fy[ind];
            rzk = fz[indk]-fz[ind];
            
            fx[ind] += radius*(rx*rxi+ry*rxj+rz*rxk);
            fy[ind] += radius*(rx*ryi+ry*ryj+rz*ryk);
            fz[ind] += radius*(rx*rzi+ry*rzj+rz*rzk);
          }
    }
      
    // Build array
    PyObject* tpl = K_ARRAY::buildArray(*f, varString, im, jm, km);
    RELEASESHAREDS(array, f);  
    return tpl;
  }
  else if (res == 2)
  {
    posx = K_ARRAY::isCoordinateXPresent(varString);
    posy = K_ARRAY::isCoordinateYPresent(varString);
    posz = K_ARRAY::isCoordinateZPresent(varString);
    if (posx == -1 || posy == -1 || posz == -1)
    {
      RELEASESHAREDU(array, f, cn);
      PyErr_SetString(PyExc_TypeError,
                      "perturbate: can't find coordinates in array.");
      return NULL;
    }
    posx++; posy++; posz++;
    
    E_Float* fx = f->begin(posx);
    E_Float* fy = f->begin(posy);
    E_Float* fz = f->begin(posz);
    
    // Random perturbation
    E_LONG idum = -1;

    E_Int npts = f->getSize();
    vector< vector<E_Int> > cVN(npts);
    if (strcmp(eltType, "NGON") == 0) K_CONNECT::connectNG2VNbrs(*cn, cVN);
    else K_CONNECT::connectEV2VNbrs(*cn, cVN);
    E_Float lx, ly, lz;

    for (E_Int i = 0; i < npts; i++)
    {
      rx = K_NOISE::stdRand(&idum);
      ry = 0.; rz = 0.;
      if (dim > 1) ry = K_NOISE::stdRand(&idum);
      if (dim > 2) rz = K_NOISE::stdRand(&idum);
      vector<E_Int>& pt = cVN[i];
      l = 1.e6;
      for (size_t j = 0; j < pt.size(); j++)
      {
        ind = pt[j]-1;
        lx = fx[i] - fx[ind];
        ly = fy[i] - fy[ind];
        lz = fz[i] - fz[ind];
        l = std::min(l, lx*lx+ly*ly+lz*lz);
      }
      l = sqrt(l);
      fx[i] += radius*rx*l;
      fy[i] += radius*ry*l;
      fz[i] += radius*rz*l;
    }
    
    // Build array
    PyObject* tpl = K_ARRAY::buildArray(*f, varString, *cn, -1, eltType);
    RELEASESHAREDU(array, f, cn);
    return tpl;
  }
  else
  {
    PyErr_SetString(PyExc_TypeError,
                    "perturbate: invalid array.");
    return NULL;
  }
}

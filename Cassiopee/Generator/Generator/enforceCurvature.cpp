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
# include <vector>

using namespace std;
using namespace K_CONST;
using namespace K_FUNC;
using namespace K_FLD;

// ============================================================================
/* Enforce curvature of a 1D curve in a distribution */
// ============================================================================
PyObject* K_GENERATOR::enforceCurvature(PyObject* self, PyObject* args)
{
  PyObject* arrayD; PyObject* arrayC;
  E_Float power;
  if (!PYPARSETUPLE_(args, OO_ R_, &arrayD, &arrayC, &power)) return NULL;
  
  // Check array
  E_Int im, jm, km, imd, jmd, kmd;
  FldArrayF* fd; FldArrayF* fc;
  FldArrayI* cnd; FldArrayI* cnc;
  char *varStringd; char *varStringc;
  char* eltTypec; char* eltTyped;
  FldArrayF coord; 
  vector<E_Int> posc; vector<E_Int> posd;
  
  E_Int resc = K_ARRAY::getFromArray3(arrayC, varStringc, fc, im, jm, km,
                                      cnc, eltTypec); 
  E_Int resd = K_ARRAY::getFromArray3(arrayD, varStringd, fd, imd, jmd, kmd,
                                      cnd, eltTyped);
  
  FldArrayF s, cm, xd, xrk, crk;
  E_Int posxc, posyc, poszc;
  E_Int posxd, posyd, poszd;
  if (resc == 1 && resd == 1)
  {
    posxc = K_ARRAY::isCoordinateXPresent(varStringc);
    posyc = K_ARRAY::isCoordinateYPresent(varStringc);
    poszc = K_ARRAY::isCoordinateZPresent(varStringc);
    posxd = K_ARRAY::isCoordinateXPresent(varStringd);
    posyd = K_ARRAY::isCoordinateYPresent(varStringd);
    poszd = K_ARRAY::isCoordinateZPresent(varStringd);

    if (posxc == -1 || posyc == -1 || poszc == -1 ||
        posxd == -1 || posyd == -1 || poszd == -1)
    {
      RELEASESHAREDS(arrayC, fc);
      RELEASESHAREDS(arrayD, fd);
      PyErr_SetString(PyExc_TypeError,
                      "enforceCurvature: coordinates not found.");
      return NULL;
    }
    
    posxc++; posyc++; poszc++;
    posxd++; posyd++; poszd++;

    // Distribution data
    E_Int N = imd;

    // Compute curvilinear absciss of curve
    s.malloc(im);
    E_Float l;
    s[0] = 0.;

    E_Float* xc = fc->begin(posxc);
    E_Float* yc = fc->begin(posyc);
    E_Float* zc = fc->begin(poszc);
    E_Float xx, yy, zz;
    for (E_Int i = 0; i < im-1; i++)
    {
      xx = xc[i+1] - xc[i];
      yy = yc[i+1] - yc[i];
      zz = zc[i+1] - zc[i];
      l =  xx*xx + yy*yy + zz*zz;
      s[i+1] = s[i] + sqrt(l);
    }   
    E_Float ltot = s[im-1];
    for (E_Int i = 0; i < im; i++) s[i] = s[i] / ltot;

    // Compute angle
    FldArrayF angle(im);
    E_Float dirVect[3];
    dirVect[0] = 0; dirVect[1] = 0; dirVect[2] = 1;
    K_COMPGEOM::compCurvatureAngle(im, xc, yc, zc, dirVect, angle);

    // Calcul des split points
    E_Float splitAngle = 60.; //  degrees
    list<E_Int> splitIndex;
    for (E_Int i = 1; i < im-1; i++) 
    {
      if (angle[i] < 180-splitAngle) splitIndex.push_back(i);
      if (angle[i] > 180+splitAngle) splitIndex.push_back(i);
    }

    // Compute curvature of curve
    FldArrayF curv(im);
    K_COMPGEOM::compCurvature(im, xc, yc, zc, curv);
    //for (E_Int i = 0; i < im; i++) printf("%d %f\n", i, 1./curv[i]);

    E_Float* c = curv.begin();
    for (E_Int i = 0; i < im; i++) c[i] = K_FUNC::E_abs(c[i]);
    for (E_Int i = 0; i < im; i++) c[i] = K_FUNC::E_max(c[i], 1.e-6);
    for (E_Int i = 0; i < im; i++) c[i] = K_FUNC::E_min(c[i], 1.e6);
    for (E_Int i = 0; i < im; i++) c[i] = 1./c[i];
    for (E_Int i = 0; i < im; i++) c[i] = pow(c[i], power);

    // Lissage
    FldArrayF curv2(im); curv2 = curv;
    E_Float* c2 = curv2.begin();
    for (E_Int nit = 0; nit < 5; nit++)
    {
      for (E_Int i = 1; i < im-1; i++)
      {
        c2[i] = c2[i] + 0.1*(c[i-1]-2*c[i]+c[i+1]);
      }
      curv = curv2;
    }

    // Normalisation de la courbure
    E_Float sum = 0.;
    for (E_Int i = 0; i < im-1; i++)
    {
      sum = sum + 0.5*(c[i]+c[i+1])*(s[i+1]-s[i]);
    }
    for (E_Int i = 0; i < im; i++) c[i] = c[i]/(sum*(N-1));
    
    // Min max courbure
    E_Float cmin = 1.e6;
    for (E_Int i = 0; i < im; i++) cmin = K_FUNC::E_min(cmin, c[i]);
    E_Float cmax = -1.e6;
    for (E_Int i = 0; i < im; i++) cmax = K_FUNC::E_max(cmax, c[i]);

    // Minimum step: on autorise un facteur 2 au min par rapport au regulier
    //E_Float hmin = (1./(N-1.))/2.;
    // Maximum step: on autorise un facteur 2 au max
    E_Float hmax = (1./(N-1.))*2.;

    // Redistribute points following curvature curv
    E_Float alphaMin = 0.1;
    E_Float alphaMax = 40.;
    E_Float alpha = 0.5*(alphaMin+alphaMax);
    E_Int is;
    E_Float cr;
    E_Int it = 0;
    xd.malloc(N);
    xrk.malloc(4);
    crk.malloc(4);
    
    while (it < 100)
    {
      xd[0] = 0.;
      is = 0;
      for (E_Int i = 0; i < N-1; i++)
      {
        is = 0;
        fa(xd[i], s, &is);
        xrk[0] = xd[i] + alpha*c[is]*0.5;
        crk[0] = c[is];
        fa(xrk[0], s, &is);
        xrk[1] = xd[i] + alpha*c[is]*0.5;
        crk[1] = c[is];
        fa(xrk[1], s, &is);
        xrk[2] = xd[i] + alpha*c[is]*0.5;
        crk[2] = c[is];
        fa( xrk[2], s, &is);
        crk[3] = c[is];
        cr = alpha/6.*(crk[0]+2*crk[1]+2*crk[2]+crk[3]);
        cr = K_FUNC::E_min(cr, hmax);
        xd[i+1] = xd[i] + cr;
      }
      //printf("%f %f\n", xd[N-1], alpha);
      if (xd[N-1] > 1.+1.e-6)
      {
        alphaMax = 0.5*(alphaMin+alphaMax);
        alpha = 0.5*(alphaMin+alphaMax);
      }
      else if (xd[N-1] < 1.-1.e-6)
      {
        alphaMin = 0.5*(alphaMin+alphaMax);
        alpha = 0.5*(alphaMin+alphaMax);
      }
      else break;
      it++;
    }
    
    //printf("upper bound=%f\n", xd[N-1]);
    
    // A mapping between 0 and 1 is enforced
    E_Float a = 1. / xd[N-1];
    E_Float b = 0.;
    for (E_Int i = 0; i < N; i++) xd[i] = a*xd[i]+b;
    xd[N-1] = 1.;

    // Distribution modification
    E_Int Njmd = N * jmd;
    for (E_Int i = 0; i < N; i++)
      for (E_Int j = 0; j < jmd; j++)
        for (E_Int k = 0; k < kmd; k++)
        {
          E_Int ind = i + j*N + k*Njmd;
          (*fd)(ind, posxd) = xd[i];
        }

    // Enforce des split points
    for (list<E_Int>::iterator it = splitIndex.begin();
         it != splitIndex.end(); it++)
    {
      FldArrayF* out = new FldArrayF();
      E_Float split = s[*it];
      E_Int is = 0;
      fa(split, *fd, &is);
      //E_Float eh = alpha*c[*it];
      E_Float eh = K_FUNC::E_min((*fd)[is+1]-(*fd)[is], (*fd)[is]-(*fd)[is-1]);
      E_Int supp = K_FUNC::E_max(N/10, 10);
      E_Int add = 0;
      E_Int niout, njout, nkout;
      enforceCommon("enforceX", varStringd, N, jmd, kmd, 
		    *fd, split-0.5*eh, eh, supp, add, 
		    *out, niout, njout, nkout, 0);
      (*fd) = (*out);
      delete out;
      N = niout;
    }

    // Enforce loop point (if any)
    E_Float pt = (xc[0]-xc[im-1])*(xc[0]-xc[im-1])+
      (yc[0]-yc[im-1])*(yc[0]-yc[im-1])+
      (zc[0]-zc[im-1])*(zc[0]-zc[im-1]);
    if (pt < 1.e-10)
    {
      FldArrayF* out = new FldArrayF();
      E_Float eh = K_FUNC::E_min((*fd)[1]-(*fd)[0], (*fd)[N-1]-(*fd)[N-2]);
      E_Int supp = K_FUNC::E_max(N/10, 10);
      E_Int add = 0;
      E_Int niout, njout, nkout;
      enforceCommon( "enforcePlusX", varStringd, N, jmd, kmd, 
                     *fd, 0., eh, supp, add, 
                     *out, niout, njout, nkout, 0 );
      (*fd) = (*out);
      N = niout;
      enforceCommon( "enforceMoinsX", varStringd, N, jmd, kmd, 
                     *fd, 1., eh, supp, add, 
                     *out, niout, njout, nkout, 0 );
      (*fd) = (*out);
      N = niout;
      delete out;
    }

    // Build array
    PyObject* tpl = K_ARRAY::buildArray3(*fd, varStringd, N, jmd, kmd);
    RELEASESHAREDS(arrayC, fc);
    RELEASESHAREDS(arrayD, fd);
    return tpl;
  }
  else if (resc == 2 || resd == 2)
  {
    RELEASESHAREDB(resc, arrayC, fc, cnc);
    RELEASESHAREDB(resd, arrayD, fd, cnd);
      
    PyErr_SetString(PyExc_TypeError,
                    "enforceCurvature: not used for unstructured arrays.");
    return NULL;
  }
  else 
  {
    PyErr_SetString(PyExc_TypeError,
                    "enforceCurvature: invalid type of array.");
    return NULL;
  }
}

//=============================================================================
// Curvilinear absciss search
// x: point coordinate we are looking for
// s: curvilinear absciss
// i is start for search
// return the index of x in i
//=============================================================================
void K_GENERATOR::fa(E_Float x, K_FLD::FldArrayF& s, E_Int* i)
{
  E_Int size = s.getSize();
  while (*i < size && s[*i] < x) *i = *i+1;
  *i = *i-1;
  if (*i < 0) *i = 0;
}

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

#include "generator.h"
using namespace K_FLD;

//============================================================================
/* Maillage en chapeau cousu (avec une couture au milieu) */
//============================================================================
PyObject* K_GENERATOR::stitchedHat(PyObject* self, PyObject* args)
{
  PyObject* array; E_Float epsilon;
  if (!PYPARSETUPLE_(args, O_ R_, &array, &epsilon)) return NULL;
 
  // Check array
  E_Int im, jm, km;
  FldArrayF* f; FldArrayI* cn;
  char* varString; char* eltType;
  
  E_Int res = 
    K_ARRAY::getFromArray3(array, varString, f, im, jm, km, cn, eltType);

  if (res != 1 && res != 2)
  {
    PyErr_SetString(PyExc_TypeError,
      "stitchedHat: invalid array.");
    return NULL;
  }

  if (res != 1)
  {
    RELEASESHAREDU(array, f, cn);
    PyErr_SetString(PyExc_TypeError,
                    "stitchedHat: array must be a structured array.");
    return NULL;
  }
  
  // traitement
  if (im < 1 || jm < 1 || km != 1)
  {
    delete f;
    PyErr_SetString(PyExc_TypeError,
    "stitchedHat: array must be a i or a ij-array.");
    return NULL;
  }
  if (f->getNfld() != 3) 
  {
    printf("Warning: stitchedHat: only coordinates are considered.\n");
  }
  E_Int posx = K_ARRAY::isCoordinateXPresent(varString);
  E_Int posy = K_ARRAY::isCoordinateYPresent(varString);
  E_Int posz = K_ARRAY::isCoordinateZPresent(varString);    
  if (posx == -1 || posy == -1 || posz == -1)
  {
    delete f;
    PyErr_SetString(PyExc_TypeError,
                    "stitchedHat: can't find coordinates in array.");
    return NULL;
  }
  
  posx++; posy++; posz++;
  E_Int api = f->getApi();
  E_Int npts = im*jm;
  FldArrayF* sol = new FldArrayF(2*npts,3);
  E_Float* xt = sol->begin(1);
  E_Float* yt = sol->begin(2);
  E_Float* zt = sol->begin(3);
  E_Float* x0 = f->begin(posx);
  E_Float* y0 = f->begin(posy);
  E_Float* z0 = f->begin(posz);
  E_Float xb, yb ,zb, d1, d2;
  E_Float px, py, pz;
  E_Float* nx = new E_Float[npts];
  E_Float* ny = new E_Float[npts];
  E_Float* nz = new E_Float[npts];
  E_Int* found = new E_Int[npts];
  E_Int* hanging = new E_Int[npts];
  E_Int indnear;

  // Dimensionnement de eps (pour la recherche)
  E_Float dd = 0.;
  for (E_Int i = 0; i < npts; i++)
  {
    for (E_Int j = 0; j < npts; j++)
      dd = 
        K_FUNC::E_max(dd, (x0[i]-x0[j])*(x0[i]-x0[j])+ 
                      (y0[i]-y0[j])*(y0[i]-y0[j]) + 
                      (z0[i]-z0[j])*(z0[i]-z0[j]) );
  }
  dd = sqrt(dd);
  E_Int Nd = 3000;
  E_Float eps = K_FUNC::E_min(dd/Nd, epsilon);
  Nd = E_Int(dd/eps);

  // Calcul des normales
  normals(npts, x0, y0, z0, nx, ny, nz);

  for (E_Int i = 0; i < npts; i++)
  { found[i] = 0; hanging[i] = 0; }

  // Calcul de l'axe median
  for (E_Int i = 0; i < npts; i++)
  {
    xt[i] = x0[i]; yt[i] = y0[i]; zt[i] = z0[i];

    if (found[i] == 0)
    {
      xb = xt[i]; yb = yt[i]; zb = zt[i];

      for (E_Int e = 0; e < Nd; e++)
      {
        px = x0[i] + e*eps*nx[i];
        py = y0[i] + e*eps*ny[i];
        pz = z0[i] + e*eps*nz[i];
        d1 = (x0[i]-px)*(x0[i]-px)+(y0[i]-py)*(y0[i]-py)+
          (z0[i]-pz)*(z0[i]-pz);
        d2 = dist2Curve(npts, x0, y0, z0, found, px, py, pz, i, indnear);
        d1 = sqrt(d1); d2 = sqrt(d2);
        if ( K_FUNC::fEqualZero(d1-d2, epsilon) == true )
        {
          xb = px; yb = py; zb = pz;
          if (found[indnear] == 0)
          {
            xt[indnear+npts] = xb; 
            yt[indnear+npts] = yb; 
            zt[indnear+npts] = zb;
            found[indnear] = 1;
          }
          else hanging[i] = 1;
          found[i] = 1;
          break;
        }
      }
      //if (i == 0) printf("%f %f %f %d\n", xb,yb,zb,found[i]);
      xt[i+npts] = xb; yt[i+npts] = yb; zt[i+npts] = zb;
    }
  }

  //for (E_Int i = 0; i < npts; i++)
  //  if (found[i] == 0) printf("point %d has not been projected.\n", i);
  
  delete [] nx; delete [] ny; delete [] nz; delete [] found;

  // Projection des pts orphelins
  E_Int orph;
  E_Float d;
  E_Float dmin = 1.e6; E_Int indmin = -1;
  for (E_Int i = 0; i < npts; i++)
  {
    orph = hanging[i]; dmin = 1.e6;
    for (E_Int j = 0; j < i; j++)
    {
      d = (xt[npts+i] - xt[npts+j])*(xt[npts+i] - xt[npts+j])+
          (yt[npts+i] - yt[npts+j])*(yt[npts+i] - yt[npts+j])+
          (zt[npts+i] - zt[npts+j])*(zt[npts+i] - zt[npts+j]);
      if (d < dmin && hanging[j] == 0) {dmin = d; indmin = j;}
      //printf("%f %f\n", d, eps);
      //if (d < eps*eps) {orph = 1; break;}
    } 
    for (E_Int j = i+1; j < npts; j++)
    {
      d = (xt[npts+i] - xt[npts+j])*(xt[npts+i] - xt[npts+j])+
          (yt[npts+i] - yt[npts+j])*(yt[npts+i] - yt[npts+j])+
          (zt[npts+i] - zt[npts+j])*(zt[npts+i] - zt[npts+j]);
      if (d < dmin && hanging[j] == 0) {dmin = d; indmin = j;}
      //printf("%f %f\n", d, eps);
      //if (d < eps*eps) {orph = 1; break;}
    }
    if (orph == 1)
    {
      //printf("pt %d est orphelin, nearest %d \n", i, indmin);
      xt[i+npts] = xt[indmin+npts];
      yt[i+npts] = yt[indmin+npts];
      zt[i+npts] = zt[indmin+npts];
    }
  }

  E_Int im2 = im; E_Int jm2 = jm; E_Int km2 = 2;
  if (jm == 1) {jm2 = 2; km2 = 1;}
  PyObject* tpl = K_ARRAY::buildArray3(*sol, "x,y,z", im2, jm2, km2, api);
  delete sol;
  RELEASESHAREDS(array, f);  
  delete [] hanging;
  return tpl;
}

//=============================================================================
void K_GENERATOR::normals(E_Int n, E_Float* x, E_Float* y, E_Float* z,
                          E_Float* nx, E_Float *ny, E_Float* nz)
{
  E_Float p1, p2, p3, p4, p5, p6, n1, n2, n3, r, ri;
  E_Float n1ref, n2ref, n3ref, alpha, beta, nrefp;
  E_Int ind = 0;

  // Premiere normale de reference valide
  ind = 1;
  r = 0;

  while (r <= 1.e-9)
  {
    p1 = x[ind+1]-x[ind];
    p2 = y[ind+1]-y[ind];
    p3 = z[ind+1]-z[ind];
    ri = 1./sqrt(p1*p1+p2*p2+p3*p3);
    p1 = p1*ri; p2 = p2*ri; p3 = p3*ri;
    p4 = x[ind-1]-x[ind];
    p5 = y[ind-1]-y[ind];
    p6 = z[ind-1]-z[ind];
    ri = 1./sqrt(p4*p4+p5*p5+p6*p6);
    p4 = p4*ri; p5 = p5*ri; p6 = p6*ri;
    n1ref = p1+p4; n2ref = p2+p5; n3ref = p3+p6;
    r = sqrt(n1ref*n1ref + n2ref*n2ref + n3ref*n3ref);
    if (r > 1.e-9)
    {
      ri = 1./r;
      n1ref = n1ref * ri; n2ref = n2ref * ri; n3ref = n3ref * ri;
      //printf("normals ref : %d, %g %g %g\n", ind, n1ref, n2ref, n3ref);
      break;
    }
    ind++;
  }

  for (E_Int i = 1; i < n-1; i++)
  {
    p1 = x[i+1]-x[i];
    p2 = y[i+1]-y[i];
    p3 = z[i+1]-z[i];
    ri = 1./sqrt(p1*p1+p2*p2+p3*p3);
    p1 = p1*ri; p2 = p2*ri; p3 = p3*ri;
    p4 = x[i-1]-x[i];
    p5 = y[i-1]-y[i];
    p6 = z[i-1]-z[i];
    ri = 1./sqrt(p4*p4+p5*p5+p6*p6);
    p4 = p4*ri; p5 = p5*ri; p6 = p6*ri;
    n1 = p1+p4; n2 = p2+p5; n3 = p3+p6;
    r = n1*n1 + n2*n2 + n3*n3;
   
    if (r < 1.e-6) // alignes...
    {
      nrefp = n1ref*p1 + n2ref*p2 + n3ref*p3;
      beta = 1./sqrt(1.-(nrefp*nrefp));
      alpha = -beta*nrefp;
      n1 = alpha*p1+beta*n1ref; 
      n2 = alpha*p2+beta*n2ref;
      n3 = alpha*p3+beta*n3ref;
      r = n1*n1 + n2*n2 + n3*n3;
      //printf("check : %f %f\n", p1*p1+p2*p2+p3*p3, n1ref*n1ref+n2ref*n2ref+n3ref*n3ref);
      //printf("must be 1 : %f\n", r);
    }
    else 
    {
      ri = 1./sqrt(r);
      n1 = n1*ri; n2 = n2*ri; n3 = n3*ri;
      //printf("std calculus must be 1 : %f\n", n1*n1+n2*n2+n3*n3);
    }
    if (n1*n1ref + n2*n2ref + n3*n3ref < 0)
    {n1 = -n1; n2 = -n2; n3 = -n3;}
    n1ref = n1; n2ref = n2; n3ref = n3;

    nx[i] = n1;
    ny[i] = n2;
    nz[i] = n3;
  }
  
  p1 = x[1]-x[0];
  p2 = y[1]-y[0];
  p3 = z[1]-z[0];
  ri = 1./sqrt(p1*p1+p2*p2+p3*p3);
  p1 = p1*ri; p2 = p2*ri; p3 = p3*ri;
  p4 = x[n-2]-x[0];
  p5 = y[n-2]-y[0];
  p6 = z[n-2]-z[0];
  ri = 1./sqrt(p4*p4+p5*p5+p6*p6);
  p4 = p4*ri; p5 = p5*ri; p6 = p6*ri;
  n1 = p1+p4; n2 = p2+p5; n3 = p3+p6;
  r = n1*n1 + n2*n2 + n3*n3;

  if (r < 1.e-6) // alignes...
  {
    nrefp = n1ref*p1 + n2ref*p2 + n3ref*p3;
    beta = 1./sqrt(1.-nrefp*nrefp);
    alpha = -beta*nrefp;
    n1 = alpha*p1+beta*n1ref; 
    n2 = alpha*p2+beta*n2ref; 
    n3 = alpha*p3+beta*n3ref;
    r = n1*n1 + n2*n2 + n3*n3; // must be 1
  }
  else 
  {
    ri = 1./sqrt(r);
    n1 = n1*ri; n2 = n2*ri; n3 = n3*ri;
  }
  if (n1*n1ref + n2*n2ref + n3*n3ref < 0)
  {n1 = -n1; n2 = -n2; n3 = -n3;}
  n1ref = n1; n2ref = n2; n3ref = n3;

  //printf("%f %f %f \n", n1,n2,n3);
  nx[0] = n1;
  ny[0] = n2;
  nz[0] = n3;
  nx[n-1] = n1;
  ny[n-1] = n2;
  nz[n-1] = n3;
}

//=============================================================================
// Calcul le point le plus proche d'un point (px,py,pz) sur un courbe
// (x,y,z) mais different de ind
// Retourne la distance au carre
//=============================================================================
E_Float K_GENERATOR::dist2Curve(E_Int n, E_Float* x, E_Float* y, E_Float* z, E_Int* found,
                                E_Float px, E_Float py, E_Float pz, E_Int ind,
                                E_Int& nearest)
{
  E_Float d = 1.e6;
  nearest = -1;
  E_Float dp;

  if (ind == 0)
  {
    // On doit eviter 0 et n-1
    for (E_Int i = 1; i < n-1; i++)
    {
      dp = (x[i]-px)*(x[i]-px) + (y[i]-py)*(y[i]-py) +
        (z[i]-pz)*(z[i]-pz);
      if (dp < d) {d = dp; nearest = i;}
    }
  }
  else
  {
    // On doit eviter ind
    for (E_Int i = 0; i < ind; i++)
    {
      dp = (x[i]-px)*(x[i]-px) + (y[i]-py)*(y[i]-py) +
        (z[i]-pz)*(z[i]-pz);
      if (dp < d) {d = dp; nearest = i;}
    }
    for (E_Int i = ind+1; i < n; i++)
    {
      dp = (x[i]-px)*(x[i]-px) + (y[i]-py)*(y[i]-py) +
        (z[i]-pz)*(z[i]-pz);
      if (dp < d) {d = dp; nearest = i;}
    }
  }
  return d;
}

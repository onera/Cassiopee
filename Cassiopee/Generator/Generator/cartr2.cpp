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

// Grille cartesienne avec facteurs d'expansion

#include "generator.h"

using namespace std;
using namespace K_FLD;
using namespace K_FUNC; 

// ============================================================================
/* Create a cartesian mesh of nixnjxnk points 
   IN: x0, y0, z0: origine de la grille
   IN: hi, hj, hk: pas initial de la grille 
   IN: ri, rj, rk: facteur d'expansion dans chaque direction
   IN: xf, yf, zf: fin de la grille

   OUT: array definissant le maillage cree. 
   On cherche le nombre de points pour créer le maillage, on rectifie les 
   facteurs d'expansion pour respecter la condition de fin de grille */
// ============================================================================

// sans double h
E_Float f1(E_Float Xo, E_Float Xf, E_Float H, E_Float N, E_Float R)
{
  return pow(R,N-1) + (1.-R)*fabs(Xf-Xo)/H - 1.;
}
E_Int getN1(E_Float Xo, E_Float Xf, E_Float H, E_Float R)
{
  E_Float res = (1.-R)*fabs(Xf-Xo)/H - 1.;
  res = log(-res) / log(R);
  E_Int N = rint(res);
  return N+1;
}
E_Float f1prime(E_Float Xo, E_Float Xf, E_Float H, E_Float N, E_Float R)
{
  return (N-1)*pow(R,N-2) - fabs(Xf-Xo)/H;
}

// double h left
E_Float f2(E_Float Xo, E_Float Xf, E_Float H, E_Float N, E_Float R)
{
  return pow(R,N-2) + (1. - R)*(-1 + fabs(Xf-Xo)/H) - 1.;
}
E_Int getN2(E_Float Xo, E_Float Xf, E_Float H, E_Float R)
{
  E_Float res = (1. - R)/H*(-H + fabs(Xf-Xo)) - 1.;
  res = log(-res) / log(R);
  E_Int N = rint(res);
  return N+2;
}
E_Float f2prime(E_Float Xo, E_Float Xf, E_Float H, E_Float N, E_Float R)
{
  return (N-2)*pow(R,N-3) - (-1 + fabs(Xf-Xo)/H);
}

// double h right
E_Float f3(E_Float Xo, E_Float Xf, E_Float H, E_Float N, E_Float R)
{
  return pow(R,N-3)*(2.*R-1) + (1. - R)*fabs(Xf-Xo)/H - 1.;
}
E_Int getN3(E_Float Xo, E_Float Xf, E_Float H, E_Float R)
{
  E_Float res = ( (1. - R)*(fabs(Xf-Xo)/H) - 1.)/(2.*R-1);
  res = log(-res) / log(R);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 
  E_Int N = rint(res);
  return N+3;
}
E_Float f3prime(E_Float Xo, E_Float Xf, E_Float H, E_Float N, E_Float R)
{
  return (N-3.)*(2.*R-1)*pow(R,N-4)+2.*pow(R,N-3)-fabs(Xf-Xo)/H;
}

// double h des deux cotes
E_Float f4(E_Float Xo, E_Float Xf, E_Float H, E_Float N, E_Float R)
{
  return pow(R,N-4)*(2.*R-1.) + R*(1.-fabs(Xf-Xo)/H) + fabs(Xf-Xo)/H -2.;
}
E_Int getN4(E_Float Xo, E_Float Xf, E_Float H, E_Float R)
{
  E_Float res = (-R* ( (fabs(Xf-Xo)/H) - 1.) + fabs(Xf-Xo)/H -2. )/(2.*R-1.);
  res = log(-res) / log(R);
  E_Int N = rint(res);
  return N+4;
}
E_Float f4prime(E_Float Xo, E_Float Xf, E_Float H, E_Float N, E_Float R)
{
  return (N-4.)*(2.*R-1)*pow(R,N-5) +2.*pow(R,N-4)+(1-fabs(Xf-Xo)/H);
}

E_Float Newton(E_Float Xo, E_Float Xf, E_Float H, E_Int N, E_Float Rini, 
               E_Float (*f)(E_Float, E_Float, E_Float, E_Float, E_Float), 
               E_Float (*fprime)(E_Float, E_Float, E_Float, E_Float, E_Float))
{
  E_Float res = Rini;
  //printf("%f \n ",res); fflush(stdout);
  E_Int Nit = 0;
  while (fabs(f(Xo,Xf,H,N,res)) >= 0.000001 && Nit<30)
  {
    res -= f(Xo,Xf,H,N,res) / fprime(Xo,Xf,H,N,res);
    // printf("%f \n ",res) ; fflush(stdout);
    // printf("%f \n ", f(Xo,Xf,H,N,res) ) ; fflush(stdout);
    Nit += 1;
  } 
  // printf("r final methode = %f \n ",res) ; fflush(stdout);
  return res;
}

PyObject* K_GENERATOR::cartr2(PyObject* self, PyObject* args)
{
  E_Float xo, yo, zo;
  E_Float hi, hj, hk;
  E_Float riinput, rjinput, rkinput;
  E_Float xf, yf, zf;
  E_Int doubleLefti, doubleLeftj, doubleLeftk, doubleRighti, doubleRightj, doubleRightk;
  E_Int skeleton = 0;
  E_Int api = 1;
  if (!PYPARSETUPLE_(args, TRRR_ TRRR_ TRRR_ TRRR_ TIII_ TIII_ II_, 
                    &xo, &yo, &zo, &hi, &hj, &hk, &riinput, &rjinput, &rkinput, 
                    &xf, &yf, &zf, 
                    &doubleLefti, &doubleLeftj, &doubleLeftk, 
                    &doubleRighti, &doubleRightj, &doubleRightk,
                    &api, &skeleton))
  {
    return NULL;
  }

  if ( ( (fabs(xf-xo) /hi ) * (riinput -1) + 1 < 0) || ( (fabs(yf-yo) /hj ) * (rjinput -1) +1 < 0 ) ||( (fabs(zf-zo) /hk ) * (rkinput -1) +1 < 0 ) )
  {
    PyErr_SetString(PyExc_ValueError, 
                    "Can not generate mesh.\n Condition not met: (Xf-Xo) /H ) * (R-1) +1 < 0.");
    return NULL;
  }

  //if ((K_FUNC::fEqual(riinput, 1.0) == true) && (K_FUNC::fEqual( fmod((xf-xo), hi),0) == false) )
  //{
  //printf("Remainder of %f / %f is %f\n", xf-xo, hi, fmod((xf-xo),hi));
  //riinput += 0.001;
  //printf("Warning: condition on r not met.\n ");
  //printf("Warning: ri set to %f\n", riinput);
  //}
  
  //if ((K_FUNC::fEqual(rjinput,1.0) == true) && (K_FUNC::fEqual( fmod(yf-yo,hj), 0)==false) )
  //{
  //  rjinput += 0.001;
  //printf("Warning: condition on r not met.\n");
  //printf("Warning: rj set to %f\n", rjinput);
  //}

  //if ((K_FUNC::fEqual(rkinput,1.0) == true) && (K_FUNC::fEqual(fmod(zf-zo,hk), 0)==false) )
  //{
  //  rkinput += 0.001;
  //printf("Warning: condition on r not met.\n");
  //printf("Warning: rk set to %f.\n" , rkinput);
  //}

  //E_Int (*getN)(E_Float Xo, E_Float Xf, E_Float H, E_Float R);
  E_Float (*f)(E_Float Xo, E_Float Xf, E_Float H, E_Float N, E_Float R);
  E_Float (*fprime)(E_Float Xo, E_Float Xf, E_Float H, E_Float N, E_Float R);

  E_Float ri, rj, rk;
  E_Int ni, nj, nk;
  if (K_FUNC::fEqual(riinput, 1.0) == true)
  {
    ni = int(std::abs(xf - xo) / hi)+1;
    if (ni > 1) hi = std::abs(xf - xo)/(ni-1);
    else hi = 1.;
    ri = 1.;
    //if (xf > xo) printf("end: %f %f\n", xf, xo+hi*(ni-1));
    //else printf("end: %f %f\n", xf, xo-hi*(ni-1));
    //printf("h %f -> %f\n", hio, hi);
  }
  else
  {
    if (doubleLefti == 0 && doubleRighti == 0)
    {
      ni = getN1(xo, xf, hi, riinput);
      f = &f1;
      fprime = &f1prime;
    }
    else if (doubleLefti == 1 && doubleRighti == 0)
    {
      ni = getN2(xo, xf, hi, riinput);
      f = &f2;
      fprime = &f2prime;
    }
    else if (doubleLefti == 0 && doubleRighti == 1) 
    {
      ni = getN3(xo, xf, hi, riinput);
      f = &f3;
      fprime = &f3prime;
    }
    else
    {
      ni = getN4(xo, xf, hi, riinput);
      f = &f4;
      fprime = &f4prime;
    }

    ri = Newton(xo,xf,hi,ni,riinput,f,fprime);
    //if (fEqual(ri, 1.) == false) ni += 1;
    //if (xf > xo) printf("end: %f %f\n", xf, xo+hi*((-1. + pow(ri ,ni-1) ) / (-1. + ri)));
    //else printf("end: %f %f\n", xf, xo-hi*((-1. + pow(ri ,ni-1) ) / (-1. + ri)));
    //printf("r %f -> %f\n", riinput, ri);
  } 

  if (K_FUNC::fEqual(rjinput, 1.) == true)
  {
    nj = int(std::abs(yf - yo) / hj)+1;
    if (nj > 1) hj = std::abs(yf - yo)/(nj-1);
    else hj = 1.;
    rj = 1.;
  }
  else
  {
    if (doubleLeftj == 0 && doubleRightj == 0)
    {
      nj = getN1(yo, yf, hj, rjinput);
      f = &f1;
      fprime = &f1prime;
    }
    else if (doubleLeftj == 1 && doubleRightj == 0)
    {
      nj = getN2(yo, yf, hj, rjinput);
      f = &f2;
      fprime = &f2prime;
    }
    else if (doubleLeftj == 0 && doubleRightj == 1) 
    {
      nj = getN3(yo, yf, hj, rjinput);
      f = &f3;
      fprime = &f3prime;
    }
    else
    {
      nj = getN4(yo, yf, hj, rjinput);
      f = &f4;
      fprime = &f4prime;
    }

    rj = Newton(yo,yf,hj,nj,rjinput,f,fprime);
    //if (fEqual(rj, 1.) == false) nj += 1;
  }
  
  if (K_FUNC::fEqual(rkinput, 1.0) == true)
  {
    nk = int(std::abs(zf - zo) / hk)+1;
    if (nk > 1) hk = std::abs(zf - zo)/(nk-1);
    else hk = 1;
    rk = 1.;
  }
  else
  {
    if (doubleLeftk == 0 && doubleRightk == 0)
    {
      nk = getN1(zo, zf, hk, rkinput);
      f = &f1;
      fprime = &f1prime;
    }
    else if (doubleLeftk == 1 && doubleRightk == 0)
    {
      nk = getN2(zo, zf, hk, rkinput);
      f = &f2;
      fprime = &f2prime;
    }
    else if (doubleLeftk == 0 && doubleRightk == 1) 
    {
      nk = getN3(zo, zf, hk, rkinput);
      f = &f3;
      fprime = &f3prime;
    }
    else
    {
      nk = getN4(zo, zf, hk, rkinput);
      f = &f4;
      fprime = &f4prime;
    }

    rk = Newton(zo,zf,hk,nk,rkinput,f,fprime);
    //if (fEqual(rk, 1.) == false) nk += 1;
  }

  if (K_FUNC::fEqual(xo,xf)==true)
  {
    ni = 1; hi = 1.;
  }
  if (K_FUNC::fEqual(yo,yf)==true)
  {
    nj = 1; hj = 1.;
  }
  if (K_FUNC::fEqual(zo,zf)==true)
  {
    nk = 1; hk = 1.;
  }

  //printf("Info: number of mesh cell = %i (i), %i (j), %i (k)\n ",ni,nj,nk); fflush(stdout);
  if (skeleton == 1) return Py_BuildValue("llldddddd", ni, nj, nk, ri, rj, rk, hi, hj, hk);

  E_Int i, j, k, ind;
  // Create cartesian mesh
  PyObject* tpl;
  tpl = K_ARRAY::buildArray3(3, "x,y,z", ni, nj, nk, api);
  
  K_FLD::FldArrayF* f2; K_FLD::FldArrayI* c;
  char* varString; char* eltType;
  K_ARRAY::getFromArray3(tpl, varString, f2, ni, nj, nk, c, eltType);

  E_Int nij = ni*nj;
  E_Int nijk = ni*nj*nk;

  E_Float* xt = f2->begin(1);
  E_Float* yt = f2->begin(2);
  E_Float* zt = f2->begin(3);

  if (K_FUNC::fEqual(ri, 1.0) == true)
  {
    #pragma omp parallel for default(shared) private(k,j,i,ind)
    for (ind = 0; ind < nijk; ind++)
    {
      k = ind/nij;
      j = (ind-k*nij)/ni;
      i = ind-j*ni-k*nij;
      xt[ind] = xo + hi * i; 
    }   
  }   
  else if (doubleLefti == 0 && doubleRighti == 0)
  {
    #pragma omp parallel for default(shared) private(k,j,i,ind)
    for (ind = 0; ind < nijk; ind++)
    {
      k = ind/nij;
      j = (ind-k*nij)/ni;
      i = ind-j*ni-k*nij;
      xt[ind] = xo + hi * ( ( -1. + pow(ri, i) ) / (-1. + ri) ); 
    }
  }
  else if (doubleLefti == 1 && doubleRighti == 0)
  {
    #pragma omp parallel for default(shared) private(k,j,i,ind)
    for (ind = 0; ind < nijk; ind++)
    {
      k = ind/nij;
      j = (ind-k*nij)/ni;
      i = ind-j*ni-k*nij;
      if (i == 0) xt[ind] = xo;
      else xt[ind] = xo + hi + hi * ((-1. + pow(ri, (i-1))) / (-1. + ri));
    }
  }
  else if (doubleLefti == 0 && doubleRighti == 1)
  {
    #pragma omp parallel for default(shared) private(k,j,i,ind)
    for (ind = 0; ind < nijk; ind++)
    {
      k = ind/nij;
      j = (ind-k*nij)/ni;
      i = ind-j*ni-k*nij;
      if (i == ni-1) xt[ind] = xo + hi*pow(ri,ni-3) + hi*((-1.+pow(ri,ni-2))/(-1.+ri));
      else xt[ind] = xo + hi*((-1.+pow(ri,i)) / (-1.+ri));
    }
  }
  else // all double 
  {
     #pragma omp parallel for default(shared) private(k,j,i,ind)
    for (ind = 0; ind < nijk; ind++)
    {
      k = ind/nij;
      j = (ind-k*nij)/ni;
      i = ind-j*ni-k*nij;
      if (i == 0) xt[ind] = xo;
      else if (i == ni-1) xt[ind] = xo + hi + hi*pow(ri,ni-4) + hi*((-1.+pow(ri,ni-3))/(-1.+ri));
      else xt[ind] = xo + hi + hi*((-1.+pow(ri,(i-1))) / (-1.+ri));
    } 
  }

  if (K_FUNC::fEqual(rj, 1.0) == true)
  {
    #pragma omp parallel for default(shared) private(k,j,i,ind)
    for (ind = 0; ind < nijk; ind++)
    {
      k = ind/nij;
      j = (ind-k*nij)/ni;
      i = ind-j*ni-k*nij;
      yt[ind] = yo + hj * j; 
    }   
  }   
  else if (doubleLeftj == 0 && doubleRightj == 0)
  {
    #pragma omp parallel for default(shared) private(k,j,i,ind)
    for (ind = 0; ind < nijk; ind++)
    {
      k = ind/nij;
      j = (ind-k*nij)/ni;
      i = ind-j*ni-k*nij;
      yt[ind] = yo + hj * ( ( -1. + pow(rj, j) ) / (-1. + rj) ); 
    }
  }
  else if (doubleLeftj == 1 && doubleRightj == 0)
  {
    #pragma omp parallel for default(shared) private(k,j,i,ind)
    for (ind = 0; ind < nijk; ind++)
    {
      k = ind/nij;
      j = (ind-k*nij)/ni;
      i = ind-j*ni-k*nij;
      if (j == 0) yt[ind] = yo;
      else yt[ind] = yo + hj + hj * ((-1. + pow(rj, (j-1))) / (-1. + rj));
    }
  }
  else if (doubleLeftj == 0 && doubleRightj == 1)
  {
    #pragma omp parallel for default(shared) private(k,j,i,ind)
    for (ind = 0; ind < nijk; ind++)
    {
      k = ind/nij;
      j = (ind-k*nij)/ni;
      i = ind-j*ni-k*nij;
      if (j == nj-1) yt[ind] = yo + hj*pow(rj,nj-3) + hj*((-1.+pow(rj,nj-2))/(-1.+rj));
      else yt[ind] = yo + hj*((-1.+pow(rj,j)) / (-1.+rj));
    }
  }
  else // all double 
  {
     #pragma omp parallel for default(shared) private(k,j,i,ind)
    for (ind = 0; ind < nijk; ind++)
    {
      k = ind/nij;
      j = (ind-k*nij)/ni;
      i = ind-j*ni-k*nij;
      if (j == 0) yt[ind] = yo;
      else if (j == nj-1) yt[ind] = yo + +hj + hj*pow(rj,nj-4) + hj*((-1.+pow(rj,nj-3))/(-1.+rj));
      else yt[ind] = yo + hj + hj*((-1.+pow(rj,(j-1))) / (-1.+rj));
    } 
  }

  if (K_FUNC::fEqual(rk, 1.0) == true)
  {
    #pragma omp parallel for default(shared) private(k,j,i,ind)
    for (ind = 0; ind < nijk; ind++)
    {
      k = ind/nij;
      j = (ind-k*nij)/ni;
      i = ind-j*ni-k*nij;
      zt[ind] = zo + hk * k; 
    }   
  }   
  else if (doubleLeftk == 0 && doubleRightk == 0)
  {
    #pragma omp parallel for default(shared) private(k,j,i,ind)
    for (ind = 0; ind < nijk; ind++)
    {
      k = ind/nij;
      j = (ind-k*nij)/ni;
      i = ind-j*ni-k*nij;
      zt[ind] = zo + hk * ( ( -1. + pow(rk, k) ) / (-1. + rk) ); 
    }
  }
  else if (doubleLeftk == 1 && doubleRightk == 0)
  {
    #pragma omp parallel for default(shared) private(k,j,i,ind)
    for (ind = 0; ind < nijk; ind++)
    {
      k = ind/nij;
      j = (ind-k*nij)/ni;
      i = ind-j*ni-k*nij;
      if (k == 0) zt[ind] = zo;
      else zt[ind] = zo + hk + hk * ((-1. + pow(rk, (k-1))) / (-1. + rk));
    }
  }
  else if (doubleLeftk == 0 && doubleRightk == 1)
  {
    #pragma omp parallel for default(shared) private(k,j,i,ind)
    for (ind = 0; ind < nijk; ind++)
    {
      k = ind/nij;
      j = (ind-k*nij)/ni;
      i = ind-j*ni-k*nij;
      if (k == nk-1) zt[ind] = zo + hk*pow(rk,nk-3) + hk*((-1.+pow(rk,nk-2))/(-1.+rk));
      else zt[ind] = zo + hk*((-1.+pow(rk,k)) / (-1.+rk));
    }
  }
  else // all double 
  {
     #pragma omp parallel for default(shared) private(k,j,i,ind)
    for (ind = 0; ind < nijk; ind++)
    {
      k = ind/nij;
      j = (ind-k*nij)/ni;
      i = ind-j*ni-k*nij;
      if (k == 0) zt[ind] = zo;
      else if (k == nk-1) zt[ind] = zo + hk + hk*pow(rk,nk-4) + hk*((-1.+pow(rk,nk-3))/(-1.+rk));
      else zt[ind] = zo + hk + hk*((-1.+pow(rk,(k-1))) / (-1.+rk));
    } 
  }
  
  // Return array
  RELEASESHAREDS(tpl, f2);
  return tpl;
}

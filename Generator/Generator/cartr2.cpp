/*    
    Copyright 2013-2022 Onera.

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
   On cherche le nombre de points pour créer le maillage, on rectifie les facteurs d'expansion pour respecter la condition de fin de grille*/
// ============================================================================
float f(float Xo, float Xf, float H, float N, float R)
{
  return pow(R,N) - R * ( fabs(Xf-Xo) / H ) + ( fabs(Xf-Xo) / H ) - 1.;
}

float fderiv(float Xo, float Xf, float H, float N, float R)
{
  return N * pow(R,N-1) - ( fabs(Xf-Xo) / H );
}

float NewtonApproche(float Xo, float Xf, float H, int N, float Rini)
{
  E_Float res=Rini;
  //printf("%f \n ",res); fflush(stdout);
  E_Int Nit=0;
  while (fabs(f(Xo,Xf,H,N,res)) >= 0.000001 && Nit<30)
  {
    res -= f(Xo,Xf,H,N,res) / fderiv(Xo,Xf,H,N,res);
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
  E_Int skeleton = 0;
  E_Int api = 1;
  if (!PYPARSETUPLE(args, 
                    "(ddd)(ddd)(ddd)(ddd)ll", "(ddd)(ddd)(ddd)(ddd)ii", 
                    "(fff)(fff)(fff)(fff)ll", "(fff)(fff)(fff)(fff)ii",
                    &xo, &yo, &zo, &hi, &hj, &hk, &riinput, &rjinput, &rkinput, &xf, &yf, &zf, &api))
  {
    return NULL;
  }

  if ( ( (fabs(xf-xo) /hi ) * (riinput -1) + 1 < 0) || ( (fabs(yf-yo) /hj ) * (rjinput -1) +1 < 0 ) ||( (fabs(zf-zo) /hk ) * (rkinput -1) +1 < 0 ) )
  {
    PyErr_SetString(PyExc_ValueError, 
                    "Can not generate mesh.\n Condition not met: (Xf-Xo) /H ) * (R-1) +1 < 0.");
    return NULL;
  }

  if ((K_FUNC::fEqual(riinput, 1.0) == true) && (K_FUNC::fEqual( fmod((xf-xo), hi),0) == false) )
  {
    //printf("Remainder of %f / %f is %f\n", xf-xo, hi, fmod((xf-xo),hi));
    riinput += 0.001;
    //printf("Warning: condition on r not met.\n ");
    //printf("Warning: ri set to %f\n", riinput);
  }
  
  if ((K_FUNC::fEqual(rjinput,1.0) == true) && (K_FUNC::fEqual( fmod(yf-yo,hj), 0)==false) )
  {
    rjinput += 0.001;
    //printf("Warning: condition on r not met.\n");
    //printf("Warning: rj set to %f\n", rjinput);
  }

  if ((K_FUNC::fEqual(rkinput,1.0) == true) && (K_FUNC::fEqual(fmod(zf-zo,hk), 0)==false) )
  {
    rkinput += 0.001;
    //printf("Warning: condition on r not met.\n");
    //printf("Warning: rk set to %f.\n" , rkinput);
  }
  E_Float niapp, njapp, nkapp;
  if (K_FUNC::fEqual(riinput, 1.0) == true)
  {
    niapp = fabs(xf - xo) / hi;
  }
  else
  {
    niapp = log( (fabs(xf - xo) / hi) * (riinput - 1) +1) / log(riinput);
  } 

  if (K_FUNC::fEqual(rjinput, 1.) == true)
  {
    njapp = fabs(yf - yo) / hj;
  }
  else
  {
    njapp = log( (fabs(yf - yo) / hj) * (rjinput - 1) +1) / log(rjinput);
  }
  
  if (K_FUNC::fEqual(rkinput, 1.0) == true)
  {
    nkapp = fabs(zf - zo) / hk;
  }
  else
  {  
    nkapp = log( (fabs(zf - zo) / hk) * (rkinput - 1) +1) / log(rkinput);
  }

  if (K_FUNC::fEqual(xo,xf)==true)
  {
    niapp=0;
  }
  if (K_FUNC::fEqual(yo,yf)==true)
  {
    njapp=0;
  }
  if (K_FUNC::fEqual(zo,zf)==true)
  {
    nkapp=0;
  }
  E_Int ni = floor(niapp) + 1;
  E_Int nj = floor(njapp) + 1;
  E_Int nk = floor(nkapp) + 1;
  
  // printf("ni partie entiere = %i \n ",ni) ; fflush(stdout);
  // printf("nj partie entiere= %i \n ",nj) ; fflush(stdout);
  // printf("nk partie entiere= %i \n ",nk) ; fflush(stdout);

  E_Float ri = NewtonApproche(xo,xf,hi,ni,riinput);
  E_Float rj = NewtonApproche(yo,yf,hj,nj,rjinput);
  E_Float rk = NewtonApproche(zo,zf,hk,nk,rkinput);
  
  //printf("Info: final ri=%f\n ",ri) ; fflush(stdout);
  //printf("Info: final rj=%f\n ",rj) ; fflush(stdout);
  //printf("Info: final rk=%f\n ",rk) ; fflush(stdout);
  
  if (fEqual(ri, 1.) == false)
  {
    ni += 1;
  }
  
  if (fEqual(rj, 1.) == false)
  {
    nj += 1;
  } 

  if (fEqual(rk, 1.) == false)
  {
    nk += 1;
  } 

  //printf("Info: number of mesh cell = %i (i), %i (j), %i (k)\n ",ni,nj,nk); fflush(stdout);
  if (skeleton == 1) return Py_BuildValue("lll", ni, nj, nk);

  E_Int i, j, k, ind;
  // Create cartesian mesh
  PyObject* tpl;
  tpl = K_ARRAY::buildArray2(3, "x,y,z", ni, nj, nk, api);
  
  K_FLD::FldArrayF* f; K_FLD::FldArrayI* c;
  char* varString; char* eltType;
  K_ARRAY::getFromArray2(tpl, varString, f, ni, nj, nk, c, eltType);

  E_Int nij = ni*nj;
  E_Int nijk = ni*nj*nk;

  E_Float* xt = f->begin(1);
  E_Float* yt = f->begin(2);
  E_Float* zt = f->begin(3);

  if (K_FUNC::fEqual(ri, 1.0) == false)
  {
    #pragma omp parallel for default(shared) private(k,j,i,ind)
    for (ind = 0; ind < nijk; ind++)
    {
      k = ind/nij;
      j = (ind-k*nij)/ni;
      i = ind-j*ni-k*nij;
      if (xf<xo)
      {
        xt[ind] = xo - hi * ( ( -1. + pow(ri ,i) ) / (-1. + ri) );
      }
      else
      {
        xt[ind] = xo + hi * ( ( -1. + pow(ri ,i) ) / (-1. + ri) );
      } 
    }
  }
  else
  {
    #pragma omp parallel for default(shared) private(k,j,i,ind)
    for (ind = 0; ind < nijk; ind++)
    {
      k = ind/nij;
      j = (ind-k*nij)/ni;
      i = ind-j*ni-k*nij;
      if (xf<xo)
      {
        xt[ind] = xo - hi * i;
      } 
      else
      {
        xt[ind] = xo + hi * i;
      } 
    }    
  }   
    
  if (fEqual(rj, 1.0) == false)
  {
    #pragma omp parallel for default(shared) private(k,j,i,ind)
    for (ind = 0; ind < nijk; ind++)
    {
      k = ind/nij;
      j = (ind-k*nij)/ni;
      i = ind-j*ni-k*nij;
      if (yf<yo)
      {
        yt[ind] = yo - hj * ( ( -1. + pow(rj, j) ) / (-1. + rj) );
      } 
      else
      {
        yt[ind] = yo + hj * ( ( -1. + pow(rj, j) ) / (-1. + rj) );
      }
    }
  }
  else
  {
    #pragma omp parallel for default(shared) private(k,j,i,ind)
    for (ind = 0; ind < nijk; ind++)
    {
      k = ind/nij;
      j = (ind-k*nij)/ni;
      i = ind-j*ni-k*nij;
      if (yf<yo)
      {
        yt[ind] = yo - hj * j;
      }
      else
      {
        yt[ind] = yo + hj * j;
      }      
    } 
  } 

  if (fEqual(rk, 1.0) == false)
  { 
    #pragma omp parallel for default(shared) private(k,j,i,ind)
    for (ind = 0; ind < nijk; ind++)
    {
      k = ind/nij;
      j = (ind-k*nij)/ni;
      i = ind-j*ni-k*nij;
      if (zf<zo)
      {
        zt[ind] = zo - hk * ( ( -1. + pow(rk ,k) ) / (-1. + rk) );
      }
      else
      {
        zt[ind] = zo + hk * ( ( -1. + pow(rk ,k) ) / (-1. + rk) );
      }  
    }    
  }
  else

  {
    #pragma omp parallel for default(shared) private(k,j,i,ind)
    for (ind = 0; ind < nijk; ind++)
    {
      k = ind/nij;
      j = (ind-k*nij)/ni;
      i = ind-j*ni-k*nij;
      if (zf < zo)
      {
        zt[ind] = zo - hk * k;
      } 
      else
      {
        zt[ind] = zo + hk * k;
      }
    }   
  }

  // Return array
  RELEASESHAREDS(tpl, f);
  return tpl;
}

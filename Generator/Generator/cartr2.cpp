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
  return pow(R,N) - R * ( (Xf-Xo) / H ) + ( (Xf-Xo) / H ) -1 ;
}

float fderiv(float Xo, float Xf, float H, float N, float R)
{
  return N * pow(R,N-1) - ( (Xf-Xo) / H ) ;
}



float NewtonApproche(float Xo, float Xf, float H, int N, float Rini)
{
  E_Float res=Rini;
  // printf("%f \n ",res) ; fflush(stdout);
  E_Int Nit=0;
  while( fabs(f(Xo,Xf,H,N,res)) >= 0.00001 && Nit<20)
  {
    res -=  f(Xo,Xf,H,N,res) / fderiv(Xo,Xf,H,N,res) ;
    // printf("%f \n ",res) ; fflush(stdout);
    // printf("%f \n ", f(Xo,Xf,H,N,res) ) ; fflush(stdout);
    Nit+=1;
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
  E_Int api = 1;
  if (!PYPARSETUPLE(args, 
                    "(ddd)(ddd)(ddd)(ddd)l", "(ddd)(ddd)(ddd)(ddd)i", 
                    "(fff)(fff)(fff)(fff)l", "(fff)(fff)(fff)(fff)i",
                    &xo, &yo, &zo, &hi, &hj, &hk, &riinput, &rjinput, &rkinput, &xf, &yf, &zf, &api))
  {
    return NULL;
  }

  if ( ( (xf-xo) /hi ) * (riinput -1) +1 < 0 || ( (yf-yo) /hj ) * (rjinput -1) +1 < 0 ||( (zf-zo) /hk ) * (rkinput -1) +1 < 0 )
  {
    PyErr_SetString(PyExc_ValueError, 
                    "Maillage impossible, veuillez changer les paramètres d'entrée \n Conditon à respecter: (Xf-Xo) /H ) * (R -1) +1 < 0");
    return NULL;
  }

  if ( (riinput == 1.0) && (fmod(xf-xo,hi)  !=0)  )
  {
    riinput+=0.1;
    printf("Condition sur r impossible à satisfaire, la valeur va être changé\n ");
    printf("nouveau ri = %f\n", riinput);
  }
  if ( (rjinput == 1.0) && ( fmod(yf-yo,hj) !=0) )
  {
    rjinput+=0.1;
    printf("Condition sur r impossible à satisfaire, la valeur va être changé\n");
    printf("nouveau rj = %f\n", rjinput);
  }
  if ( (rkinput == 1.0) && ( fmod(zf-zo,hk) !=0) )
  {
    rkinput+=0.1;
    printf("Condition sur r impossible à satisfaire, la valeur va être changé\n");
    printf("nouveau rk = %f\n" , rkinput);
  }
  E_Float niapp, njapp, nkapp;
  if (riinput == 1.0 )
  {
    niapp = (xf - xo) / hi ;
  }
  else
  {
    niapp = log( ((xf - xo) / hi) * (riinput - 1) +1) / log(riinput) ;
  } 

  if (rjinput == 1.0)
  {
    njapp = (yf - yo) / hj ;
  }
  else
  {
    njapp = log( ((yf - yo) / hj) * (rjinput - 1) +1) / log(rjinput) ;
  }
  
  if (rkinput == 1.0)
  {
    nkapp = (zf - zo) / hk ;
  }
  else
  {  
    nkapp = log( ((zf - zo) / hk) * (rkinput - 1) +1) / log(rkinput) ;
  }

  // printf(" ni approche = %f \n ",niapp) ; fflush(stdout);
  // printf(" nj approche = %f \n ",njapp) ; fflush(stdout);
  // printf(" nk approche = %f \n ",nkapp) ; fflush(stdout);

  
  E_Int ni = floor(niapp) + 1 ;
  E_Int nj = floor(njapp) + 1 ;
  E_Int nk = floor(nkapp) + 1 ;
  
  // printf("ni partie entiere = %i \n ",ni) ; fflush(stdout);
  // printf("nj partie entiere= %i \n ",nj) ; fflush(stdout);
  // printf("nk partie entiere= %i \n ",nk) ; fflush(stdout);

  E_Float ri = NewtonApproche(xo,xf,hi,ni,riinput);
  E_Float rj = NewtonApproche(yo,yf,hj,nj,rjinput);
  E_Float rk = NewtonApproche(zo,zf,hk,nk,rkinput);
  
  printf("ri Final =%f \n ",ri) ; fflush(stdout);
  printf("rj Final =%f \n ",rj) ; fflush(stdout);
  printf("rk Final =%f \n ",rk) ; fflush(stdout);
  
  if (ri != 1.)
  {
    ni+=1 ;
  } 
  
  if (rj != 1.)
  {
    nj+=1 ;
  } 

  if (rk != 1.)
  {
    nk+=1 ;
  } 

  printf("nombre mailles en i = %i \n ",ni) ; fflush(stdout);
  printf("nombre mailles en j = %i \n ",nj) ; fflush(stdout);
  printf("nombre mailles en k = %i \n ",nk) ; fflush(stdout);

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
  

  if (ri != 1.0) 
  {
    #pragma omp parallel for default(shared) private(k,j,i,ind)
    for (ind = 0; ind < nijk; ind++)
    {
      k = ind/nij;
      j = (ind-k*nij)/ni;
      i = ind-j*ni-k*nij;
      xt[ind] = xo + hi * ( ( -1 + pow(ri ,i) ) / (-1 + ri) );
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
      xt[ind] = xo + hi * i ;
    } 
        
    }   
    
  if (rj != 1.0)
    {
      #pragma omp parallel for default(shared) private(k,j,i,ind)
      for (ind = 0; ind < nijk; ind++)
      {
        k = ind/nij;
        j = (ind-k*nij)/ni;
        i = ind-j*ni-k*nij;
        yt[ind] = yo + hj * ( ( -1 + pow(rj ,j) ) / (-1 + rj) );
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
        yt[ind] = yo + hj * j ;
      }
        
    } 

    if (rk != 1.0)
    { 
      #pragma omp parallel for default(shared) private(k,j,i,ind)
      for (ind = 0; ind < nijk; ind++)
      {
        k = ind/nij;
        j = (ind-k*nij)/ni;
        i = ind-j*ni-k*nij;
        zt[ind] = zo + hk * ( ( -1 + pow(rk ,k) ) / (-1 + rk) );
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
        zt[ind] = zo + hk * k ;
      }
        
    }

  // Return array
  RELEASESHAREDS(tpl, f);
  return tpl;
}

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

# include "CompGeom/compGeom.h"
# include "parallel.h"

using namespace K_FLD;
extern "C"
{
  void k6onedmap_(const E_Int& ni,
                  const E_Float* x, const E_Float* y, const E_Float* z,
                  const E_Int& no,
                  const E_Float* distrib,
                  E_Float* xo, E_Float* yo, E_Float* zo,
                  E_Float* s, E_Float* dx, E_Float* dy, E_Float* dz);
}

//===========================================================================
/* spline a partir des points de controle 
   IN: im: nbre de pts de controles
   IN: xt, yt, zt: coordonnees des pts de controle
   IN: ordern: ordre de la spline
   IN: N: nbre de pts voulus sur la spline
   OUT: PF: coord. de la spline.
*/
//===========================================================================
void K_COMPGEOM::spline(E_Int im, E_Int ordern, E_Int N, 
                        E_Float* xt, E_Float* yt, E_Float* zt, 
                        FldArrayF& PF)
{
  E_Int Nbpoints = N;

  E_Int n = im-1; 
  E_Int c = ordern; 
  E_Int m = n+c; 

  E_Float* PF1 = PF.begin(1);
  E_Float* PF2 = PF.begin(2);
  E_Float* PF3 = PF.begin(3);

  if (im == 2)
  {
    E_Float delta = 1./(N-1.);
    E_Float dx12 = delta * (xt[1]-xt[0]);
    E_Float dy12 = delta * (yt[1]-yt[0]);
    E_Float dz12 = delta * (zt[1]-zt[0]);
    
#pragma omp parallel default(shared) if (N > __MIN_SIZE_MEAN__)
    {
#pragma omp for
      for (E_Int i = 0; i < N; i++)
      {
        PF1[i] = xt[0] + i*dx12;
        PF2[i] = yt[0] + i*dy12;
        PF3[i] = zt[0] + i*dz12;
      }
    } 
  }

  else
  {
    // Vecteur de noeuds
    FldArrayF x(m+1); x.setAllValuesAtNull();
    
    for (E_Int i = 1; i <= m;i++)
    {
      if( i < c || i > n+1 ) x[i] = x[i-1];      
      else x[i] = x[i-1]+1.;
    }    

    E_Float pas = x[m]/(Nbpoints-1);

#pragma omp parallel default(shared) if (Nbpoints > __MIN_SIZE_MEAN__)
    {
      E_Float tmp1;
      E_Float tmp2;
      E_Float tmp3;
      E_Float t = 0.;
#pragma omp for
      for (E_Int i1 = 0; i1 < Nbpoints; i1++)
      {
        t = i1*pas;
        
        tmp1 = 0.;
        tmp2 = 0.;
        tmp3 = 0.;
        
        evalSpline(t, x, n, c, 
                   xt, yt, zt, tmp1, tmp2, tmp3);
        
        PF1[i1] = tmp1;
        PF2[i1] = tmp2;
        PF3[i1] = tmp3;
        
      }
    }
  }
}

//=============================================================================
void K_COMPGEOM::spline2D(E_Int im, E_Int jm, E_Int ordern,
                          E_Int N, E_Int orderm, E_Int M, 
                          E_Float* xt, E_Float* yt, E_Float* zt, 
                          FldArrayF& PF)
{
  E_Int n = im-1; 
  E_Int m = jm-1;
  E_Int p = n+ordern; 
  E_Int q = m+orderm;

  PF.setAllValuesAtNull();

  FldArrayF x(p+1);
  x.setAllValuesAtNull();
  
  for (E_Int i = 1; i <= p;i++)
  {
    if( i < ordern || i > n+1 )
      x[i] = x[i-1];
    
    else x[i] = x[i-1]+1.;
  }
  
  FldArrayF y(q+1);
  y.setAllValuesAtNull();
  
  for (E_Int i = 1; i <= q;i++)
  {
    if( i < orderm || i > m+1 )
      y[i] = y[i-1];
    
    else y[i] = y[i-1]+1.;    
  }

  E_Float pasn = x[p]/(N-1);
  E_Float pasm = y[q]/(M-1);

 
  E_Float* PF1 = PF.begin(1);
  E_Float* PF2 = PF.begin(2);
  E_Float* PF3 = PF.begin(3);

#pragma omp parallel default(shared) if (N > __MIN_SIZE_MEAN__)
  {
  E_Int ind, i, j;

  E_Float tmp1 = 0.;
  E_Float tmp2 = 0.;
  E_Float tmp3 = 0.;  
  E_Float u = 0.;
  E_Float v = 0.;

  FldArrayF P(m+1,3);
  E_Float* P1 = P.begin(1);
  E_Float* P2 = P.begin(2);
  E_Float* P3 = P.begin(3);

#pragma omp for
    for (i = 0; i < N; i++)
    {
      u= i*pasn;

      P.setAllValuesAtNull();
      for (j = 0; j < m+1; j++)
      {
        evalSpline(u, x, n, ordern, 
                   xt+im*j, yt+im*j, zt+im*j, tmp1, tmp2, tmp3);
        P1[j] = P1[j] + tmp1;
        P2[j] = P2[j] + tmp2;
        P3[j] = P3[j] + tmp3;
      }
      v = 0.;
      
      for (j = 0; j < M; j++)
      {
        v = j*pasm;
        
        evalSpline(v, y, m, orderm, P.begin(1), P.begin(2), P.begin(3),
                   tmp1, tmp2, tmp3);
        
        ind = i + j * N;
        PF1[ind] = PF1[ind] + tmp1;
        PF2[ind] = PF2[ind] + tmp2;
        PF3[ind] = PF3[ind] + tmp3;
      }
    }
  }
}

//===========================================================================
/* Creation et evaluation de la spline
   IN : t : parametre ou la spline est evaluee
   IN : x : vecteurs des noeuds parametres
   IN : n : nombre de pts de controle - 1
   IN : c : ordre de la spline
   OUT : xo, yo, zo : coord. du point de la spline correspondant a t
*/
//===========================================================================
void K_COMPGEOM::evalSpline(E_Float t, FldArrayF& x, E_Int n, E_Int c, 
                            E_Float* xt, E_Float* yt, E_Float* zt,
                            E_Float& xo, E_Float& yo, E_Float& zo)
{
  
  E_Int m = n+c;
  E_Int d = c-1;
  E_Float alpha,alpha1;
  E_Int ind, i;

  if ( K_FUNC::fEqualZero(t - x[0]) == true )
  {
    xo = xt[0]; yo = yt[0]; zo = zt[0];
    return;
  }
  if ( K_FUNC::fEqualZero(t - x[m]) == true) 
  {
    xo = xt[n]; yo = yt[n]; zo = zt[n];
    return;
  }
  
  // trouver k tq t in [xk,xk+1[ 
  E_Int k = -1;
  for (i = 0; i <= m-1; i++)
  {     
    if ( t >= x[i] && t < x[i+1] )
      k = i;
  }
  
  FldArrayF xPk(d+1,d+1);
  FldArrayF yPk(d+1,d+1);
  FldArrayF zPk(d+1,d+1);
  
  E_Float* xPk1 = xPk.begin(1);
  E_Float* yPk1 = yPk.begin(1);
  E_Float* zPk1 = zPk.begin(1);

  for (i = k-d; i <= k; i++)
  {
    ind = i-k+d;
    xPk1[ind] = xt[i];
    yPk1[ind] = yt[i];
    zPk1[ind] = zt[i];
  }

  for (E_Int j = 1; j <= d; j++)
  {
    E_Float* xPkj = xPk.begin(j);
    E_Float* yPkj = yPk.begin(j);
    E_Float* zPkj = zPk.begin(j);

    E_Float* xPkj1 = xPk.begin(j+1);
    E_Float* yPkj1 = yPk.begin(j+1);
    E_Float* zPkj1 = zPk.begin(j+1);
    
    for (i = k; i >= k-d+j; i--)
    {
      alpha = (t - x[i])/(x[i+c-j] - x[i]);
      alpha1 = 1.-alpha;
      ind = i-k+d;
      xPkj1[ind] = alpha1 * xPkj[ind-1] + alpha*xPkj[ind];
      yPkj1[ind] = alpha1 * yPkj[ind-1] + alpha*yPkj[ind];
      zPkj1[ind] = alpha1 * zPkj[ind-1] + alpha*zPkj[ind];
    }
  }
 
  E_Float* xPkd1 = xPk.begin(d+1);
  E_Float* yPkd1 = yPk.begin(d+1);
  E_Float* zPkd1 = zPk.begin(d+1);
  xo = xPkd1[k-k+d];
  yo = yPkd1[k-k+d];
  zo = zPkd1[k-k+d];
}

//=============================================================================
// Calcul Ni,p(t)
// IN : p : ordre de la spline
//=============================================================================
E_Float N(E_Int n, E_Int p, E_Float t, E_Float* x)
{
  E_Int m = n+p;
  E_Int d = p-1;
  E_Float alpha, alpha1, x0;
  E_Int ind, i;
  if ( K_FUNC::fEqualZero(t - x[0]) == true )
  {
    return 0.;
  }
  if ( K_FUNC::fEqualZero(t - x[m]) == true) 
  {
    return 0.;
  }
  
  // trouver k tq t in [xk,xk+1[ 
  E_Int k = -1;
  for (i = 0; i <= m-1; i++)
  {     
    if ( t >= x[i] && t < x[i+1] )
      k = i;
  }
  
  FldArrayF xPk(d+1,d+1);
  E_Float* xPk1 = xPk.begin(1);

  for (i = k-d; i <= k; i++)
  {
    ind = i-k+d;
    xPk1[ind] = 1.;
  }

  for (E_Int j = 1; j <= d; j++)
  {
    E_Float* xPkj = xPk.begin(j);
    E_Float* xPkj1 = xPk.begin(j+1);
    
    for (i = k; i >= k-d+j; i--)
    {
      alpha = (t - x[i])/(x[i+p-j] - x[i]);
      alpha1 = 1. - alpha;
      ind = i-k+d;
      xPkj1[ind] = alpha1 * xPkj[ind-1] + alpha*xPkj[ind];
    }
  }
 
  E_Float* xPkd1 = xPk.begin(d+1);
  x0 = xPkd1[k-k+d];
  return x0;
}

//===========================================================================
/* 
   Calcul la courbe spline a partir des points de controle
   avec une discretisation reguliere.
   IN : n : nb de pts de controle
   IN : ordern : ordre de la spline
   IN : N : nombre de pts sur la courbe finale(si density>0, N est ignore)
   IN : density : npts par unite de longueur de la bezier resultante
   IN : xt, yt, zt : pts de controle de la courbe
   OUT : coord : tableau des coordonnees des pts de la courbe de Bezier 
*/
//===========================================================================
void K_COMPGEOM::regularSpline(E_Int n, E_Int ordern, E_Int N, E_Float density,
                               E_Float* xt, E_Float* yt, E_Float* zt,
                               FldArrayF& coord)
{
  // approx de la longueur de la spline par la longueur des pts de controle
  E_Float len = 0.;
#pragma omp parallel default(shared) if (n > __MIN_SIZE_MEAN__)
  {
  E_Int i1;
  E_Float dx, dy, dz;
#pragma omp for reduction(+:len)
    for (E_Int i = 1; i < n; i++)
    {
      i1 = i-1;
      dx = xt[i] - xt[i1];
      dy = yt[i] - yt[i1];
      dz = zt[i] - zt[i1];
      len = len + sqrt(dx*dx+dy*dy+dz*dz);
    }
  }  
  E_Int npts;
  if (density > 0) npts = E_Int(len*density)+1;
  else npts = N;
  E_Int npts0 = 10*npts;
  FldArrayF coord0(npts0, 3);
  coord0.setAllValuesAtNull();
  E_Float* coordx0 = coord0.begin(1);
  E_Float* coordy0 = coord0.begin(2);
  E_Float* coordz0 = coord0.begin(3);
  spline(n, ordern, npts0, xt, yt, zt, coord0);

  // vraie longueur de la spline
  len = 0.;
#pragma omp parallel default(shared) if (npts0 > __MIN_SIZE_MEAN__)
  {
    E_Int i1;
    E_Float dx, dy, dz;
#pragma omp for reduction(+:len)
    for (E_Int i = 1; i < npts0; i++)
    {
      i1 = i-1;
      dx = coordx0[i] - coordx0[i1];
      dy = coordy0[i] - coordy0[i1];
      dz = coordz0[i] - coordz0[i1];
      len = len + sqrt(dx*dx+dy*dy+dz*dz);
    }
  }
  if (density > 0) npts = E_Int(len*density)+1;
  else npts = N;
 
  // spline reguliere
  coord.malloc(npts, 3);
  E_Float* coordx = coord.begin(1);
  E_Float* coordy = coord.begin(2);
  E_Float* coordz = coord.begin(3);

  FldArrayF sp(npts0);
  FldArrayF dxp(npts0);
  FldArrayF dyp(npts0);
  FldArrayF dzp(npts0);
  FldArrayF fd(npts);
  E_Float* fdx = fd.begin(1);
  E_Float delta = 1./(npts-1.);
#pragma omp parallel for default(shared) if (npts > __MIN_SIZE_MEAN__)
  for (E_Int i = 0; i < npts; i++) fdx[i] = delta*i;

  k6onedmap_(npts0, coordx0, coordy0, coordz0,
             npts, fd.begin(1),
             coordx, coordy, coordz,
             sp.begin(), dxp.begin(), dyp.begin(), dzp.begin());
}

//=============================================================================
/* 
   Calcul la surface spline a partir des points de controle
   avec une discretisation reguliere.
   IN: n,m: nb de pts de controle
   IN: orderm, orderm: ordre de la spline
   IN: N,M: nbre de pts sur la surface (si density<0)
   IN: density: npts par unite de longueur de la bezier resultante
   IN: xt, yt, zt: pts de controle de la courbe
   OUT: coord: tableau des coordonnees des pts de la courbe de Bezier 
*/
//=============================================================================
void K_COMPGEOM::regularSpline2D(E_Int n, E_Int m, E_Int ordern, E_Int N, 
                                 E_Int orderm, E_Int M, 
                                 E_Float density, 
                                 E_Float* xt, E_Float* yt, E_Float* zt,
                                 FldArrayF& coord, E_Int& niout, E_Int& njout)
{
  // approx de la longueur de la spline par la longueur des pts de controle
  E_Float leni = 0.;
  E_Float lenj = 0.;

#pragma omp parallel default(shared) if (m*n > __MIN_SIZE_MEAN__)
  {
    E_Int ind, ind1, j1;
    E_Float len=0.;
    E_Float dx, dy, dz;
    E_Float leni_private=0., lenj_private=0.;
#pragma omp for nowait
    for (E_Int j = 0; j < m; j++)
    { 
      len = 0.;
      for (E_Int i = 1; i < n; i++)
      {
        ind = i+j*n; ind1 = i-1+j*n;
        dx = xt[ind] - xt[ind1];
        dy = yt[ind] - yt[ind1];
        dz = zt[ind] - zt[ind1];
        len = len + sqrt(dx*dx+dy*dy+dz*dz);
      }
      leni_private = K_FUNC::E_max(leni_private, len);
    }

#pragma omp for nowait
    for (E_Int i = 0; i < n; i++)
    {
      len = 0.;
      for (E_Int j = 1; j < m; j++)
      {
        j1 = j-1;
        ind = i+j*n; ind1 = i+j1*n;
        dx = xt[ind] - xt[ind1];
        dy = yt[ind] - yt[ind1];
        dz = zt[ind] - zt[ind1];
        len = len + sqrt(dx*dx+dy*dy+dz*dz);
      }
      lenj_private = K_FUNC::E_max(lenj_private, len);
    }
#pragma omp critical
    {
      leni = K_FUNC::E_max(leni, leni_private);
      lenj = K_FUNC::E_max(lenj, lenj_private);
    }
  }
  E_Int nptsi, nptsj;
  if (density > 0) 
  { nptsi = E_Int(leni*density)+1; nptsj = E_Int(lenj*density)+1; }
  else { nptsi = N; nptsj = M; }
  E_Int nptsi0 = 10*nptsi;
  E_Int nptsj0 = 10*nptsj;

  FldArrayF coord0(nptsi0*nptsj0, 3); coord0.setAllValuesAtNull();
  E_Float* coordx0 = coord0.begin(1);
  E_Float* coordy0 = coord0.begin(2);
  E_Float* coordz0 = coord0.begin(3);
  spline2D(n, m, ordern, nptsi0, orderm, nptsj0, xt, yt, zt, coord0);

  // Remaille suivant i
  // vraie longueur de la spline suivant i
  E_Float len = 0.;
#pragma omp parallel default(shared) if (nptsi0 > __MIN_SIZE_MEAN__)
  {
    E_Float dx, dy, dz;
    E_Int ind, ind0;    
#pragma omp for reduction(+:len)
    for (E_Int i = 1; i < nptsi0; i++)
    {
      ind = i; ind0 = i-1;
      dx = coordx0[ind] - coordx0[ind0];
      dy = coordy0[ind] - coordy0[ind0];
      dz = coordz0[ind] - coordz0[ind0];
      len = len + sqrt(dx*dx+dy*dy+dz*dz);
    }
  }

  if (density > 0) nptsi = E_Int(len*density)+1;
  else nptsi = N;

  FldArrayF sp(nptsi0);
  FldArrayF dxp(nptsi0);
  FldArrayF dyp(nptsi0);
  FldArrayF dzp(nptsi0);
  FldArrayF fd(nptsi);
  E_Float* fdx = fd.begin(1);
  E_Float delta = 1./(nptsi-1.);
#pragma omp parallel for default(shared) if (nptsi > __MIN_SIZE_MEAN__)
  for (E_Int i = 0; i < nptsi; i++) fdx[i] = delta*i;
  FldArrayF coord1(nptsi*nptsj0, 3); coord1.setAllValuesAtNull();
  E_Float* coordx1 = coord1.begin(1);
  E_Float* coordy1 = coord1.begin(2);
  E_Float* coordz1 = coord1.begin(3);

  for (E_Int j = 0; j < nptsj0; j++)
  {
    k6onedmap_(nptsi0, coordx0+j*nptsi0, coordy0+j*nptsi0, coordz0+j*nptsi0,
               nptsi, fd.begin(1),
               coordx1+j*nptsi, coordy1+j*nptsi, coordz1+j*nptsi,
               sp.begin(), dxp.begin(), dyp.begin(), dzp.begin());
  }

  // Remaille en j
  // vraie longueur de la spline suivant j
  len = 0.;
#pragma omp parallel default(shared) if (nptsj0 > __MIN_SIZE_MEAN__)
  {
    E_Float dx, dy, dz;
    E_Int ind, ind0;
#pragma omp for reduction(+:len)
    for (E_Int j = 1; j < nptsj0; j++)
    {
      ind = j*nptsi0; ind0 = (j-1)*nptsi0;
      dx = coordx0[ind] - coordx0[ind0];
      dy = coordy0[ind] - coordy0[ind0];
      dz = coordz0[ind] - coordz0[ind0];
      len = len + sqrt(dx*dx+dy*dy+dz*dz);
    }
  }

  if (density > 0) nptsj = E_Int(len*density)+1;
  else nptsj = M;

  sp.malloc(nptsj0);
  dxp.malloc(nptsj0);
  dyp.malloc(nptsj0);
  dzp.malloc(nptsj0);
  fd.malloc(nptsj);
  fdx = fd.begin(1);
  delta = 1./(nptsj-1.);
#pragma omp parallel for default(shared) if (nptsj > __MIN_SIZE_MEAN__)
  for (E_Int j = 0; j < nptsj; j++) fdx[j] = delta*j;
  coord.malloc(nptsi*nptsj, 3); coord.setAllValuesAtNull();

  E_Float* coordx = coord.begin(1);
  E_Float* coordy = coord.begin(2);
  E_Float* coordz = coord.begin(3);

  FldArrayF t0(nptsj0, 3);
  E_Float* t0x = t0.begin(1);
  E_Float* t0y = t0.begin(2);
  E_Float* t0z = t0.begin(3);
  FldArrayF t(nptsj, 3);
  E_Float* tx = t.begin(1);
  E_Float* ty = t.begin(2);
  E_Float* tz = t.begin(3);

  for (E_Int i = 0; i < nptsi; i++)
  {
    for (E_Int j = 0; j < nptsj0; j++) 
    { t0x[j] = coordx1[i+j*nptsi];
      t0y[j] = coordy1[i+j*nptsi];
      t0z[j] = coordz1[i+j*nptsi]; }
    k6onedmap_(nptsj0, t0x, t0y, t0z,
               nptsj, fd.begin(1),
               tx, ty, tz,
               sp.begin(), dxp.begin(), dyp.begin(), dzp.begin());
    for (E_Int j = 0; j < nptsj; j++) 
    { coordx[i+j*nptsi] = tx[j]; 
      coordy[i+j*nptsi] = ty[j]; 
      coordz[i+j*nptsi] = tz[j]; }
  }
  niout = nptsi; njout = nptsj;
}

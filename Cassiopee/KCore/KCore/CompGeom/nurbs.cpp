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
# include <stdio.h>
# include "parallel.h"
# include "CompGeom/compGeom.h"

using namespace K_FLD;

//===========================================================================
/* nurbs a partir des points de controle */
//===========================================================================
void K_COMPGEOM::nurbs(E_Int im, E_Int ordern, E_Int N, 
                       E_Float* xt, E_Float* yt, E_Float* zt, E_Float* W, 
                       FldArrayF& coord)
{
  E_Int Nbpoints = N;
  E_Int n = im-1; 
  E_Int c = ordern; 
  E_Int m = n+c; 
  E_Int d = c-1;

  E_Float* coordx = coord.begin(1);
  E_Float* coordy = coord.begin(2);
  E_Float* coordz = coord.begin(3);
  
  if (im == 2)
  {
    E_Float delta = 1./(N-1.);
    E_Float dx12 = delta * (xt[1]-xt[0]);
    E_Float dy12 = delta * (yt[1]-yt[0]);
    E_Float dz12 = delta * (zt[1]-zt[0]);
    
    #pragma omp parallel
    {
      #pragma omp for
      for (E_Int i = 0; i < N; i++)
      {
        coordx[i] = xt[0] + i*dx12;
        coordy[i] = yt[0] + i*dy12;
        coordz[i] = zt[0] + i*dz12;
      } 
    }
  }
  else
  {
    FldArrayF x(m+1);
    x.setAllValuesAtNull();
    
    for (E_Int i = 1; i <= m;i++)
    {
      if (i < c || i > n+1) x[i] = x[i-1];
      else x[i] = x[i-1] + 1.;
    }

    E_Float pas = x[m]/(Nbpoints-1);
    
    #pragma omp parallel
    {
      FldArrayF B(m+1,d+1);
      E_Float* Bd = B.begin(d+1);
      E_Float tmp1, tmp2, tmp3, tmp4;
      E_Float t = 0.;
      #pragma omp for
      for (E_Int i1 = 0; i1 < Nbpoints; i1++)
      {
        t = i1*pas;        
        tmp1 = 0.; tmp2 = 0.; tmp3 = 0.; tmp4 = 0.;
        evalNurbs(t, x, n, c, B);
        
        for(E_Int i = 0; i <= n; i++)
        {
          tmp1 = tmp1 + Bd[i] * xt[i]*W[i];
          tmp2 = tmp2 + Bd[i] * yt[i]*W[i];
          tmp3 = tmp3 + Bd[i] * zt[i]*W[i];
          tmp4 = tmp4 + Bd[i] * W[i];
        }
        coordx[i1] = tmp1/tmp4;
        coordy[i1] = tmp2/tmp4;
        coordz[i1] = tmp3/tmp4; 
      }
    }
  }
}

//=============================================================================
// Creation et evaluation de la  surface 2D nurbs
//=============================================================================
void K_COMPGEOM::nurbs2D(E_Int im, E_Int jm, E_Int ordern,
                         E_Int N, E_Int orderm, E_Int M, 
                         E_Float* xt, E_Float* yt, E_Float* zt, E_Float* W,
                         FldArrayF& coord)
{ 
  E_Int n = im-1;
  E_Int m = jm-1; 

  E_Int p = n+ordern;
  E_Int q = m+orderm;
  E_Float* coordx = coord.begin(1);
  E_Float* coordy = coord.begin(2);
  E_Float* coordz = coord.begin(3);
  
  FldArrayF P(m+1,3);
  coord.setAllValuesAtNull();
  
  FldArrayF x(p+1);
  x.setAllValuesAtNull();
  
  for (E_Int i = 1; i <= p; i++)
  {
    if (i < ordern || i > n+1) x[i] = x[i-1];
    else x[i] = x[i-1]+1.;  
  }   
  
  FldArrayF y(q+1);
  y.setAllValuesAtNull();
  
  for (E_Int j = 1; j <= q;j++)
  {
    if (j < orderm || j > m+1) y[j] = y[j-1];
    else y[j] = y[j-1]+1.;    
  }
  
  E_Float pasn = x[p]/(N-1);
  E_Float pasm = y[q]/(M-1);

  #pragma omp parallel 
  {
    FldArrayF Bn(p+1,ordern);
    FldArrayF Bm(q+1,orderm);
    E_Float* Bmm = Bm.begin(orderm);
    E_Float* Bnn = Bn.begin(ordern);   
    E_Float* Bnm = Bn.begin(orderm);   
    E_Float u = 0.;
    E_Float v = 0.;
    E_Float tmpx = 0.;
    E_Float tmpy = 0.; 
    E_Float tmpz = 0.;
    E_Float norm = 0.;
    #pragma omp for
    for (E_Int j = 0; j < M; j++)
    { 
      v = j*pasm;
      
      evalNurbs(v,y,m,orderm,Bm);
      
      u = 0.;
      
      for (E_Int i = 0; i < N; i++)
      {
        u = i*pasn;
        evalNurbs(u, x, n, ordern, Bn);
        
        for (E_Int jjj = 0; jjj <= m; jjj++)
          for (E_Int iii = 0; iii <= n; iii++)
          {
            E_Int ind = iii + jjj * im;
            tmpx = tmpx + Bmm[jjj]*Bnn[iii]*W[ind]*xt[ind];
            tmpy = tmpy + Bmm[jjj]*Bnn[iii]*W[ind]*yt[ind];
            tmpz = tmpz + Bmm[jjj]*Bnn[iii]*W[ind]*zt[ind];
          }
        
        for (E_Int jj = 0; jj <= m; jj++)
          for (E_Int ii = 0; ii <= n; ii++)
            norm = norm + Bnm[ii]*Bmm[jj]*W[ii+im*jj];
        
        coordx[i+N*j] = tmpx/norm;    
        coordy[i+N*j] = tmpy/norm; 
        coordz[i+N*j] = tmpz/norm;
        tmpx = 0.; tmpy = 0.; tmpz = 0.; norm = 0.;
      }
    }
  }
}

//=============================================================================
// Creation et evaluation de la nurbs
//=============================================================================
void K_COMPGEOM::evalNurbs(E_Float t, FldArrayF& x, 
                           E_Int n, E_Int c, 
                           FldArrayF& B)
{
  E_Int m = n+c;
  E_Int d = c-1;
  E_Float C1,C2,C1num,C1den,C2num,C2den;

  B.setAllValuesAtNull();
  E_Float* B1 = B.begin(1);

  for (E_Int j = 0; j < d; j++)
  {  
    if (K_FUNC::fEqualZero(t - x[0], 1.e-6) == true)
      B1[j] = 1.;
  }

  for (E_Int j = d; j < m; j++)
  {
    if ( (t >= x[j] ) && (t < x[j+1]) ) 
      B1[j] = 1.;
  }
  
  for (E_Int j = n; j <= m; j++)
  {  
    if (K_FUNC::fEqualZero(t - x[m], 1.e-6) == true)
      B1[j] = 1.;
  }

  // Calcul des N de facon recursives
  for (E_Int k = 2; k <= d+1; k++)
  {
    E_Float* Bk = B.begin(k);
    E_Float* Bk1 = B.begin(k-1);
    for (E_Int i = 0 ; i <= n ; i++) 
    {
      C1num = t - x[i];
      C1den = x[i+k-1] - x[i];
      C2num = x[i+k] - t;
      C2den = x[i+k] - x[i+1];
      
      if (K_FUNC::fEqualZero(C1den) == true ) C1 = 0;
      else C1 = C1num/C1den;
    
      if (K_FUNC::fEqualZero(C2den) == true ) C2 = 0;
      else C2 = C2num/C2den;    
    
      Bk[i] = C1 * Bk1[i] + C2 * Bk1[i+1];
    }
  }
}

//===========================================================================
/* 
   Calcul la courbe nurbs a partir des points de controle
   avec une discretisation reguliere.
   IN : n : nb de pts de controle
   IN : ordern : ordre de la spline
   IN : N : nombre de pts sur la courbe finale(si density>0, N est ignore)
   IN : density : npts par unite de longueur de la bezier resultante
   IN : xt, yt, zt : pts de controle de la courbe
   OUT : coord : tableau des coordonnees des pts de la courbe de Bezier 
*/
//===========================================================================
void K_COMPGEOM::regularNurbs(E_Int n, E_Int ordern, E_Int N, E_Float density,
                              E_Float* xt, E_Float* yt, E_Float* zt,
                              E_Float* W,
                              FldArrayF& coord)
{
  // approx de la longueur de la nurbs par la longueur des pts de controle
  E_Float len = 0.;
  E_Int nthreads = __NUMTHREADS__;
  E_Float* lp = new E_Float [nthreads];

  #pragma omp parallel
  {
    E_Float dx, dy, dz;
    E_Int i1;
    E_Int it = __CURRENT_THREAD__;
    lp[it] = 0.;
    #pragma omp for
    for (E_Int i = 1; i < n; i++)
    {
      i1 = i-1;
      dx = xt[i] - xt[i1];
      dy = yt[i] - yt[i1];
      dz = zt[i] - zt[i1];
      lp[it] += sqrt(dx*dx+dy*dy+dz*dz);
    }
  }
  for (E_Int it = 0; it < nthreads; it++)
  {
    len += lp[it];
  }

  E_Int npts;
  if (density > 0.) npts = E_Int(len*density)+1;
  else npts = N;
  E_Int npts0 = 10*npts;
  FldArrayF coord0(npts0, 3);
  coord0.setAllValuesAtNull();
  E_Float* coordx0 = coord0.begin(1);
  E_Float* coordy0 = coord0.begin(2);
  E_Float* coordz0 = coord0.begin(3);
  nurbs(n, ordern, npts0, xt, yt, zt, W, coord0);

  // vraie longueur de la nurbs
  len = 0.;
  #pragma omp parallel
  {
    E_Float dx, dy, dz;
    E_Int i1;
    E_Int it = __CURRENT_THREAD__;
    lp[it] = 0.;
    #pragma omp for reduction(+:len)
    for (E_Int i = 1; i < npts0; i++)
    {
      i1 = i-1;
      dx = coordx0[i] - coordx0[i1];
      dy = coordy0[i] - coordy0[i1];
      dz = coordz0[i] - coordz0[i1];
      lp[it] += sqrt(dx*dx+dy*dy+dz*dz);
    }
  }
  for (E_Int it = 0; it < nthreads; it++)
  {
    len += lp[it];
  }
  if (density > 0.) npts = E_Int(len*density)+1;
  else npts = N;
 
  // nurbs reguliere
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
  #pragma omp parallel for
  for (E_Int i = 0; i < npts; i++) fdx[i] = delta*i;

  K_COMPGEOM::onedmap(npts0, coordx0, coordy0, coordz0,
		      npts, fd.begin(1),
		      coordx, coordy, coordz,
		      sp.begin(), dxp.begin(), dyp.begin(), dzp.begin());
  delete [] lp;
}

//=============================================================================
/* 
   Calcul la surface nurbs a partir des points de controle
   avec une discretisation reguliere.
   IN: n,m: nb de pts de controle
   IN: orderm, orderm: ordre de la spline
   IN: N,M: nbre de pts sur la surface (si density<0)
   IN: density: npts par unite de longueur de la bezier resultante
   IN: xt, yt, zt: pts de controle de la courbe
   OUT: coord: tableau des coordonnees des pts de la courbe de Bezier 
*/
//=============================================================================
void K_COMPGEOM::regularNurbs2D(E_Int n, E_Int m, E_Int ordern, E_Int N, 
                                E_Int orderm, E_Int M, 
                                E_Float density, 
                                E_Float* xt, E_Float* yt, E_Float* zt, 
                                E_Float* W,
                                FldArrayF& coord, E_Int& niout, E_Int& njout)
{
  // approx de la longueur de la spline par la longueur des pts de controle
  E_Float leni = 0.;
  E_Float lenj = 0.;

  #pragma omp parallel
  {
    E_Int ind, ind1, j1;
    E_Float len = 0.;
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
  if (density > 0.) 
  { nptsi = E_Int(leni*density)+1; nptsj = E_Int(lenj*density)+1; }
  else { nptsi = N; nptsj = M; }
  E_Int nptsi0 = 10*nptsi;
  E_Int nptsj0 = 10*nptsj;

  FldArrayF coord0(nptsi0*nptsj0, 3); coord0.setAllValuesAtNull();
  E_Float* coordx0 = coord0.begin(1);
  E_Float* coordy0 = coord0.begin(2);
  E_Float* coordz0 = coord0.begin(3);
  nurbs2D(n, m, ordern, nptsi0, orderm, nptsj0, xt, yt, zt, W, coord0);

  E_Int nthreads = __NUMTHREADS__;
  E_Float* lp = new E_Float [nthreads];

  // Remaille suivant i
  // vraie longueur de la spline suivant i
  E_Float len = 0.;
  #pragma omp parallel
  {
    E_Float dx, dy, dz;
    E_Int ind, ind0;
    E_Int it = __CURRENT_THREAD__;
    lp[it] = 0.;
    #pragma omp for
    for (E_Int i = 1; i < nptsi0; i++)
    {
      ind = i; ind0 = i-1;
      dx = coordx0[ind] - coordx0[ind0];
      dy = coordy0[ind] - coordy0[ind0];
      dz = coordz0[ind] - coordz0[ind0];
      lp[it] += sqrt(dx*dx+dy*dy+dz*dz);
    }
  }
  for (E_Int it = 0; it < nthreads; it++)
  {
    len += lp[it];
  }

  if (density > 0.) nptsi = E_Int(len*density)+1;
  else nptsi = N;

  FldArrayF sp(nptsi0);
  FldArrayF dxp(nptsi0);
  FldArrayF dyp(nptsi0);
  FldArrayF dzp(nptsi0);
  FldArrayF fd(nptsi);
  E_Float* fdx = fd.begin(1);
  E_Float delta = 1./(nptsi-1.);
  #pragma omp parallel for
  for (E_Int i = 0; i < nptsi; i++) fdx[i] = delta*i;
  FldArrayF coord1(nptsi*nptsj0, 3); coord1.setAllValuesAtNull();
  E_Float* coordx1 = coord1.begin(1);
  E_Float* coordy1 = coord1.begin(2);
  E_Float* coordz1 = coord1.begin(3);

  for (E_Int j = 0; j < nptsj0; j++)
  {
    K_COMPGEOM::onedmap(nptsi0, coordx0+j*nptsi0, coordy0+j*nptsi0, coordz0+j*nptsi0,
			nptsi, fd.begin(1),
			coordx1+j*nptsi, coordy1+j*nptsi, coordz1+j*nptsi,
			sp.begin(), dxp.begin(), dyp.begin(), dzp.begin());
  }

  // Remaille en j
  // vraie longueur de la spline suivant j
  len = 0.;
  #pragma omp parallel
  {
    E_Float dx, dy, dz;
    E_Int ind, ind0;
    E_Int it = __CURRENT_THREAD__;
    lp[it] = 0.;

    #pragma omp for
    for (E_Int j = 1; j < nptsj0; j++)
    {
      ind = j*nptsi0; ind0 = (j-1)*nptsi0;
      dx = coordx0[ind] - coordx0[ind0];
      dy = coordy0[ind] - coordy0[ind0];
      dz = coordz0[ind] - coordz0[ind0];
      lp[it] += sqrt(dx*dx+dy*dy+dz*dz);
    }
  }
  for (E_Int it = 0; it < nthreads; it++)
  {
    len += lp[it];
  }

  if (density > 0.) nptsj = E_Int(len*density)+1;
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
    { 
      t0x[j] = coordx1[i+j*nptsi];
      t0y[j] = coordy1[i+j*nptsi];
      t0z[j] = coordz1[i+j*nptsi]; 
    }
    K_COMPGEOM::onedmap(nptsj0, t0x, t0y, t0z,
			nptsj, fd.begin(1),
			tx, ty, tz,
			sp.begin(), dxp.begin(), dyp.begin(), dzp.begin());
    for (E_Int j = 0; j < nptsj; j++) 
    { 
      coordx[i+j*nptsi] = tx[j]; 
      coordy[i+j*nptsi] = ty[j]; 
      coordz[i+j*nptsi] = tz[j]; 
    }
  }
  niout = nptsi; njout = nptsj;
  delete [] lp;
}

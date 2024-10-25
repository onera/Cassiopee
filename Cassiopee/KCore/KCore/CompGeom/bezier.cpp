/*    
    Copyright 2013-2024 Onera.

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
# include <math.h>

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
//=============================================================================
/* 
   Calcul la courbe de bezier a partir des points de controle.
   IN: n: nb de pts de controle
   IN: N: nb de pts totaux dans la courbe de Bezier finale
   IN: xt, yt, zt: pts de controle de la courbe
   OUT: coord: tableau des coordonnees des pts de la courbe de Bezier 
        resultante (doit etre deja alloue a Nx3)
*/
//=============================================================================
void K_COMPGEOM::bezier(E_Int n, E_Int N, 
                        E_Float* xt, E_Float* yt, E_Float* zt,
                        FldArrayF& coord)
{
  n = n-1;
  E_Float pas = 1./(N-1);
  coord.setAllValuesAtNull();
  E_Float* coordx = coord.begin(1);
  E_Float* coordy = coord.begin(2);
  E_Float* coordz = coord.begin(3);

#pragma omp parallel for default(shared)
  for (E_Int k = 0; k < N; k++)
  {
    E_Float t = k*pas;
    for (E_Int i = 0; i <= n ; i++)
    {
      E_Float B = Bernstein(i, n, t);
      coordx[k] = coordx[k]+B*xt[i];
      coordy[k] = coordy[k]+B*yt[i];
      coordz[k] = coordz[k]+B*zt[i];
    }
  }
}

//=============================================================================
/* 
   Calcul la courbe de bezier a partir des points de controle.
   IN: nxm: nb de pts de controle
   IN: NxM: nb de pts totaux dans la courbe de Bezier finale
   IN: xt, yt, zt: pts de controle de la courbe
   OUT: coord: tableau des coordonnees des pts de la courbe de Bezier 
        resultante (doit etre deja alloue a NxMx3)
*/
//===========================================================================
void K_COMPGEOM::bezier2D(E_Int n, E_Int m, E_Int N, E_Int M, 
                          E_Float* xt, E_Float* yt, E_Float* zt,
                          FldArrayF& coord)
{
  n = n-1; m = m-1;
  
  E_Float pasN = 1./(N-1);
  E_Float pasM = 1./(M-1);
  coord.setAllValuesAtNull();
  E_Float* coordx = coord.begin(1);
  E_Float* coordy = coord.begin(2);
  E_Float* coordz = coord.begin(3);
  E_Int np1 = n+1;
#pragma omp parallel default(shared) if (M*N > __MIN_SIZE_MEAN__)
  {
    E_Float B, B1;
    E_Float u, v;
    E_Int k, l;
    E_Int ind1, ind2;
#pragma omp for
    for (E_Int kl = 0; kl < M*N; kl++)
    {
      k = int(kl/N); l = kl%N;
      v = k*pasM;
      u = l*pasN;
      ind2 = l + k*N;
      for (E_Int j = 0; j <= m ; j++)
      { 
        B1 = Bernstein(j, m, u);
        for (E_Int i = 0; i <= n ; i++)
        {
          B = B1*Bernstein(i, n, v);
          ind1 = i + j*np1;
          coordx[ind2] = coordx[ind2] + B*xt[ind1];
          coordy[ind2] = coordy[ind2] + B*yt[ind1];
          coordz[ind2] = coordz[ind2] + B*zt[ind1];
        }
      }
    }
  }
}

//=============================================================================
E_Int K_COMPGEOM::factorielle(E_Int N) 
{
  if (N <= 1) return 1;
  else return(N*K_COMPGEOM::factorielle(N-1));
}

E_Int K_COMPGEOM::combinations(E_Int n, E_Int i)
{
    if (i > n) return 0;
    if (i == 0 or i == n) return 1;
    if (i > n - i) i = n - i; // symmetry
    E_Int comb = 1;
    for (E_Int j = 0; j < i; j++) comb = comb*(n - j)/(j + 1);
    return comb;
}

//=============================================================================
E_Float K_COMPGEOM::Bernstein(E_Int i, E_Int n, E_Float t) 
{
  return combinations(n, i) * pow(t, i) * pow(1.-t, n-i);
}

//===========================================================================
/* 
   Calcul la courbe de bezier a partir des points de controle
   avec une discretisation reguliere.
   IN: n: nb de pts de controle
   IN: N: nombre de pts sur la courbe finale(si density>0, N est ignore)
   IN: density: npts par unite de longueur de la bezier resultante
   IN: xt, yt, zt: pts de controle de la courbe
   OUT: coord: tableau des coordonnees des pts de la courbe de Bezier 
*/
//===========================================================================
void K_COMPGEOM::regularBezier(E_Int n, E_Int N, E_Float density,
                               E_Float* xt, E_Float* yt, E_Float* zt,
                               FldArrayF& coord)
{
  // approx de la longueur de la bezier par la longueur des pts de controle
  E_Float len = 0.;
#pragma omp parallel default(shared) if (n > __MIN_SIZE_MEAN__)
  {
    E_Float dx, dy, dz;
    E_Int i1;
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
  bezier(n, npts0, xt, yt, zt, coord0);

  // vraie longueur de la bezier
  len = 0.;
#pragma omp parallel default(shared) if (npts0 > __MIN_SIZE_MEAN__)
  {
    E_Float dx, dy, dz;
    E_Int i1;
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
 
  // bezier reguliere
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
   Calcul la surface de bezier a partir des points de controle
   avec une discretisation reguliere.
   IN: n: nb de pts de controle
   IN: N,M: nbre de pts sur la surface (si density<0)
   IN: density: npts par unite de longueur de la bezier resultante
   IN: xt, yt, zt: pts de controle de la courbe
   OUT: coord: tableau des coordonnees des pts de la courbe de Bezier 
*/
//=============================================================================
void K_COMPGEOM::regularBezier2D(E_Int n, E_Int m, E_Int N, E_Int M, 
                                 E_Float density, 
                                 E_Float* xt, E_Float* yt, E_Float* zt,
                                 FldArrayF& coord, E_Int& niout, E_Int& njout)
{
  // approx de la longueur de la bezier par la longueur des pts de controle
  E_Float len = 0.;
  E_Float leni=0.;
  E_Float lenj=0.;
  E_Float dx, dy, dz;
  E_Int ind, ind1, j1;
#pragma omp parallel default(shared) if (m*n > __MIN_SIZE_MEAN__)
  {
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
        len += sqrt(dx*dx+dy*dy+dz*dz);
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
        len += sqrt(dx*dx+dy*dy+dz*dz);
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
  bezier2D(n, m, nptsi0, nptsj0, xt, yt, zt, coord0);

  // Remaille suivant i
  // vraie longueur de la bezier suivant i
  len = 0.;
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
  // vraie longueur de la bezier suivant j
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

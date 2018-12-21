/*    
    Copyright 2013-2019 Onera.

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
#include "StructZone.h"
#include "Functions.h"
#include <math.h>

//=============================================================================
/* 
   Calcul les normales pour les interfaces exterieures d'une zone structuree.
   Pour chaque noeud, on calcul la moyenne des normales des triangles
   adjacents au noeud.
   Allocate zone->surf.
   Warning: les normales ne sont pas necessairement exterieures.
*/
//=============================================================================
void StructZone::compNorm()
{
  int nij = ni*nj;
  int nik = ni*nk;
  int njk = nj*nk;
  int ni1 = ni-1;
  int nj1 = nj-1;
  int nk1 = nk-1;
  int nim = MIN(ni, 2);
  int njm = MIN(nj, 2);
  int nkm = MIN(nk, 2);

  // Nombre d'interfaces
  int nbElti = njk*nim;
  int nbEltj = nik*njm;
  int nbEltk = nij*nkm;
  int nbElt = nbElti+nbEltj+nbEltk;

  surf = new float [nbElt*3];
  float* surfx = surf;
  float* surfy = surfx + nbElti;
  float* surfz = surfy + nbElti;

  float* surfx2 = surf + 3*nbElti;
  float* surfy2 = surfx2 + nbEltj;
  float* surfz2 = surfy2 + nbEltj;

  float* surfx3 = surf + 3*nbElti + 3*nbEltj;
  float* surfy3 = surfx3 + nbEltk;
  float* surfz3 = surfy3 + nbEltk;

//#pragma omp parallel default(shared)
  {
    int n1, n2, n3, n4, n5, i, j, k, p, indp;
    double x0, x1, x2, x3, x4;
    double y0, y1, y2, y3, y4;
    double z0, z1, z2, z3, z4;
    double vx, vy, vz;
    double vx1, vy1, vz1;
    double vx2, vy2, vz2;
    double vx3, vy3, vz3;
    double vx4, vy4, vz4;
    float normi;

// Interfaces en i
//#pragma omp for
  for (int ind = 0; ind < nbElti; ind++)
  {
    p = ind/njk; // plane
    indp = ind - p*njk;
    k = indp/nj;
    j = (indp-k*nj);
    i = p*ni1;
    n1 = i+j*ni+k*nij;
    n2 = j < nj1 ? n1+ni : n1;
    n3 = j > 0 ? n1-ni : n1;
    n4 = k < nk1 ? n1+nij : n1;
    n5 = k > 0 ? n1-nij : n1;

    x0 = x[n1]; y0 = y[n1]; z0 = z[n1];
    x1 = x[n2]; y1 = y[n2]; z1 = z[n2];
    x2 = x[n3]; y2 = y[n3]; z2 = z[n3];
    x3 = x[n4]; y3 = y[n4]; z3 = z[n4];
    x4 = x[n5]; y4 = y[n5]; z4 = z[n5];

    // Normal to triangle
    vx1 = (y1-y0)*(z3-z0)-(z1-z0)*(y3-y0);
    vy1 = (z1-z0)*(x3-x0)-(x1-x0)*(z3-z0);
    vz1 = (x1-x0)*(y3-y0)-(y1-y0)*(x3-x0);

    vx2 = (y3-y0)*(z2-z0)-(z3-z0)*(y2-y0);
    vy2 = (z3-z0)*(x2-x0)-(x3-x0)*(z2-z0);
    vz2 = (x3-x0)*(y2-y0)-(y3-y0)*(x2-x0);
    
    vx3 = (y2-y0)*(z4-z0)-(z2-z0)*(y4-y0);
    vy3 = (z2-z0)*(x4-x0)-(x2-x0)*(z4-z0);
    vz3 = (x2-x0)*(y4-y0)-(y2-y0)*(x4-x0);

    vx4 = (y4-y0)*(z1-z0)-(z4-z0)*(y1-y0);
    vy4 = (z4-z0)*(x1-x0)-(x4-x0)*(z1-z0);
    vz4 = (x4-x0)*(y1-y0)-(y4-y0)*(x1-x0);
        
    vx = (vx1 + vx2 + vx3 + vx4);
    vy = (vy1 + vy2 + vy3 + vy4);
    vz = (vz1 + vz2 + vz3 + vz4);

    normi = invsqrt((float)(vx*vx+vy*vy+vz*vz));
    
    surfx[ind] = vx*normi;
    surfy[ind] = vy*normi;
    surfz[ind] = vz*normi;
  }

  // Interfaces en j
//#pragma omp for
  for (int ind = 0; ind < nbEltj; ind++)
  {
    p = ind/(ni*nk); // plane
    indp = ind - p*nik;
    k = indp/ni;
    j = p*nj1;
    i = indp-k*ni;
    n1 = i+j*ni+k*nij;
    n2 = i < ni1 ? n1+1 : n1;
    n3 = i > 0 ? n1-1 : n1;
    n4 = k < nk1 ? n1+nij : n1;
    n5 = k > 0 ? n1-nij : n1;

    x0 = x[n1]; y0 = y[n1]; z0 = z[n1];
    x1 = x[n2]; y1 = y[n2]; z1 = z[n2];
    x2 = x[n3]; y2 = y[n3]; z2 = z[n3];
    x3 = x[n4]; y3 = y[n4]; z3 = z[n4];
    x4 = x[n5]; y4 = y[n5]; z4 = z[n5];
        
    // Normal to triangle
    vx1 = (y1-y0)*(z3-z0)-(z1-z0)*(y3-y0);
    vy1 = (z1-z0)*(x3-x0)-(x1-x0)*(z3-z0);
    vz1 = (x1-x0)*(y3-y0)-(y1-y0)*(x3-x0);

    vx2 = (y3-y0)*(z2-z0)-(z3-z0)*(y2-y0);
    vy2 = (z3-z0)*(x2-x0)-(x3-x0)*(z2-z0);
    vz2 = (x3-x0)*(y2-y0)-(y3-y0)*(x2-x0);
        
    vx3 = (y2-y0)*(z4-z0)-(z2-z0)*(y4-y0);
    vy3 = (z2-z0)*(x4-x0)-(x2-x0)*(z4-z0);
    vz3 = (x2-x0)*(y4-y0)-(y2-y0)*(x4-x0);
    
    vx4 = (y4-y0)*(z1-z0)-(z4-z0)*(y1-y0);
    vy4 = (z4-z0)*(x1-x0)-(x4-x0)*(z1-z0);
    vz4 = (x4-x0)*(y1-y0)-(y4-y0)*(x1-x0);
    
    vx = (vx1 + vx2 + vx3 + vx4);
    vy = (vy1 + vy2 + vy3 + vy4);
    vz = (vz1 + vz2 + vz3 + vz4);
    
    normi = invsqrt((float)(vx*vx+vy*vy+vz*vz));
    
    surfx2[ind] = vx*normi;
    surfy2[ind] = vy*normi;
    surfz2[ind] = vz*normi;
  }

  // Interfaces en k
//#pragma omp for
  for (int ind = 0; ind < nbEltk; ind++)
  {
    p = ind/nij; // plane
    indp = ind - p*nij;
    k = p*nk1;
    j = indp/ni;
    i = indp-j*ni;
    n1 = i+j*ni+k*nij;
    n2 = i < ni1 ? n1+1 : n1;
    n3 = i > 0 ? n1-1 : n1;
    n4 = j < nj1 ? n1+ni : n1;
    n5 = j > 0 ? n1-ni : n1;
        
    x0 = x[n1]; y0 = y[n1]; z0 = z[n1];
    x1 = x[n2]; y1 = y[n2]; z1 = z[n2];
    x2 = x[n3]; y2 = y[n3]; z2 = z[n3];
    x3 = x[n4]; y3 = y[n4]; z3 = z[n4];
    x4 = x[n5]; y4 = y[n5]; z4 = z[n5];
    
    // Normal to triangle
    vx1 = (y1-y0)*(z3-z0)-(z1-z0)*(y3-y0);
    vy1 = (z1-z0)*(x3-x0)-(x1-x0)*(z3-z0);
    vz1 = (x1-x0)*(y3-y0)-(y1-y0)*(x3-x0);

    vx2 = (y3-y0)*(z2-z0)-(z3-z0)*(y2-y0);
    vy2 = (z3-z0)*(x2-x0)-(x3-x0)*(z2-z0);
    vz2 = (x3-x0)*(y2-y0)-(y3-y0)*(x2-x0);
        
    vx3 = (y2-y0)*(z4-z0)-(z2-z0)*(y4-y0);
    vy3 = (z2-z0)*(x4-x0)-(x2-x0)*(z4-z0);
    vz3 = (x2-x0)*(y4-y0)-(y2-y0)*(x4-x0);

    vx4 = (y4-y0)*(z1-z0)-(z4-z0)*(y1-y0);
    vy4 = (z4-z0)*(x1-x0)-(x4-x0)*(z1-z0);
    vz4 = (x4-x0)*(y1-y0)-(y4-y0)*(x1-x0);
        
    vx = (vx1 + vx2 + vx3 + vx4);
    vy = (vy1 + vy2 + vy3 + vy4);
    vz = (vz1 + vz2 + vz3 + vz4);
    
    normi = invsqrt((float)(vx*vx+vy*vy+vz*vz));
    
    surfx3[ind] = vx*normi;
    surfy3[ind] = vy*normi;
    surfz3[ind] = vz*normi;
  }
  }
}

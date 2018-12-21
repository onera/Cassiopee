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
#include "UnstructZone.h"
#include "Functions.h"
#include <math.h>                      
#include <stdio.h>

#define NORMTRI x0 = x[n1]; y0 = y[n1]; z0 = z[n1];     \
  x1 = x[n2]; y1 = y[n2]; z1 = z[n2];                   \
  x2 = x[n3]; y2 = y[n3]; z2 = z[n3];                   \
  vx = (y1-y0)*(z2-z0)-(z1-z0)*(y2-y0);                 \
  vy = (z1-z0)*(x2-x0)-(x1-x0)*(z2-z0);               \
  vz = (x1-x0)*(y2-y0)-(y1-y0)*(x2-x0);               \
  surfx[n1] += vx;                                     \
  surfy[n1] += vy;                                     \
  surfz[n1] += vz;                                     \
  surfx[n2] += vx;                                     \
  surfy[n2] += vy;                                     \
  surfz[n2] += vz;                                     \
  surfx[n3] += vx;                                     \
  surfy[n3] += vy;                                     \
  surfz[n3] += vz; 

#define NORMTRI2 x0 = x[n1]; y0 = y[n1]; z0 = z[n1];     \
  x1 = x[n2]; y1 = y[n2]; z1 = z[n2];                   \
  x2 = x[n3]; y2 = y[n3]; z2 = z[n3];                   \
  vx = (y1-y0)*(z2-z0)-(z1-z0)*(y2-y0);                 \
  vy = (z1-z0)*(x2-x0)-(x1-x0)*(z2-z0);                 \
  vz = (x1-x0)*(y2-y0)-(y1-y0)*(x2-x0);                 \
  normi = invsqrt((float)(vx*vx+vy*vy+vz*vz));          \
  vx = normi*vx; vy = normi*vy; vz = normi*vz;

#define NORMTRI3 x0 = xc; y0 = yc; z0 = zc;             \
  x1 = x[n1]; y1 = y[n1]; z1 = z[n1];                   \
  x2 = x[n2]; y2 = y[n2]; z2 = z[n2];                   \
  vx = (y1-y0)*(z2-z0)-(z1-z0)*(y2-y0);                 \
  vy = (z1-z0)*(x2-x0)-(x1-x0)*(z2-z0);                 \
  vz = (x1-x0)*(y2-y0)-(y1-y0)*(x2-x0);                 \
  normi = invsqrt((float)(vx*vx+vy*vy+vz*vz));          \
  vx = normi*vx; vy = normi*vy; vz = normi*vz;

#define NORMQUAD x0 = x[n1]; y0 = y[n1]; z0 = z[n1];    \
  x1 = x[n2]; y1 = y[n2]; z1 = z[n2];                   \
  x2 = x[n3]; y2 = y[n3]; z2 = z[n3];                   \
  x3 = x[n4]; y3 = y[n4]; z3 = z[n4];                   \
  vx1 = (y1-y0)*(z2-z0)-(z1-z0)*(y2-y0);                \
  vy1 = (z1-z0)*(x2-x0)-(x1-x0)*(z2-z0);                \
  vz1 = (x1-x0)*(y2-y0)-(y1-y0)*(x2-x0);                \
  vx2 = (y2-y0)*(z3-z0)-(z2-z0)*(y3-y0);                \
  vy2 = (z2-z0)*(x3-x0)-(x2-x0)*(z3-z0);                \
  vz2 = (x2-x0)*(y3-y0)-(y2-y0)*(x3-x0);                \
  vx = vx1+vx2;                                         \
  vy = vy1+vy2;                                         \
  vz = vz1+vz2;                                         \
  surfx[n1] += vx;                                       \
  surfy[n1] += vy;                                       \
  surfz[n1] += vz;                                       \
  surfx[n2] += vx;                                       \
  surfy[n2] += vy;                                       \
  surfz[n2] += vz;                                       \
  surfx[n3] += vx;                                       \
  surfy[n3] += vy;                                       \
  surfz[n3] += vz;                                       \
  surfx[n4] += vx;                                       \
  surfy[n4] += vy;                                       \
  surfz[n4] += vz;

#define NORMQUAD2 x0 = x[n1]; y0 = y[n1]; z0 = z[n1];    \
  x1 = x[n2]; y1 = y[n2]; z1 = z[n2];                   \
  x2 = x[n3]; y2 = y[n3]; z2 = z[n3];                   \
  x3 = x[n4]; y3 = y[n4]; z3 = z[n4];                   \
  vx1 = (y1-y0)*(z2-z0)-(z1-z0)*(y2-y0);                \
  vy1 = (z1-z0)*(x2-x0)-(x1-x0)*(z2-z0);                \
  vz1 = (x1-x0)*(y2-y0)-(y1-y0)*(x2-x0);                \
  vx2 = (y2-y0)*(z3-z0)-(z2-z0)*(y3-y0);                \
  vy2 = (z2-z0)*(x3-x0)-(x2-x0)*(z3-z0);                \
  vz2 = (x2-x0)*(y3-y0)-(y2-y0)*(x3-x0);                \
  vx = vx1+vx2;                                         \
  vy = vy1+vy2;                                         \
  vz = vz1+vz2;                                         \
  normi = invsqrt((float)(vx*vx+vy*vy+vz*vz));          \
  vx = normi*vx; vy = normi*vy; vz = normi*vz;

//=============================================================================
/* 
   Calcul les normales pour chaque noeud d'une zone non-structuree.
   Allocate zone->surf.
   Warning: les normales ne sont pas necessairement exterieures.
*/
//=============================================================================
void UnstructZone::compNorm()
{
  int i, n1, n2, n3, n4;
  double x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3;
  double vx1, vy1, vz1, vx2, vy2, vz2, vx, vy, vz;
  float normi;

  int ne2 = 2*ne;
  int ne3 = 3*ne;

  // TRI: stocke en vertex + lissage gouraud
  if (eltType == 2)
  {
    surf = new float[np * 3];
    float* surfx = surf;
    float* surfy = surfx + np;
    float* surfz = surfy + np;

//#pragma omp parallel default(shared) private(i,n1,n2,n3,x0,y0,z0,x1,y1,z1,x2,y2,z2,vx,vy,vz,normi)
    {
//#pragma omp for
      for (i = 0; i < 3*np; i++) surf[i] = 0.;
      
//#pragma omp for
      for (i = 0; i < ne; i++)
      {
        n1 = connect[i]-1;
        n2 = connect[i+ne]-1;
        n3 = connect[i+ne2]-1;
        NORMTRI;
      }
//#pragma omp for
      for (i = 0; i < np; i++)
      {
        vx = surfx[i]; vy = surfy[i]; vz = surfz[i];
        normi = invsqrt((float)(vx*vx+vy*vy+vz*vz));
        surfx[i] *= normi;
        surfy[i] *= normi;
        surfz[i] *= normi;
      }
    }
  }

  // QUADS: stocke en vertex + gouraud
  else if (eltType == 3)
  {
    surf = new float[np * 3];
    float* surfx = surf;
    float* surfy = surfx + np;
    float* surfz = surfy + np;

//#pragma omp parallel default(shared) private(i,n1,n2,n3,n4,x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3,vx1,vy1,vz1,vx2,vy2,vz2,vx,vy,vz,normi)
    {
//#pragma omp for
      for (i = 0; i < 3*np; i++) surf[i] = 0.;
    
//#pragma omp for
      for (i = 0; i < ne; i++)
      {
        n1 = connect[i]-1;
        n2 = connect[i+ne]-1;
        n3 = connect[i+ne2]-1;
        n4 = connect[i+ne3]-1;
        NORMQUAD;
      }

//#pragma omp for
      for (i = 0; i < np; i++)
      {
        vx = surfx[i]; vy = surfy[i]; vz = surfz[i];
        normi = invsqrt((float)(vx*vx+vy*vy+vz*vz));
        surfx[i] *= normi;
        surfy[i] *= normi;
        surfz[i] *= normi;
      }
    }
  }

  // TETRA: stocke par facette sans lissage
  if (eltType == 4)
  {
    surf = new float[ne * 4 * 3];

    float* surfx = surf;
    float* surfy = surfx + 4*ne;
    float* surfz = surfy + 4*ne;
    
//#pragma omp parallel for default(shared) private(i,n1,n2,n3,x0,y0,z0,x1,y1,z1,x2,y2,z2,vx,vy,vz,normi)
    for (i = 0; i < ne; i++)
    {
      n1 = connect[i]-1;
      n2 = connect[i+ne]-1;
      n3 = connect[i+ne2]-1;
      NORMTRI2;
      surfx[4*i] = vx; surfy[4*i] = vy; surfz[4*i] = vz;
      n1 = connect[i]-1;
      n2 = connect[i+ne]-1;
      n3 = connect[i+ne*3]-1;
      NORMTRI2;
      surfx[4*i+1] = vx; surfy[4*i+1] = vy; surfz[4*i+1] = vz;
      n1 = connect[i+ne]-1;
      n2 = connect[i+ne2]-1;
      n3 = connect[i+ne*3]-1;
      NORMTRI2;
      surfx[4*i+2] = vx; surfy[4*i+2] = vy; surfz[4*i+2] = vz;
      n1 = connect[i+ne2]-1;
      n2 = connect[i]-1;
      n3 = connect[i+ne*3]-1;
      NORMTRI2;
      surfx[4*i+3] = vx; surfy[4*i+3] = vy; surfz[4*i+3] = vz;
    }
  }

  // PENTA: stocke par facette sans lissage
  if (eltType == 5)
  {
    surf = new float[ne * 5 * 3];

    float* surfx = surf;
    float* surfy = surfx + 5*ne;
    float* surfz = surfy + 5*ne;
  
//#pragma omp parallel for default(shared) private(i,n1,n2,n3,n4,x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3,vx1,vy1,vz1,vx2,vy2,vz2,vx,vy,vz,normi)  
    for (i = 0; i < ne; i++)
    {
      n1 = connect[i]-1;
      n2 = connect[i+ne]-1;
      n3 = connect[i+ne2]-1;
      NORMTRI2;
      surfx[5*i] = vx; surfy[5*i] = vy; surfz[5*i] = vz;
      n1 = connect[i+ne*3]-1;
      n2 = connect[i+ne*4]-1;
      n3 = connect[i+ne*5]-1;
      NORMTRI2;
      surfx[5*i+1] = vx; surfy[5*i+1] = vy; surfz[5*i+1] = vz;
      n1 = connect[i]-1;
      n2 = connect[i+ne]-1;
      n3 = connect[i+ne*4]-1;
      n4 = connect[i+ne*3]-1;
      NORMQUAD2;
      surfx[5*i+2] = vx; surfy[5*i+2] = vy; surfz[5*i+2] = vz;
      n1 = connect[i+ne]-1;
      n2 = connect[i+ne*2]-1;
      n3 = connect[i+ne*5]-1;
      n4 = connect[i+ne*4]-1;
      NORMQUAD2;
      surfx[5*i+3] = vx; surfy[5*i+3] = vy; surfz[5*i+3] = vz;
      n1 = connect[i]-1;
      n2 = connect[i+ne*2]-1;
      n3 = connect[i+ne*5]-1;
      n4 = connect[i+ne*3]-1;
      NORMQUAD2;
      surfx[5*i+4] = vx; surfy[5*i+4] = vy; surfz[5*i+4] = vz;
    }
  }

  // PYRA: stocke par facette sans lissage
  if (eltType == 6)
  {
    surf = new float[ne * 5 * 3];

    float* surfx = surf;
    float* surfy = surfx + 5*ne;
    float* surfz = surfy + 5*ne;
    
//#pragma omp parallel for default(shared) private(i,n1,n2,n3,n4,x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3,vx1,vy1,vz1,vx2,vy2,vz2,vx,vy,vz,normi) 
    for (i = 0; i < ne; i++)
    {
      n1 = connect[i]-1;
      n2 = connect[i+ne]-1;
      n3 = connect[i+ne*4]-1;
      NORMTRI2;
      surfx[5*i] = vx; surfy[5*i] = vy; surfz[5*i] = vz;
      n1 = connect[i+ne]-1;
      n2 = connect[i+ne2]-1;
      n3 = connect[i+ne*4]-1;
      NORMTRI2;
      surfx[5*i+1] = vx; surfy[5*i+1] = vy; surfz[5*i+1] = vz;
      n1 = connect[i+ne2]-1;
      n2 = connect[i+ne*3]-1;
      n3 = connect[i+ne*4]-1;
      NORMTRI2;
      surfx[5*i+2] = vx; surfy[5*i+2] = vy; surfz[5*i+2] = vz;
      n1 = connect[i+ne*3]-1;
      n2 = connect[i]-1;
      n3 = connect[i+ne*4]-1;
      NORMTRI2;
      surfx[5*i+3] = vx; surfy[5*i+3] = vy; surfz[5*i+3] = vz;
      n1 = connect[i]-1;
      n2 = connect[i+ne]-1;
      n3 = connect[i+ne*2]-1;
      n4 = connect[i+ne*3]-1;
      NORMQUAD2;
      surfx[5*i+4] = vx; surfy[5*i+4] = vy; surfz[5*i+4] = vz;
    }
  }

  // HEXA
  else if (eltType == 7)
  {
    surf = new float[ne * 6 * 3];

    float* surfx = surf;
    float* surfy = surfx + 6*ne;
    float* surfz = surfy + 6*ne;

//#pragma omp parallel for default(shared) private(i,n1,n2,n3,n4,x0,y0,z0,x1,y1,z1,x2,y2,z2,x3,y3,z3,vx1,vy1,vz1,vx2,vy2,vz2,vx,vy,vz,normi)   
    for (i = 0; i < ne; i++)
    {
      n1 = connect[i]-1;
      n2 = connect[i+ne]-1;
      n3 = connect[i+ne2]-1;
      n4 = connect[i+ne3]-1;
      NORMQUAD2;
      surfx[6*i] = vx; surfy[6*i] = vy; surfz[6*i] = vz;
      n1 = connect[i+ne]-1;
      n2 = connect[i+5*ne]-1;
      n3 = connect[i+6*ne]-1;
      n4 = connect[i+ne2]-1;
      NORMQUAD2;
      surfx[6*i+1] = vx; surfy[6*i+1] = vy; surfz[6*i+1] = vz;
      n1 = connect[i+5*ne]-1;
      n2 = connect[i+4*ne]-1;
      n3 = connect[i+7*ne]-1;
      n4 = connect[i+6*ne]-1;
      NORMQUAD2;
      surfx[6*i+2] = vx; surfy[6*i+2] = vy; surfz[6*i+2] = vz;
      n1 = connect[i]-1;
      n2 = connect[i+4*ne]-1;
      n3 = connect[i+7*ne]-1;
      n4 = connect[i+ne3]-1;
      NORMQUAD2;
      surfx[6*i+3] = vx; surfy[6*i+3] = vy; surfz[6*i+3] = vz;
      n1 = connect[i]-1;
      n2 = connect[i+ne]-1;
      n3 = connect[i+5*ne]-1;
      n4 = connect[i+4*ne]-1;
      NORMQUAD2;
      surfx[6*i+4] = vx; surfy[6*i+4] = vy; surfz[6*i+4] = vz;
      n1 = connect[i+ne3]-1;
      n2 = connect[i+ne2]-1;
      n3 = connect[i+6*ne]-1;
      n4 = connect[i+7*ne]-1;
      NORMQUAD2;
      surfx[6*i+5] = vx; surfy[6*i+5] = vy; surfz[6*i+5] = vz;
    }
  }
  // NGON
  else if (eltType == 10)
  {
    int nf = connect[0];
    surf = new float[nf * 3];

    float* surfx = surf;
    float* surfy = surfx + nf;
    float* surfz = surfy + nf;
    int c = 2;
    int nd, l;
    float inv, vxm, vym, vzm, xc, yc, zc;

    for (i = 0; i < nf; i++)
    {
      nd = connect[c]; // nbre de noeuds de la face

      // traitement pour elements 3D
      // Face barycenter
      xc = 0.; yc = 0.; zc = 0.;
      for (l = 0; l < nd; l++)
      {
        n1 = connect[c+l+1]-1;
        xc += x[n1]; yc += y[n1]; zc += z[n1];
      }
      inv = 1./nd;
      xc = xc*inv; yc = yc*inv; zc = zc*inv;
      // Mean face normal (very approx)
      vxm = 0.; vym = 0.; vzm = 0.;
      for (l = 0; l < nd-1; l++)
      {
        n1 = connect[c+l+1]-1; n2 = connect[c+l+2]-1;
        NORMTRI3;
        vxm += vx; vym += vy; vzm += vz;
      }
      n1 = connect[c+nd]-1; n2 = connect[c+1]-1;
      NORMTRI3;
      vxm += vx; vym += vy; vzm += vz;
      
      vxm *= inv; vym *= inv; vzm *= inv;
      
      surfx[i] = vxm; surfy[i] = vym; surfz[i] = vzm;

      c += nd+1;
    }

    // Correction pour elements 2D (prise en compte d'un seul triangle)
    for (i = 0; i < nelts2D; i++)
    {
      int elt = posElts2D[i];
      int* ptrelt = &connect[elt];
      int nf = ptrelt[0];
      int face;
      int* ptrface;

      // barycentre de l'element
      xc = 0.; yc = 0.; zc = 0.;
      for (int j = 0; j < nf; j++)
      {
        face = ptrelt[j+1]-1;
        ptrface = &connect[posFaces[face]];
        n1 = ptrface[1]-1; n2 = ptrface[2]-1;
        xc += x[n1]; yc += y[n1]; zc += z[n1];
        xc += x[n2]; yc += y[n2]; zc += z[n2];
      }
      inv = 1./(2.*nf);
      xc = xc*inv; yc = yc*inv; zc = zc*inv;

      // face 1
      face = ptrelt[1]-1;
      ptrface = &connect[posFaces[face]];
      n1 = ptrface[1]-1; n2 = ptrface[2]-1;
      // face 2 
      face = ptrelt[2]-1;
      ptrface = &connect[posFaces[face]];
      n3 = ptrface[1]-1; n4 = ptrface[2]-1;

      if (n1 == n3 || n1 == n4)
      { n3 = n1; n1 = n2; n2 = n3; } // swap n1 et n2 - signe normale

      x0 = x[n1]; y0 = y[n1]; z0 = z[n1];   
      x1 = x[n2]; y1 = y[n2]; z1 = z[n2];                 
      x2 = xc; y2 = yc; z2 = zc;                  
      vx = (y1-y0)*(z2-z0)-(z1-z0)*(y2-y0);                
      vy = (z1-z0)*(x2-x0)-(x1-x0)*(z2-z0);              
      vz = (x1-x0)*(y2-y0)-(y1-y0)*(x2-x0);              
      normi = invsqrt((float)(vx*vx+vy*vy+vz*vz));
      vx = normi*vx; vy = normi*vy; vz = normi*vz;
      face = ptrelt[1]-1; // stocke sur la premiere face
      surfx[face] = vx; surfy[face] = vy; surfz[face] = vz;
    }
  }
}

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
#include "UnstructZone.h"
#include "Functions.h"
#include <string.h>
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
   Calcul les normales pour chaque noeud ou facette d'une zone non-structuree.
   Si nx,ny,nz present en noeuds, on les choisit pour TRI et QUAD.
   Allocate zone->surf.
   Warning: les normales ne sont pas necessairement exterieures.
*/
//=============================================================================
void UnstructZone::compNorm()
{

  // Detecte si les normales sont presentes aux noeuds dans les champs
  E_Float* pnx = NULL; E_Float* pny = NULL; E_Float* pnz = NULL;
  if (eltType[0] == 2 || eltType[0] == 3) // TRI or QUAD car normale en vertex
  {
    for (E_Int nv = 0; nv < nfield; nv++)
    {
      if (strcmp(varnames[nv], "_nx_") == 0) pnx = f[nv];
      if (strcmp(varnames[nv], "_ny_") == 0) pny = f[nv];
      if (strcmp(varnames[nv], "_nz_") == 0) pnz = f[nv];
    }
    if (pnx != NULL && pny != NULL && pnz != NULL)
    {
      float* surfp = new float[np * 3];
      surf.push_back(surfp);
      for (E_Int i = 0; i < np; i++) surfp[i] = pnx[i];
      for (E_Int i = 0; i < np; i++) surfp[i+np] = pny[i];
      for (E_Int i = 0; i < np; i++) surfp[i+2*np] = pnz[i];
      return;
    }
  }

  // Calcule les normales
  E_Int i, n1, n2, n3, n4;
  double x0, y0, z0, x1, y1, z1, x2, y2, z2, x3, y3, z3;
  double vx1, vy1, vz1, vx2, vy2, vz2, vx, vy, vz;
  float normi;

  E_Int eltType0 = eltType[0];

  // TRI or QUAD: stocke en vertex + lissage gouraud
  if (eltType0 == 2 || eltType0 == 3)
  {
    float* surfp = new float[np * 3];
    surf.push_back(surfp);
    float* surfx = surfp;
    float* surfy = surfx + np;
    float* surfz = surfy + np;

    for (i = 0; i < 3*np; i++) surfp[i] = 0.;

    for (size_t nc = 0; nc < connect.size(); nc++)
    { 
      E_Int ne1 = nec[nc];
      E_Int ne2 = 2*ne1;
      E_Int ne3 = 3*ne1;
      E_Int* zconnect = connect[nc];

      if (eltType[nc] == 2)
      {
        for (i = 0; i < ne1; i++)
        {
          n1 = zconnect[i]-1;
          n2 = zconnect[i+ne1]-1;
          n3 = zconnect[i+ne2]-1;
          NORMTRI; // acumulate normals
        }
      }
      else if (eltType[nc] == 3)
      {
        for (i = 0; i < ne1; i++)
        {
          n1 = zconnect[i]-1;
          n2 = zconnect[i+ne1]-1;
          n3 = zconnect[i+ne2]-1;
          n4 = zconnect[i+ne3]-1;
          NORMQUAD;
        }
      }
    }

    for (i = 0; i < np; i++) // normalize
    {
      vx = surfx[i]; vy = surfy[i]; vz = surfz[i];
      normi = invsqrt((float)(vx*vx+vy*vy+vz*vz));
      surfx[i] *= normi;
      surfy[i] *= normi;
      surfz[i] *= normi;
    }
  }

  // pour les autres elements, les normales sont par facettes
  if (eltType0 == 4 || eltType0 == 5 || eltType0 == 6 || eltType0 == 7)
  {
    for (size_t nc = 0; nc < connect.size(); nc++)
    {
      E_Int eltType1 = eltType[nc];

      // TETRA: stocke par facette sans lissage
      if (eltType1 == 4)
      {
        E_Int* zconnect = connect[nc];
        E_Int ne1 = nec[nc];
        E_Int ne2 = 2*ne1;
        E_Int ne3 = 3*ne1;
        float* surfp = new float[ne1 * 4 * 3];
        surf.push_back(surfp);
        float* surfx = surfp;
        float* surfy = surfx + 4*ne1;
        float* surfz = surfy + 4*ne1;
    
        for (i = 0; i < ne1; i++)
        {
          n1 = zconnect[i]-1;
          n2 = zconnect[i+ne1]-1;
          n3 = zconnect[i+ne2]-1;
          NORMTRI2;
          surfx[4*i] = vx; surfy[4*i] = vy; surfz[4*i] = vz;
          n1 = zconnect[i]-1;
          n2 = zconnect[i+ne1]-1;
          n3 = zconnect[i+ne3]-1;
          NORMTRI2;
          surfx[4*i+1] = vx; surfy[4*i+1] = vy; surfz[4*i+1] = vz;
          n1 = zconnect[i+ne1]-1;
          n2 = zconnect[i+ne2]-1;
          n3 = zconnect[i+ne3]-1;
          NORMTRI2;
          surfx[4*i+2] = vx; surfy[4*i+2] = vy; surfz[4*i+2] = vz;
          n1 = zconnect[i+ne2]-1;
          n2 = zconnect[i]-1;
          n3 = zconnect[i+ne3]-1;
          NORMTRI2;
          surfx[4*i+3] = vx; surfy[4*i+3] = vy; surfz[4*i+3] = vz;
        }
      }
      // HEXA
      else if (eltType1 == 7)
      {
        E_Int* zconnect = connect[nc];
        E_Int ne1 = nec[nc];
        E_Int ne2 = 2*ne1;
        E_Int ne3 = 3*ne1;
        float* surfp = new float[ne1 * 6 * 3];
        surf.push_back(surfp);

        float* surfx = surfp;
        float* surfy = surfx + 6*ne1;
        float* surfz = surfy + 6*ne1;

        for (i = 0; i < ne1; i++)
        {
          n1 = zconnect[i]-1;
          n2 = zconnect[i+ne1]-1;
          n3 = zconnect[i+ne2]-1;
          n4 = zconnect[i+ne3]-1;
          NORMQUAD2;
          surfx[6*i] = vx; surfy[6*i] = vy; surfz[6*i] = vz;
          n1 = zconnect[i+ne1]-1;
          n2 = zconnect[i+5*ne1]-1;
          n3 = zconnect[i+6*ne1]-1;
          n4 = zconnect[i+ne2]-1;
          NORMQUAD2;
          surfx[6*i+1] = vx; surfy[6*i+1] = vy; surfz[6*i+1] = vz;
          n1 = zconnect[i+5*ne1]-1;
          n2 = zconnect[i+4*ne1]-1;
          n3 = zconnect[i+7*ne1]-1;
          n4 = zconnect[i+6*ne1]-1;
          NORMQUAD2;
          surfx[6*i+2] = vx; surfy[6*i+2] = vy; surfz[6*i+2] = vz;
          n1 = zconnect[i]-1;
          n2 = zconnect[i+4*ne1]-1;
          n3 = zconnect[i+7*ne1]-1;
          n4 = zconnect[i+ne3]-1;
          NORMQUAD2;
          surfx[6*i+3] = vx; surfy[6*i+3] = vy; surfz[6*i+3] = vz;
          n1 = zconnect[i]-1;
          n2 = zconnect[i+ne1]-1;
          n3 = zconnect[i+5*ne1]-1;
          n4 = zconnect[i+4*ne1]-1;
          NORMQUAD2;
          surfx[6*i+4] = vx; surfy[6*i+4] = vy; surfz[6*i+4] = vz;
          n1 = zconnect[i+ne3]-1;
          n2 = zconnect[i+ne2]-1;
          n3 = zconnect[i+6*ne1]-1;
          n4 = zconnect[i+7*ne1]-1;
          NORMQUAD2;
          surfx[6*i+5] = vx; surfy[6*i+5] = vy; surfz[6*i+5] = vz;
        }        
      }
      // PENTA: stocke par facette sans lissage
      else if (eltType1 == 5)
      {
        E_Int* zconnect = connect[nc];
        E_Int ne1 = nec[nc];
        E_Int ne2 = 2*ne1;
        E_Int ne3 = 3*ne1;
        float* surfp = new float[ne1 * 5 * 3];
        surf.push_back(surfp);
        float* surfx = surfp;
        float* surfy = surfx + 5*ne1;
        float* surfz = surfy + 5*ne1;
  
        for (i = 0; i < ne1; i++)
        {
          n1 = zconnect[i]-1;
          n2 = zconnect[i+ne1]-1;
          n3 = zconnect[i+ne2]-1;
          NORMTRI2;
          surfx[5*i] = vx; surfy[5*i] = vy; surfz[5*i] = vz;
          n1 = zconnect[i+ne3]-1;
          n2 = zconnect[i+ne1*4]-1;
          n3 = zconnect[i+ne1*5]-1;
          NORMTRI2;
          surfx[5*i+1] = vx; surfy[5*i+1] = vy; surfz[5*i+1] = vz;
          n1 = zconnect[i]-1;
          n2 = zconnect[i+ne1]-1;
          n3 = zconnect[i+ne1*4]-1;
          n4 = zconnect[i+ne3]-1;
          NORMQUAD2;
          surfx[5*i+2] = vx; surfy[5*i+2] = vy; surfz[5*i+2] = vz;
          n1 = zconnect[i+ne1]-1;
          n2 = zconnect[i+ne2]-1;
          n3 = zconnect[i+ne1*5]-1;
          n4 = zconnect[i+ne1*4]-1;
          NORMQUAD2;
          surfx[5*i+3] = vx; surfy[5*i+3] = vy; surfz[5*i+3] = vz;
          n1 = zconnect[i]-1;
          n2 = zconnect[i+ne2]-1;
          n3 = zconnect[i+ne1*5]-1;
          n4 = zconnect[i+ne3]-1;
          NORMQUAD2;
          surfx[5*i+4] = vx; surfy[5*i+4] = vy; surfz[5*i+4] = vz;
        }
      }
      // PYRA: stocke par facette sans lissage
      else if (eltType1 == 6)
      {
        E_Int* zconnect = connect[nc];
        E_Int ne1 = nec[nc];
        E_Int ne2 = 2*ne1;
        E_Int ne3 = 3*ne1;
        float* surfp = new float[ne1 * 5 * 3];
        surf.push_back(surfp);

        float* surfx = surfp;
        float* surfy = surfx + 5*ne1;
        float* surfz = surfy + 5*ne1;
    
        for (i = 0; i < ne1; i++)
        {
          n1 = zconnect[i]-1;
          n2 = zconnect[i+ne1]-1;
          n3 = zconnect[i+ne1*4]-1;
          NORMTRI2;
          surfx[5*i] = vx; surfy[5*i] = vy; surfz[5*i] = vz;
          n1 = zconnect[i+ne1]-1;
          n2 = zconnect[i+ne2]-1;
          n3 = zconnect[i+ne1*4]-1;
          NORMTRI2;
          surfx[5*i+1] = vx; surfy[5*i+1] = vy; surfz[5*i+1] = vz;
          n1 = zconnect[i+ne2]-1;
          n2 = zconnect[i+ne3]-1;
          n3 = zconnect[i+ne1*4]-1;
          NORMTRI2;
          surfx[5*i+2] = vx; surfy[5*i+2] = vy; surfz[5*i+2] = vz;
          n1 = zconnect[i+ne3]-1;
          n2 = zconnect[i]-1;
          n3 = zconnect[i+ne1*4]-1;
          NORMTRI2;
          surfx[5*i+3] = vx; surfy[5*i+3] = vy; surfz[5*i+3] = vz;
          n1 = zconnect[i]-1;
          n2 = zconnect[i+ne1]-1;
          n3 = zconnect[i+ne2]-1;
          n4 = zconnect[i+ne3]-1;
          NORMQUAD2;
          surfx[5*i+4] = vx; surfy[5*i+4] = vy; surfz[5*i+4] = vz;
        }
      }
    }

  }

  // NGON
  else if (eltType0 == 10)
  {
    E_Int* zconnect = connect[0];
    E_Int nf = zconnect[0];
    float* surfp = new float[nf * 3];
    surf.push_back(surfp);

    float* surfx = surfp;
    float* surfy = surfx + nf;
    float* surfz = surfy + nf;
    E_Int c = 2;
    E_Int nd, l;
    float inv, vxm, vym, vzm, xc, yc, zc;

    for (i = 0; i < nf; i++)
    {
      nd = zconnect[c]; // nbre de noeuds de la face

      // traitement pour elements 3D
      // Face barycenter
      xc = 0.; yc = 0.; zc = 0.;
      for (l = 0; l < nd; l++)
      {
        n1 = zconnect[c+l+1]-1;
        xc += x[n1]; yc += y[n1]; zc += z[n1];
      }
      inv = 1./nd;
      xc = xc*inv; yc = yc*inv; zc = zc*inv;
      // Mean face normal (very approx)
      vxm = 0.; vym = 0.; vzm = 0.;
      for (l = 0; l < nd-1; l++)
      {
        n1 = zconnect[c+l+1]-1; n2 = zconnect[c+l+2]-1;
        NORMTRI3;
        vxm += vx; vym += vy; vzm += vz;
      }
      n1 = zconnect[c+nd]-1; n2 = zconnect[c+1]-1;
      NORMTRI3;
      vxm += vx; vym += vy; vzm += vz;
      
      vxm *= inv; vym *= inv; vzm *= inv;
      
      surfx[i] = vxm; surfy[i] = vym; surfz[i] = vzm;

      c += nd+1;
    }

    // Correction pour elements 2D (prise en compte d'un seul triangle)
    E_Int face, elt;
    E_Int* ptrface;
    for (i = 0; i < nelts2D; i++)
    {
      elt = posElts2D[i];
      E_Int* ptrelt = &zconnect[elt];
      nf = ptrelt[0];

      // barycentre de l'element
      xc = 0.; yc = 0.; zc = 0.;
      for (E_Int j = 0; j < nf; j++)
      {
        face = ptrelt[j+1]-1;
        ptrface = &zconnect[posFaces[face]];
        n1 = ptrface[1]-1; n2 = ptrface[2]-1;
        xc += x[n1]; yc += y[n1]; zc += z[n1];
        xc += x[n2]; yc += y[n2]; zc += z[n2];
      }
      inv = 1./(2.*nf);
      xc = xc*inv; yc = yc*inv; zc = zc*inv;

      // face 1
      face = ptrelt[1]-1;
      ptrface = &zconnect[posFaces[face]];
      n1 = ptrface[1]-1; n2 = ptrface[2]-1;
      // face 2 
      face = ptrelt[2]-1;
      ptrface = &zconnect[posFaces[face]];
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

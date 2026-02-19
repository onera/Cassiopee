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
#define PLOTTRI glNormal3f(surfx[n1], surfy[n1], surfz[n1]); \
  glVertex3d(x[n1], y[n1], z[n1]);                        \
  glNormal3f(surfx[n2], surfy[n2], surfz[n2]);            \
  glVertex3d(x[n2], y[n2], z[n2]);                        \
  glNormal3f(surfx[n3], surfy[n3], surfz[n3]);            \
  glVertex3d(x[n3], y[n3], z[n3]);

// Plot un triangle avec une normale par facette (pas de lissage)
#define PLOTTRI2 glNormal3f(surfx[f], surfy[f], surfz[f]); \
  glVertex3d(x[n1], y[n1], z[n1]);                         \
  glVertex3d(x[n2], y[n2], z[n2]);                         \
  glVertex3d(x[n3], y[n3], z[n3]);

// Fait un PLOTTRI avec blanking
#define PLOTTRIB ret1 = _pref.blanking->f(this, n1, zonep->blank, zonet);   \
  ret2 = _pref.blanking->f(this, n2, zonep->blank, zonet);                  \
  ret3 = _pref.blanking->f(this, n3, zonep->blank, zonet);                  \
  if (ret1*ret2*ret3 != 0) { PLOTTRI; }

// Fait un PLOTTRI2 avec blanking
#define PLOTTRI2B ret1 = _pref.blanking->f(this, n1, zonep->blank, zonet); \
  ret2 = _pref.blanking->f(this, n2, zonep->blank, zonet);                 \
  ret3 = _pref.blanking->f(this, n3, zonep->blank, zonet);                 \
  if (ret1*ret2*ret3 != 0) { PLOTTRI2; }

// Plot un quad avec les normales aux vertex (lissage gouraud)
#define PLOTQUAD glNormal3f(surfx[n1], surfy[n1], surfz[n1]); \
  glVertex3d(x[n1], y[n1], z[n1]);                            \
  glNormal3f(surfx[n2], surfy[n2], surfz[n2]);            \
  glVertex3d(x[n2], y[n2], z[n2]);                        \
  glNormal3f(surfx[n3], surfy[n3], surfz[n3]);            \
  glVertex3d(x[n3], y[n3], z[n3]);                        \
  glNormal3f(surfx[n4], surfy[n4], surfz[n4]);            \
  glVertex3d(x[n4], y[n4], z[n4]);

// Plot un quad avec une normale par facette (pas de lissage)
#define PLOTQUAD2 glNormal3f(surfx[f], surfy[f], surfz[f]);     \
  glVertex3d(x[n1], y[n1], z[n1]);                              \
  glVertex3d(x[n2], y[n2], z[n2]);                              \
  glVertex3d(x[n3], y[n3], z[n3]);                              \
  glVertex3d(x[n4], y[n4], z[n4]);

// Plot un quad avec des coords de texture
#define PLOTQUAD3 glNormal3f(surfx[n1], surfy[n1], surfz[n1]); \
  glTexCoord3f(0,0,0);                                         \
  glVertex3d(x[n1], y[n1], z[n1]);                             \
  glNormal3f(surfx[n2], surfy[n2], surfz[n2]);            \
  glTexCoord3f(1,0,0);                                    \
  glVertex3d(x[n2], y[n2], z[n2]);                        \
  glNormal3f(surfx[n3], surfy[n3], surfz[n3]);            \
  glTexCoord3f(1,1,0);                                    \
  glVertex3d(x[n3], y[n3], z[n3]);                        \
  glNormal3f(surfx[n4], surfy[n4], surfz[n4]);            \
  glTexCoord3f(0,1,0);                                    \
  glVertex3d(x[n4], y[n4], z[n4]);

// Plot un hexa avec des coords de texture
#define PLOTHEXA glNormal3f(surfx[n1], surfy[n1], surfz[n1]);   \
  glTexCoord3f(0,0,0);                                          \
  glVertex3d(x[n1], y[n1], z[n1]);                              \
  glNormal3f(surfx[n4], surfy[n4], surfz[n4]);                  \
  glTexCoord3f(0,1,0);                                          \
  glVertex3d(x[n4], y[n4], z[n4]);                              \
  glNormal3f(surfx[n3], surfy[n3], surfz[n3]);                  \
  glTexCoord3f(1,1,0);                                          \
  glVertex3d(x[n3], y[n3], z[n3]);                              \
  glNormal3f(surfx[n2], surfy[n2], surfz[n2]);                  \
  glTexCoord3f(1,0,0);                                          \
  glVertex3d(x[n2], y[n2], z[n2]);                        \
  glNormal3f(surfx[n5], surfy[n5], surfz[n5]);            \
  glTexCoord3f(0,0,1);                                    \
  glVertex3d(x[n5], y[n5], z[n5]);                        \
  glNormal3f(surfx[n6], surfy[n6], surfz[n6]);            \
  glTexCoord3f(1,0,1);                                    \
  glVertex3d(x[n6], y[n6], z[n6]);                        \
  glNormal3f(surfx[n7], surfy[n7], surfz[n7]);            \
  glTexCoord3f(1,1,1);                                    \
  glVertex3d(x[n7], y[n7], z[n7]);                        \
  glNormal3f(surfx[n8], surfy[n8], surfz[n8]);            \
  glTexCoord3f(0,1,1);                                    \
  glVertex3d(x[n8], y[n8], z[n8]);                        \
  glNormal3f(surfx[n1], surfy[n1], surfz[n1]);            \
  glTexCoord3f(0,0,0);                                    \
  glVertex3d(x[n1], y[n1], z[n1]);                        \
  glNormal3f(surfx[n2], surfy[n2], surfz[n2]);            \
  glTexCoord3f(1,0,0);                                    \
  glVertex3d(x[n2], y[n2], z[n2]);                        \
  glNormal3f(surfx[n6], surfy[n6], surfz[n6]);            \
  glTexCoord3f(1,0,1);                                    \
  glVertex3d(x[n6], y[n6], z[n6]);                        \
  glNormal3f(surfx[n5], surfy[n5], surfz[n5]);            \
  glTexCoord3f(0,0,1);                                    \
  glVertex3d(x[n5], y[n5], z[n5]);                        \
  glNormal3f(surfx[n4], surfy[n4], surfz[n4]);            \
  glTexCoord3f(0,1,0);                                    \
  glVertex3d(x[n4], y[n4], z[n4]);                        \
  glNormal3f(surfx[n8], surfy[n8], surfz[n8]);            \
  glTexCoord3f(0,1,1);                                    \
  glVertex3d(x[n8], y[n8], z[n8]);                        \
  glNormal3f(surfx[n7], surfy[n7], surfz[n7]);            \
  glTexCoord3f(1,1,1);                                    \
  glVertex3d(x[n7], y[n7], z[n7]);                        \
  glNormal3f(surfx[n3], surfy[n3], surfz[n3]);            \
  glTexCoord3f(1,1,0);                                    \
  glVertex3d(x[n3], y[n3], z[n3]);                        \
  glNormal3f(surfx[n1], surfy[n1], surfz[n1]);            \
  glTexCoord3f(0,0,0);                                    \
  glVertex3d(x[n1], y[n1], z[n1]);                        \
  glNormal3f(surfx[n5], surfy[n5], surfz[n5]);            \
  glTexCoord3f(0,0,1);                                    \
  glVertex3d(x[n5], y[n5], z[n5]);                        \
  glNormal3f(surfx[n8], surfy[n8], surfz[n8]);            \
  glTexCoord3f(0,1,1);                                    \
  glVertex3d(x[n8], y[n8], z[n8]);                        \
  glNormal3f(surfx[n4], surfy[n4], surfz[n4]);            \
  glTexCoord3f(0,1,0);                                    \
  glVertex3d(x[n4], y[n4], z[n4]);                        \
  glNormal3f(surfx[n2], surfy[n2], surfz[n2]);            \
  glTexCoord3f(1,0,0);                                    \
  glVertex3d(x[n2], y[n2], z[n2]);                        \
  glNormal3f(surfx[n3], surfy[n3], surfz[n3]);            \
  glTexCoord3f(1,1,0);                                    \
  glVertex3d(x[n3], y[n3], z[n3]);                        \
  glNormal3f(surfx[n7], surfy[n7], surfz[n7]);            \
  glTexCoord3f(1,1,1);                                    \
  glVertex3d(x[n7], y[n7], z[n7]);                        \
  glNormal3f(surfx[n6], surfy[n6], surfz[n6]);            \
  glTexCoord3f(1,0,1);                                    \
  glVertex3d(x[n6], y[n6], z[n6]);

// Fait un plotquad avec blanking
#define PLOTQUADB ret1 = _pref.blanking->f(this, n1, zonep->blank, zonet);   \
  ret2 = _pref.blanking->f(this, n2, zonep->blank, zonet);               \
  ret3 = _pref.blanking->f(this, n3, zonep->blank, zonet);               \
  ret4 = _pref.blanking->f(this, n4, zonep->blank, zonet);               \
  if (ret1*ret2*ret3*ret4 != 0) { PLOTQUAD; }

// Fait un plotquad2 avec blanking
#define PLOTQUAD2B ret1 = _pref.blanking->f(this, n1, zonep->blank, zonet); \
  ret2 = _pref.blanking->f(this, n2, zonep->blank, zonet);               \
  ret3 = _pref.blanking->f(this, n3, zonep->blank, zonet);               \
  ret4 = _pref.blanking->f(this, n4, zonep->blank, zonet);               \
  if (ret1*ret2*ret3*ret4 != 0) { PLOTQUAD2; }

  double* x = zonep->x;
  double* y = zonep->y;
  double* z = zonep->z;

  E_Int np = zonep->np;

  for (size_t nc = 0; nc < zonep->connect.size(); nc++) {

  E_Int ne = zonep->nec[nc];
  E_Int ne2 = 2*ne;
  E_Int ne3 = 3*ne;
  E_Int eltType = zonep->eltType[nc];
  E_Int* connect = zonep->connect[nc];
  
  if (eltType == 2) // TRI
  {
    float* surfp = zonep->surf[0];
    float* surfx = surfp;
    float* surfy = surfx + np;
    float* surfz = surfy + np;
    
    glBegin(GL_TRIANGLES);      
    
    if (zonep->blank == 0)
    {
      for (i = 0; i < ne; i++)
      {
        n1 = connect[i]-1;
        n2 = connect[i+ne]-1;
        n3 = connect[i+ne2]-1;
        PLOTTRI;
      }
    }
    else
    {
      // With blanking
      for (i = 0; i < ne; i++)
      {
        n1 = connect[i]-1;
        n2 = connect[i+ne]-1;
        n3 = connect[i+ne2]-1;
        PLOTTRIB;
      }
    }
    glEnd();
  }
  else if (eltType == 3) // QUADS
  {
    float* surfp = zonep->surf[0];
    float* surfx = surfp;
    float* surfy = surfx + np;
    float* surfz = surfy + np;
    
    glBegin(GL_QUADS);      
    
    if (zonep->blank == 0)
    {
      for (i = 0; i < ne; i++)
      {
        n1 = connect[i]-1;
        n2 = connect[i+ne]-1;
        n3 = connect[i+ne2]-1;
        n4 = connect[i+ne3]-1;
        PLOTQUAD;
      }
    }
    else
    {
      // With blanking
      for (i = 0; i < ne; i++)
      {
        n1 = connect[i]-1;
        n2 = connect[i+ne]-1;
        n3 = connect[i+ne2]-1;
        n4 = connect[i+ne3]-1;
        PLOTQUADB;
      }
    }
    glEnd();
  }
  else if (eltType == 4) // TETRA
  {
    float* surfp = zonep->surf[nc];
    float* surfx = surfp;
    float* surfy = surfx + 4*ne;
    float* surfz = surfy + 4*ne;
    
    glBegin(GL_TRIANGLES);
    
    if (zonep->blank == 0)
    {
      for (i = 0; i < ne; i++)
      {
        n1 = connect[i]-1;
        n2 = connect[i+ne]-1;
        n3 = connect[i+ne2]-1;
        f = 4*i;
        PLOTTRI2;
        n1 = connect[i]-1;
        n2 = connect[i+ne]-1;
        n3 = connect[i+ne*3]-1;
        f = 4*i+1;
        PLOTTRI2;
        n1 = connect[i+ne]-1;
        n2 = connect[i+ne2]-1;
        n3 = connect[i+ne*3]-1;
        f = 4*i+2;
        PLOTTRI2;
        n1 = connect[i+ne2]-1;
        n2 = connect[i]-1;
        n3 = connect[i+ne*3]-1;
        f = 4*i+3;
        PLOTTRI2;
      }
    }
    else
    {
      // With blanking
      for (i = 0; i < ne; i++)
      {
        n1 = connect[i]-1;
        n2 = connect[i+ne]-1;
        n3 = connect[i+ne2]-1;
        f = 4*i;
        PLOTTRI2B;
        n1 = connect[i]-1;
        n2 = connect[i+ne]-1;
        n3 = connect[i+ne*3]-1;
        f = 4*i+1;
        PLOTTRI2B;
        n1 = connect[i+ne]-1;
        n2 = connect[i+ne2]-1;
        n3 = connect[i+ne*3]-1;
        f = 4*i+2;
        PLOTTRI2B;
        n1 = connect[i+ne2]-1;
        n2 = connect[i]-1;
        n3 = connect[i+ne*3]-1;
        f = 4*i+3;
        PLOTTRI2B;
      }
    }
    glEnd();
  }
  else if (eltType == 5) // PENTA
  {
    float* surfp = zonep->surf[nc];
    float* surfx = surfp;
    float* surfy = surfx + 5*ne;
    float* surfz = surfy + 5*ne;
    
    glBegin(GL_TRIANGLES);
    if (zonep->blank == 0)
    {
      for (i = 0; i < ne; i++)
      {
        n1 = connect[i]-1;
        n2 = connect[i+ne]-1;
        n3 = connect[i+ne2]-1;
        f = 5*i;
        PLOTTRI2;
        n1 = connect[i+ne*3]-1;
        n2 = connect[i+ne*4]-1;
        n3 = connect[i+ne*5]-1;
        f = 5*i + 1;
        PLOTTRI2;
      }
    }
    else
    {
      // With blanking
      for (i = 0; i < ne; i++)
      {
        n1 = connect[i]-1;
        n2 = connect[i+ne]-1;
        n3 = connect[i+ne2]-1;
        f = 5*i;
        PLOTTRI2B;
        n1 = connect[i+ne*3]-1;
        n2 = connect[i+ne*4]-1;
        n3 = connect[i+ne*5]-1;
        f = 5*i+1;
        PLOTTRI2B;
      }
    }
    glEnd();

    glBegin(GL_QUADS);
    if (zonep->blank == 0)
    {
      for (i = 0; i < ne; i++)
      {
        n1 = connect[i]-1;
        n2 = connect[i+ne]-1;
        n3 = connect[i+ne*4]-1;
        n4 = connect[i+ne*3]-1;
        f = 5*i+2;
        PLOTQUAD2;
        n1 = connect[i+ne]-1;
        n2 = connect[i+ne*2]-1;
        n3 = connect[i+ne*5]-1;
        n4 = connect[i+ne*4]-1;
        f = 5*i+3;
        PLOTQUAD2;
        n1 = connect[i]-1;
        n2 = connect[i+ne2]-1;
        n3 = connect[i+ne*5]-1;
        n4 = connect[i+ne*3]-1;
        f = 5*i+4;
        PLOTQUAD2;
      }
    }
    else
    {
      // With blanking
      for (i = 0; i < ne; i++)
      {
        n1 = connect[i]-1;
        n2 = connect[i+ne]-1;
        n3 = connect[i+ne*4]-1;
        n4 = connect[i+ne*3]-1;
        f = 5*i+2;
        PLOTQUAD2B;
        n1 = connect[i+ne]-1;
        n2 = connect[i+ne*2]-1;
        n3 = connect[i+ne*5]-1;
        n4 = connect[i+ne*4]-1;
        f = 5*i+3;
        PLOTQUAD2B;
        n1 = connect[i]-1;
        n2 = connect[i+ne2]-1;
        n3 = connect[i+ne*5]-1;
        n4 = connect[i+ne*3]-1;
        f = 5*i+4;
        PLOTQUAD2B;
      }
    }
    glEnd();
  }
  else if (eltType == 6) // PYRA
  {
    float* surfp = zonep->surf[nc];
    float* surfx = surfp;
    float* surfy = surfx + 5*ne;
    float* surfz = surfy + 5*ne;
    
    glBegin(GL_TRIANGLES);
    if (zonep->blank == 0)
    {
      for (i = 0; i < ne; i++)
      {
        n1 = connect[i]-1;
        n2 = connect[i+ne]-1;
        n3 = connect[i+ne*4]-1;
        f = 5*i;
        PLOTTRI2;
        n1 = connect[i+ne]-1;
        n2 = connect[i+ne2]-1;
        n3 = connect[i+ne*4]-1;
        f = 5*i + 1;
        PLOTTRI2;
        n1 = connect[i+ne2]-1;
        n2 = connect[i+ne*3]-1;
        n3 = connect[i+ne*4]-1;
        f = 5*i + 2;
        PLOTTRI2;
        n1 = connect[i+ne*3]-1;
        n2 = connect[i]-1;
        n3 = connect[i+ne*4]-1;
        f = 5*i + 3;
        PLOTTRI2;
      }
    }
    else
    {
      // With blanking
      for (i = 0; i < ne; i++)
      {
        n1 = connect[i]-1;
        n2 = connect[i+ne]-1;
        n3 = connect[i+ne*4]-1;
        f = 5*i;
        PLOTTRI2B;
        n1 = connect[i+ne]-1;
        n2 = connect[i+ne2]-1;
        n3 = connect[i+ne*4]-1;
        f = 5*i + 1;
        PLOTTRI2B;
        n1 = connect[i+ne2]-1;
        n2 = connect[i+ne*3]-1;
        n3 = connect[i+ne*4]-1;
        f = 5*i + 2;
        PLOTTRI2B;
        n1 = connect[i+ne*3]-1;
        n2 = connect[i]-1;
        n3 = connect[i+ne*4]-1;
        f = 5*i + 3;
        PLOTTRI2B;
      }
    }
    glEnd();

    glBegin(GL_QUADS);
    if (zonep->blank == 0)
    {
      for (i = 0; i < ne; i++)
      {
        n1 = connect[i]-1;
        n2 = connect[i+ne]-1;
        n3 = connect[i+ne*2]-1;
        n4 = connect[i+ne*3]-1;
        f = 5*i + 4;
        PLOTQUAD2;
      }
    }
    else
    {
      // With blanking
      for (i = 0; i < ne; i++)
      {
        n1 = connect[i]-1;
        n2 = connect[i+ne]-1;
        n3 = connect[i+ne*2]-1;
        n4 = connect[i+ne*3]-1;
        f = 5*i + 4;
        PLOTQUAD2B;
      }
    }
    glEnd();
  }
  else if (eltType == 7 && zonep->material != 6) // HEXA
  {
    float* surfp = zonep->surf[nc];
    float* surfx = surfp;
    float* surfy = surfx + 6*ne;
    float* surfz = surfy + 6*ne;
    
    glBegin(GL_QUADS); 
    
    if (zonep->blank == 0)
    {
      for (i = 0; i < ne; i++)
      {
        n1 = connect[i]-1;
        n2 = connect[i+ne]-1;
        n3 = connect[i+ne2]-1;
        n4 = connect[i+ne3]-1;
        f = 6*i;
        PLOTQUAD2;
        n1 = connect[i+ne]-1;
        n2 = connect[i+5*ne]-1;
        n3 = connect[i+6*ne]-1;
        n4 = connect[i+ne2]-1;
        f = 6*i+1;
        PLOTQUAD2;
        n1 = connect[i+5*ne]-1;
        n2 = connect[i+4*ne]-1;
        n3 = connect[i+7*ne]-1;
        n4 = connect[i+6*ne]-1;
        f = 6*i+2;
        PLOTQUAD2;
        n1 = connect[i]-1;
        n2 = connect[i+4*ne]-1;
        n3 = connect[i+7*ne]-1;
        n4 = connect[i+3*ne]-1;
        f = 6*i+3;
        PLOTQUAD2;
        n1 = connect[i]-1;
        n2 = connect[i+ne]-1;
        n3 = connect[i+5*ne]-1;
        n4 = connect[i+4*ne]-1;
        f = 6*i+4;
        PLOTQUAD2;
        n1 = connect[i+3*ne]-1;
        n2 = connect[i+2*ne]-1;
        n3 = connect[i+6*ne]-1;
        n4 = connect[i+7*ne]-1;
        f = 6*i+5;
        PLOTQUAD2;
      }
    }
    else
    {
      // With blanking
      for (i = 0; i < ne; i++)
      {
        n1 = connect[i]-1;
        n2 = connect[i+ne]-1;
        n3 = connect[i+ne2]-1;
        n4 = connect[i+ne3]-1;
        f = 6*i;
        PLOTQUAD2B;
        n1 = connect[i+ne]-1;
        n2 = connect[i+5*ne]-1;
        n3 = connect[i+6*ne]-1;
        n4 = connect[i+ne2]-1;
        f = 6*i+1;
        PLOTQUAD2B;
        n1 = connect[i+5*ne]-1;
        n2 = connect[i+4*ne]-1;
        n3 = connect[i+7*ne]-1;
        n4 = connect[i+6*ne]-1;
        f = 6*i+2;
        PLOTQUAD2B;
        n1 = connect[i]-1;
        n2 = connect[i+4*ne]-1;
        n3 = connect[i+7*ne]-1;
        n4 = connect[i+3*ne]-1;
        f = 6*i+3;
        PLOTQUAD2B;
        n1 = connect[i]-1;
        n2 = connect[i+ne]-1;
        n3 = connect[i+5*ne]-1;
        n4 = connect[i+4*ne]-1;
        f = 6*i+4;
        PLOTQUAD2B;
        n1 = connect[i+3*ne]-1;
        n2 = connect[i+2*ne]-1;
        n3 = connect[i+6*ne]-1;
        n4 = connect[i+7*ne]-1;
        f = 6*i+5;
        PLOTQUAD2B;
      }
    }
    glEnd();
  }
  else if (eltType == 7 && zonep->material == 6) // HEXA for SMOKE
  {
    float* surfp = zonep->surf[nc];
    float* surfx = surfp;
    float* surfy = surfx + np;
    float* surfz = surfy + np;
    
    glBegin(GL_QUADS);
    for (i = 0; i < ne; i++)
    {
      n1 = connect[i]-1;
      n2 = connect[i+ne]-1;
      n3 = connect[i+ne2]-1;
      n4 = connect[i+ne3]-1;
      n5 = connect[i+4*ne]-1;
      n6 = connect[i+5*ne]-1;
      n7 = connect[i+6*ne]-1;
      n8 = connect[i+7*ne]-1;
      PLOTHEXA;
    }
    glEnd();
  }
  else if (eltType == 10) // NGON
  {
    E_Int nf = connect[0];
    E_Int l, nd, c;
    float* surfp = zonep->surf[nc];
    float* surfx = surfp;
    float* surfy = surfx + nf;
    float* surfz = surfy + nf;
    E_Int next, elt, drawn, blank;
    E_Int prev, j, first, face;
    E_Int* ptrelt;
    E_Int* ptrface;
    E_Int na1, na2, nb1, nb2;

    if (zonep->blank == 0)
    {
      // Faces des elements 3D
      if (zonep->ne != zonep->nelts2D)
      {
        c = 2;
        glBegin(GL_TRIANGLES);
        for (i = 0; i < nf; i++)
        {
          nd = connect[c]; 
          if (nd == 3) // TRI
          {
            f = i;
            n1 = connect[c+1]-1;
            n2 = connect[c+2]-1;
            n3 = connect[c+3]-1;
            PLOTTRI2;
          }
          c += nd+1;
        }
        glEnd();

        c = 2;
        glBegin(GL_QUADS);
        for (i = 0; i < nf; i++)
        {
          nd = connect[c]; 
          if (nd == 4) // QUAD
          {
            f = i;
            n1 = connect[c+1]-1;
            n2 = connect[c+2]-1;
            n3 = connect[c+3]-1;
            n4 = connect[c+4]-1;
            PLOTQUAD2;
          }
          c += nd+1;
        }
        glEnd();

        c = 2;
        for (i = 0; i < nf; i++)
        {
          nd = connect[c]; 
          if (nd > 4) // Elt 3D
          {
            glNormal3f(surfx[i], surfy[i], surfz[i]);
            glBegin(GL_POLYGON);
            for (l = 0; l < nd; l++)
            {
              n1 = connect[c+l+1]-1;
              glVertex3d(x[n1], y[n1], z[n1]);
            }
            glEnd();
          }
          c += nd+1;
        }
      }

      // Elements 2D
      if (zonep->ne == zonep->nelts2D)
      {
        surfx = surfp;
        surfy = surfx + np;
        surfz = surfy + np;

        // Find TRIs
        for (i = 0; i < zonep->nelts2D; i++)
        {
          glBegin(GL_TRIANGLES);
          for (i = 0; i < zonep->nelts2D; i++)
          {
            elt = zonep->posElts2D[i];
            ptrelt = &connect[elt];
            nf = ptrelt[0];
            if (nf == 3)
            {
              face = ptrelt[1]-1;
              ptrface = &connect[zonep->posFaces[face]];
              n1 = ptrface[1]-1; n2 = ptrface[2]-1;
              face = ptrelt[2]-1;
              ptrface = &connect[zonep->posFaces[face]];
              na1 = ptrface[1]-1; na2 = ptrface[2]-1;
              n3 = na1;
              if (na1 == n1) n3 = na2;
              else if (na1 == n2) n3 = na2;
              PLOTTRI;
            }
          }
          glEnd();
        }
        // Find QUADS
        for (i = 0; i < zonep->nelts2D; i++)
        {
          glBegin(GL_QUADS);
          for (i = 0; i < zonep->nelts2D; i++)
          {
            elt = zonep->posElts2D[i];
            ptrelt = &connect[elt];
            nf = ptrelt[0];
            if (nf == 4)
            {
              face = ptrelt[1]-1;
              ptrface = &connect[zonep->posFaces[face]];
              n1 = ptrface[1]-1; n2 = ptrface[2]-1;
              face = ptrelt[2]-1;
              ptrface = &connect[zonep->posFaces[face]];
              na1 = ptrface[1]-1; na2 = ptrface[2]-1;
              face = ptrelt[3]-1;
              ptrface = &connect[zonep->posFaces[face]];
              nb1 = ptrface[1]-1; nb2 = ptrface[2]-1;
              if (na1 == n1) 
              { 
                n4 = na2;
                if (nb1 != n1 && nb1 != n2 && nb1 != n4) n3 = nb1;
                else n3 = nb2;
              }
              else if (na1 == n2) 
              {
                n3 = na2;
                if (nb1 != n1 && nb1 != n2 && nb1 != n3) n4 = nb1;
                else n4 = nb2;
              }
              else if (na2 == n1) 
              {
                n4 = na1;
                if (nb1 != n1 && nb1 != n2 && nb1 != n4) n3 = nb1;
                else n3 = nb2;
              }
              else if (na2 == n2) 
              {
                n3 = na1;
                if (nb1 != n1 && nb1 != n2 && nb1 != n3) n4 = nb1;
                else n4 = nb2;
              }
              else if (nb1 == n1)
              {
                n4 = nb2;
                if (na1 != n1 && na1 != n2 && na1 != n4) n3 = na1;
                else n3 = na2;
              }
              else if (nb1 == n2) 
              {
                n3 = nb2;
                if (na1 != n1 && na1 != n2 && na1 != n3) n4 = na1;
                else n4 = na2;
              }
              else if (nb2 == n1) 
              {
                n4 = nb1;
                if (na1 != n1 && na1 != n2 && na1 != n4) n3 = na1;
                else n3 = na2;
              }
              else // if (nb2 == n2) 
              {
                n3 = nb1;
                if (na1 != n1 && na1 != n2 && na1 != n3) n4 = na1;
                else n4 = na2;
              }
              PLOTQUAD;
            }
          }
          glEnd();
        }

        for (i = 0; i < zonep->nelts2D; i++)
        {
          elt = zonep->posElts2D[i];
          ptrelt = &connect[elt];
          nf = ptrelt[0];
          if (nf == 3 || nf == 4) continue;
          
          glBegin(GL_POLYGON);
          drawn = 0;
          face = ptrelt[1]-1;
          //glNormal3f(surfx[face], surfy[face], surfz[face]);
          ptrface = &connect[zonep->posFaces[face]];
          n1 = ptrface[1]-1; n2 = ptrface[2]-1;

          face = ptrelt[2]-1;
          ptrface = &connect[zonep->posFaces[face]];
          n3 = ptrface[1]-1; n4 = ptrface[2]-1;

          if (n2 == n3 || n2 == n4)
          { first = n1; prev = n1; next = n2; }
          else { first = n2; prev = n2; next = n1; }

          glNormal3f(surfx[prev], surfy[prev], surfz[prev]);
          glVertex3d(x[prev], y[prev], z[prev]);

          glNormal3f(surfx[next], surfy[next], surfz[next]);
          glVertex3d(x[next], y[next], z[next]);
          drawn++;
        
          // Cherche
          while (drawn < nf)
          {
            for (j = 2; j <= nf; j++)
            {
              face = ptrelt[j]-1;
              ptrface = &connect[zonep->posFaces[face]];
              n1 = ptrface[1]-1; n2 = ptrface[2]-1;
              if (n1 == next && n2 != prev)
              { 
                glNormal3f(surfx[n2], surfy[n2], surfz[n2]);
                glVertex3d(x[n2], y[n2], z[n2]); prev = n1; next = n2; drawn++; break; 
              }
              else if (n2 == next && n1 != prev)
              {  
                glNormal3f(surfx[n1], surfy[n1], surfz[n1]);
                glVertex3d(x[n1], y[n1], z[n1]); prev = n2; next = n1; drawn++; break; 
              }
            }
            if (j == nf+1) drawn++; // pour eviter les boucles infinies
          }
          if (next != first) 
          {
            glNormal3f(surfx[first], surfy[first], surfz[first]);
            glVertex3d(x[first], y[first], z[first]); // force close
          }
          glEnd();
        }
      }
    }
    else // blanking
    {
      // Faces des elements 3D
      if (zonep->ne != zonep->nelts2D)
      {
        c = 2;
        glBegin(GL_TRIANGLES);
        for (i = 0; i < nf; i++)
        {
          nd = connect[c]; 
          if (nd == 3) // TRI
          {
            f = i;
            n1 = connect[c+1]-1;
            n2 = connect[c+2]-1;
            n3 = connect[c+3]-1;
            PLOTTRI2B;
          }
          c += nd+1;
        }
        glEnd();

        c = 2;
        glBegin(GL_QUADS);
        for (i = 0; i < nf; i++)
        {
          nd = connect[c]; 
          if (nd == 4) // QUAD
          {
            f = i;
            n1 = connect[c+1]-1;
            n2 = connect[c+2]-1;
            n3 = connect[c+3]-1;
            n4 = connect[c+4]-1;
            PLOTQUAD2B;
          }
          c += nd+1;
        }
        glEnd();

        c = 2;
        for (i = 0; i < nf; i++)
        {
          nd = connect[c]; 
          if (nd > 4) // elt 3D
          {
            blank = 0;
            for (l = 0; l < nd; l++)
            {
              n1 = connect[c+l+1]-1;
              if (_pref.blanking->f(this, n1, zonep->blank, zonet) == 0)
              { blank = 1; break; }
            }
            if (blank == 0)
            {
              glNormal3f(surfx[i], surfy[i], surfz[i]);
              glBegin(GL_POLYGON);
              for (l = 0; l < nd; l++)
              {
                n1 = connect[c+l+1]-1;
                glVertex3d(x[n1], y[n1], z[n1]);
              }
              glEnd();
            }
          }
          c += nd+1;
        }
      }

      // Elements 2D
      if (zonep->ne == zonep->nelts2D)
      {
        surfx = surfp;
        surfy = surfx + np;
        surfz = surfy + np;
        
        // Find TRIs
        for (i = 0; i < zonep->nelts2D; i++)
        {
          glBegin(GL_TRIANGLES);
          for (i = 0; i < zonep->nelts2D; i++)
          {
            elt = zonep->posElts2D[i];
            ptrelt = &connect[elt];
            nf = ptrelt[0];
            if (nf == 3)
            {
              face = ptrelt[1]-1;
              ptrface = &connect[zonep->posFaces[face]];
              n1 = ptrface[1]-1; n2 = ptrface[2]-1;
              face = ptrelt[2]-1;
              ptrface = &connect[zonep->posFaces[face]];
              na1 = ptrface[1]-1; na2 = ptrface[2]-1;
              n3 = na1;
              if (na1 == n1) n3 = na2;
              else if (na1 == n2) n3 = na2;
              PLOTTRIB;
            }
          }
          glEnd();
        }

        // Find QUADS
        for (i = 0; i < zonep->nelts2D; i++)
        {
          glBegin(GL_QUADS);
          for (i = 0; i < zonep->nelts2D; i++)
          {
            elt = zonep->posElts2D[i];
            ptrelt = &connect[elt];
            nf = ptrelt[0];
            if (nf == 4)
            {
              face = ptrelt[1]-1;
              ptrface = &connect[zonep->posFaces[face]];
              n1 = ptrface[1]-1; n2 = ptrface[2]-1;
              face = ptrelt[2]-1;
              ptrface = &connect[zonep->posFaces[face]];
              na1 = ptrface[1]-1; na2 = ptrface[2]-1;
              face = ptrelt[3]-1;
              ptrface = &connect[zonep->posFaces[face]];
              nb1 = ptrface[1]-1; nb2 = ptrface[2]-1;
              if (na1 == n1) 
              { 
                n4 = na2;
                if (nb1 != n1 && nb1 != n2 && nb1 != n4) n3 = nb1;
                else n3 = nb2;
              }
              else if (na1 == n2) 
              {
                n3 = na2;
                if (nb1 != n1 && nb1 != n2 && nb1 != n3) n4 = nb1;
                else n4 = nb2;
              }
              else if (na2 == n1) 
              {
                n4 = na1;
                if (nb1 != n1 && nb1 != n2 && nb1 != n4) n3 = nb1;
                else n3 = nb2;
              }
              else if (na2 == n2) 
              {
                n3 = na1;
                if (nb1 != n1 && nb1 != n2 && nb1 != n3) n4 = nb1;
                else n4 = nb2;
              }
              else if (nb1 == n1)
              {
                n4 = nb2;
                if (na1 != n1 && na1 != n2 && na1 != n4) n3 = na1;
                else n3 = na2;
              }
              else if (nb1 == n2) 
              {
                n3 = nb2;
                if (na1 != n1 && na1 != n2 && na1 != n3) n4 = na1;
                else n4 = na2;
              }
              else if (nb2 == n1) 
              {
                n4 = nb1;
                if (na1 != n1 && na1 != n2 && na1 != n4) n3 = na1;
                else n3 = na2;
              }
              else // if (nb2 == n2) 
              {
                n3 = nb1;
                if (na1 != n1 && na1 != n2 && na1 != n3) n4 = na1;
                else n4 = na2;
              }
              PLOTQUADB;
            }
          }
          glEnd();
        }

        for (i = 0; i < zonep->nelts2D; i++)
        {
          elt = zonep->posElts2D[i];
          ptrelt = &connect[elt];
          nf = ptrelt[0];
          if (nf == 3 || nf == 4) continue;
          
          blank = 0;
          for (E_Int j = 1; j <= nf; j++)
          {
            face = ptrelt[1]-1;
            ptrface = &connect[zonep->posFaces[face]];
            n1 = ptrface[1]-1;
            n2 = ptrface[2]-1;
            if (_pref.blanking->f(this, n1, zonep->blank, zonet) == 0)
            { blank = 1; break; }
            if (_pref.blanking->f(this, n2, zonep->blank, zonet) == 0)
            { blank = 1; break; }
          }
          if (blank == 0)
          {
            glBegin(GL_POLYGON);
            elt = zonep->posElts2D[i];
            ptrelt = &connect[elt];
            nf = ptrelt[0];
            drawn = 0;
          
            face = ptrelt[1]-1;
            //glNormal3f(surfx[face], surfy[face], surfz[face]);

            ptrface = &connect[zonep->posFaces[face]];
            n1 = ptrface[1]-1; first = n1;
            n2 = ptrface[2]-1;
            glNormal3f(surfx[n1], surfy[n1], surfz[n1]);
            glVertex3d(x[n1], y[n1], z[n1]);
            glNormal3f(surfx[n2], surfy[n2], surfz[n2]);
            glVertex3d(x[n2], y[n2], z[n2]);
            prev = n1; next = n2;
            drawn++;
        
            // Cherche
            while (drawn < nf)
            {
              for (j = 2; j <= nf; j++)
              {
                face = ptrelt[j]-1;
                ptrface = &connect[zonep->posFaces[face]];
                n1 = ptrface[1]-1;
                n2 = ptrface[2]-1;
                if (n1 == next && n2 != prev)
                { 
                  glNormal3f(surfx[n2], surfy[n2], surfz[n2]);
                  glVertex3d(x[n2], y[n2], z[n2]); prev = n1; next = n2; drawn++; break; 
                }
                else if (n2 == next && n1 != prev)
                { 
                  glNormal3f(surfx[n1], surfy[n1], surfz[n1]);
                  glVertex3d(x[n1], y[n1], z[n1]); prev = n2; next = n1; drawn++; break; 
                }
              }
              if (j == nf+1) drawn++; // pour eviter les boucles infinies
            }
            if (next != first) 
            { 
              glNormal3f(surfx[first], surfy[first], surfz[first]);
              glVertex3d(x[first], y[first], z[first]); // force close
            }
            glEnd();
          }
        }
      }
    }
  }

  // Pour les BAR
  if (eltType == 1)
  {
    glPolygonOffset(-1.,-10.); // force offset
    glBegin(GL_LINES);
    if (zonep->blank == 0)
    {
      for (i = 0; i < ne; i++)
      {
        n1 = connect[i]-1;
        n2 = connect[i+ne]-1;
        glVertex3d(x[n1], y[n1], z[n1]);
        glVertex3d(x[n2], y[n2], z[n2]);
      }
    }
    else
    {
      for (i = 0; i < ne; i++)
      {
        n1 = connect[i]-1;
        n2 = connect[i+ne]-1;
        ret1 = _pref.blanking->f(this, n1, zonep->blank, zonet);
        ret2 = _pref.blanking->f(this, n2, zonep->blank, zonet);
        
        if (ret1*ret2 != 0)
        {
          glVertex3d(x[n1], y[n1], z[n1]);
          glVertex3d(x[n2], y[n2], z[n2]);
        }
      }
    }
    glEnd();
  }

  // Pour les NGONS 1D
  if (eltType == 10 && zonep->nelts1D > 0)
  {
    E_Int elt, face1, face2, posface1, posface2;
    glBegin(GL_LINES);
    if (zonep->blank == 0)
    {
      for (i = 0; i < zonep->nelts1D; i++)
      {
        elt = zonep->posElts1D[i];
        E_Int* ptrelt = &connect[elt];
        face1 = ptrelt[1]-1;
        face2 = ptrelt[2]-1;
        posface1 = zonep->posFaces[face1];
        posface2 = zonep->posFaces[face2];
        n1 = connect[posface1+1]-1;
        n2 = connect[posface2+1]-1;
        glVertex3d(x[n1], y[n1], z[n1]);
        glVertex3d(x[n2], y[n2], z[n2]);
      }
    }
    else
    {
      for (i = 0; i < zonep->nelts1D; i++)
      {
        elt = zonep->posElts1D[i];
        E_Int* ptrelt = &connect[elt];
        face1 = ptrelt[1]-1;
        face2 = ptrelt[2]-1;
        posface1 = zonep->posFaces[face1];
        posface2 = zonep->posFaces[face2];
        n1 = connect[posface1+1]-1;
        n2 = connect[posface2+1]-1;
        ret1 = _pref.blanking->f(this, n1, zonep->blank, zonet);
        ret2 = _pref.blanking->f(this, n2, zonep->blank, zonet);
        if (ret1*ret2 != 0)
        {
          glVertex3d(x[n1], y[n1], z[n1]);
          glVertex3d(x[n2], y[n2], z[n2]);
        }
      }
    }
    glEnd();
  }

  } // connects
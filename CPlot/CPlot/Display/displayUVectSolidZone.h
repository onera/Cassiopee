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
// change this to draw 2 triangles instead of one quad
//#define GL_QUADS_ARE GL_QUADS/GL_TRIANGLES
//#define PLOTQUAD PLOTQUADQ/PLOTSQUADT
//#define PLOTQUAD2 PLOTQUADQ2/PLOTQUADT2

#define PLOTTRI r = f1[n1]*deltai+0.5;                              \
  g = f2[n1]*deltai+0.5;                                            \
  b = f3[n1]*deltai+0.5;                                            \
  glColor3f(r, g, b);                                               \
  glNormal3f(surfx[n1], surfy[n1], surfz[n1]);                      \
  glVertex3d(x[n1], y[n1], z[n1]);                                  \
  r = f1[n2]*deltai+0.5;                                            \
  g = f2[n2]*deltai+0.5;                                            \
  b = f3[n2]*deltai+0.5;                                            \
  glColor3f(r, g, b);                                               \
  glNormal3f(surfx[n2], surfy[n2], surfz[n2]);                      \
  glVertex3d(x[n2], y[n2], z[n2]);                                  \
  r = f1[n3]*deltai+0.5;                                            \
  g = f2[n3]*deltai+0.5;                                            \
  b = f3[n3]*deltai+0.5;                                            \
  glColor3f(r, g, b);                                               \
  glNormal3f(surfx[n3], surfy[n3], surfz[n3]);                      \
  glVertex3d(x[n3], y[n3], z[n3]);
                              
#define PLOTTRI2 r = f1[n1]*deltai+0.5;                             \
  g = f2[n1]*deltai+0.5;                                            \
  b = f3[n1]*deltai+0.5;                                            \
  glColor3f(r, g, b);                                               \
  glNormal3f(surfx[ff], surfy[ff], surfz[ff]);                      \
  glVertex3d(x[n1], y[n1], z[n1]);                                  \
  r = f1[n2]*deltai+0.5;                                            \
  g = f2[n2]*deltai+0.5;                                            \
  b = f3[n2]*deltai+0.5;                                            \
  glColor3f(r, g, b);                                               \
  glVertex3d(x[n2], y[n2], z[n2]);                                  \
  r = f1[n3]*deltai+0.5;                                            \
  g = f2[n3]*deltai+0.5;                                            \
  b = f3[n3]*deltai+0.5;                                            \
  glColor3f(r, g, b);                                               \
  glVertex3d(x[n3], y[n3], z[n3]); 

#define PLOTTRIB ret1 = _pref.blanking->f(this, n1, zonep->blank, zonet); \
  ret2 = _pref.blanking->f(this, n2, zonep->blank, zonet);               \
  ret3 = _pref.blanking->f(this, n3, zonep->blank, zonet);               \
  if (ret1*ret2*ret3 != 0) { PLOTTRI; }

#define PLOTTRI2B ret1 = _pref.blanking->f(this, n1, zonep->blank, zonet); \
  ret2 = _pref.blanking->f(this, n2, zonep->blank, zonet);               \
  ret3 = _pref.blanking->f(this, n3, zonep->blank, zonet);               \
  if (ret1*ret2*ret3 != 0) { PLOTTRI2; }

#define PLOTQUADQ r = f1[n1]*deltai+0.5;                             \
  g = f2[n1]*deltai+0.5;                                            \
  b = f3[n1]*deltai+0.5;                                            \
  glColor3f(r, g, b);                                               \
  glNormal3f(surfx[n1], surfy[n1], surfz[n1]);                      \
  glVertex3d(x[n1], y[n1], z[n1]);                                  \
  r = f1[n2]*deltai+0.5;                                            \
  g = f2[n2]*deltai+0.5;                                            \
  b = f3[n2]*deltai+0.5;                                            \
  glColor3f(r, g, b);                                               \
  glNormal3f(surfx[n2], surfy[n2], surfz[n2]);                      \
  glVertex3d(x[n2], y[n2], z[n2]);                                  \
  r = f1[n3]*deltai+0.5;                                            \
  g = f2[n3]*deltai+0.5;                                            \
  b = f3[n3]*deltai+0.5;                                            \
  glColor3f(r, g, b);                                               \
  glNormal3f(surfx[n3], surfy[n3], surfz[n3]);                      \
  glVertex3d(x[n3], y[n3], z[n3]);                                  \
  r = f1[n4]*deltai+0.5;                                            \
  g = f2[n4]*deltai+0.5;                                            \
  b = f3[n4]*deltai+0.5;                                            \
  glColor3f(r, g, b);                                               \
  glNormal3f(surfx[n4], surfy[n4], surfz[n4]);                      \
  glVertex3d(x[n4], y[n4], z[n4]);                                    

#define PLOTQUADT r = f1[n1]*deltai+0.5;                            \
  g = f2[n1]*deltai+0.5;                                            \
  b = f3[n1]*deltai+0.5;                                            \
  glColor3f(r, g, b);                                               \
  glNormal3f(surfx[n1], surfy[n1], surfz[n1]);                      \
  glVertex3d(x[n1], y[n1], z[n1]);                                  \
  r = f1[n2]*deltai+0.5;                                            \
  g = f2[n2]*deltai+0.5;                                            \
  b = f3[n2]*deltai+0.5;                                            \
  glColor3f(r, g, b);                                               \
  glNormal3f(surfx[n2], surfy[n2], surfz[n2]);                      \
  glVertex3d(x[n2], y[n2], z[n2]);                                  \
  r = f1[n4]*deltai+0.5;                                            \
  g = f2[n4]*deltai+0.5;                                            \
  b = f3[n4]*deltai+0.5;                                            \
  glColor3f(r, g, b);                                                   \
  glNormal3f(surfx[n4], surfy[n4], surfz[n4]);                          \
  glVertex3d(x[n4], y[n4], z[n4]);                                      \
  r = f1[n2]*deltai+0.5;                             \
  g = f2[n2]*deltai+0.5;                                            \
  b = f3[n2]*deltai+0.5;                                            \
  glColor3f(r, g, b);                                                   \
  glNormal3f(surfx[n2], surfy[n2], surfz[n2]);                          \
  glVertex3d(x[n2], y[n2], z[n2]);                                      \
  r = f1[n3]*deltai+0.5;                                            \
  g = f2[n3]*deltai+0.5;                                            \
  b = f3[n3]*deltai+0.5;                                            \
  glColor3f(r, g, b);                                                   \
  glNormal3f(surfx[n3], surfy[n3], surfz[n3]);                          \
  glVertex3d(x[n3], y[n3], z[n3]);                                      \
  r = f1[n4]*deltai+0.5;                                            \
  g = f2[n4]*deltai+0.5;                                            \
  b = f3[n4]*deltai+0.5;                                            \
  glColor3f(r, g, b);                                                   \
  glNormal3f(surfx[n4], surfy[n4], surfz[n4]);                          \
  glVertex3d(x[n4], y[n4], z[n4]);

#define PLOTQUADQ2 r = f1[n1]*deltai+0.5;                            \
  g = f2[n1]*deltai+0.5;                                            \
  b = f3[n1]*deltai+0.5;                                            \
  glColor3f(r, g, b);                                                   \
  glNormal3f(surfx[ff], surfy[ff], surfz[ff]);                          \
  glVertex3d(x[n1], y[n1], z[n1]);                                      \
  r = f1[n2]*deltai+0.5;                                            \
  g = f2[n2]*deltai+0.5;                                            \
  b = f3[n2]*deltai+0.5;                                            \
  glColor3f(r, g, b);                                                   \
  glVertex3d(x[n2], y[n2], z[n2]);                                      \
  r = f1[n3]*deltai+0.5;                                            \
  g = f2[n3]*deltai+0.5;                                            \
  b = f3[n3]*deltai+0.5;                                            \
  glColor3f(r, g, b);                                                   \
  glVertex3d(x[n3], y[n3], z[n3]);                                      \
  r = f1[n4]*deltai+0.5;                                            \
  g = f2[n4]*deltai+0.5;                                            \
  b = f3[n4]*deltai+0.5;                                            \
  glColor3f(r, g, b);                                                   \
  glVertex3d(x[n4], y[n4], z[n4]);

#define PLOTQUADT2 r = f1[n1]*deltai+0.5;                            \
  g = f2[n1]*deltai+0.5;                                            \
  b = f3[n1]*deltai+0.5;                                            \
  glColor3f(r, g, b);                                                   \
  glNormal3f(surfx[ff], surfy[ff], surfz[ff]);                          \
  glVertex3d(x[n1], y[n1], z[n1]);                                      \
  r = f1[n2]*deltai+0.5;                                            \
  g = f2[n2]*deltai+0.5;                                            \
  b = f3[n2]*deltai+0.5;                                            \
  glColor3f(r, g, b);                                                   \
  glVertex3d(x[n2], y[n2], z[n2]);                                      \
  r = f1[n4]*deltai+0.5;                                            \
  g = f2[n4]*deltai+0.5;                                            \
  b = f3[n4]*deltai+0.5;                                            \
  glColor3f(r, g, b);                                                   \
  glVertex3d(x[n4], y[n4], z[n4]);                                      \
  r = f1[n2]*deltai+0.5;                            \
  g = f2[n2]*deltai+0.5;                                            \
  b = f3[n2]*deltai+0.5;                                            \
  glColor3f(r, g, b);                                                   \
  glNormal3f(surfx[ff], surfy[ff], surfz[ff]);                          \
  glVertex3d(x[n2], y[n2], z[n2]);                                      \
  r = f1[n3]*deltai+0.5;                                            \
  g = f2[n3]*deltai+0.5;                                            \
  b = f3[n3]*deltai+0.5;                                            \
  glColor3f(r, g, b);                                                   \
  glVertex3d(x[n3], y[n3], z[n3]);                                      \
  r = f1[n4]*deltai+0.5;                                            \
  g = f2[n4]*deltai+0.5;                                            \
  b = f3[n4]*deltai+0.5;                                            \
  glColor3f(r, g, b);                                                   \
  glVertex3d(x[n4], y[n4], z[n4]);

#define PLOTQUADB ret1 = _pref.blanking->f(this, n1, zonep->blank, zonet); \
  ret2 = _pref.blanking->f(this, n2, zonep->blank, zonet);               \
  ret3 = _pref.blanking->f(this, n3, zonep->blank, zonet);               \
  ret4 = _pref.blanking->f(this, n4, zonep->blank, zonet);               \
  if (ret1*ret2*ret3*ret4 != 0) { PLOTQUAD; }

#define PLOTQUAD2B ret1 = _pref.blanking->f(this, n1, zonep->blank, zonet); \
  ret2 = _pref.blanking->f(this, n2, zonep->blank, zonet);               \
  ret3 = _pref.blanking->f(this, n3, zonep->blank, zonet);               \
  ret4 = _pref.blanking->f(this, n4, zonep->blank, zonet);               \
  if (ret1*ret2*ret3*ret4 != 0) { PLOTQUAD2; }

#define PLOTNGON(n) r = f1[n]*deltai+0.5;                           \
  g = f2[n]*deltai+0.5;                                             \
  b = f3[n]*deltai+0.5;                                             \
  glColor3f(r, g, b);                                                   \
  glVertex3d(x[n], y[n], z[n]);
  
  double deltai1 = MAX(ABS(fmax1), 1.e-6);
  deltai1 = MAX(deltai1, ABS(fmin1));
  double deltai2 = MAX(ABS(fmax2), 1.e-6);
  deltai2 = MAX(deltai2, ABS(fmin2));
  double deltai3 = MAX(ABS(fmax3), 1.e-6);
  deltai2 = MAX(deltai3, ABS(fmin3));

  double deltai = deltai1;
  deltai = MAX(deltai, deltai2);
  deltai = MAX(deltai, deltai3);
  deltai = 0.5/deltai;

  // Colormap
  void (*getrgb)(Data* data, double, float*, float*, float*);
  getrgb = _pref.colorMap->f;

  int ne = zonep->ne;
  int ne2 = 2*ne; int ne3 = 3*ne;
  int ne4 = 4*ne; int ne5 = 5*ne;
  int np = zonep->np;

  double* x = zonep->x; double* y = zonep->y; double* z = zonep->z;
  int* connect = zonep->connect;

  if (zonep->eltType == 2) // TRI
  {
    float* surfx = zonep->surf;
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
  else if (zonep->eltType == 3) // QUAD
  {
    float* surfx = zonep->surf;
    float* surfy = surfx + np;
    float* surfz = surfy + np;
    
    glBegin(GL_QUADS_ARE);
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
  else if (zonep->eltType == 4) // TETRA
  {
    float* surfx = zonep->surf;
    float* surfy = surfx + ne4;
    float* surfz = surfy + ne4;
    
    glBegin(GL_TRIANGLES);
    if (zonep->blank == 0)
    {
      for (i = 0; i < ne; i++)
      {
        n1 = connect[i]-1;
        n2 = connect[i+ne]-1;
        n3 = connect[i+ne2]-1;
        ff = 4*i;
        PLOTTRI2;
        n1 = connect[i]-1;
        n2 = connect[i+ne]-1;
        n3 = connect[i+ne*3]-1;
        ff = 4*i+1;
        PLOTTRI2;
        n1 = connect[i+ne]-1;
        n2 = connect[i+ne2]-1;
        n3 = connect[i+ne*3]-1;
        ff = 4*i + 2;
        PLOTTRI2;
        n1 = connect[i+ne2]-1;
        n2 = connect[i]-1;
        n3 = connect[i+ne*3]-1;
        ff = 4*i+3;
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
        ff = 4*i;
        PLOTTRI2B;
        n1 = connect[i]-1;
        n2 = connect[i+ne]-1;
        n3 = connect[i+ne*3]-1;
        ff = 4*i+1;
        PLOTTRI2B;
        n1 = connect[i+ne]-1;
        n2 = connect[i+ne2]-1;
        n3 = connect[i+ne*3]-1;
        ff = 4*i+2;
        PLOTTRI2B;
        n1 = connect[i+ne2]-1;
        n2 = connect[i]-1;
        n3 = connect[i+ne*3]-1;
        ff = 4*i+3;
        PLOTTRI2B;
      }
    }
    glEnd();
  }
  else if (zonep->eltType == 5) // PENTA
  {
    float* surfx = zonep->surf;
    float* surfy = surfx + ne5;
    float* surfz = surfy + ne5;
    
    glBegin(GL_TRIANGLES);
    if (zonep->blank == 0)
    {
      for (i = 0; i < ne; i++)
      {
        n1 = connect[i]-1;
        n2 = connect[i+ne]-1;
        n3 = connect[i+ne2]-1;
        ff = 5*i;
        PLOTTRI2;
        n1 = connect[i+ne*3]-1;
        n2 = connect[i+ne*4]-1;
        n3 = connect[i+ne*5]-1;
        ff = 5*i+1;
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
        ff = 5*i;
        PLOTTRI2B;
        n1 = connect[i+ne*3]-1;
        n2 = connect[i+ne*4]-1;
        n3 = connect[i+ne*5]-1;
        ff = 5*i+1;
        PLOTTRI2B;
      }
    }
    glEnd();
    
    glBegin(GL_QUADS_ARE);
    if (zonep->blank == 0)
    {
      for (i = 0; i < ne; i++)
      {
        n1 = connect[i]-1;
        n2 = connect[i+ne]-1;
        n3 = connect[i+ne*4]-1;
        n4 = connect[i+ne*3]-1;
        ff = 5*i+2;
        PLOTQUAD2;
        n1 = connect[i+ne]-1;
        n2 = connect[i+ne*2]-1;
        n3 = connect[i+ne*5]-1;
        n4 = connect[i+ne*4]-1;
        ff = 5*i+3;
        PLOTQUAD2;
        n1 = connect[i]-1;
        n2 = connect[i+ne2]-1;
        n3 = connect[i+ne*5]-1;
        n4 = connect[i+ne*3]-1;
        ff = 5*i+4;
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
        ff = 5*i+2;
        PLOTQUAD2B;
        n1 = connect[i+ne]-1;
        n2 = connect[i+ne*2]-1;
        n3 = connect[i+ne*5]-1;
        n4 = connect[i+ne*4]-1;
        ff = 5*i+3;
        PLOTQUAD2B;
        n1 = connect[i]-1;
        n2 = connect[i+ne2]-1;
        n3 = connect[i+ne*5]-1;
        n4 = connect[i+ne*3]-1;
        ff = 5*i+4;
        PLOTQUAD2B;
      }
    }
    glEnd();
  } 
  else if (zonep->eltType == 6) // PYRA
  {
    float* surfx = zonep->surf;
    float* surfy = surfx + ne5;
    float* surfz = surfy + ne5;
    
    glBegin(GL_TRIANGLES);
    if (zonep->blank == 0)
    {
      for (i = 0; i < ne; i++)
      {
        n1 = connect[i]-1;
        n2 = connect[i+ne]-1;
        n3 = connect[i+ne*4]-1;
        ff = 5*i;
        PLOTTRI2;
        n1 = connect[i+ne]-1;
        n2 = connect[i+ne*2]-1;
        n3 = connect[i+ne*4]-1;
        ff = 5*i+1;
        PLOTTRI2;
        n1 = connect[i+ne2]-1;
        n2 = connect[i+ne*3]-1;
        n3 = connect[i+ne*4]-1;
        ff = 5*i+2;
        PLOTTRI2;
        n1 = connect[i+ne*3]-1;
        n2 = connect[i]-1;
        n3 = connect[i+ne*4]-1;
        ff = 5*i+3;
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
        ff = 5*i;
        PLOTTRI2B;
        n1 = connect[i+ne]-1;
        n2 = connect[i+ne*2]-1;
        n3 = connect[i+ne*4]-1;
        ff = 5*i+1;
        PLOTTRI2B;
        n1 = connect[i+ne2]-1;
        n2 = connect[i+ne*3]-1;
        n3 = connect[i+ne*4]-1;
        ff = 5*i+2;
        PLOTTRI2B;
        n1 = connect[i+ne*3]-1;
        n2 = connect[i]-1;
        n3 = connect[i+ne*4]-1;
        ff = 5*i+3;
        PLOTTRI2B;
      }
    }
    glEnd();
    
    glBegin(GL_QUADS_ARE);
    if (zonep->blank == 0)
    {
      for (i = 0; i < ne; i++)
      {
        n1 = connect[i]-1;
        n2 = connect[i+ne]-1;
        n3 = connect[i+ne*2]-1;
        n4 = connect[i+ne*3]-1;
        ff = 5*i+4;
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
        ff = 5*i+4;
        PLOTQUAD2B;
      }
    }
    glEnd();
  } 
  else if (zonep->eltType == 7) // HEXA
  {
    float* surfx = zonep->surf;
    float* surfy = surfx + 6*ne;
    float* surfz = surfy + 6*ne;
    
    glBegin(GL_QUADS_ARE);      
    if (zonep->blank == 0)
    {
      for (i = 0; i < ne; i++)
      {
        n1 = connect[i]-1;
        n2 = connect[i+ne]-1;
        n3 = connect[i+ne2]-1;
        n4 = connect[i+ne3]-1;
        ff = 6*i;
        PLOTQUAD2;
        n1 = connect[i+ne]-1;
        n2 = connect[i+ne5]-1;
        n3 = connect[i+6*ne]-1;
        n4 = connect[i+ne2]-1;
        ff = 6*i+1;
        PLOTQUAD2;
        n1 = connect[i+ne5]-1;
        n2 = connect[i+ne4]-1;
        n3 = connect[i+7*ne]-1;
        n4 = connect[i+6*ne]-1;
        ff = 6*i+2;
        PLOTQUAD2;
        n1 = connect[i]-1;
        n2 = connect[i+ne4]-1;
        n3 = connect[i+7*ne]-1;
        n4 = connect[i+3*ne]-1;
        ff = 6*i+3;
        PLOTQUAD2;
        n1 = connect[i]-1;
        n2 = connect[i+ne]-1;
        n3 = connect[i+ne5]-1;
        n4 = connect[i+ne4]-1;
        ff = 6*i+4;
        PLOTQUAD2;
        n1 = connect[i+3*ne]-1;
        n2 = connect[i+2*ne]-1;
        n3 = connect[i+6*ne]-1;
        n4 = connect[i+7*ne]-1;
        ff = 6*i+5;
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
        ff = 6*i;
        PLOTQUAD2B;
        n1 = connect[i+ne]-1;
        n2 = connect[i+ne5]-1;
        n3 = connect[i+6*ne]-1;
        n4 = connect[i+ne2]-1;
        ff = 6*i+1;
        PLOTQUAD2B;
        n1 = connect[i+ne5]-1;
        n2 = connect[i+ne4]-1;
        n3 = connect[i+7*ne]-1;
        n4 = connect[i+6*ne]-1;
        ff = 6*i+2;
        PLOTQUAD2B;
        n1 = connect[i]-1;
        n2 = connect[i+ne4]-1;
        n3 = connect[i+7*ne]-1;
        n4 = connect[i+3*ne]-1;
        ff = 6*i+3;
        PLOTQUAD2B;
        n1 = connect[i]-1;
        n2 = connect[i+ne]-1;
        n3 = connect[i+ne5]-1;
        n4 = connect[i+ne4]-1;
        ff = 6*i+4;
        PLOTQUAD2B;
        n1 = connect[i+3*ne]-1;
        n2 = connect[i+2*ne]-1;
        n3 = connect[i+6*ne]-1;
        n4 = connect[i+7*ne]-1;
        ff = 6*i+5;
        PLOTQUAD2B;
      }
    }
    glEnd();
  }
  else if (zonep->eltType == 10) // NGON
  {
    int nf = connect[0];
    int l, nd, c;
    float* surfx = zonep->surf;
    float* surfy = surfx + nf;
    float* surfz = surfy + nf;
    int next, prev;
    
    if (zonep->blank == 0)
    {
      // Faces des elements 3D
      c = 2;
      glBegin(GL_TRIANGLES);
      for (i = 0; i < nf; i++)
      {
        nd = connect[c]; 
        if (nd == 3) // TRI
        {
          n1 = connect[c+1]-1;
          n2 = connect[c+2]-1;
          n3 = connect[c+3]-1;
          PLOTTRI2;
        }
        c += nd+1;
      }
      glEnd();

      c = 2;
      glBegin(GL_QUADS_ARE); 
      for (i = 0; i < nf; i++)
      {
        nd = connect[c];
        if (nd == 4) // QUAD
        {
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
        nd = connect[c]; // nbre de noeuds de la face
        if (nd > 4) // elt 3D
        {
          glNormal3f(surfx[i], surfy[i], surfz[i]);
          glBegin(GL_POLYGON);
          for (l = 0; l < nd; l++)
          {
            n1 = connect[c+l+1]-1;
            PLOTNGON(n1);
          }
          glEnd();
        }
        c += nd+1;
      }

      // Elements 2D
      for (i = 0; i < zonep->nelts2D; i++)
      {
        glBegin(GL_POLYGON);
        int elt = zonep->posElts2D[i];
        int* ptrelt = &connect[elt];
        int nf = ptrelt[0];
        int drawn = 0;
        int j, first;

        int face = ptrelt[1]-1;
        glNormal3f(surfx[face], surfy[face], surfz[face]);
        int* ptrface = &connect[zonep->posFaces[face]];
        n1 = ptrface[1]-1; first = n1;
        n2 = ptrface[2]-1;
        PLOTNGON(n1);
        PLOTNGON(n2);
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
            { PLOTNGON(n2);
              prev = n1; next = n2; drawn++; break; }
            else if (n2 == next && n1 != prev)
            { 
              PLOTNGON(n1);
              prev = n2; next = n1; drawn++; break; }
          }
          if (j == nf+1) drawn++; // pour eviter les boucles infinies
        }
        if (next != first) glVertex3d(x[first], y[first], z[first]); // force close
        glEnd();
      }
    }
    else // blanking
    {
      // Faces des elements 3D
      c = 2;
      glBegin(GL_TRIANGLES);
      for (i = 0; i < nf; i++)
      {
        nd = connect[c]; 
        if (nd == 3) // TRI
        {
          n1 = connect[c+1]-1;
          n2 = connect[c+2]-1;
          n3 = connect[c+3]-1;
          PLOTTRI2B;
        }
        c += nd+1;
      }
      glEnd();

      c = 2;
      glBegin(GL_QUADS_ARE);
      for (i = 0; i < nf; i++)
      {
        nd = connect[c];
        if (nd == 4) // QUAD
        {
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
        nd = connect[c]; // nbre de noeuds de la face
        if (nd > 2) // elt 3D
        {
          int blank = 0;
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
              PLOTNGON(n1);
            }
            glEnd();
          }
        }
        c += nd+1;
      }

      // Elements 2D
      for (i = 0; i < zonep->nelts2D; i++)
      {
        int elt = zonep->posElts2D[i];
        int* ptrelt = &connect[elt];
        int nf = ptrelt[0];

        int blank = 0;
        for (int j = 1; j <= nf; j++)
        {
          int face = ptrelt[1]-1;
          int* ptrface = &connect[zonep->posFaces[face]];
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
          int elt = zonep->posElts2D[i];
          int* ptrelt = &connect[elt];
          int nf = ptrelt[0];
          int drawn = 0;
          int j, first;

          int face = ptrelt[1]-1;
          glNormal3f(surfx[face], surfy[face], surfz[face]);
          int* ptrface = &connect[zonep->posFaces[face]];
          n1 = ptrface[1]-1; first = n1;
          n2 = ptrface[2]-1;
          PLOTNGON(n1); 
          PLOTNGON(n2);
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
                PLOTNGON(n2);
                prev = n1; next = n2; drawn++; break; }
              else if (n2 == next && n1 != prev)
              { 
                PLOTNGON(n1);
                prev = n2; next = n1; drawn++; break; }
            }
            if (j == nf+1) drawn++; // pour eviter les boucles infinies
          }
          if (next != first) glVertex3d(x[first], y[first], z[first]); // force close
          glEnd();
        }
      }
    }
  }

  // Pour les BAR
  if (zonep->eltType == 1)
  {
    glLineWidth(3.);
    glPolygonOffset(-1.,-10.); // force offset
    glBegin(GL_LINES);
    if (zonep->blank == 0)
    {
      for (i = 0; i < ne; i++)
      {
        n1 = connect[i]-1;
        n2 = connect[i+ne]-1;
        r = (f1[n1]-fmin1)*deltai;         
        g = (f2[n1]-fmin2)*deltai;
        b = (f3[n1]-fmin3)*deltai; 
        glColor3f(r, g, b+offb); 
        glVertex3d(x[n1], y[n1], z[n1]);
        r = (f1[n2]-fmin1)*deltai;         
        g = (f2[n2]-fmin2)*deltai;
        b = (f3[n2]-fmin3)*deltai; 
        glColor3f(r, g, b+offb); 
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
          r = (f1[n1]-fmin1)*deltai;         
          g = (f2[n1]-fmin2)*deltai;
          b = (f3[n1]-fmin3)*deltai; 
          glColor3f(r, g, b+offb); 
          glVertex3d(x[n1], y[n1], z[n1]);
          r = (f1[n2]-fmin1)*deltai;         
          g = (f2[n2]-fmin2)*deltai;
          b = (f3[n2]-fmin3)*deltai; 
          glColor3f(r, g, b+offb); 
          glVertex3d(x[n2], y[n2], z[n2]);
        }
      }
    }
    glEnd();
    glLineWidth(1.);
  }

  // Pour les NGONS 1D
  if (zonep->eltType == 10 && zonep->nelts1D > 0)
  {
    glLineWidth(3.);
    glBegin(GL_LINES);
    if (zonep->blank == 0)
    {
      for (i = 0; i < zonep->nelts1D; i++)
      {
        int elt = zonep->posElts1D[i];
        int* ptrelt = &connect[elt];
        int face1 = ptrelt[1]-1;
        int face2 = ptrelt[2]-1;
        int posface1 = zonep->posFaces[face1];
        int posface2 = zonep->posFaces[face2];
        n1 = connect[posface1+1]-1;
        n2 = connect[posface2+1]-1;
        r = (f1[n1]-fmin1)*deltai;         
        g = (f2[n1]-fmin2)*deltai;
        b = (f3[n1]-fmin3)*deltai; 
        glColor3f(r, g, b+offb); 
        glVertex3d(x[n1], y[n1], z[n1]);
        r = (f1[n2]-fmin1)*deltai;         
        g = (f2[n2]-fmin2)*deltai;
        b = (f3[n2]-fmin3)*deltai; 
        glColor3f(r, g, b+offb); 
        glVertex3d(x[n2], y[n2], z[n2]);
      }
    }
    else
    {
      for (i = 0; i < zonep->nelts1D; i++)
      {
        int elt = zonep->posElts1D[i];
        int* ptrelt = &connect[elt];
        int face1 = ptrelt[1]-1;
        int face2 = ptrelt[2]-1;
        int posface1 = zonep->posFaces[face1];
        int posface2 = zonep->posFaces[face2];
        n1 = connect[posface1+1]-1;
        n2 = connect[posface2+1]-1;
        ret1 = _pref.blanking->f(this, n1, zonep->blank, zonet);
        ret2 = _pref.blanking->f(this, n2, zonep->blank, zonet);
        if (ret1*ret2 != 0)
        {
          r = (f1[n1]-fmin1)*deltai;         
          g = (f2[n1]-fmin2)*deltai;
          b = (f3[n1]-fmin3)*deltai; 
          glColor3f(r, g, b+offb); 
          glVertex3d(x[n1], y[n1], z[n1]);
          r = (f1[n2]-fmin1)*deltai;         
          g = (f2[n2]-fmin2)*deltai;
          b = (f3[n2]-fmin3)*deltai; 
          glColor3f(r, g, b+offb);
          glVertex3d(x[n2], y[n2], z[n2]);
        }
      }
    }
    glEnd();
    glLineWidth(1.);
  }

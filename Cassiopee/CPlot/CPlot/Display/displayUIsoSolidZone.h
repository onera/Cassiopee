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
#define MAXT2 \
  g = fabs(f[n2]-f[n1])*deltai;                                         \
  b = fabs(f[n3]-f[n1])*deltai;                                         \
  g = g * invsqrt((x[n2]-x[n1])*(x[n2]-x[n1])+(y[n2]-y[n1])*(y[n2]-y[n1])+(z[n2]-z[n1])*(z[n2]-z[n1])); \
  b = b * invsqrt((x[n3]-x[n1])*(x[n3]-x[n1])+(y[n3]-y[n1])*(y[n3]-y[n1])+(z[n3]-z[n1])*(z[n3]-z[n1])); \
  g = MAX(g, b);
#define MAXQ2 \
  g = fabs(f[n2]-f[n1])*deltai;                                         \
  b = fabs(f[n3]-f[n1])*deltai;                                         \
  a = fabs(f[n4]-f[n1])*deltai;                                         \
  g = g * invsqrt((x[n2]-x[n1])*(x[n2]-x[n1])+(y[n2]-y[n1])*(y[n2]-y[n1])+(z[n2]-z[n1])*(z[n2]-z[n1])); \
  b = b * invsqrt((x[n3]-x[n1])*(x[n3]-x[n1])+(y[n3]-y[n1])*(y[n3]-y[n1])+(z[n3]-z[n1])*(z[n3]-z[n1])); \
  a = a * invsqrt((x[n4]-x[n1])*(x[n4]-x[n1])+(y[n4]-y[n1])*(y[n4]-y[n1])+(z[n4]-z[n1])*(z[n4]-z[n1])); \
  g = MAX(g, b);                                                        \
  g = MAX(g, a);

#ifndef __SHADERS__
#define PLOTTRI getrgb(this, (f[n1]-fmin)*deltai, &r, &g, &b);          \
  glColor3f(r, g, b+offb);                                              \
  glNormal3f(surfx[n1], surfy[n1], surfz[n1]);                          \
  glVertex3d(x[n1], y[n1], z[n1]);                                      \
  getrgb(this, (f[n2]-fmin)*deltai, &r, &g, &b);                        \
  glColor3f(r, g, b+offb);                                              \
  glNormal3f(surfx[n2], surfy[n2], surfz[n2]);                          \
  glVertex3d(x[n2], y[n2], z[n2]);                                      \
  getrgb(this, (f[n3]-fmin)*deltai, &r, &g, &b);                        \
  glColor3f(r, g, b+offb);                                              \
  glNormal3f(surfx[n3], surfy[n3], surfz[n3]);                          \
  glVertex3d(x[n3], y[n3], z[n3]);                                   

#define PLOTTRI2 getrgb(this, (f[n1]-fmin)*deltai, &r, &g, &b);         \
  glColor3f(r, g, b+offb);                                              \
  glNormal3f(surfx[ff], surfy[ff], surfz[ff]);                          \
  glVertex3d(x[n1], y[n1], z[n1]);                                      \
  getrgb(this, (f[n2]-fmin)*deltai, &r, &g, &b);                        \
  glColor3f(r, g, b+offb);                                              \
  glVertex3d(x[n2], y[n2], z[n2]);                                      \
  getrgb(this, (f[n3]-fmin)*deltai, &r, &g, &b);                        \
  glColor3f(r, g, b+offb);                                              \
  glVertex3d(x[n3], y[n3], z[n3]);                
#else
#define PLOTTRI r = (f[n1]-fmin)*deltai;                                \
  glColor3f(r, 0.f, 0.f);                                               \
  glNormal3f(surfx[n1], surfy[n1], surfz[n1]);                          \
  glVertex3d(x[n1], y[n1], z[n1]);                                      \
  r = (f[n2]-fmin)*deltai;                                              \
  glColor3f(r, 0.f, 0.f);                                               \
  glNormal3f(surfx[n2], surfy[n2], surfz[n2]);                          \
  glVertex3d(x[n2], y[n2], z[n2]);                                      \
  r = (f[n3]-fmin)*deltai;                                              \
  glColor3f(r, 0.f, 0.f);                                               \
  glNormal3f(surfx[n3], surfy[n3], surfz[n3]);                          \
  glVertex3d(x[n3], y[n3], z[n3]);                                   

#define PLOTTRI2 r =(f[n1]-fmin)*deltai;                                \
  glColor3f(r, 0.f, 0.f);                                               \
  glNormal3f(surfx[ff], surfy[ff], surfz[ff]);                          \
  glVertex3d(x[n1], y[n1], z[n1]);                                      \
  r = (f[n2]-fmin)*deltai;                                              \
  glColor3f(r, 0.f, 0.f);                                               \
  glVertex3d(x[n2], y[n2], z[n2]);                                      \
  r = (f[n3]-fmin)*deltai;                                              \
  glColor3f(r, 0.f, 0.f);                                               \
  glVertex3d(x[n3], y[n3], z[n3]); 
#endif

#define PLOTTRIB ret1 = _pref.blanking->f(this, n1, zonep->blank, zonet); \
  ret2 = _pref.blanking->f(this, n2, zonep->blank, zonet);               \
  ret3 = _pref.blanking->f(this, n3, zonep->blank, zonet);               \
  if (ret1*ret2*ret3 != 0) { PLOTTRI; }

#define PLOTTRI2B ret1 = _pref.blanking->f(this, n1, zonep->blank, zonet); \
  ret2 = _pref.blanking->f(this, n2, zonep->blank, zonet);               \
  ret3 = _pref.blanking->f(this, n3, zonep->blank, zonet);               \
  if (ret1*ret2*ret3 != 0) { PLOTTRI2; }

#ifndef __SHADERS__
#define PLOTQUAD getrgb(this, (f[n1]-fmin)*deltai, &r, &g, &b);         \
  glColor3f(r, g, b+offb);                                              \
  glNormal3f(surfx[n1], surfy[n1], surfz[n1]);                          \
  glVertex3d(x[n1], y[n1], z[n1]);                                      \
  getrgb(this, (f[n2]-fmin)*deltai, &r, &g, &b);                        \
  glColor3f(r, g, b+offb);                                              \
  glNormal3f(surfx[n2], surfy[n2], surfz[n2]);                          \
  glVertex3d(x[n2], y[n2], z[n2]);                                      \
  getrgb(this, (f[n3]-fmin)*deltai, &r, &g, &b);                        \
  glColor3f(r, g, b+offb);                                              \
  glNormal3f(surfx[n3], surfy[n3], surfz[n3]);                          \
  glVertex3d(x[n3], y[n3], z[n3]);                                      \
  getrgb(this, (f[n4]-fmin)*deltai, &r, &g, &b);                        \
  glColor3f(r, g, b+offb);                                              \
  glNormal3f(surfx[n4], surfy[n4], surfz[n4]);                          \
  glVertex3d(x[n4], y[n4], z[n4]);                                    

#define PLOTQUAD2 getrgb(this, (f[n1]-fmin)*deltai, &r, &g, &b);        \
  glColor3f(r, g, b+offb);                                              \
  glNormal3f(surfx[ff], surfy[ff], surfz[ff]);                          \
  glVertex3d(x[n1], y[n1], z[n1]);                                      \
  getrgb(this, (f[n2]-fmin)*deltai, &r, &g, &b);                        \
  glColor3f(r, g, b+offb);                                              \
  glVertex3d(x[n2], y[n2], z[n2]);                                      \
  getrgb(this, (f[n3]-fmin)*deltai, &r, &g, &b);                        \
  glColor3f(r, g, b+offb);                                              \
  glVertex3d(x[n3], y[n3], z[n3]);                                      \
  getrgb(this, (f[n4]-fmin)*deltai, &r, &g, &b);                        \
  glColor3f(r, g, b+offb);                                              \
  glVertex3d(x[n4], y[n4], z[n4]);
#else
#define PLOTQUAD r = (f[n1]-fmin)*deltai;                               \
  glColor3f(r, 0.f, 0.f);                                               \
  glNormal3f(surfx[n1], surfy[n1], surfz[n1]);                          \
  glVertex3d(x[n1], y[n1], z[n1]);                                      \
  r = (f[n2]-fmin)*deltai;                                              \
  glColor3f(r, 0.f, 0.f);                                               \
  glNormal3f(surfx[n2], surfy[n2], surfz[n2]);                          \
  glVertex3d(x[n2], y[n2], z[n2]);                                      \
  r = (f[n3]-fmin)*deltai;                                              \
  glColor3f(r, 0.f, 0.f);                                               \
  glNormal3f(surfx[n3], surfy[n3], surfz[n3]);                          \
  glVertex3d(x[n3], y[n3], z[n3]);                                      \
  r = (f[n4]-fmin)*deltai;                                              \
  glColor3f(r, 0.f, 0.f);                                               \
  glNormal3f(surfx[n4], surfy[n4], surfz[n4]);                          \
  glVertex3d(x[n4], y[n4], z[n4]);                                    

#define PLOTQUAD2 r = (f[n1]-fmin)*deltai;                              \
  glColor3f(r, 0.f, 0.f);                                               \
  glNormal3f(surfx[ff], surfy[ff], surfz[ff]);                          \
  glVertex3d(x[n1], y[n1], z[n1]);                                      \
  r = (f[n2]-fmin)*deltai;                                              \
  glColor3f(r, 0.f, 0.f);                                               \
  glVertex3d(x[n2], y[n2], z[n2]);                                      \
  r = (f[n3]-fmin)*deltai;                                              \
  glColor3f(r, 0.f, 0.f);                                               \
  glVertex3d(x[n3], y[n3], z[n3]);                                      \
  r = (f[n4]-fmin)*deltai;                                              \
  glColor3f(r, 0.f, 0.f);                                               \
  glVertex3d(x[n4], y[n4], z[n4]);
#endif

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

#ifndef __SHADERS__
#define PLOTNGON(n) getrgb(this, (f[n]-fmin)*deltai, &r, &g, &b);       \
  glColor3f(r, g, b+offb);                                              \
  glVertex3d(x[n], y[n], z[n]); 
#define PLOTNGON2(n) getrgb(this, (f[n]-fmin)*deltai, &r, &g, &b);      \
  glColor3f(r, g, b+offb);                                              \
  glNormal3f(surfx[n], surfy[n], surfz[n]);                             \
  glVertex3d(x[n], y[n], z[n]);
  #else
#define PLOTNGON(n) r = (f[n]-fmin)*deltai;                             \
  glColor3f(r, 0., 0.);                                                 \
  glVertex3d(x[n], y[n], z[n]);
#define PLOTNGON2(n) r = (f[n]-fmin)*deltai;                            \
  glColor3f(r, 0., 0.);                                                 \
  glNormal3f(surfx[n], surfy[n], surfz[n]);                             \
  glVertex3d(x[n], y[n], z[n]);
 #endif

  double fmin, fmax;
  fmax = maxf[nofield]; fmin = minf[nofield];
  double deltai = MAX(fmax-fmin, ISOCUTOFF);
  deltai = 1./deltai;

  // Colormap
#ifndef __SHADERS__
  void (*getrgb)(Data* data, double, float*, float*, float*);
  getrgb = _pref.colorMap->f;
  double r2, r3, r4, a;
  #endif

  double* f = zonep->f[nofield];

  E_Int np = zonep->np;
  double* x = zonep->x; 
  double* y = zonep->y; 
  double* z = zonep->z;
  
  for (size_t nc = 0; nc < zonep->connect.size(); nc++) {

  E_Int eltType = zonep->eltType[nc];
  E_Int* connect = zonep->connect[nc];
  E_Int ne = zonep->nec[nc];

  E_Int ne2 = 2*ne; E_Int ne3 = 3*ne;
  E_Int ne4 = 4*ne; E_Int ne5 = 5*ne;
  
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
  else if (eltType == 3) // QUAD
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
  else if (eltType == 5) // PENTA
  {
    float* surfp = zonep->surf[nc];
    float* surfx = surfp;
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
    
    glBegin(GL_QUADS);
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
  else if (eltType == 6) // PYRA
  {
    float* surfp = zonep->surf[nc];
    float* surfx = surfp;
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
    
    glBegin(GL_QUADS);
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
  else if (eltType == 7) // HEXA
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
  else if (eltType == 10) // NGON
  {
    E_Int nf = connect[nc];
    E_Int l, nd, c, j;
    float* surfp = zonep->surf[nc];
    float* surfx = surfp;
    float* surfy = surfx + nf;
    float* surfz = surfy + nf;
    E_Int prev, next, first, face;
    E_Int elt, drawn, blank;
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
            n1 = connect[c+1]-1;
            n2 = connect[c+2]-1;
            n3 = connect[c+3]-1;
            glBegin(GL_TRIANGLES);
            glNormal3f(surfx[i], surfy[i], surfz[i]);
            PLOTNGON(n1); PLOTNGON(n2); PLOTNGON(n3);
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
            n1 = connect[c+1]-1;
            n2 = connect[c+2]-1;
            n3 = connect[c+3]-1;
            n4 = connect[c+4]-1;
            glNormal3f(surfx[i], surfy[i], surfz[i]);
            PLOTNGON(n1); PLOTNGON(n2); PLOTNGON(n3); PLOTNGON(n4);    
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
            glBegin(GL_POLYGON);
            glNormal3f(surfx[i], surfy[i], surfz[i]);
            for (l = 0; l < nd; l++)
            {
              n1 = connect[c+l+1]-1;
              PLOTNGON(n1);
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
          drawn = 0;
          if (nf == 3 || nf == 4) continue;
          glBegin(GL_POLYGON);
          
          face = ptrelt[1]-1;
          //glNormal3f(surfx[face], surfy[face], surfz[face]);
          ptrface = &connect[zonep->posFaces[face]];
          n1 = ptrface[1]-1; first = n1;
          n2 = ptrface[2]-1;
          PLOTNGON2(n1);
          PLOTNGON2(n2);
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
              { PLOTNGON2(n2);
                prev = n1; next = n2; drawn++; break; }
              else if (n2 == next && n1 != prev)
              { PLOTNGON2(n1);
                prev = n2; next = n1; drawn++; break; }
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
            blank = 0;
            for (l = 0; l < nd; l++)
            {
              n1 = connect[c+l+1]-1;
              if (_pref.blanking->f(this, n1, zonep->blank, zonet) == 0)
              { blank = 1; break; }
            }
            if (blank == 0)
            {
              n1 = connect[c+1]-1;
              n2 = connect[c+2]-1;
              n3 = connect[c+3]-1;
              glNormal3f(surfx[i], surfy[i], surfz[i]);
              PLOTNGON(n1); PLOTNGON(n2); PLOTNGON(n3);
            }
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
            blank = 0;
            for (l = 0; l < nd; l++)
            {
              n1 = connect[c+l+1]-1;
              if (_pref.blanking->f(this, n1, zonep->blank, zonet) == 0)
              { blank = 1; break; }
            }
            if (blank == 0)
            {
              n1 = connect[c+1]-1;
              n2 = connect[c+2]-1;
              n3 = connect[c+3]-1;
              n4 = connect[c+4]-1;
              glNormal3f(surfx[i], surfy[i], surfz[i]);
              PLOTNGON(n1); PLOTNGON(n2); PLOTNGON(n3); PLOTNGON(n4);
            }
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
              glBegin(GL_POLYGON);
              glNormal3f(surfx[i], surfy[i], surfz[i]);
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
            PLOTNGON2(n1);
            PLOTNGON2(n2);
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
                { PLOTNGON2(n2);
                  prev = n1; next = n2; drawn++; break; }
                else if (n2 == next && n1 != prev)
                { 
                  PLOTNGON2(n1);
                  prev = n2; next = n1; drawn++; break; }
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
    glLineWidth(3.);
    glPolygonOffset(-1.,-10.); // force offset
    glBegin(GL_LINES);
    if (zonep->blank == 0)
    {
      for (i = 0; i < ne; i++)
      {
        n1 = connect[i]-1;
        n2 = connect[i+ne]-1;
        #ifndef __SHADERS__
          getrgb(this, (f[n1]-fmin)*deltai, &r, &g, &b);
          glColor3f(r, g, b+offb);
          glVertex3d(x[n1], y[n1], z[n1]);
          getrgb(this, (f[n2]-fmin)*deltai, &r, &g, &b);
          glColor3f(r, g, b+offb); 
          glVertex3d(x[n2], y[n2], z[n2]);
        #else
          r =(f[n1]-fmin)*deltai;
          glColor3f(r, 0.f, 0.f);
          glVertex3d(x[n1], y[n1], z[n1]);
          r =(f[n2]-fmin)*deltai;
          glColor3f(r, 0.f, 0.f);
          glVertex3d(x[n2], y[n2], z[n2]);
        #endif
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
          #ifndef __SHADERS__
            getrgb(this, (f[n1]-fmin)*deltai, &r, &g, &b);
            glColor3f(r, g, b+offb); 
            glVertex3d(x[n1], y[n1], z[n1]);
            getrgb(this, (f[n2]-fmin)*deltai, &r, &g, &b);
            glColor3f(r, g, b+offb); 
            glVertex3d(x[n2], y[n2], z[n2]);
          #else
            r =(f[n1]-fmin)*deltai;
            glColor3f(r, 0.f, 0.f);
            glVertex3d(x[n1], y[n1], z[n1]);
            r =(f[n2]-fmin)*deltai;
            glColor3f(r, 0.f, 0.f);
            glVertex3d(x[n2], y[n2], z[n2]);
          #endif
        }
      }
    }
    glEnd();
    glLineWidth(1.);
  }

  // Pour les NGONS 1D
  if (eltType == 10 && zonep->nelts1D > 0)
  {
    glLineWidth(3.);
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
        #ifndef __SHADERS__
          getrgb(this, (f[n1]-fmin)*deltai, &r, &g, &b);
          glColor3f(r, g, b+offb);
          glVertex3d(x[n1], y[n1], z[n1]);
          getrgb(this, (f[n2]-fmin)*deltai, &r, &g, &b);
          glColor3f(r, g, b+offb); 
          glVertex3d(x[n2], y[n2], z[n2]);
        #else
          r =(f[n1]-fmin)*deltai;
          glColor3f(r, 0.f, 0.f);
          glVertex3d(x[n1], y[n1], z[n1]);
          r =(f[n2]-fmin)*deltai;
          glColor3f(r, 0.f, 0.f);
          glVertex3d(x[n2], y[n2], z[n2]);
        #endif
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
          #ifndef __SHADERS__
            getrgb(this, (f[n1]-fmin)*deltai, &r, &g, &b);
            glColor3f(r, g, b+offb); 
            glVertex3d(x[n1], y[n1], z[n1]);
            getrgb(this, (f[n2]-fmin)*deltai, &r, &g, &b);
            glColor3f(r, g, b+offb); 
            glVertex3d(x[n2], y[n2], z[n2]);
          #else
            r =(f[n1]-fmin)*deltai;
            glColor3f(r, 0.f, 0.f);
            glVertex3d(x[n1], y[n1], z[n1]);
            r =(f[n2]-fmin)*deltai;
            glColor3f(r, 0.f, 0.f);
            glVertex3d(x[n2], y[n2], z[n2]);
          #endif
        }
      }
    }
    glEnd();
    glLineWidth(1.);
  }

  } // connects
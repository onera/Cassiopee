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
#ifndef __SHADERS__
#define PLOT getrgb(this, (f[n1]-fmin)*deltai, &r, &g, &b);             \
  glColor3f(r, g, b+offb);                                              \
  glNormal3f(surfx[n1s], surfy[n1s], surfz[n1s]);                       \
  glVertex3d(x[n1], y[n1], z[n1]);                                      \
  getrgb(this, (f[n2]-fmin)*deltai, &r, &g, &b);                        \
  glColor3f(r, g, b+offb);                                              \
  glNormal3f(surfx[n2s], surfy[n2s], surfz[n2s]);                       \
  glVertex3d(x[n2], y[n2], z[n2]);                                      \
  getrgb(this, (f[n4]-fmin)*deltai, &r, &g, &b);                        \
  glColor3f(r, g, b+offb);                                              \
  glNormal3f(surfx[n4s], surfy[n4s], surfz[n4s]);                       \
  glVertex3d(x[n4], y[n4], z[n4]);                                      \
  getrgb(this, (f[n3]-fmin)*deltai, &r, &g, &b);                        \
  glColor3f(r, g, b+offb);                                              \
  glNormal3f(surfx[n3s], surfy[n3s], surfz[n3s]);                       \
  glVertex3d(x[n3], y[n3], z[n3]);
#else // pour les shaders f est envoye dans r
#define PLOTG r = (f[n1]-fmin)*deltai;                                   \
  glColor3f(r, 0., 0.);                                                 \
  glNormal3f(surfx[n1s], surfy[n1s], surfz[n1s]);                       \
  glVertex3d(x[n1], y[n1], z[n1]);                                      \
  r = (f[n2]-fmin)*deltai;                                              \
  glColor3f(r, 0., 0.);                                                 \
  glNormal3f(surfx[n2s], surfy[n2s], surfz[n2s]);                       \
  glVertex3d(x[n2], y[n2], z[n2]);                                      \
  r = (f[n4]-fmin)*deltai;                                              \
  glColor3f(r, 0., 0.);                                                 \
  glNormal3f(surfx[n4s], surfy[n4s], surfz[n4s]);                       \
  glVertex3d(x[n4], y[n4], z[n4]);                                      \
  r = (f[n3]-fmin)*deltai;                                              \
  glColor3f(r, 0., 0.);                                                 \
  glNormal3f(surfx[n3s], surfy[n3s], surfz[n3s]);                       \
  glVertex3d(x[n3], y[n3], z[n3]);

#define MAXN2                                                           \
  g = fabs(f[n2]-f[n1])*deltai;                                         \
  b = fabs(f[n3]-f[n1])*deltai;                                         \
  a = fabs(f[n4]-f[n1])*deltai;                                         \
  g = g * invsqrt((x[n2]-x[n1])*(x[n2]-x[n1])+(y[n2]-y[n1])*(y[n2]-y[n1])+(z[n2]-z[n1])*(z[n2]-z[n1])); \
  b = b * invsqrt((x[n3]-x[n1])*(x[n3]-x[n1])+(y[n3]-y[n1])*(y[n3]-y[n1])+(z[n3]-z[n1])*(z[n3]-z[n1])); \
  a = a * invsqrt((x[n4]-x[n1])*(x[n4]-x[n1])+(y[n4]-y[n1])*(y[n4]-y[n1])+(z[n4]-z[n1])*(z[n4]-z[n1])); \
  g = MAX(g, b);							\
  g = MAX(g, a);							\
  b = 0.;

#define PLOT r = (f[n1]-fmin)*deltai;                                   \
  glColor3f(r, 0.f, 0.f);                                               \
  glNormal3f(surfx[n1s], surfy[n1s], surfz[n1s]);                       \
  glVertex3d(x[n1], y[n1], z[n1]);                                      \
  r = (f[n2]-fmin)*deltai;                                              \
  glColor3f(r, 0.f, 0.f);                                               \
  glNormal3f(surfx[n2s], surfy[n2s], surfz[n2s]);                       \
  glVertex3d(x[n2], y[n2], z[n2]);                                      \
  r = (f[n4]-fmin)*deltai;                                              \
  glColor3f(r, 0.f, 0.f);                                               \
  glNormal3f(surfx[n4s], surfy[n4s], surfz[n4s]);                       \
  glVertex3d(x[n4], y[n4], z[n4]);                                      \
  r = (f[n3]-fmin)*deltai;                                              \
  glColor3f(r, 0.f, 0.f);                                               \
  glNormal3f(surfx[n3s], surfy[n3s], surfz[n3s]);                       \
  glVertex3d(x[n3], y[n3], z[n3]);
#endif

#define PLOTB ret1 = _pref.blanking->f(this, n1, zonep->blank, zone);   \
  ret2 = _pref.blanking->f(this, n2, zonep->blank, zone);               \
  ret3 = _pref.blanking->f(this, n3, zonep->blank, zone);               \
  ret4 = _pref.blanking->f(this, n4, zonep->blank, zone);               \
  if (ret1*ret2*ret3*ret4 != 0) { PLOT; }    

  E_Int i, j, k, n1, n2, n3, n4, n1s, n2s, n3s, n4s;
  float r, g, b;
#ifndef __SHADERS__
  float a;
#endif
  E_Int ret1, ret2, ret3, ret4;
  //double r2, r3, r4;

  // scalar field
  double* f = zonep->f[nofield];
  double fmin, fmax;
  fmax = maxf[nofield]; fmin = minf[nofield];
  double deltai = MAX(fmax-fmin, ISOCUTOFF);
  //if (fmax-fmin < 1.e-10) { deltai = 1.; fmax = fmin+1.; }
  //else deltai = fmax-fmin;
  deltai = 1./deltai;

  // Colormap
  void (*getrgb)(Data* data, double, float*, float*, float*);
  getrgb = _pref.colorMap->f;

  // Grid dimensions
  E_Int ni = zonep->ni;
  E_Int nj = zonep->nj;
  E_Int nk = zonep->nk;
  if (ptrState->dim == 2) nk = 1;
  E_Int nij = ni*nj;

  E_Int nim = MIN(ni, E_Int(2));
  E_Int njm = MIN(nj, E_Int(2));
  E_Int nkm = MIN(nk, E_Int(2));
  E_Int nbElti = nj*nk*nim;
  E_Int nbEltj = ni*nk*njm;
  E_Int nbEltk = nij*nkm;

  E_Int nistepj = ni*stepj;
  E_Int nijstepk = nij*stepk;
  E_Int njstepk = nj*stepk;
  E_Int nistepk = ni*stepk;
  E_Int inci = (nim-1)*nj*nk;
  E_Int incj = (njm-1)*ni*nk;
  E_Int inck = (nkm-1)*nij;

  double* x = zonep->x; 
  double* y = zonep->y; 
  double* z = zonep->z;
  float* surfp = zonep->surf[0];
  float* surfx = surfp;
  float* surfy = surfx + nbElti;
  float* surfz = surfy + nbElti;

  g = 0; b = 0;

  glBegin(GL_QUADS);  

  // I Plane
  i = zonep->iPlane;
  if (i >= 0) // display single i plane
  {
    if (zonep->blank == 0)
    {
      // No blanking
      for (k = 0; k < nk-stepk; k += stepk)
        for (j = 0; j < nj-stepj; j += stepj)
        {
          n1s = j+k*nj;
          n2s = n1s+stepj;
          n3s = n1s+njstepk;
          n4s = n3s+stepj;
          n1 = i+j*ni+k*nij;
          n2 = n1+nistepj;
          n3 = n1+nijstepk;
          n4 = n3+nistepj;
          PLOT;
        }
    }
    else
    {
      // With blanking
      for (k = 0; k < nk-stepk; k += stepk)
        for (j = 0; j < nj-stepj; j += stepj)
        {
          n1s = j+k*nj;
          n2s = n1s+stepj;
          n3s = n1s+njstepk;
          n4s = n3s+stepj;
          n1 = i+j*ni+k*nij;
          n2 = n1+nistepj;
          n3 = n1+nijstepk;
          n4 = n3+nistepj;
          PLOTB;
        }
    }
  }
  else if (i == -1) // display both planes
  {
    if (zonep->blank == 0)
    {
      // No blanking
      for (k = 0; k < nk-stepk; k += stepk)
        for (j = 0; j < nj-stepj; j += stepj)
        {
          n1s = j+k*nj;
          n2s = n1s+stepj;
          n3s = n1s+njstepk;
          n4s = n3s+stepj;
          n1 = j*ni+k*nij;
          n2 = n1+nistepj;
          n3 = n1+nijstepk;
          n4 = n3+nistepj;
          PLOT;
        }
      for (k = 0; k < nk-stepk; k += stepk)
        for (j = 0; j < nj-stepj; j += stepj)
        {
          n1s = inci+j+k*nj;
          n2s = n1s+stepj;
          n3s = n1s+njstepk;
          n4s = n3s+stepj;
          n1 = ni-1+j*ni+k*nij;
          n2 = n1+nistepj;
          n3 = n1+nijstepk;
          n4 = n3+nistepj;
          PLOT;
        }
    }
    else
    {
      // With blanking
      for (k = 0; k < nk-stepk; k += stepk)
        for (j = 0; j < nj-stepj; j += stepj)
        {
          n1s = j+k*nj;
          n2s = n1s+stepj;
          n3s = n1s+njstepk;
          n4s = n3s+stepj;
          n1 = j*ni+k*nij;
          n2 = n1+nistepj;
          n3 = n1+nijstepk;
          n4 = n3+nistepj;
          PLOTB;
        }
      for (k = 0; k < nk-stepk; k += stepk)
        for (j = 0; j < nj-stepj; j += stepj)
        {
          n1s = inci+j+k*nj;
          n2s = n1s+stepj;
          n3s = n1s+njstepk;
          n4s = n3s+stepj;
          n1 = ni-1+j*ni+k*nij;
          n2 = n1+nistepj;
          n3 = n1+nijstepk;
          n4 = n3+nistepj;
          PLOTB;
        }
    }
  }

  surfp = zonep->surf[0];
  surfx = surfp + 3*nbElti;
  surfy = surfx + nbEltj;
  surfz = surfy + nbEltj;

  // J Plane
  j = zonep->jPlane;
  if (j >= 0)
  {
    if (zonep->blank == 0)
    {
      // No blanking
      for (k = 0; k < nk-stepk; k += stepk)
        for (i = 0; i < ni-stepi; i += stepi)
        {
          n1s = i+k*ni;
          n2s = n1s+stepi;
          n3s = n1s+nistepk;
          n4s = n3s+stepi;
          n1 = i+j*ni+k*nij;
          n2 = n1+stepi;
          n3 = n1+nijstepk;
          n4 = n3+stepi;
          PLOT;
        }
    }
    else
    {
      // With blanking
      for (k = 0; k < nk-stepk; k += stepk)
        for (i = 0; i < ni-stepi; i += stepi)
        {
          n1s = i+k*ni;
          n2s = n1s+stepi;
          n3s = n1s+nistepk;
          n4s = n3s+stepi;
          n1 = i+j*ni+k*nij;
          n2 = n1+stepi;
          n3 = n1+nijstepk;
          n4 = n3+stepi;
          PLOTB;
        }
    }
  }
  else if (j == -1)
  {
    if (zonep->blank == 0)
    {
      // No blanking
      for (k = 0; k < nk-stepk; k += stepk)
        for (i = 0; i < ni-stepi; i += stepi)
        {
          n1s = i+k*ni;
          n2s = n1s+stepi;
          n3s = n1s+nistepk;
          n4s = n3s+stepi;
          n1 = i+k*nij;
          n2 = n1+stepi;
          n3 = n1+nijstepk;
          n4 = n3+stepi;
          PLOT;
        }
      for (k = 0; k < nk-stepk; k += stepk)
        for (i = 0; i < ni-stepi; i += stepi)
        {
          n1s = i+k*ni+incj;
          n2s = n1s+stepi;
          n3s = n1s+nistepk;
          n4s = n3s+stepi;
          n1 = i+(nj-1)*ni+k*nij;
          n2 = n1+stepi;
          n3 = n1+nijstepk;
          n4 = n3+stepi;
          PLOT;
        }
    }
    else
    {
      // With blanking
      for (k = 0; k < nk-stepk; k += stepk)
        for (i = 0; i < ni-stepi; i += stepi)
        {
          n1s = i+k*ni;
          n2s = n1s+stepi;
          n3s = n1s+nistepk;
          n4s = n3s+stepi;
          n1 = i+k*nij;
          n2 = n1+stepi;
          n3 = n1+nijstepk;
          n4 = n3+stepi;
          PLOTB;
        }
      for (k = 0; k < nk-stepk; k += stepk)
        for (i = 0; i < ni-stepi; i += stepi)
        {
          n1s = i+k*ni+incj;
          n2s = n1s+stepi;
          n3s = n1s+nistepk;
          n4s = n3s+stepi;
          n1 = i+(nj-1)*ni+k*nij;
          n2 = n1+stepi;
          n3 = n1+nijstepk;
          n4 = n3+stepi;
          PLOTB;
        }
    }
  }
  
  surfp = zonep->surf[0];
  surfx = surfp + 3*nbElti + 3*nbEltj;
  surfy = surfx + nbEltk;
  surfz = surfy + nbEltk;

  // K Plane
  k = zonep->kPlane;
  if (k >= 0)
  {
    if (zonep->blank == 0)
    {
      // No blanking
      for (j = 0; j < nj-stepj; j += stepj)
        for (i = 0; i < ni-stepi; i += stepi)
        {
          n1s = i+j*ni;
          n2s = n1s+stepi;
          n3s = n1s+nistepj;
          n4s = n3s+stepi;
          n1 = i+j*ni+k*nij;
          n2 = n1+stepi;
          n3 = n1+nistepj;
          n4 = n3+stepi;
          PLOT;
        }
    }
    else
    {
      // With blanking
      for (j = 0; j < nj-stepj; j += stepj)
        for (i = 0; i < ni-stepi; i += stepi)
        {
          n1s = i+j*ni;
          n2s = n1s+stepi;
          n3s = n1s+nistepj;
          n4s = n3s+stepi;
          n1 = i+j*ni+k*nij;
          n2 = n1+stepi;
          n3 = n1+nistepj;
          n4 = n3+stepi;
          PLOTB;
        }
    }
  }
  else
  {
    if (zonep->blank == 0)
    {
      // No blanking
      for (j = 0; j < nj-stepj; j += stepj)
        for (i = 0; i < ni-stepi; i += stepi)
        {
          n1s = i+j*ni;
          n2s = n1s+stepi;
          n3s = n1s+nistepj;
          n4s = n3s+stepi;
          n1 = i+j*ni;
          n2 = n1+stepi;
          n3 = n1+nistepj;
          n4 = n3+stepi;
          PLOT;
        }
      for (j = 0; j < nj-stepj; j += stepj)
        for (i = 0; i < ni-stepi; i += stepi)
        {
          n1s = i+j*ni+inck;
          n2s = n1s+stepi;
          n3s = n1s+nistepj;
          n4s = n3s+stepi;
          n1 = i+j*ni+(nk-1)*nij;
          n2 = n1+stepi;
          n3 = n1+nistepj;
          n4 = n3+stepi;
          PLOT;
        }
    }
    else if (k == -1)
    {
      // With blanking
      for (j = 0; j < nj-stepj; j += stepj)
        for (i = 0; i < ni-stepi; i += stepi)
        {
          n1s = i+j*ni;
          n2s = n1s+stepi;
          n3s = n1s+nistepj;
          n4s = n3s+stepi;
          n1 = i+j*ni;
          n2 = n1+stepi;
          n3 = n1+nistepj;
          n4 = n3+stepi;
          PLOTB;
        }
      for (j = 0; j < nj-stepj; j += stepj)
        for (i = 0; i < ni-stepi; i += stepi)
        {
          n1s = i+j*ni+inck;
          n2s = n1s+stepi;
          n3s = n1s+nistepj;
          n4s = n3s+stepi;
          n1 = i+j*ni+(nk-1)*nij;
          n2 = n1+stepi;
          n3 = n1+nistepj;
          n4 = n3+stepi;
          PLOTB;
        }
    }
  }
  glEnd();

  // Pour les lignes
  if (ni*nj == 1 || ni*nk == 1 || nj*nk == 1)
  {
    glLineWidth(3.);
    glPolygonOffset(-1.,-10.); // force offset
    glBegin(GL_LINES);
    E_Int nie, nje, nke;
    nie = ni; nje = nj; nke = nk;
    if (ni*nj == 1) nke = nke-1;
    if (ni*nk == 1) nje = nje-1;
    if (nj*nk == 1) nie = nie-1;
    if (zonep->blank == 0)
    {
      // No blanking
      for (k = 0; k < nke; k++)
        for (j = 0; j < nje; j++)
          for (i = 0; i < nie; i++)
          {
            n1 = i+j*ni+k*nij;
            n2 = n1+1;
            #ifndef __SHADERS__
              getrgb(this, (f[n1]-fmin)*deltai, &r, &g, &b);
              glColor3f(r, g, b+offb);
              glVertex3d(x[n1], y[n1], z[n1]);
              getrgb(this, (f[n1]-fmin)*deltai, &r, &g, &b);
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
      for (k = 0; k < nke; k++)
        for (j = 0; j < nje; j++)
          for (i = 0; i < nie; i++)
          {
            n1 = i+j*ni+k*nij;
            n2 = n1+1;
            ret1 = _pref.blanking->f(this, n1, zonep->blank, zone);
            ret2 = _pref.blanking->f(this, n2, zonep->blank, zone);
            if (ret1*ret2 != 0)
            { 
              #ifndef __SHADERS__
                getrgb(this, (f[n1]-fmin)*deltai, &r, &g, &b);
                glColor3f(r, g, b+offb);
                glVertex3d(x[n1], y[n1], z[n1]);
                getrgb(this, (f[n1]-fmin)*deltai, &r, &g, &b);
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

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
#define PLOT13 glNormal3f(surfx[n1s], surfy[n1s], surfz[n1s]); \
  glVertex3d(x[n1], y[n1], z[n1]);                             \
  glNormal3f(surfx[n3s], surfy[n3s], surfz[n3s]);  \
  glVertex3d(x[n3], y[n3], z[n3]);

// Plot 2 vertex n2 et n4
#define PLOT24 glNormal3f(surfx[n2s], surfy[n2s], surfz[n2s]); \
  glVertex3d(x[n2], y[n2], z[n2]);                             \
  glNormal3f(surfx[n4s], surfy[n4s], surfz[n4s]);       \
  glVertex3d(x[n4], y[n4], z[n4]);

// Plot a quad n1, n2, n3, n4
#define PLOT glNormal3f(surfx[n1s], surfy[n1s], surfz[n1s]); \
  glVertex3d(x[n1], y[n1], z[n1]);                        \
  glNormal3f(surfx[n2s], surfy[n2s], surfz[n2s]);         \
  glVertex3d(x[n2], y[n2], z[n2]);                        \
  glNormal3f(surfx[n4s], surfy[n4s], surfz[n4s]);         \
  glVertex3d(x[n4], y[n4], z[n4]);                        \
  glNormal3f(surfx[n3s], surfy[n3s], surfz[n3s]);         \
  glVertex3d(x[n3], y[n3], z[n3]);

// Blanking de n1, n3
#define PLOTB13 ret1 = _pref.blanking->f(this, n1, zonep->blank, zone); \
  ret3 = _pref.blanking->f(this, n3, zonep->blank, zone); \
  ret13 = ret1*ret3;

// Blanking de n2, n4 + PLOT
#define PLOTB24 ret2 = _pref.blanking->f(this, n2, zonep->blank, zone); \
  ret4 = _pref.blanking->f(this, n4, zonep->blank, zone); \
  ret24 = ret2*ret4; \
  if (ret13*ret24 != 0) { PLOT; } \
  n1 = n2; n3 = n4; ret13 = ret24; n1s = n2s; n3s = n4s;

// Plot a quad with blanking
#define PLOTB ret1 = _pref.blanking->f(this, n1, zonep->blank, zone);   \
  ret2 = _pref.blanking->f(this, n2, zonep->blank, zone);               \
  ret3 = _pref.blanking->f(this, n3, zonep->blank, zone);               \
  ret4 = _pref.blanking->f(this, n4, zonep->blank, zone);               \
  if (ret1*ret2*ret3*ret4 != 0) { PLOT; }

// Plot un hexa avec des coords de texture
#define PLOTHEXA glTexCoord3f(0,0,0);                           \
  glVertex3d(x[n1], y[n1], z[n1]);                              \
  glTexCoord3f(0,1,0);                                          \
  glVertex3d(x[n4], y[n4], z[n4]);                              \
  glTexCoord3f(1,1,0);                                          \
  glVertex3d(x[n3], y[n3], z[n3]);                              \
  glTexCoord3f(1,0,0);                                          \
  glVertex3d(x[n2], y[n2], z[n2]);                              \
  glTexCoord3f(0,0,1);                                          \
  glVertex3d(x[n5], y[n5], z[n5]);                              \
  glTexCoord3f(1,0,1);                                          \
  glVertex3d(x[n6], y[n6], z[n6]);                              \
  glTexCoord3f(1,1,1);                                          \
  glVertex3d(x[n7], y[n7], z[n7]);                              \
  glTexCoord3f(0,1,1);                                          \
  glVertex3d(x[n8], y[n8], z[n8]);                              \
  glTexCoord3f(0,0,0);                                          \
  glVertex3d(x[n1], y[n1], z[n1]);                              \
  glTexCoord3f(1,0,0);                                          \
  glVertex3d(x[n2], y[n2], z[n2]);                              \
  glTexCoord3f(1,0,1);                                          \
  glVertex3d(x[n6], y[n6], z[n6]);                              \
  glTexCoord3f(0,0,1);                                          \
  glVertex3d(x[n5], y[n5], z[n5]);                              \
  glTexCoord3f(0,1,0);                                          \
  glVertex3d(x[n4], y[n4], z[n4]);                              \
  glTexCoord3f(0,1,1);                                          \
  glVertex3d(x[n8], y[n8], z[n8]);                              \
  glTexCoord3f(1,1,1);                                          \
  glVertex3d(x[n7], y[n7], z[n7]);                              \
  glTexCoord3f(1,1,0);                                          \
  glVertex3d(x[n3], y[n3], z[n3]);                              \
  glTexCoord3f(0,0,0);                                          \
  glVertex3d(x[n1], y[n1], z[n1]);                              \
  glTexCoord3f(0,0,1);                                          \
  glVertex3d(x[n5], y[n5], z[n5]);                              \
  glTexCoord3f(0,1,1);                                          \
  glVertex3d(x[n8], y[n8], z[n8]);                              \
  glTexCoord3f(0,1,0);                                          \
  glVertex3d(x[n4], y[n4], z[n4]);                              \
  glTexCoord3f(1,0,0);                                          \
  glVertex3d(x[n2], y[n2], z[n2]);                              \
  glTexCoord3f(1,1,0);                                          \
  glVertex3d(x[n3], y[n3], z[n3]);                              \
  glTexCoord3f(1,1,1);                                          \
  glVertex3d(x[n7], y[n7], z[n7]);                              \
  glTexCoord3f(1,0,1);                                          \
  glVertex3d(x[n6], y[n6], z[n6]);

  // Grid dimensions
  E_Int ni = zonep->ni;
  E_Int nj = zonep->nj;
  E_Int nk = zonep->nk;
  E_Int n1s, n2s, n3s, n4s;
  if (ptrState->dim == 2) nk = 1;
  E_Int nij = ni*nj;
  E_Int nistepj = ni*stepj;
  E_Int nijstepk = nij*stepk;
  
  E_Int nim = MIN(ni, E_Int(2));
  E_Int njm = MIN(nj, E_Int(2));
  E_Int nkm = MIN(nk, E_Int(2));
  E_Int nbElti = nj*nk*nim;
  E_Int nbEltj = ni*nk*njm;
  E_Int nbEltk = ni*nj*nkm;

  double* x = zonep->x;
  double* y = zonep->y;
  double* z = zonep->z;

  float* surfp = zonep->surf[0];
  float* surfx = surfp;
  float* surfy = surfx + nbElti;
  float* surfz = surfy + nbElti;
  
  // Traitement Smoke
  if (zonep->material == 6) // Cube smoke
  {
    glBegin(GL_QUADS);
    n1 = 0; n2 = 1; n3 = 3; n4 = 2; n5 = 4; n6= 5; n7 = 7; n8 = 6;
    PLOTHEXA;
    glEnd();
    return;
  }
  
  if (zonep->blank == 0) // no blanking
  {
    // I Plane
    i = zonep->iPlane;
    if (i >= 0) // single plane
    {
      for (k = 0; k < nk-stepk; k += stepk)
      {
        n1s = k*nj;
        n1 = i+k*nij;
        n3s = n1s+nj*stepk;
        n3 = n1+nijstepk;
        glBegin(GL_QUAD_STRIP);
        PLOT13;
        for (j = stepj; j < nj; j += stepj)
        {
          n2s = n1s+j;
          n2 = n1+j*ni;
          n4s = n2s+nj*stepk;
          n4 = n2+nijstepk;
          PLOT24;
        }
        glEnd();
      }
    }
    else if (i == -1) // both planes
    {
      for (k = 0; k < nk-stepk; k += stepk)
      {
        n1s = k*nj;
        n1 = k*nij;
        n3s = n1s+nj*stepk;
        n3 = n1+nijstepk;
        glBegin(GL_QUAD_STRIP);
        PLOT13;
        for (j = stepj; j < nj; j += stepj)
        {
          n2s = n1s+j;
          n2 = n1+j*ni;
          n4s = n2s+nj*stepk;
          n4 = n2+nijstepk;
          PLOT24;
        }
        glEnd();
      }

      for (k = 0; k < nk-stepk; k += stepk)
      {
        n1s = (nim-1)*nj*nk+k*nj;
        n1 = (ni-1)+k*nij;
        n3s = n1s+nj*stepk;
        n3 = n1+nijstepk;
        glBegin(GL_QUAD_STRIP);
        PLOT13;
        for (j = stepj; j < nj; j += stepj)
        {
          n2s = n1s+j;
          n2 = n1+j*ni;
          n4s = n2s+nj*stepk;
          n4 = n2+nijstepk;
          PLOT24;
        }
        glEnd();
      }
    }
    
    surfx = surfp + 3*nbElti;
    surfy = surfx + nbEltj;
    surfz = surfy + nbEltj;
  
    // J Plane
    j = zonep->jPlane;
    if (j >= 0)
    {
      for (k = 0; k < nk-stepk; k += stepk)
      {
        n1s = k*ni;
        n1 = j*ni+k*nij;
        n3s = n1s+ni*stepk;
        n3 = n1+nijstepk;
        glBegin(GL_QUAD_STRIP);
        PLOT13;
        for (i = stepi; i < ni; i += stepi)
        {
          n2s = n1s+i;
          n2 = n1+i;
          n4s = n2s+ni*stepk;
          n4 = n2+nijstepk;
          PLOT24;
        }
        glEnd();
      }
    }
    else if (j == -1)
    {
      for (k = 0; k < nk-stepk; k += stepk)
      {
        n1s = k*ni;
        n1 = k*nij;
        n3s = n1s+ni*stepk;
        n3 = n1+nijstepk;
        glBegin(GL_QUAD_STRIP);
        PLOT13;
        for (i = stepi; i < ni; i += stepi)
        {
          n2s = n1s+i;
          n2 = n1+i;
          n4s = n2s+ni*stepk;
          n4 = n2+nijstepk;
          PLOT24;
        }
        glEnd();
      }
      for (k = 0; k < nk-stepk; k += stepk)
      {
        n1s = (njm-1)*ni*nk+k*ni;
        n1 = (nj-1)*ni + k*nij;
        n3s = n1s+ni*stepk;
        n3 = n1+nijstepk;
        glBegin(GL_QUAD_STRIP);
        PLOT13;
        for (i = stepi; i < ni; i += stepi)
        {
          n2s = n1s+i;
          n2 = n1+i;
          n4s = n2s+ni*stepk;
          n4 = n2+nijstepk;
          PLOT24;
        }
        glEnd();
      }
    }
   
    surfx = surfp + 3*nbElti + 3*nbEltj;
    surfy = surfx + nbEltk;
    surfz = surfy + nbEltk;

    // K Plane
    k = _szones[zone]->kPlane;
    if (k >= 0)
    {
      for (j = 0; j < nj-stepj; j += stepj)
      {
        n1s = j*ni;
        n1 = j*ni+k*nij;
        n3s = n1s+nistepj;
        n3 = n1+nistepj;
        glBegin(GL_QUAD_STRIP);
        PLOT13;
        for (i = stepi; i < ni; i += stepi)
        {
          n2s = n1s+i;
          n2 = n1+i;
          n4s = n2s+nistepj;
          n4 = n2+nistepj;
          PLOT24;
        }
        glEnd();
      }
    }
    else if (k == -1)
    {
      for (j = 0; j < nj-stepj; j += stepj)
      {
        n1s = j*ni;
        n1 = j*ni;
        n3s = n1s+nistepj;
        n3 = n1+nistepj;
        glBegin(GL_QUAD_STRIP);
        PLOT13;
        for (i = stepi; i < ni; i += stepi)
        {
          n2s = n1s+i;
          n2 = n1+i;
          n4s = n2s+nistepj;
          n4 = n2+nistepj;
          PLOT24;
        }
        glEnd();
      }
      for (j = 0; j < nj-stepj; j += stepj)
      {
        n1s = (nkm-1)*ni*nj+j*ni;
        n1 = j*ni+(nk-1)*nij;
        n3s = n1s+nistepj;
        n3 = n1+nistepj;
        glBegin(GL_QUAD_STRIP);
        PLOT13;
        for (i = stepi; i < ni; i += stepi)
        {
          n2s = n1s+i;
          n2 = n1+i;
          n4s = n2s+nistepj;
          n4 = n2+nistepj;
          PLOT24;
        }
        glEnd();
      }
    }
    
  }
  else // blanking --
  {
    glBegin(GL_QUADS);

    // I Plane
    i = zonep->iPlane;
    if (i >= 0) // single plane
    {
      for (k = 0; k < nk-stepk; k += stepk)
      {
        n1s = k*nj;
        n1 = i+k*nij;
        n3s = n1s+nj*stepk;
        n3 = n1+nijstepk;
        PLOTB13;
        for (j = stepj; j < nj; j += stepj)
        {
          n2s = n1s+stepj;
          n2 = n1+stepj*ni;
          n4s = n2s+nj*stepk;
          n4 = n2+nijstepk;
          PLOTB24;
        }
      }
    }
    else if (i == -1) // both planes
    {
      for (k = 0; k < nk-stepk; k += stepk)
      {
        n1s = k*nj;
        n1 = k*nij;
        n3s = n1s+nj*stepk;
        n3 = n1+nijstepk;
        PLOTB13;
        for (j = stepj; j < nj; j += stepj)
        {
          n2s = n1s+stepj;
          n2 = n1+stepj*ni;
          n4s = n2s+nj*stepk;
          n4 = n2+nijstepk;
          PLOTB24;
        }
      }
      for (k = 0; k < nk-stepk; k += stepk)
      {
        n1s = (nim-1)*nj*nk+k*nj;
        n1 = (ni-1)+k*nij;
        n3s = n1s+nj*stepk;
        n3 = n1+nijstepk;
        PLOTB13;
        for (j = stepj; j < nj; j += stepj)
        {
          n2s = n1s+stepj;
          n2 = n1+stepj*ni;
          n4s = n2s+nj*stepk;
          n4 = n2+nijstepk;
          PLOTB24;
        }
      }
    }
  
    surfx = surfp + 3*nbElti;
    surfy = surfx + nbEltj;
    surfz = surfy + nbEltj;
    
    // J Plane
    j = zonep->jPlane;
    if (j >= 0)
    {
      for (k = 0; k < nk-stepk; k += stepk)
      {
        n1s = k*ni;
        n1 = j*ni+k*nij;
        n3s = n1s+ni*stepk;
        n3 = n1+nijstepk;
        PLOTB13;
        for (i = stepi; i < ni; i += stepi)
        {
          n2s = n1s+stepi;
          n2 = n1+stepi;
          n4s = n2s+ni*stepk;
          n4 = n2+nijstepk;
          PLOTB24;
        }
      }
    }
    else if (j == -1)
    {
      for (k = 0; k < nk-stepk; k += stepk)
      {
        n1s = k*ni;
        n1 = k*nij;
        n3s = n1s+ni*stepk;
        n3 = n1+nijstepk;
        PLOTB13;
        for (i = stepi; i < ni; i += stepi)
        {
          n2s = n1s+stepi;
          n2 = n1+stepi;
          n4s = n2s+ni*stepk; 
          n4 = n2+nijstepk;
          PLOTB24;
        }
      }
      for (k = 0; k < nk-stepk; k += stepk)
      {
        n1s = (njm-1)*ni*nk+k*ni;
        n1 = (nj-1)*ni+k*nij;
        n3s = n1s+ni*stepk;
        n3 = n1+nijstepk;
        PLOTB13;
        for (i = stepi; i < ni; i += stepi)
        {
          n2s = n1s+stepi;
          n2 = n1+stepi;
          n4s = n2s+ni*stepk;
          n4 = n2+nijstepk;
          PLOTB24;
        }
      }
    }
  
    surfx = surfp + 3*nbElti + 3*nbEltj;
    surfy = surfx + nbEltk;
    surfz = surfy + nbEltk;

    // K Plane
    k = _szones[zone]->kPlane;
    if (k >= 0)
    {
      for (j = 0; j < nj-stepj; j += stepj)
      {
        n1s = j*ni;
        n1 = j*ni+k*nij;
        n3s = n1s+nistepj;
        n3 = n1+nistepj;
        PLOTB13;
        for (i = stepi; i < ni; i += stepi)
        {
          n2s = n1s+stepi;
          n2 = n1+stepi;
          n4s = n2s+nistepj;
          n4 = n2+nistepj;
          PLOTB24;
        }
      }
    }
    else if (k == -1)
    {
      for (j = 0; j < nj-stepj; j += stepj)
      {
        n1s = j*ni;
        n1 = j*ni;
        n3s = n1s+nistepj;
        n3 = n1+nistepj;
        PLOTB13;
        for (i = stepi; i < ni; i += stepi)
        {
          n2s = n1s+stepi;
          n2 = n1+stepi;
          n4s = n2s+nistepj;
          n4 = n2+nistepj;
          PLOTB24;
        }
      }
      for (j = 0; j < nj-stepj; j += stepj)
      {
        n1s = (nkm-1)*ni*nj+j*ni;
        n1 = j*ni+(nk-1)*nij;
        n3s = n1s+nistepj;
        n3 = n1+nistepj;
        PLOTB13;
        for (i = stepi; i < ni; i += stepi)
        {
          n2s = n1s+stepi;
          n2 = n1+stepi;
          n4s = n2s+nistepj;
          n4 = n2+nistepj;
          PLOTB24;
        }
      }
    }
    glEnd();
  }
  
  // Pour les lignes
  if (ni*nj == 1 || ni*nk == 1 || nj*nk == 1)
  {
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
            glVertex3d(x[n1], y[n1], z[n1]);
            glVertex3d(x[n2], y[n2], z[n2]);
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
              glVertex3d(x[n1], y[n1], z[n1]);
              glVertex3d(x[n2], y[n2], z[n2]);
            }
          }
    }
    glEnd();
  }

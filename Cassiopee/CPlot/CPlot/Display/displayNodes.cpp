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
#include "Data.h"

#define PLOTBILLBOARD \
    xi = x[i]; yi = y[i]; zi = z[i]; \
    if (radiusField >= 0) d = zonep->f[radiusField][i]; \
    else { \
    dx = xi - xcam; dy = yi - ycam; dz = zi - zcam; \
    dist = dx*dx + dy*dy + dz*dz; \
    d = sqrt(dist)*dref; } \
    pru0 = d*(right[0] + up[0]); \
    pru1 = d*(right[1] + up[1]); \
    pru2 = d*(right[2] + up[2]); \
    mru0 = d*(right[0] - up[0]); \
    mru1 = d*(right[1] - up[1]); \
    mru2 = d*(right[2] - up[2]); \
    pt1[0] = xi - pru0; \
    pt1[1] = yi - pru1; \
    pt1[2] = zi - pru2; \
    pt2[0] = xi + mru0; \
    pt2[1] = yi + mru1; \
    pt2[2] = zi + mru2; \
    pt3[0] = xi + pru0; \
    pt3[1] = yi + pru1; \
    pt3[2] = zi + pru2; \
    pt4[0] = xi - mru0; \
    pt4[1] = yi - mru1; \
    pt4[2] = zi - mru2; \
    glTexCoord2f(0.0, 0.0); glVertex3dv(pt1); \
    glTexCoord2f(1.0, 0.0); glVertex3dv(pt2); \
    glTexCoord2f(1.0, 1.0); glVertex3dv(pt3); \
    glTexCoord2f(0.0, 1.0); glVertex3dv(pt4);

//=============================================================================
/*
  Display les array nodes d'une zone non structuree.
*/
//=============================================================================
void Data::displayNodes()
{
  double pt1[3]; double pt2[3]; double pt3[3]; double pt4[3];
  double xi, yi, zi;

  double viewMatrix[16];
  glGetDoublev(GL_MODELVIEW_MATRIX, viewMatrix);

  double right[3];
  right[0] = viewMatrix[0];
  right[1] = viewMatrix[4];
  right[2] = viewMatrix[8];
  double up[3];
  up[0] = viewMatrix[1];
  up[1] = viewMatrix[5];
  up[2] = viewMatrix[9];
  
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);  
  //glEnable(GL_TEXTURE_2D);
  glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
  
  double xcam = _view.xcam;
  double ycam = _view.ycam;
  double zcam = _view.zcam;
  ptrState->alpha = 1.;
  float color1[3];
  float color2[3];
  float r, g, b;
  double dref = 0.005;
  double d, dx, dy, dz, dist, alphaSav;
  double pru0, pru1, pru2, mru0, mru1, mru2;
  void (*getrgb)(Data* data, double, float*, float*, float*);
  getrgb = _plugins.zoneColorMap->f;

  E_Int zone = 0;
  while (zone < _numberOfUnstructZones)
  {
    UnstructZone* zonep = _uzones[zone];
    E_Int eltType = zonep->eltType[0];

    // if zone is activated and nodes and in frustum
    if (zonep->active == 1 && eltType == 0 && isInFrustum(zonep, _view) == 1)
    {
      if (ptrState->mode == RENDER && zonep->material == 9) // billboard
      { 
        alphaSav = ptrState->alpha;
        if (zonep->blending != -1.) ptrState->alpha = zonep->blending;
        displayBillBoards(zonep, zone);
        //glDisable(GL_CULL_FACE); 
        glDepthMask(GL_TRUE); 
        ptrState->alpha = alphaSav;
      }
      else // Display les noeuds sous forme de texture node
      {
#ifdef __SHADERS__
        //if (_shaders.currentShader() != 0)
        //{
        //  glActiveTexture(GL_TEXTURE0);
        //  glBindTexture(GL_TEXTURE_2D, _texNodes);
        //  _shaders.activate((short unsigned int)0);
        //}
#endif
        // look for radius field (if any)
        char** v = zonep->varnames;
        E_Int nf = zonep->nfield;
        E_Int radiusField = -1;
        for (E_Int i = 0; i < nf; i++)
        {
            if (strcmp(v[i], "radius") == 0) { radiusField = i; break; }
        }
      
        // Color
        switch (ptrState->meshStyle)
        {
          case 0:
            // Couleur uniforme blanche (defaut)
            color1[0] = 0.5; color1[1] = 0.5; color1[2] = 1.;
            color2[0] = 0.1; color2[1] = 0.1; color2[2] = 1;
            break;
          
          case 1:
            getrgb(this, zone*1./_numberOfUnstructZones, &r, &g, &b); 
            color1[0] = r; color1[1] = g; color1[2] = b;
            if (b > 0.8 && r < 0.2 && g < 0.2)
            { color1[0] = r; color1[1] = b; color1[2] = g; }
            color2[0] = 0.1; color2[1] = 0.1; color2[2] = 1;
            break;

          case 2:
            getrgb(this, zone*1./_numberOfUnstructZones, &r, &g, &b); 
            color1[0] = g; color1[1] = r;  color1[2] = b;
            if (r > 0.8 && g < 0.2 && b < 0.2) 
            {color1[0] = g; color1[1] = r; color1[2] = b;}
            color2[0] = 0.3; color2[1] = 0.3; color2[2] = 1.;
            break;

          default:
            color1[0] = 0.95; color1[1] = 0.95; color1[2] = 1.;
            color2[0] = 0.1;  color2[1] = 0.1;  color2[2] = 1;
            break;
        }
      
        // Ecrasement si renderTag
        if (zonep->colorR != -1.)
        {color1[0] = zonep->colorR; 
          color1[1] = zonep->colorG; 
          color1[2] = zonep->colorB;}
      
#include "selection.h"

        // scalar field
        if (ptrState->mode == SCALARFIELD)
        {
          float r, g, b, offb, blend;
          int nofield = ptrState->scalarField;
          double* f = zonep->f[nofield];
          double fmin, fmax;
          if (_niso[nofield] == -1) // max field values
          { fmax = maxf[nofield]; fmin = minf[nofield]; }
          else { fmax = _isoMax[nofield]; fmin = _isoMin[nofield]; }
          double deltai = MAX(fmax-fmin, ISOCUTOFF);
          deltai = 1./deltai;
      
          offb = 0.;
          if (zonep->selected == 1 && zonep->active == 1)
          {
            offb = 0.2; // blue color offset for active zone
            blend = ptrState->alpha;
          }
          else if (zonep->active == 0) blend = 0.2;
          else blend = ptrState->alpha; // active=1

          glBegin(GL_QUADS);  
          double *x = zonep->x;
          double *y = zonep->y;
          double *z = zonep->z;
          for (int i = 0; i < zonep->np; i++)
          {
            getrgb(this, (f[i]-fmin)*deltai, &r, &g, &b);
            glColor4f(r, g, b+offb, blend);
            PLOTBILLBOARD;
          }
          glEnd();
        }
        else if (ptrState->mode == VECTORFIELD)
        {
          float offb, blend;
          E_Int nofield1 = ptrState->vectorField1;
          E_Int nofield2 = ptrState->vectorField2;
          E_Int nofield3 = ptrState->vectorField3;
          
          double* pr = zonep->f[nofield1];
          double* pg = zonep->f[nofield2];
          double* pb = zonep->f[nofield3];
          
          /*
          double fmin, fmax;
          if (_niso[nofield] == -1) // max field values
          { fmax = maxf[nofield]; fmin = minf[nofield]; }
          else { fmax = _isoMax[nofield]; fmin = _isoMin[nofield]; }
          double deltai = MAX(fmax-fmin, ISOCUTOFF);
          deltai = 1./deltai;
          */

          offb = 0.;
          if (zonep->selected == 1 && zonep->active == 1)
          {
            offb = 0.2; // blue color offset for active zone
            blend = ptrState->alpha;
          }
          else if (zonep->active == 0) blend = 0.2;
          else blend = ptrState->alpha; // active=1

          glBegin(GL_QUADS);  
          double *x = zonep->x;
          double *y = zonep->y;
          double *z = zonep->z;
          for (E_Int i = 0; i < zonep->np; i++)
          {
            glColor4f(pr[i], pg[i], pb[i]+offb, blend);
            PLOTBILLBOARD;
          }
          glEnd();
        }
        else // no field
        {
          glBegin(GL_QUADS);  
          double *x = zonep->x;
          double *y = zonep->y;
          double *z = zonep->z;
          for (E_Int i = 0; i < zonep->np; i++)
          {
            PLOTBILLBOARD;
          }
          glEnd();
        }
      } // std
    }
    zone++;
  }

#ifdef __SHADERS__
  _shaders.activate((short unsigned int)0);
#endif
  glActiveTexture(GL_TEXTURE0);
  glDisable(GL_TEXTURE_2D);
  glDisable(GL_BLEND);
}

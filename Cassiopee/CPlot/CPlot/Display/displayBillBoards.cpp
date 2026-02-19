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
#include "../Data.h"

//=============================================================================
/*
  Display les noeuds d'un maillage sous forme de billboards.
  (mode RENDER).
  On suppose que les zones billboard sont apres les zones opaques!
*/
//=============================================================================
void Data::displayBillBoards(Zone* zonep, E_Int zone)
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
  
  double xcam = _view.xcam;
  double ycam = _view.ycam;
  double zcam = _view.zcam;
  float color1[3];
  float color2[3];
  float r, g, b;
  double dref = 0.02 * zonep->shaderParam1;
  double d, dx, dy, dz, dist;
  double pru0, pru1, pru2, mru0, mru1, mru2;
  void (*getrgb)(Data* data, double, float*, float*, float*);
  getrgb = _plugins.colorMap->next->f;

  // Color
  switch (ptrState->meshStyle)
  {
    case 0:
      // Couleur uniforme blanche (defaut)
      color1[0] = 0.95; color1[1] = 0.95; color1[2] = 1.;
      color2[0] = 0.1;  color2[1] = 0.1;  color2[2] = 1;
      break;
      
    case 1:
      getrgb(this, zone*1./_numberOfUnstructZones, &r, &g, &b); 
      color1[0] = r; color1[1] = g; color1[2] = b;
      if (b > 0.9 && r < 0.1 && g < 0.1) color1[2] = 1.-b;
      color2[0] = 0.1;  color2[1] = 0.1;  color2[2] = 1;
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

  if (zonep->selected == 1)
    glColor4f(color2[0], color2[1], color2[2], ptrState->alpha);
  else
    glColor4f(color1[0], color1[1], color1[2], ptrState->alpha);

#ifdef __SHADERS__
  if (zonep->colorR < -1.5 && zonep->shaderParam2 < 0.1) // iso
  {
    E_Int shader = _shaders.shader_id(shader::iso_sphere_billboarding);
    ptrState->billBoardNi = 1;
    ptrState->billBoardNj = 1;
    ptrState->billBoardWidth = 1;
    ptrState->billBoardHeight = 1;
    glActiveTexture(GL_TEXTURE1);
    if (_texColormap == 0) createColormapTexture();
    fillColormapTexture((int)_pref.colorMap->varName[0]-48);
    glBindTexture(GL_TEXTURE_1D, _texColormap);
    if (_shaders.currentShader() != shader) _shaders.activate((short unsigned int)shader);
    _shaders[shader]->setUniform("colormap", (int)1);
    int nofield = -int(zonep->colorR)-2;
    if (_niso[nofield] == -1)
    {
        _shaders[shader]->setUniform("niso", (float)ptrState->niso);
        _shaders[shader]->setUniform("alpha", (float)1.);
        _shaders[shader]->setUniform("beta", (float)0.);
        _shaders[shader]->setUniform("amin", (float)0.);
        _shaders[shader]->setUniform("amax", (float)1.);
    }
    else
    {
        double rmin, rmax, alpha, beta;
        double deltai = MAX(maxf[nofield]-minf[nofield], ISOCUTOFF);
        rmin = (_isoMin[nofield] -minf[nofield])/deltai;
        rmax = (_isoMax[nofield] -minf[nofield])/deltai;
        alpha = 1./MAX(rmax-rmin, ISOCUTOFF); beta = -rmin*alpha;
        _shaders[shader]->setUniform("niso", (float)_niso[nofield]);
        _shaders[shader]->setUniform("alpha", (float)alpha);
        _shaders[shader]->setUniform("beta", (float)beta);
        double amin = (_isoAlphaMin[nofield] - minf[nofield])/deltai;
        double amax = (_isoAlphaMax[nofield] - minf[nofield])/deltai;
        _shaders[shader]->setUniform("amin", (float)amin);
        _shaders[shader]->setUniform("amax", (float)amax);
    }
  }
  else if (zonep->selected == 1 && zonep->active == 1) triggerShader(*zonep, zonep->material, 0., color2);
  else triggerShader(*zonep, zonep->material, 0., color1);
#endif

  srand(1); // init aleatoire

  // Scale base sur la bb globale
  dx = xmin - xcam; dy = ymin - ycam; dz = zmin - zcam;
  dist = dx*dx + dy*dy + dz*dz;
  dx = xmax - xcam; dy = ymax - ycam; dz = zmax - zcam;
  dist = K_FUNC::E_max(dist, dx*dx + dy*dy + dz*dz);
  if (ptrState->billBoardSize <= 1.e-6) d = sqrt(dist)*dref;
  else d = ptrState->billBoardSize*dref;
  //printf("%g %g\n", dref, ptrState->billBoardSize);

  // Billboard size
  E_Int bw = (E_Int)(ptrState->billBoardWidth);
  E_Int bh = (E_Int)(ptrState->billBoardHeight);
  E_Int Ni = (E_Int)(ptrState->billBoardNi);
  E_Int Nj = (E_Int)(ptrState->billBoardNj);
  float rt = (bh*Ni*1./(bw*Nj));
  //printf("Width=%d %d - %g\n",bw,bh,rt);
  //rt = 1.; // DBX

  glBegin(GL_QUADS);  
  double *x = zonep->x;
  double *y = zonep->y;
  double *z = zonep->z;

  pru0 = d*(right[0] + rt*up[0]);
  pru1 = d*(right[1] + rt*up[1]);
  pru2 = d*(right[2] + rt*up[2]);
  mru0 = d*(right[0] - rt*up[0]);
  mru1 = d*(right[1] - rt*up[1]);
  mru2 = d*(right[2] - rt*up[2]);
  int npts = zonep->npts;

  // compute cam distance
  double* di = new double [npts];
  double* ran = new double [npts];
  double dmin = 1.e30, dmax = -1.e30;
  double randConst = 1./99.;

  // look for billBoard field (field to choose image in multi-image billboards)
  char** v = zonep->varnames;
  E_Int nf = zonep->nfield;
  for (E_Int i = 0; i < nf; i++)
  {
    if (strcmp(v[i], "TBB__") == 0) ptrState->billBoardT = i;
  }
  //ptrState->billBoardT = -1;

  // look for radius field (tell the size of billboards)
  v = zonep->varnames;
  nf = zonep->nfield;
  E_Int radiusField = -1;
  for (E_Int i = 0; i < nf; i++)
  {
    if (strcmp(v[i], "radius") == 0) { radiusField = i; break; }
  }
  
  // Compute ran field (choose image in billboard)
  // Compute di distance to camera field
  for (E_Int i = 0; i < npts; i++)
  {
    xi = x[i]; yi = y[i]; zi = z[i];
    di[i] = sqrt((xcam-xi)*(xcam-xi)+(ycam-yi)*(ycam-yi)+(zcam-zi)*(zcam-zi));
    dmin = std::min(di[i], dmin);
    dmax = std::max(di[i], dmax);
    if (ptrState->billBoardT == -1) ran[i] = rand()%100*randConst; // entre 0 et 1
    else ran[i] = zonep->f[ptrState->billBoardT][i]; // doit etre entre 0 et 1
  }

  E_Int NSplit = 25;
  dmin = dmin-1.e-6;
  dmax = dmax+1.e-6;
  double delta;
  delta = (dmax-dmin)/(NSplit-1.);
  double range, ranged;
  E_Int e1,e2;

  // Field for iso (if needed)
  E_Int nofield=0; double* f=NULL; float deltai=1.; float fmin=0.;
  if (zonep->colorR < -1.5) // iso 
  {
    nofield = -int(zonep->colorR)-2;
    f = zonep->f[nofield];
    deltai = 1./MAX(maxf[nofield]-minf[nofield], ISOCUTOFF);
    fmin = minf[nofield];
  }

  // render
  for (E_Int n = NSplit-1; n >= 0; n--)
  {
    range = dmin+delta*n;
    ranged = range+delta;

    if (zonep->shaderParam2 < 0.1) // sphere shader
    {
      for (E_Int i = 0; i < npts; i++)
      {
        if (di[i] > range && di[i] <= ranged)
        {
          if (zonep->colorR < -1.5) // for iso color
          { 
            glColor4f((f[i]-fmin)*deltai,0.,0.,zonep->blending);
          }

          xi = x[i]; yi = y[i]; zi = z[i];
          
          if (radiusField >= 0)
          {
            d = zonep->f[radiusField][i];
            pru0 = d*(right[0] + rt*up[0]);
            pru1 = d*(right[1] + rt*up[1]);
            pru2 = d*(right[2] + rt*up[2]);
            mru0 = d*(right[0] - rt*up[0]);
            mru1 = d*(right[1] - rt*up[1]);
            mru2 = d*(right[2] - rt*up[2]);
          }

          pt1[0] = xi - pru0;
          pt1[1] = yi - pru1;
          pt1[2] = zi - pru2;
        
          pt2[0] = xi + mru0;
          pt2[1] = yi + mru1;
          pt2[2] = zi + mru2;
        
          pt3[0] = xi + pru0;
          pt3[1] = yi + pru1;
          pt3[2] = zi + pru2;
        
          pt4[0] = xi - mru0;
          pt4[1] = yi - mru1;
          pt4[2] = zi - mru2;
    
          glTexCoord2f(0.0, 0.0); glVertex3dv(pt1);
          glTexCoord2f(1.0, 0.0); glVertex3dv(pt2);
          glTexCoord2f(1.0, 1.0); glVertex3dv(pt3);
          glTexCoord2f(0.0, 1.0); glVertex3dv(pt4);
        }
      }
    }
    else // billboard shader texture
    {
      for (E_Int i = 0; i < npts; i++)
      {
        if (di[i] > range && di[i] <= ranged)
        {
          xi = x[i]; yi = y[i]; zi = z[i];
          
          pt1[0] = xi - pru0;
          pt1[1] = yi - pru1;
          pt1[2] = zi - pru2;
          
          pt2[0] = xi + mru0;
          pt2[1] = yi + mru1;
          pt2[2] = zi + mru2;
          
          pt3[0] = xi + pru0;
          pt3[1] = yi + pru1;
          pt3[2] = zi + pru2;
          
          pt4[0] = xi - mru0;
          pt4[1] = yi - mru1;
          pt4[2] = zi - mru2;
          
          // Le champ T est transmis dans r avec la couleur encodee
          e1 = E_Int(ran[i]*255.); // encode entre 0 et 255
          e2 = E_Int(color1[0]*255.); // encode entre 0 et 255
          //printf("encode %d %d %f\n", e1,e2,(e1+2*e2)/768.);
          glColor4f((e1+255*e2)/65280., color1[1], color1[2], ptrState->alpha);
        
          //printf("%f %f %f - %f %f %f\n", pt1[0], pt1[1], pt1[2],
          //       pt2[0], pt2[1], pt2[2]);
          glTexCoord2f(0.0, 0.0); glVertex3dv(pt1);
          glTexCoord2f(1.0, 0.0); glVertex3dv(pt2);
          glTexCoord2f(1.0, 1.0); glVertex3dv(pt3);
          glTexCoord2f(0.0, 1.0); glVertex3dv(pt4);
        }
      }
    }
  }
  glEnd();
  delete [] di; delete [] ran;
}

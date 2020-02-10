/*
    Copyright 2013-2020 Onera.

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

#define PLOTNODE xi = x[i]; yi = y[i]; zi = z[i];               \
  dx = xi - xcam; dy = yi - ycam; dz = zi - zcam;               \
  dist = dx*dx + dy*dy + dz*dz;                                 \
  d = sqrt(dist)*dref;                                          \
  pru0 = d*(right[0] + up[0]);                                  \
  pru1 = d*(right[1] + up[1]);                                  \
  pru2 = d*(right[2] + up[2]);                                  \
  mru0 = d*(right[0] - up[0]);                                  \
  mru1 = d*(right[1] - up[1]);                                  \
  mru2 = d*(right[2] - up[2]);                                  \
  pt1[0] = xi - pru0;                                           \
  pt1[1] = yi - pru1;                                           \
  pt1[2] = zi - pru2;                                           \
  pt2[0] = xi + mru0;                                           \
  pt2[1] = yi + mru1;                                           \
  pt2[2] = zi + mru2;                                           \
  pt3[0] = xi + pru0;                                           \
  pt3[1] = yi + pru1;                                           \
  pt3[2] = zi + pru2;                                           \
  pt4[0] = xi - mru0;                                           \
  pt4[1] = yi - mru1;                                           \
  pt4[2] = zi - mru2;                                           \
  glVertex3dv(pt1); glVertex3dv(pt2);                           \
  glVertex3dv(pt3); glVertex3dv(pt4);

//=============================================================================
/*
  Display une zone en mesh.
  IN: zonep: pointeur sur la zone a afficher
  IN: zone: le no de zone non structure
  IN: zonet: le no de la zone dans liste globale des zones
*/
//=============================================================================
void Data::displayUMeshZone(UnstructZone* zonep, int zone, int zonet)
{
  if ( zonep->_is_high_order == true )
  {
    displayUMeshZone_ho(zonep, zone, zonet);
    return;
  }
  int i, n1, n2, ret1, ret2, ret;

  // Style colors
  float color1[3]; float color2[3];

  // Colormap
  float r, g, b;
  void (*getrgb)(Data* data, double, float*, float*, float*);
  getrgb = _plugins.colorMap->next->f;

  // For node rendering (1D zones)
  double dref = 0.004;
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
  double dx, dy, dz, dist, d;
  double pru0, pru1, pru2, mru0, mru1, mru2;
  double pt1[3]; double pt2[3]; double pt3[3]; double pt4[3];
  
  E_Float nz = 1./_numberOfUnstructZones;
#include "meshStyles.h"

  if (zonep->eltType == 1 || zonep->eltType == 0 || (zonep->eltType == 10 && zonep->nelts1D > 0)) 
  { glLineWidth(3.); color2[0] = 0.1; color2[1] = 0.1; color2[2] = 1.; }
      
#include "selection.h"
    
    /*if ( zonep->_is_high_order == true )
    {
      int ishader = 0;
      if ( zonep->eltType == UnstructZone::TRI )
        ishader = 1;  // OK, element de type Tri_6
      if ( zonep->eltType == UnstructZone::QUAD )
        ishader = 2;  // OK, element de type Quad_8
      if ( not this->_shaders.has_tesselation() ) {
        this->_shaders.set_tesselation( ishader );
      }

      this->_shaders.activate( (short unsigned int)this->_shaders.shader_id(0) );
      std::cerr << "Shader id ::: " << this->_shaders.currentShader() << std::flush << std::endl;
      int t_inner = this->ptrState->inner_tesselation;
      int t_outer = this->ptrState->outer_tesselation;
      this->_shaders[ this->_shaders.currentShader() ]->setUniform( "uInner", (float)t_inner );
      this->_shaders[ this->_shaders.currentShader() ]->setUniform( "uOuter", (float)t_outer );
      glPatchParameteri( GL_PATCH_VERTICES, zonep->eltSize );
    }*/

#include "displayUMeshZone.h"

  // For BARS or NODES or 1D NGONS: display node
  if (eltType == 1 || eltType ==  0 || (eltType == 10 && zonep->nelts1D > 0))
  {
    glBegin(GL_QUADS);
    if (zonep->blank == 0)
    {
      for (i = 0; i < zonep->np; i++)
      {
        PLOTNODE;
      }
    }
    else
    {
      for (i = 0; i < zonep->np; i++)
      {
        ret = _pref.blanking->f(this, i, zonep->blank, zone);
        if (ret != 0)
        {
          PLOTNODE;
        }
      }
    }
    glEnd();
  }
  glLineWidth(1.);
}

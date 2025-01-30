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
#include "../Data.h"

//=============================================================================
/*
  Display une zone en solid ou en material (mode SOLID ou RENDER)
*/
//=============================================================================
void Data::displayUSolidHOZone(UnstructZone *zonep, E_Int zone, E_Int zonet)
{
  E_Int i, n1, n2, n3, n4, n5, n6, n7, n8;
  E_Int ret1, ret2, ret3, ret4, f;

  // Style
  float color1[ 3 ];
  float color2[ 3 ];

  // Colormap
  float r, g, b;
  void ( *getrgb )( Data * data, double, float *, float *, float * );
  getrgb = _plugins.zoneColorMap->f;

  E_Float nz = 1. / _numberOfUnstructZones;
  E_Int eltType0 = zonep->eltType[0];
  bool is1D = ( (eltType0 == 1) | (eltType0 == 10 && zonep->nelts1D > 0) );

  #include "solidStyles.h"
  #include "selection.h"

  // scale
  E_Float s = MAX(zonep->xmax - zonep->xmin, zonep->ymax - zonep->ymin);
  s = MAX(s, zonep->zmax - zonep->zmin);
  s = 100. / (s + 1.e-12);

  // Only for textured rendering, we use vect display =======================
  if (ptrState->mode == RENDER && zonep->material == 14 && zonep->texu != NULL)  // Textured rendering
  {
    // Sans doute également à modifier pour le high order ?
    triggerShader(*zonep, zonep->material, s, color1);

    E_Int ff = 0;
    double offb = 0.;
#undef PLOT
    double* f1 = zonep->texu;
    double* f2 = zonep->texv;
    double* f3 = zonep->texw;
    double fmin1, fmax1, fmin2, fmax2, fmin3, fmax3;
    fmax1 = 0.; fmin1 = 1.;
    fmax2 = 0.; fmin2 = 1.;
    fmax3 = 0.; fmin3 = 1.;
#define GL_QUADS_ARE GL_QUADS
#define PLOTQUAD PLOTQUADQ
#define PLOTQUAD2 PLOTQUADQ2
#include "displayUVectSolidZone.h"
    glLineWidth( 1. );
    return;
  }
  // END Textured rendering ===========================================
  // Activation du shader de tesselation :
  int ishader = 0;
  if (eltType0 == UnstructZone::TRI)
    ishader = 1;  // OK, element de type Tri_6, TRI_9, TRI_10, TRI_12 or TRI_15
  if (eltType0 == UnstructZone::QUAD)
    ishader = 2;  // OK, element de type Quad_8 ou Quad_9, QUAD_12, QUAD_16, QUAD_25
  this->_shaders.set_tesselation(ishader);
  if (ptrState->mode == RENDER)
  {
    if (zonep->selected == 1 && zonep->active == 1)
      triggerShader(*zonep, zonep->material, s, color2);
    else
      triggerShader(*zonep, zonep->material, s, color1);
  } 
  else 
  {
    if (zonep->selected == 1 && zonep->active == 1)
      triggerShader(*zonep, 0, s, color2);
    else
      triggerShader(*zonep, 0, s, color1);
  }
  // Pour eviter de tracer le low order sans faire expres :-)
  unsigned short idShader = this->_shaders.currentShader();
  int t_inner = this->ptrState->inner_tesselation;
  int t_outer = this->ptrState->outer_tesselation;
  this->_shaders[ idShader ]->setUniform( "uInner", (float)t_inner );
  this->_shaders[ idShader ]->setUniform( "uOuter", (float)t_outer );
  this->_shaders[ idShader ]->setUniform( "patch_size", (int)zonep->eltSize[0] );

  glPatchParameteri(GL_PATCH_VERTICES, zonep->eltSize[0]);
  glBegin(GL_PATCHES);

  E_Int ind_elt, ind;
  E_Int ne = zonep->nec[0];
  E_Int stride = ne;
  E_Int* connect = zonep->connect[0];
  double* x = zonep->x;
  double* y = zonep->y;
  double* z = zonep->z;
    
  for (E_Int ielts = 0; ielts < ne; ++ielts) 
  {
    ind_elt = ielts;
    for (E_Int inode = 0; inode < zonep->eltSize[0]; inode++) 
    {
      ind = connect[ind_elt + inode * stride]-1;
      glVertex3f( (float)x[ind], (float)y[ind], (float)z[ind] );
    }
  }
  glEnd();
  glLineWidth(1.);
  return;
}

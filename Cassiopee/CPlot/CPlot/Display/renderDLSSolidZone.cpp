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
#include "../DataDL.h"
#include "../ZoneImplDL.h"

//=============================================================================
/*
  Display une zone en solid ou en material.
*/
//=============================================================================
void DataDL::renderGPUSSolidZone(StructZone* zonep, E_Int zone)
{
  // Style
  float color1[3]; float color2[3];

  // Colormap
  float r, g, b;
  void (*getrgb)(Data* data, double, float*, float*, float*);
  getrgb = _plugins.zoneColorMap->f;

  E_Float nz = 1./_numberOfStructZones;

  bool is1D = ((zonep->ni*zonep->nj == 1) | (zonep->ni*zonep->nk == 1) | (zonep->nj*zonep->nk == 1));
 
  #include "solidStyles.h"
  #include "selection.h"

  if (is1D == true && ptrState->mode == MESH) return; // already drawn in mesh

  // scale
  E_Float s = MAX(zonep->xmax-zonep->xmin, zonep->ymax-zonep->ymin);
  s = MAX(s, zonep->zmax-zonep->zmin);
  s = 100./(s+1.e-12);

  // Only for textured rendering, we use vect display =======================
  if (ptrState->mode == RENDER && zonep->material == 14 && zonep->texu != NULL) // Textured rendering
  {
#ifdef __SHADERS__
      triggerShader(*zonep, zonep->material, s, color1);
#endif
      double* f1 = zonep->texu;
      double* f2 = zonep->texv;
      double* f3 = zonep->texw;
      double fmin1, fmax1, fmin2, fmax2, fmin3, fmax3;
      fmax1 = 0.; fmin1 = 1.;
      fmax2 = 0.; fmin2 = 1.;
      fmax3 = 0.; fmin3 = 1.;
  
      E_Int stepi, stepj, stepk;
      computeSteps(zonep, stepi, stepj, stepk);
      #undef PLOT
      #define GL_QUADS_ARE GL_QUADS
      #define PLOT PLOTQ
      #include "displaySVectSolidZone.h"
      glLineWidth(1.);
      return;
  }
  // END Textured rendering ============================================

#ifdef __SHADERS__
  if (ptrState->mode == RENDER)
  {
    if (zonep->selected == 1 && zonep->active == 1) 
      triggerShader(*zonep, zonep->material, s, color2);
    else triggerShader(*zonep, zonep->material, s, color1);
  }
  else
  {
    if (zonep->selected == 1 && zonep->active == 1) 
      triggerShader(*zonep, 0, s, color2);
    else triggerShader(*zonep, 0, s, color1);
  }
#endif

  ZoneImplDL* zImpl = static_cast<ZoneImplDL*>(zonep->ptr_impl);
  glCallList(zImpl->_DLsolid);
  glLineWidth(1.);
}

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
   Display billboards.
*/
//=============================================================================
void Data::displayAllBillBoards()
{
  E_Int zone;

  // Enable blending
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  
  zone = 0;
  while (zone < _numberOfStructZones)
  {
    StructZone* zonep = _szones[zone];

    // if zone is active and in frustum
    if ((zonep->active == 1 || 
         (zonep->active == 0 && ptrState->ghostifyDeactivatedZones == 1))
        && isInFrustum(zonep, _view) == 1)
    {
      double alphaSav = ptrState->alpha;
      if (ptrState->mode == RENDER && zonep->blending != -1.)
      {
        ptrState->alpha = zonep->blending;
        if (ptrState->alpha < 0.9999) 
        { /*glEnable(GL_CULL_FACE);*/ glDepthMask(GL_FALSE); }
      }
      if (ptrState->mode == RENDER && zonep->material == 9) // Billboarding
      {
        displayBillBoards((Zone*)zonep, zone);
      }
      ptrState->alpha = alphaSav;
      glDisable(GL_CULL_FACE); glDepthMask(GL_TRUE);
    }
    zone++;
  }

  zone = 0;
  while (zone < _numberOfUnstructZones)
  {
    UnstructZone* zonep = _uzones[zone];

    // if zone is active and in frustum
    if ((zonep->active == 1 ||
         (zonep->active == 0 && ptrState->ghostifyDeactivatedZones == 1))
        && isInFrustum(zonep, _view) == 1)
    {
      double alphaSav = ptrState->alpha;
      if (ptrState->mode == RENDER && zonep->blending != -1.)
      {
        ptrState->alpha = zonep->blending;
        if (ptrState->alpha < 0.9999) 
        { /*glEnable(GL_CULL_FACE);*/ glDepthMask(GL_FALSE); }
      }

      if (ptrState->mode == RENDER && zonep->material == 9) // Billboarding
      {
        displayBillBoards((Zone*)zonep, zone);
      }
      ptrState->alpha = alphaSav;
      glDisable(GL_CULL_FACE); glDepthMask(GL_TRUE);
    }
    zone++;
  }

  noLight();
#ifdef __SHADERS__
  _shaders.activate((short unsigned int)0);
#endif
  glDisable(GL_BLEND);
  glColor4f(1., 1., 1., 1.);
}

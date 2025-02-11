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
#include "DataDL.h"
#include "ZoneImplDL.h"

// Affichage des maillages non structures

//=============================================================================
/* 
   displayUMesh: display the unstructured zones as MESH
   IN: _uzones: unstructured zones definition
   IN: _numberOfUnstructZones: number of U-Zones.
*/
//=============================================================================
void DataDL::displayUMesh()
{  
  if (_numberOfUnstructZones == 0) return;
  E_Int zone, zonet, isDL;

  // Enable blending
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  
  // Render style
  switch (ptrState->meshStyle)
  {
    case 0:
      // Couleur uniforme rouge + fond solide blanc
      glEnable(GL_POLYGON_OFFSET_FILL);
      glEnable(GL_POLYGON_OFFSET_LINE);
      glPolygonOffset(1., 1.);
      break;

    case 1:
      // Multicolor
      ptrState->alpha = 0.5;
      glEnable(GL_BLEND);
      glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
      glEnable(GL_LINE_SMOOTH);
      glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
      glDepthMask(GL_FALSE);
      break;

    case 2:
      // Multicolor avec fond solide
      glEnable(GL_POLYGON_OFFSET_FILL);
      glEnable(GL_POLYGON_OFFSET_LINE);
      glPolygonOffset(1., 1.);
      break;
      
    case 3:
      // Monocolor black + fond solide cyan
      glEnable(GL_POLYGON_OFFSET_FILL);
      glEnable(GL_POLYGON_OFFSET_LINE);
      glPolygonOffset(1., 1.);
      break;

    case 4:
      // Monocolor cyan + fond solide black
      glEnable(GL_POLYGON_OFFSET_FILL);
      glEnable(GL_POLYGON_OFFSET_LINE);
      glPolygonOffset(1., 1.);
      break;
  }

  zone = 0;
  while (zone < _numberOfUnstructZones)
  {
    zonet = zone + _numberOfStructZones;
    UnstructZone* zonep = _uzones[zone];
    ZoneImplDL* zoneImpl = static_cast<ZoneImplDL*>(zonep->ptr_impl);

    // if zone is active and in frustum
    if (zonep->active == 1 && isInFrustum(zonep, _view) == 1)
    {
#if (USEDLMESH == 1)
      isDL = zoneImpl->_DLmesh;
#else
      isDL = zoneImpl->_DLsolid;
#endif
      if (ptrState->simplifyOnDrag == 1 && ptrState->ondrag == 1) displayUBBZone(zonep);
      else
      {
        if (isDL == 0) displayUMeshZone(zonep, zone, zonet);
        else renderGPUUMeshZone(zonep, zone, zonet);
      }
    }
    zone++;
  }
  
  glDepthMask(GL_TRUE);
  glDisable(GL_POLYGON_OFFSET_FILL);
  glDisable(GL_POLYGON_OFFSET_LINE);
  glDisable(GL_LINE_SMOOTH);
  glDisable(GL_BLEND);

  ptrState->alpha = 1.;
  glColor3f(1., 1., 1.);

  float alphaSav = ptrState->alpha;
  E_Int solidStyleSav = ptrState->solidStyle;
  ptrState->alpha = 0.01;
  switch (ptrState->meshStyle)
  {
    case 0:
      ptrState->solidStyle = 2; ptrState->alpha = 1.;
      break;

    case 2:
      ptrState->solidStyle = 1; ptrState->alpha = 0.95;
      break;

    case 3:
      ptrState->solidStyle = 3; ptrState->alpha = 1.;
      break;
  
    case 4:
      ptrState->solidStyle = 4; ptrState->alpha = 1.;
      break;
  }
  displayUSolid();  
  ptrState->alpha = alphaSav;
  ptrState->solidStyle = solidStyleSav;
}

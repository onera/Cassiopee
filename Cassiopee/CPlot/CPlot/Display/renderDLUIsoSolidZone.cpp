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

//=============================================================================
/*
  Display une zone en iso solides pour un champ scalaire.
  IN: zonep: pointeur sur la zone a afficher
  IN: zonet: le no de la zone dans liste globale des zones
  IN: nofield: le no du champ a afficher.
*/
//=============================================================================
void DataDL::renderUIsoSolidZone(UnstructZone* zonep, E_Int zonet, E_Int nofield)
{
  float offb;
  double blend;
  
  ZoneImplDL* zoneImpl = static_cast<ZoneImplDL*>(zonep->ptrImpl);
  
  // Blending
  blend = 1.;
#include "selection2.h"

  E_Int eltType0 = zonep->eltType[0];
  bool is1D = false;
  if ((eltType0 == 1) || (eltType0 == 10 && zonep->nelts1D > 0)) is1D = true;
  
#ifdef __SHADERS__
  int curr = _shaders.currentShader();
  if (curr != 0) _shaders[curr]->setUniform("blend", (float)blend);
  glColor4f(0.,0.,0., blend); // pour imposer blend

  if (is1D)
  {
    if (curr != 0) _shaders[curr]->setUniform("lightOn", (int)0); // impose isoLight off on 1D meshes
  }
#endif
  
  glCallList(zoneImpl->_DLiso);

#ifdef __SHADERS__
  if (is1D)
  {
    if (ptrState->isoLight == 1 && ptrState->dim == 3)
    {
      if (curr != 0) _shaders[curr]->setUniform("lightOn", (int)1); // put back the isoLight value found in the CPlot state
    }
  }
#endif
}

//=============================================================================
/*
  Display une zone en iso solides pour un champ vectoriel.
  IN: zonep: pointeur sur la zone a afficher
  IN: zonet: le no de la zone dans liste globale des zones
  IN: nofield1, nofield2, nofield3: les no des champs a afficher.
*/
//=============================================================================
void DataDL::renderUIsoSolidZone(UnstructZone* zonep, E_Int zonet,
				 E_Int nofield1, E_Int nofield2, E_Int nofield3)
{
  float offb;
  double blend;
  
  ZoneImplDL* zoneImpl = static_cast<ZoneImplDL*>(zonep->ptrImpl);
  
  // Blending
  blend = 1.;
#include "selection2.h"

#ifdef __SHADERS__
  int curr = _shaders.currentShader();
  if (curr != 0) _shaders[curr]->setUniform("blend", (float)blend);
  glColor4f(0.,0.,0., blend); // pour imposer blend
#endif
  
  glCallList(zoneImpl->_DLiso);
}

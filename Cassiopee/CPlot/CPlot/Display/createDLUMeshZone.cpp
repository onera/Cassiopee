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
  Display une zone en mesh.
*/
//=============================================================================
void DataDL::createGPUUMeshZone(UnstructZone* zonep, E_Int zone, E_Int zonet)
{
  if (zonep->_is_high_order == true) 
  {
    createGPUUMeshZoneHO(zonep, zone, zonet);
    return;
  }
  E_Int i, n1, n2, ret1, ret2;
  ZoneImplDL* zImpl = static_cast<ZoneImplDL*>(zonep->ptr_impl);
  
  zImpl->_DLmesh = glGenLists(1);
  glNewList(zImpl->_DLmesh, GL_COMPILE);
#include "displayUMeshZone.h"
  glEndList();
}

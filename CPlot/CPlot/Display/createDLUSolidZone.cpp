/*    
    Copyright 2013-2019 Onera.

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
void DataDL::createGPUUSolidZone(UnstructZone* zonep, int zone, int zonet)
{
  int i, n1, n2, n3, n4, n5, n6, n7, n8;
  int ret1, ret2, ret3, ret4, f;
  ZoneImplDL* zImpl = static_cast<ZoneImplDL*>(zonep->ptr_impl);
  zImpl->_DLsolid = glGenLists(1);
  glNewList(zImpl->_DLsolid, GL_COMPILE);
#include "displayUSolidZone.h"
  glEndList();
}

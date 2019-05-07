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
#include <iostream>
#include "../DataDL.h"
#include "../ZoneImplDL.h"
#include "../Data.h"

#define PLOTHO                                    \
               glVertex3d(x[n1], y[n1], z[n1]);   \
               glVertex3d(x[n3], y[n3], z[n3]);   \
               glVertex3d(x[n2], y[n2], z[n2]);
#define PLOTBHO ret1 = _pref.blanking->f(this, n1, zonep->blank, zonet); \
                ret2 = _pref.blanking->f(this, n2, zonep->blank, zonet); \
                if (ret1*ret2 != 0) { PLOTHO; }

//=============================================================================
/*
  Display une zone en mesh.
*/
//=============================================================================
void DataDL::createGPUUMeshZone_ho(UnstructZone* zonep, int zone, int zonet)
{
  int i, n1, n2, n3, ret1, ret2;
  ZoneImplDL* zImpl = static_cast<ZoneImplDL*>(zonep->ptr_impl);

  zImpl->_DLmesh = glGenLists(1);
  glNewList(zImpl->_DLmesh, GL_COMPILE);
  #include "displayUMeshZone_ho.h"
  glEndList();
}

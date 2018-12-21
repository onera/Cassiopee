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
void DataDL::createGPUSMeshZone(StructZone* zonep, int zone)
{
  int i, n1, n2, j, k, plane;
  int stepi, stepj, stepk;
  int nis, njs, nks;
  int nie, nje, nke;
  int ret1, ret2;
  
  // Grid dimensions
  int ni = zonep->ni;
  int nj = zonep->nj;
  int nk = zonep->nk;
  if (ptrState->dim == 2) nk = 1;
  int nij = ni*nj;

  // steps
  stepi = 1; stepj = 1; stepk = 1;
  ZoneImplDL* zImpl = static_cast<ZoneImplDL*>(zonep->ptr_impl);

  zImpl->_DLmesh = glGenLists(1);
  glNewList(zImpl->_DLmesh, GL_COMPILE);
#include "displaySMeshZone.h"
  glEndList();
}

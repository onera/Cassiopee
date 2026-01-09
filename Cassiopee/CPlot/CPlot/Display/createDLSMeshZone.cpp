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
void DataDL::createGPUSMeshZone(StructZone* zonep, E_Int zone)
{
  E_Int i, n1, n2, j, k, plane;
  E_Int stepi, stepj, stepk;
  E_Int nis, njs, nks;
  E_Int nie, nje, nke;
  E_Int ret1, ret2;
  
  // Grid dimensions
  E_Int ni = zonep->ni;
  E_Int nj = zonep->nj;
  E_Int nk = zonep->nk;
  if (ptrState->dim == 2) nk = 1;
  E_Int nij = ni*nj;

  // steps
  stepi = 1; stepj = 1; stepk = 1;
  ZoneImplDL* zImpl = static_cast<ZoneImplDL*>(zonep->ptr_impl);

  zImpl->_DLmesh = glGenLists(1);
  glNewList(zImpl->_DLmesh, GL_COMPILE);
#include "displaySMeshZone.h"
  glEndList();
}

/*    
    Copyright 2013-2018 Onera.

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
  Cree la DL pour un champ scalaire.
  Cette DL depend de:
  - nofield
  - minf, maxf (min/max globaux de nofield)
  - la colormap (si pas de shader)
  - offb=0 (couleur de selection)

  isoMin et isoMax sont geres par le shader.
  Les steps sont forces a 1. La couleur de selection
  est forcee a 0.
*/
//=============================================================================
void DataDL::createGPUSIsoSolidZone(StructZone* zonep, int zone, int nofield)
{
  int stepi, stepj, stepk;
  float offb;
  stepi = 1; stepj = 1; stepk = 1;
  offb = 0.;
  ZoneImplDL* zImpl = static_cast<ZoneImplDL*>(zonep->ptr_impl);
  zImpl->_DLiso = glGenLists(1);
  glNewList(zImpl->_DLiso, GL_COMPILE);
#include "displaySIsoSolidZone.h"
  glEndList();
}

//=============================================================================
/*
  Cree la DL pour un champ vectoriel.
  Cette DL depend de:
  - nofield1, nofield2, nofield3
  - minf, maxf (min/max globaux de nofield1, nofield2, nofield3)
  - offb=0 (couleur de selection)

  isoMin et isoMax sont geres par le shader.
  Les steps sont forces a 1. La couleur de selection
  est forcee a 0.
*/
//=============================================================================
void DataDL::createGPUSIsoSolidZone(StructZone* zonep, int zone, int nofield1,
				    int nofield2, int nofield3)
{
  int stepi, stepj, stepk;
  float offb;
  stepi = 1; stepj = 1; stepk = 1;
  offb = 0.;
  ZoneImplDL* zImpl = static_cast<ZoneImplDL*>(zonep->ptr_impl);

  zImpl->_DLiso = glGenLists(1);
  glNewList(zImpl->_DLiso, GL_COMPILE);
#undef PLOT
#include "displaySVectSolidZone.h"
  glEndList();
}

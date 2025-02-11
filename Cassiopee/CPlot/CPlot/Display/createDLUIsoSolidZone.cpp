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

#include "../DataDL.h"

//=============================================================================
/*
  Cree la DL pour un champ scalaire.
  Cette DL depend de:
  - nofield
  - minf, maxf (min/max globaux de nofield)
  - la colormap (si pas de shader)
  - offb=0 (couleur de selection)
*/
//=============================================================================
void DataDL::createGPUUIsoSolidZone(UnstructZone* zonep, E_Int zone, E_Int zonet, 
				    E_Int nofield)
{  
  E_Int i, n1, n2, n3, n4;
  float r, g, b, offb;
  E_Int ret1, ret2, ret3, ret4, ff;
  offb = 0.;
  ZoneImplDL* zImpl = static_cast<ZoneImplDL*>(zonep->ptr_impl);

  zImpl->_DLiso = glGenLists(1);
  glNewList(zImpl->_DLiso, GL_COMPILE);
#include "displayUIsoSolidZone.h"
  glEndList();
}

//=============================================================================
/*
  Cree la DL pour un champ vectoriel.
  Cette DL depend de:
  - nofield1, nofield2, nofield3
  - minf, maxf (min/max globaux de nofield1, nofield2, nofield3)
  - offb=0 (couleur de selection)
*/
//=============================================================================
void DataDL::createGPUUIsoSolidZone(UnstructZone* zonep, E_Int zone, E_Int zonet, 
				    E_Int nofield1, E_Int nofield2, E_Int nofield3)
{  
  E_Int i, n1, n2, n3, n4;
  float r, g, b, offb;
  E_Int ret1, ret2, ret3, ret4, ff;
  offb = 0.;

  ZoneImplDL* zImpl = static_cast<ZoneImplDL*>(zonep->ptr_impl);
  zImpl->_DLiso = glGenLists(1);
  glNewList(zImpl->_DLiso, GL_COMPILE);
#undef PLOTTRI
#undef PLOTTRI2
#undef PLOTQUAD
#undef PLOTQUAD2
#undef PLOTNGON
  double* f1 = zonep->f[nofield1];
  double* f2 = zonep->f[nofield2];
  double* f3 = zonep->f[nofield3];
  double fmin1, fmax1, fmin2, fmax2, fmin3, fmax3;
  fmax1 = maxf[nofield1]; fmin1 = minf[nofield1];
  fmax2 = maxf[nofield2]; fmin2 = minf[nofield2];
  fmax3 = maxf[nofield3]; fmin3 = minf[nofield3];
#define GL_QUADS_ARE GL_TRIANGLES
#define PLOTQUAD PLOTQUADT
#define PLOTQUAD2 PLOTQUADT2
#include "displayUVectSolidZone.h"
  glEndList();
}

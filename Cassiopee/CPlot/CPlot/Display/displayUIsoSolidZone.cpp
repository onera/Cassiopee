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
#include "Data.h"

//=============================================================================
/*
  Display une zone en iso solides pour les champs scalaires.
  IN: zonep: pointeur sur la zone a afficher
  IN: zonet: le no de la zone dans liste globale des zones
  IN: nofield: le no du champ a afficher.
*/
//=============================================================================
void Data::displayUIsoSolidZone(UnstructZone* zonep, E_Int zonet, E_Int nofield)
{
  E_Int i, n1, n2, n3, n4;
  float r, /*g, b,*/ offb=0.;
  double blend;
  E_Int ret1, ret2, ret3, ret4, ff=0;

  // Blending
  blend = 1.;
#include "selection2.h"
 
#ifdef __SHADERS__
  int curr = _shaders.currentShader();
  if (curr != 0) _shaders[curr]->setUniform("blend", (float)blend);
  glColor4f(0.,0.,0., blend); // pour imposer blend
#endif

#include "displayUIsoSolidZone.h"
}

//=============================================================================
/*
  Display une zone en iso solides pour les champs vectoriels.
  IN: zonep: pointeur sur la zone a afficher
  IN: zonet: le no de la zone dans liste globale des zones
  IN: nofield1, nofield2, nofield3: les no des champs a afficher.
*/
//=============================================================================
void Data::displayUIsoSolidZone(UnstructZone* zonep, E_Int zonet,
                                E_Int nofield1, E_Int nofield2, E_Int nofield3)
{
  E_Int i, n1, n2, n3, n4;
  float r, g, b, offb=0.;
  double blend;
  E_Int ret1, ret2, ret3, ret4, ff=0;

  // Blending
  blend = 1.;
#include "selection2.h"

#ifdef __SHADERS__
  int curr = _shaders.currentShader();
  if (curr != 0) _shaders[curr]->setUniform("blend", (float)blend);
  glColor4f(0.,0.,0., blend); // pour imposer blend
#endif

#undef PLOTTRI
#undef PLOTTRI2
#undef PLOTQUAD
#undef PLOTQUAD2
#undef PLOTNGON
#undef PLOTNGON2
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
}

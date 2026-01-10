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
  IN: zone: le no de la zone dans la liste globale des zones
  IN: nofield: le no du champ a afficher.
  Cette fonction n'utilise pas la DL.
*/
//=============================================================================
void Data::displaySIsoSolidZone(StructZone* zonep, E_Int zone,
                                E_Int nofield)
{
  E_Int stepi, stepj, stepk;
  double blend; float offb;
  computeSteps(zonep, stepi, stepj, stepk);

  // Blending
  blend = 1.;
#include "selection2.h"

#ifdef __SHADERS__
  int curr = _shaders.currentShader();
  if (curr != 0) _shaders[curr]->setUniform("blend", (float)blend);
  glColor4f(0.,0.,0., blend); // pour imposer blend
#endif

#include "displaySIsoSolidZone.h"
}

//=============================================================================
/*
  Display une zone en iso solides pour les champs vectoriels.
  IN: zonep: pointeur sur la zone a afficher
  IN: zone: le no de la zone dans la liste globale des zones
  IN: nofield1, nofield2, nofield3: le no des champs a afficher.
  Cette fonction n'utilise pas la DL.
*/
//=============================================================================
void Data::displaySIsoSolidZone(StructZone* zonep, E_Int zone,
                                E_Int nofield1, E_Int nofield2, E_Int nofield3)
{
  E_Int stepi, stepj, stepk;
  double blend; float offb;
  computeSteps(zonep, stepi, stepj, stepk);

  // Blending
  blend = 1.;
#include "selection2.h"

#ifdef __SHADERS__
  int curr = _shaders.currentShader();
  if (curr != 0) _shaders[curr]->setUniform("blend", (float)blend);
  glColor4f(0.,0.,0., blend); // pour imposer blend
#endif
  
#undef PLOT
  double* f1 = zonep->f[nofield1];
  double* f2 = zonep->f[nofield2];
  double* f3 = zonep->f[nofield3];
  double fmin1, fmax1, fmin2, fmax2, fmin3, fmax3;
  fmax1 = maxf[nofield1]; fmin1 = minf[nofield1];
  fmax2 = maxf[nofield2]; fmin2 = minf[nofield2];
  fmax3 = maxf[nofield3]; fmin3 = minf[nofield3];
#define GL_QUADS_ARE GL_TRIANGLES
#define PLOT PLOTT
#include "displaySVectSolidZone.h"
  
  // Pour les lignes
  if (nij == 1 || ni*nk == 1 || nj*nk == 1)
  {
    glBegin(GL_LINES);
    E_Int nie, nje, nke;
    nie = ni; nje = nj; nke = nk;
    if (ni*nj == 1) nke = nke-1;
    if (ni*nk == 1) nje = nje-1;
    if (nj*nk == 1) nie = nie-1;
    if (zonep->blank == 0)
    {
      // No blanking
      for (k = 0; k < nke; k++)
        for (j = 0; j < nje; j++)
          for (i = 0; i < nie; i++)
          {
            n1 = i+j*ni+k*nij;
            n2 = n1+1;
            glColor4f(0, 0, 0+offb, blend);
            glVertex3d(x[n1], y[n1], z[n1]);
            glColor4f(0, 0, 0+offb, blend);
            glVertex3d(x[n2], y[n2], z[n2]);
          }
    }
    else
    {
      for (k = 0; k < nke; k++)
        for (j = 0; j < nje; j++)
          for (i = 0; i < nie; i++)
          {
            n1 = i+j*ni+k*nij;
            n2 = n1+1;
            ret1 = _pref.blanking->f(this, n1, zonep->blank, zone);
            ret2 = _pref.blanking->f(this, n2, zonep->blank, zone);
            if (ret1*ret2 != 0)
            { 
              glColor4f(0, 0, 0+offb, blend);
              glVertex3d(x[n1], y[n1], z[n1]);
              glColor4f(0, 0, 0+offb, blend);
              glVertex3d(x[n2], y[n2], z[n2]);
            }
          }
    }
    glEnd();
  }
}

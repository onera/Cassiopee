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
   Display an unstructured mesh lines.
   Only for BAR and 1D NGON
*/
//=============================================================================
void Data::displayUEdges()
{
  if (_numberOfUnstructZones == 0) return;
  if (ptrState->edgifyDeactivatedZones == 0 && ptrState->edgifyActivatedZones == 0)
    return;
  E_Int zone;

  glColor4f(0.8, 0.8, 0.8, 1.);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_LINE_SMOOTH);
  glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);

  E_Int ne, n1, n2, elt;

  zone = 0;
  while (zone < _numberOfUnstructZones)
  {
    UnstructZone* zonep = _uzones[zone];
    E_Int eltType = zonep->eltType[0];
    if ((eltType == 1) || (eltType == 10 && zonep->nelts1D > 0))
    {
      if (isInFrustum(zonep, _view) == 1)
      {
        // if zone is active and in frustum
        if ((zonep->active == 1 && ptrState->edgifyActivatedZones == 1) ||
            (zonep->active == 0 && ptrState->edgifyDeactivatedZones == 1))        
        {
          double* x = zonep->x; double* y = zonep->y; double* z = zonep->z;
          E_Int* connect = zonep->connect[0];
          if (eltType == 1) // BAR
          {
            ne = zonep->ne;
            glBegin(GL_LINES);
            for (E_Int i = 0; i < ne; i++)
            {
              n1 = connect[i]-1;
              n2 = connect[i+ne]-1;
              glVertex3d(x[n1], y[n1], z[n1]);
              glVertex3d(x[n2], y[n2], z[n2]);
            }
            glEnd();
          }
          if (eltType == 10) // NGON1D
          {
            glBegin(GL_LINES);
            for (E_Int i = 0; i < zonep->nelts1D; i++)
            {
              elt = zonep->posElts1D[i];
              E_Int* ptrelt = &connect[elt];
              E_Int face = ptrelt[1]-1; // indice de la face
              E_Int* ptrface = &connect[zonep->posFaces[face]];
              n1 = ptrface[1]-1;
              face = ptrelt[2]-1; // indice de la face
              ptrface = &connect[zonep->posFaces[face]];
              n2 = ptrface[1]-1;
              glVertex3d(x[n1], y[n1], z[n1]);
              glVertex3d(x[n2], y[n2], z[n2]);
            }
            glEnd();
          }
        }
      }
    }
    zone++;
  }
 
  glDisable(GL_LINE_SMOOTH);
  glDisable(GL_BLEND);
  ptrState->alpha = 1.; 
  glColor3f(1., 1., 1.);
}


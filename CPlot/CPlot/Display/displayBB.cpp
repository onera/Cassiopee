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

#include "../Data.h"

//============================================================================
/*
  Display bounding box of model in wire frame (stipple).
  Red box.
*/
//============================================================================
void Data::displayBB2()
{
  glColor3f(1., 0, 0);
  glLineWidth(2.);
  glLineStipple(2, 0x3F07);
  glEnable(GL_LINE_STIPPLE);
  double xmi, xma, ymi, yma, zmi, zma;

  if (ptrState->selectedZone == 0)
  {
    xmi = xmin; ymi = ymin; zmi = zmin;
    xma = xmax; yma = ymax; zma = zmax;
  }
  else
  {
    xmi = 1.e6; ymi = 1.e6; zmi = 1.e6;
    xma = -1.e6; yma = -1.e6; zma = -1.e6;
    for (int i = 0; i < _numberOfZones; i++)
    {
      Zone* zonep = _zones[i];
      if (zonep->selected == 1)
      {
        xmi = MIN(zonep->xmin, xmi);
        xma = MAX(zonep->xmax, xma);
        ymi = MIN(zonep->ymin, ymi);
        yma = MAX(zonep->ymax, yma);
        zmi = MIN(zonep->zmin, zmi);
        zma = MAX(zonep->zmax, zma);
      }
    }
  }
    
  glBegin(GL_LINES);
  glVertex3d(xmi, ymi, zmi);
  glVertex3d(xma, ymi, zmi);
  glVertex3d(xmi, ymi, zmi);
  glVertex3d(xmi, yma, zmi);
  glVertex3d(xma, ymi, zmi);
  glVertex3d(xma, yma, zmi);
  glVertex3d(xmi, yma, zmi);
  glVertex3d(xma, yma, zmi);
  glVertex3d(xmi, ymi, zma);
  glVertex3d(xma, ymi, zma);
  glVertex3d(xmi, ymi, zma);
  glVertex3d(xmi, yma, zma);
  glVertex3d(xma, ymi, zma);
  glVertex3d(xma, yma, zma);
  glVertex3d(xmi, yma, zma);
  glVertex3d(xma, yma, zma);
  glVertex3d(xmi, ymi, zmi);
  glVertex3d(xmi, ymi, zma);
  glVertex3d(xma, ymi, zmi);
  glVertex3d(xma, ymi, zma);
  glVertex3d(xmi, yma, zmi);
  glVertex3d(xmi, yma, zma);
  glVertex3d(xma, yma, zmi);
  glVertex3d(xma, yma, zma);
  glEnd();
  
  glLineWidth(1.0);
  glDisable(GL_LINE_STIPPLE);
  glColor3f(1., 1., 1.);
}
//============================================================================
/*
  Display bounding box of model.
  Ne display que les coins de la BB.
  Red box.
*/
//============================================================================
void Data::displayBB()
{
  glColor3f(1., 0, 0);
  double xmi, xma, ymi, yma, zmi, zma;
  
  if (ptrState->selectedZone == 0)
  {
    xmi = xmin; ymi = ymin; zmi = zmin;
    xma = xmax; yma = ymax; zma = zmax;
  }
  else
  {
    xmi = 1.e6; ymi = 1.e6; zmi = 1.e6;
    xma = -1.e6; yma = -1.e6; zma = -1.e6;
    for (int i = 0; i < _numberOfZones; i++)
    {
      Zone* zonep = _zones[i];
      if (zonep->selected == 1)
      {
        xmi = MIN(zonep->xmin, xmi);
        xma = MAX(zonep->xmax, xma);
        ymi = MIN(zonep->ymin, ymi);
        yma = MAX(zonep->ymax, yma);
        zmi = MIN(zonep->zmin, zmi);
        zma = MAX(zonep->zmax, zma);
      }
    }
  }

  double dx = xma - xmi;
  double dy = yma - ymi;
  double dz = zma - zmi;
  double alpha = 0.10;

  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_LINE_SMOOTH);
  glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
  glDepthMask(GL_FALSE);

  glBegin(GL_LINES);
  glVertex3d(xmi, ymi, zmi);
  glVertex3d(xmi+alpha*dx, ymi, zmi);
  glVertex3d(xmi, ymi, zmi);
  glVertex3d(xmi, ymi+alpha*dy, zmi);
  glVertex3d(xmi, ymi, zmi);
  glVertex3d(xmi, ymi, zmi+alpha*dz);

  glVertex3d(xma, ymi, zmi);
  glVertex3d(xma-alpha*dx, ymi, zmi);
  glVertex3d(xma, ymi, zmi);
  glVertex3d(xma, ymi+alpha*dy, zmi);
  glVertex3d(xma, ymi, zmi);
  glVertex3d(xma, ymi, zmi+alpha*dz);
  
  glVertex3d(xmi, yma, zmi);
  glVertex3d(xmi+alpha*dx, yma, zmi);
  glVertex3d(xmi, yma, zmi);
  glVertex3d(xmi, yma-alpha*dy, zmi);
  glVertex3d(xmi, yma, zmi);
  glVertex3d(xmi, yma, zmi+alpha*dz);

  glVertex3d(xmi, ymi, zma);
  glVertex3d(xmi+alpha*dx, ymi, zma);
  glVertex3d(xmi, ymi, zma);
  glVertex3d(xmi, ymi+alpha*dy, zma);
  glVertex3d(xmi, ymi, zma);
  glVertex3d(xmi, ymi, zma-alpha*dz);

  glVertex3d(xma, yma, zmi);
  glVertex3d(xma-alpha*dx, yma, zmi);
  glVertex3d(xma, yma, zmi);
  glVertex3d(xma, yma-alpha*dy, zmi);
  glVertex3d(xma, yma, zmi);
  glVertex3d(xma, yma, zmi+alpha*dz);

  glVertex3d(xma, ymi, zma);
  glVertex3d(xma-alpha*dx, ymi, zma);
  glVertex3d(xma, ymi, zma);
  glVertex3d(xma, ymi+alpha*dy, zma);
  glVertex3d(xma, ymi, zma);
  glVertex3d(xma, ymi, zma-alpha*dz);

  glVertex3d(xmi, yma, zma);
  glVertex3d(xmi+alpha*dx, yma, zma);
  glVertex3d(xmi, yma, zma);
  glVertex3d(xmi, yma-alpha*dy, zma);
  glVertex3d(xmi, yma, zma);
  glVertex3d(xmi, yma, zma-alpha*dz);

  glVertex3d(xma, yma, zma);
  glVertex3d(xma-alpha*dx, yma, zma);
  glVertex3d(xma, yma, zma);
  glVertex3d(xma, yma-alpha*dy, zma);
  glVertex3d(xma, yma, zma);
  glVertex3d(xma, yma, zma-alpha*dz);
  
  glEnd();
  glDepthMask(GL_TRUE);
  glDisable(GL_LINE_SMOOTH);
  glDisable(GL_BLEND);
  glColor3f(1., 1., 1.);
}

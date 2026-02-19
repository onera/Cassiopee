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
  Display la bbox d'une zone structuree 
  IN: zonep: pointeur sur la zone a afficher
*/
//=============================================================================
void Data::displaySBBZone(StructZone* zonep)
{
  glColor3f(1., 0, 0);
  glLineWidth(2.);

  double xmi = zonep->xmin;
  double ymi = zonep->ymin;
  double zmi = zonep->zmin;
  double xma = zonep->xmax;
  double yma = zonep->ymax;
  double zma = zonep->zmax;

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
  glColor3f(1., 1., 1.);
}

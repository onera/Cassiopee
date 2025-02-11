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
#include "../Data.h"
#include <math.h>

//=============================================================================
/*
  Trigger fog.
*/
//=============================================================================
void Data::fog()
{
  double d, lx;

  glEnable(GL_FOG);
  {
    //GLfloat fogColor[4] = {0.05, 0.05, 0.2, 0.5};
    GLfloat fogColor[4] = {0.5, 0.5, 0.5, 0.1};
    
    GLint fogMode = GL_LINEAR;
    glFogi (GL_FOG_MODE, fogMode);
    glFogfv (GL_FOG_COLOR, fogColor);
    glFogf (GL_FOG_DENSITY, 0.35);
    glHint (GL_FOG_HINT, GL_DONT_CARE);

    // Starting fog distance
    d = dist2BB(_view.xcam, _view.ycam, _view.zcam,
                xmin, ymin, zmin,
                xmax, ymax, zmax);
    d = sqrt(d);
    glFogf (GL_FOG_START, d);
    // Max fog at this distance
    lx = MAX( xmax-xmin, ymax-ymin);
    lx = MAX( lx, zmax-zmin );
    d = d+0.7*lx;
    glFogf (GL_FOG_END, d);
  }
}

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

//=============================================================================
// Change the mouse cursor
// cursor=0 (default), 1 (cross), 2 (wait)
// Cette fonction doit etre appelee de CPlot et non de l'interface tkInter.
//=============================================================================
void Data::setCursor(E_Int type)
{
  switch (type)
  {
    case 0:
      // INHERIT fails for some window manager
      //glutSetCursor(GLUT_CURSOR_INHERIT);
      glutSetCursor(GLUT_CURSOR_LEFT_ARROW);
      break;
    case 1: 
      glutSetCursor(GLUT_CURSOR_CROSSHAIR); 
      break;
    case 2: 
      glutSetCursor(GLUT_CURSOR_WAIT); 
      break;
    default:
      //glutSetCursor(GLUT_CURSOR_INHERIT);
      glutSetCursor(GLUT_CURSOR_LEFT_ARROW);
      break;
  }
  ptrState->cursorType = type;
}

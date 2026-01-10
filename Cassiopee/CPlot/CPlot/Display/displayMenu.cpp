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
#include "../Data.h"

// Affichage des menus sur l'ecran

//=============================================================================
// Display menu
//=============================================================================
void Data::displayMenu()
{
#define NMENU 4
  E_Int pos, posn, i, l;
  E_Int width[NMENU];

  setOrthographicProjection();
  glPushMatrix();
  glLoadIdentity();

  pos = 0;
  posn = 0;
  i = 0;
  displayDimensionMenu(&posn); width[i] = posn - pos; pos = posn; i++;
  displayVariableMenu(&posn); width[i] = posn - pos; pos = posn; i++;
  displayAxisMenu(&posn); width[i] = posn - pos; pos = posn; i++;
  displayZoneMenu(&posn); width[i] = posn - pos; pos = posn; i++;
  
  // Render the rectangle
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glColor4f(0.0, 0.0, 1., 0.3);
  glBegin(GL_QUADS);
  l = 0;
  for (i = 0; i < NMENU; i++)
  {
    glVertex3d(l+2, _view.h, 0);
    glVertex3d(l + width[i] - 2, _view.h, 0);
    glVertex3d(l + width[i] - 2, _view.h - 32, 0);
    glVertex3d(l+2, _view.h - 32, 0);
    l = l + width[i];
  }
  glEnd();
  glDisable(GL_BLEND);

  glPopMatrix();
  resetPerspectiveProjection(); 
}

//=============================================================================
// Display dimension menu
//=============================================================================
void Data::displayDimensionMenu(int* x)
{
  E_Int posx = *x + 4;
  E_Int posy = _view.h - 20;

  switch (ptrState->dim)
  {
    case 1:
      renderBitmapString(posx, posy, 0, _fontSize1, "1D");
      break;
      
    case 2:
      renderBitmapString(posx, posy, 0, _fontSize1, "2D");
      break;

    case 3:
      renderBitmapString(posx, posy, 0, _fontSisze1, "3D");
      break;
  }
  posy = _view.h - 5;
  renderBitmapString(posx, posy, 0, _fontSize1, "(m)");
      
  *x = *x + 12*2 + 4;
}

//=============================================================================
// Display variable menu
//=============================================================================
void Data::displayVariableMenu(int* x)
{
  char msg[30];
  E_Int posx = *x + 4;
  E_Int posy = _view.h - 20;

  switch (ptrState->mode)
  {
    case MESH:
      strcpy(msg, "mesh");
      break;

    case SOLID:
      strcpy(msg, "solid");
      break;
      
    case RENDER:
      strcpy(msg, "render");
      break;

    case SCALARFIELD:
      strcpy(msg, _zones[0]->varnames[ptrState->scalarField]);
      break;

    case VECTORFIELD:
      strcpy(msg, _zones[0]->varnames[ptrState->vectorField1]);
      strcat(msg, ", ");
      strcat(msg, _zones[0]->varnames[ptrState->vectorField2]);
      strcat(msg, ", ");
      strcat(msg, _zones[0]->varnames[ptrState->vectorField3]);
      break;

    default:
      strcpy(msg, "unknown mode");
      break;
  }
  renderBitmapString(posx, posy, 0, _fontSize1, msg);

  posx = posx + strlen(msg);
  posy = _view.h - 5;
  renderBitmapString(posx, posy, 0, _fontSize1, "(1)");
      
  *x = *x + 8*strlen(msg) + 4;
}

//=============================================================================
// Display axis
//=============================================================================
void Data::displayAxisMenu(int* x)
{
  char msg[30];
  E_Int posx = *x + 4;
  E_Int posy = _view.h - 20;

  if (ptrState->dim != 1)
  {
    // 2D or 3D mode
    if (ptrState->dim == 2)
    {
      switch (ptrState->var2D)
      {
        case 0:
          strcpy(msg, "(x,y)");
          break;
        case 1:
          strcpy(msg, "(x,z)");
          break;
        case 2:
          strcpy(msg, "(y,z)");
          break;
      }
    }
    else
      strcpy(msg, "(x,y,z)");
  }
  else
  {
    switch (ptrState->var1D)
    {
      case 0:
        strcpy(msg, "(x,f)");
        break;
      case 1:
        strcpy(msg, "(y,f)");
        break;
      case 2:
        strcpy(msg, "(z,f)");
        break;
      case 3:
        strcpy(msg, "(s,f)");
        break;
    }
  }
  renderBitmapString(posx, posy, 0, _fontSize1, msg);

  posx = posx + strlen(msg);
  posy = _view.h - 5;
  renderBitmapString(posx, posy, 0, _fontSize1, "(2)");
      
  *x = *x + 8*strlen(msg) + 4;
}

//=============================================================================
// Display zone
//=============================================================================
void Data::displayZoneMenu(int* x)
{
  char msg[100];
  E_Int posx = *x + 4;
  E_Int posy = _view.h - 20;
  if (ptrState->selectedZone > 0)
  {
    E_Int nz = ptrState->selectedZone-1;
    strcpy(msg, _zones[nz]->zoneName);
  }
  else
  {
    strcpy(msg, "Select zone");
  }
  renderBitmapString(posx, posy, 0, _fontSize1, msg);

  posx = posx + strlen(msg) - 5;
  posy = _view.h - 5;
  renderBitmapString(posx, posy, 0, _fontSize1, "(z-Z ; a-A)");
      
  *x = *x + 8*MAX(strlen(msg), strlen("(z-Z ; a-A)")) + 4;
}

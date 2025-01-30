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
#include "../Particles.h"

//=============================================================================
// Note on menus :
// To go in a submenu : left/right arrows or space bar
// To leave a submenu to cplot : esc
// To go up a submenu : backspace or back
// To toggle values : left/right arrows
// To check value : space bar (ex : compute a new field, display the changes)
// To exit : spacebar on exit
//=============================================================================

// Keyboard documentation text
static int keybTextSize = 20;
static const char* keybText[] = {
  "Arrows: move around",
  "o: look up",
  "p: look down",
  "f: fit view",
  "shift + arrows: straff",
  "shift + key: reverse effect of key",
  "ctrl + key: change key function",
  "m: toggle between 1D, 2D, 3D modes",
  "1: toggle primary variable",
  "2: toggle secondary variable",
  "3: toggle i,j,k mode",
  "space bar: mesh, surface display",
  "a: activate/deactivate zone",
  "z: select zone",
  "l: look for zone",
  "i,j,k: change i,j,k displayed plane",
  "c: change appearance",
  "s: toggle smoke mode",
  "r: reload files",
  "esc: display menu "
 };
// Credit text
//static const char* creditsText[] = {"CPlot (c) 2005-2016 by ONERA"};

// Different kind of values for menu
#define NONE     0
#define INT      1
#define FLOAT    2
#define PTR_VOID 3
#define PTR_INT  4
#define PTR_DOUBLE 5
#define PTR_VOID2 6
#define PTR_VOID3 7
#define TEXT_DATA 8
#define TEXT_KEYBOARD 9
// Number of items in menuItem
#define NITEM   25

static int start;
static int end = 4;
static int currentMenu = 0;
// Describe the ptrs of each menu and under-menu
static int menuSize[] = {4,6,13,15,17};
// Describe the button for each menu item of each menu
static int order[] = {12,10,11,6,     // main menu
                      1,9,            // pref menu
                      2,3,17,20,18,19,9,  // plugin menu
                      22,9,              // data set info menu
                      23,9};             // keyboard menu
static const char* submenuNames[] = {"- CPlot Menu -",
				     "- Preferences -",
				     "- Plugins -",
				     "- Data Set Info -",
				     "- Keyboard -"};

//=============================================================================
/*
   Menu.
   Display information on software state.
*/
//=============================================================================
void Data::menu()
{
  static const char* menuItem[] = {
  "Keyboard settings",                   // 0
  "Drawing accuracy : ",                 // 1
  "Blanking function : ",                // 2
  "Look-for function (l): ",             // 3
  "Add-a-variable function : ",          // 4
  "Credits",                             // 5
  "Exit",                                // 6
  "Arrows : move around",                // 7
  "Shift + arrows : straff",             // 8
  "Back",                                // 9
  "Preferences",                         // 10
  "Plugins",                             // 11
  "Data set info",                       // 12
  "Max number of particles : ",          // 13
  "Smoke particles radius : ",           // 14
  "Smoke time step : ",                  // 15
  "Smoke emission zone radius : ",       // 16
  "Selection function : ",               // 17
  "Colormap function : ",                // 18
  "Screen dump function : ",             // 19
  "Mouse click function : ",             // 20
  "Texture size i : ",                   // 21
  "Data set info text",                  // 22
  "Keyboard text",                       // 23
  "Texture size j : "                    // 24
  };

  static int menuItemType[NITEM];
  static void* menuItemValue[NITEM];

  char msg[120];
  E_Int l = 50; // current line
  E_Int i, size, o;

  // Init menu item values
  menuItemType[0] = NONE; menuItemValue[0] = NULL;
  menuItemType[1] = INT; menuItemValue[1] = &(_pref.speed);
  menuItemType[2] = PTR_INT; menuItemValue[2] = _pref.blanking;
  menuItemType[3] = PTR_VOID; menuItemValue[3] = _pref.lookFor;
  menuItemType[4] = PTR_VOID; menuItemValue[4] = _pref.addAVariable;
  menuItemType[5] = NONE; menuItemValue[5] = NULL;
  menuItemType[6] = NONE; menuItemValue[6] = NULL;
  menuItemType[7] = NONE; menuItemValue[7] = NULL;
  menuItemType[8] = NONE; menuItemValue[8] = NULL;
  menuItemType[9] = NONE; menuItemValue[9] = NULL;
  menuItemType[10] = NONE; menuItemValue[10] = NULL;
  menuItemType[11] = NONE; menuItemValue[11] = NULL;
  menuItemType[12] = NONE; menuItemValue[12] = NULL;
  menuItemType[13] = INT; menuItemValue[13] = &(_pref.maxParticles);
  menuItemType[14] = FLOAT; menuItemValue[14] = &(_pref.smokeRadius);
  menuItemType[15] = FLOAT; menuItemValue[15] = &(_pref.smokeTimeStep);
  menuItemType[16] = FLOAT; menuItemValue[16] = &(_pref.emissionRadius);
  menuItemType[17] = PTR_VOID; menuItemValue[17] = _pref.selectNextZone;
  menuItemType[18] = PTR_DOUBLE; menuItemValue[18] = _pref.colorMap;
  menuItemType[19] = PTR_VOID2; menuItemValue[19] = _pref.screenDump;
  menuItemType[20] = PTR_VOID3; menuItemValue[20] = _pref.mouseClick;
  menuItemType[21] = INT; menuItemValue[21] = &(_pref.sizeTexi);
  menuItemType[22] = TEXT_DATA; menuItemValue[22] = NULL;
  menuItemType[23] = TEXT_KEYBOARD; menuItemValue[23] = NULL;
  menuItemType[24] = INT; menuItemValue[24] = &(_pref.sizeTexj);

  // Set new functions
  ptrState->inmenu = 1;
  glutKeyboardFunc(gmenuKeyboard);
  glutSpecialFunc(gmenuArrows);

  // Disable display/reread if any
  if (ptrState->smoke == 1) glutIdleFunc(NULL);

  // Enable blending
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  ptrState->alpha = 0.3;

  // Display in alpha blending
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  ptrState->render = 0;
  display();
  ptrState->render = 1;

  // Display menu title
  E_Int sizeMax = 0;
  glColor3f(0.1, 0.1, 0.3);
  strcpy(msg, submenuNames[currentMenu]);
  size = strlen(msg);
  sizeMax = MAX(sizeMax, size);
  displayBigText((int)(_view.w/2-size*0.5*9), l, msg);
  l = l + 50;

  for (i = start; i < end; i++)
  {
    o = order[i];
    displayMenuString(o, (char*)menuItem[o], &l, &sizeMax,
                      menuItemType[o], menuItemValue[o]);
  }

  // Refresh
  //glutSwapBuffers();

  ptrState->alpha = 1.;
  glDisable(GL_BLEND);

  // Render a rectangle
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  setOrthographicProjection();
  glPushMatrix();
  glLoadIdentity();

  glColor4f(0.9, 0.9, 1., 0.4);
  glBegin(GL_QUADS);
  glVertex3d((int)(_view.w/2-sizeMax*0.5*9-10), 20, 0);
  glVertex3d((int)(_view.w/2+sizeMax*0.5*9+10), 20, 0);
  glVertex3d((int)(_view.w/2+sizeMax*0.5*9+10), 70, 0);
  glVertex3d((int)(_view.w/2-sizeMax*0.5*9-10), 70, 0);
  glEnd();

  glColor4f(0.0, 0.0, 1., 0.4);
  glBegin(GL_QUADS);
  glVertex3d((int)(_view.w/2+sizeMax*0.5*9+10), 70, 0);
  glVertex3d((int)(_view.w/2-sizeMax*0.5*9-10), 70, 0);
  glVertex3d((int)(_view.w/2-sizeMax*0.5*9-10), l, 0);
  glVertex3d((int)(_view.w/2+sizeMax*0.5*9+10), l, 0);
  glEnd();
  //glEnable(GL_DEPTH_TEST);
  glPopMatrix();
  resetPerspectiveProjection();
  glDisable(GL_BLEND);
  glutSwapBuffers();

}

//=============================================================================
// Menu Keyboards call
//=============================================================================
// Branche les differents appels pour le clavier
void gmenuKeyboard(unsigned char key, int x, int y)
{
  Data* d = Data::getInstance();
  d->menuKeyboard(key, x, y);
}
void Data::menuKeyboard(unsigned char key, int x, int y)
{
  //printf("%d\n",key);

  switch (key)
  {
    // Menu - display state information
    case 27: // esc key
      ptrState->inmenu = 0;
      glutKeyboardFunc(gkeyboard);
      glutSpecialFunc(garrows);
      glutTimerFunc(ptrState->timeStep, gtimer, 0);
      if (ptrState->smoke == 1)
        glutIdleFunc(gdisplay);
      display();
      break;

    case 127: // back key
    case 8:
      if (currentMenu == 0)
      {
        ptrState->inmenu = 0;
        glutKeyboardFunc(gkeyboard);
        glutSpecialFunc(garrows);
        glutTimerFunc(ptrState->timeStep, gtimer, 0);
        if (ptrState->smoke == 1)
          glutIdleFunc(gdisplay);
        display();
      }
      else
      {
        switch (currentMenu)
        {
          case 1:
            ptrState->currentMenuItem = 1;
            break;
          case 2:
            ptrState->currentMenuItem = 2;
            break;
          case 3:
             ptrState->currentMenuItem = 0;
             break;
          case 4:
            ptrState->currentMenuItem = 3;
          default:
            ptrState->currentMenuItem = 0;
        }
        currentMenu = 0;
        start = 0;
        end = menuSize[currentMenu];
        menu();
      }
      break;

    case 32: // space bar or enter
    case 13:
      changeMenuItemSpaceBar();
      break;
  }
}

//=============================================================================
// Branche les appels suivants les fleches
//=============================================================================
void gmenuArrows(int key, int x, int y)
{
  Data* d = Data::getInstance();
  d->menuArrows(key, x, y);
}
void Data::menuArrows(int key, int x, int y)
{
  switch (key)
  {
    case GLUT_KEY_UP: // Move through menu
      ptrState->currentMenuItem--;
      if (ptrState->currentMenuItem < start)
        ptrState->currentMenuItem = end-1;
      menu();
      break;

    case GLUT_KEY_DOWN: // Move through menu
      ptrState->currentMenuItem++;
      if (ptrState->currentMenuItem > end-1)
        ptrState->currentMenuItem = start;
      menu();
      break;

    case GLUT_KEY_LEFT: // Change menu settings or select
      changeMenuItemLeft();
      break;

    case GLUT_KEY_RIGHT:
      changeMenuItemRight();
      break;
  }
}

//=============================================================================
/*
  Display menu string
  no : the number of the menu in list
  msg : the msg string
  l : the line number
  type : 0 no value menu, 1 int value menu, 2 ptr on chain_function_void
  value menu, 3 ptr on chain_function_int menu.
  value : the menu value
*/
//=============================================================================
void Data::displayMenuString(E_Int no, char* msg, E_Int* l, E_Int* sizeMax,
                             E_Int type, void* value)
{
  static char loc[120];
  static char loc2[20];
  E_Int size;
  E_Int* vali;
  float* valf;
  struct chain_function_void* v1;
  struct chain_function_int* v2;
  struct chain_function_double* v3;
  struct chain_function_void2* v4;
  struct chain_function_void3* v5;
  E_Int n;

  if (order[ptrState->currentMenuItem] == no) glColor3f(0.8, 0.5, 0.2);
  else glColor3f(1.0, 1.0, 1.0);

  switch (type)
  {
    case NONE: // no value menu
      strcpy(loc, msg);
      break;

    case INT:
      vali = (E_Int*)(value);
      strcpy(loc, msg);
      sprintf(loc2, SF_D_, *vali);
      strcat(loc, loc2);
      break;

    case FLOAT:
      valf = (float*)(value);
      strcpy(loc, msg);
      sprintf(loc2, SF_F_, *valf);
      strcat(loc, loc2);
      break;

    case PTR_VOID:
      v1 = (struct chain_function_void*)(value);
      strcpy(loc, msg);
      if (value == NULL)
        strcat(loc, "none");
      else
        strcat(loc, v1->functionName);
      break;

    case PTR_VOID2:
      v4 = (struct chain_function_void2*)(value);
      strcpy(loc, msg);
      if (value == NULL)
        strcat(loc, "none");
      else
        strcat(loc, v4->functionName);
      break;

    case PTR_VOID3:
      v5 = (struct chain_function_void3*)(value);
      strcpy(loc, msg);
      if (value == NULL)
        strcat(loc, "none");
      else
        strcat(loc, v5->functionName);
      break;

    case PTR_INT:
      v2 = (struct chain_function_int*)(value);
      strcpy(loc, msg);
      if (value == NULL)
        strcat(loc, "none");
      else
        strcat(loc, v2->functionName);
      break;

    case PTR_DOUBLE:
      v3 = (struct chain_function_double*)(value);
      strcpy(loc, msg);
      if (value == NULL)
        strcat(loc, "none");
      else
        strcat(loc, v3->functionName);
      break;
  }

  switch (type)
  {
    case TEXT_DATA:
      // Data set info text
      n = 0;
      for (E_Int i = 0; i < _numberOfZones; i++)
        n = n + _zones[i]->npts;
      sprintf(loc, "Total number of points : " SF_D_, n);
      size = strlen(loc);
      *sizeMax = MAX(size, *sizeMax);
      displayBigText((int)(_view.w/2-size*0.5*9), *l, loc);
      *l = *l + 20;

      sprintf(loc, "Total number of blocks : " SF_D_, _numberOfZones);
      size = strlen(loc);
      *sizeMax = MAX(size, *sizeMax);
      displayBigText((int)(_view.w/2-size*0.5*9), *l, loc);
      *l = *l + 20;

      for (int i = 3; i < _zones[0]->nfield; i++)
      {
        sprintf(loc, "%s : min=%f, max=%f", _zones[0]->varnames[i],
                minf[i-3], maxf[i-3]);
        size = strlen(loc);
        *sizeMax = MAX(size, *sizeMax);
        displayBigText((int)(_view.w/2-size*0.5*9), *l, loc);
        *l = *l + 20;
      }
      *l = *l+10;
      break;

    case TEXT_KEYBOARD:
      for (int i = 0; i < keybTextSize; i++)
      {
        sprintf(loc, "%s", keybText[i]);
        size = strlen(loc);
        *sizeMax = MAX(size, *sizeMax);
        displayBigText((int)(_view.w/2-size*0.5*9), *l, loc);
        *l = *l + 20;
      }
      *l = *l+10;
      break;

    default:
      // one line items
      size = strlen(loc);
      *sizeMax = MAX(size, *sizeMax);
      displayBigText((int)(_view.w/2-size*0.5*9), *l, loc);
      *l = *l + 30;
      break;
  }
}

//=============================================================================
/*
  Change the value of the selected item.
*/
//=============================================================================
void Data::changeMenuItemRight()
{
  switch (order[ptrState->currentMenuItem])
  {
    case 0: // Keyboard settings submenu is 4
      currentMenu = 4;
      start = menuSize[currentMenu-1];
      end = menuSize[currentMenu];
      ptrState->currentMenuItem = start;
      menu();
      break;

    case 1: // Increase/decrease drawing accuracy
      _pref.speed++;
      if (_pref.speed > 5) _pref.speed = 5;
      menu();
      break;

    case 2: // Current blanking function
      if (_pref.blanking == NULL)
      {
        _pref.blanking = _plugins.blanking;
        if (_pref.blanking != NULL)
          findBlankedZones();
      }
      else
      {
        _pref.blanking = _pref.blanking->next;
        findBlankedZones();
      }
      menu();
      break;

    case 3: // Current look-for function
      if (_pref.lookFor == NULL)
        _pref.lookFor = _plugins.lookFor;
      else
        _pref.lookFor = _pref.lookFor->next;
      if (_pref.lookFor == NULL) _pref.lookFor = _plugins.lookFor;
      menu();
      break;

    case 4: // Current add-a-variable function
      if (_pref.addAVariable == NULL)
        _pref.addAVariable = _plugins.addAVariable;
      else
        _pref.addAVariable = _pref.addAVariable->next;
      menu();
      break;

    case 5: // Display credits

      break;

    case 6: // Exit menu
      ptrState->inmenu = 0;
      glutKeyboardFunc(gkeyboard);
      glutSpecialFunc(garrows);
      glutTimerFunc(ptrState->timeStep, gtimer, 0);
      if (ptrState->smoke == 1)
        glutIdleFunc(gdisplay);
      display();
      break;

    case 9: // Back from sub menu
      switch (currentMenu)
      {
        case 1:
          ptrState->currentMenuItem = 1;
          break;
        case 2:
          ptrState->currentMenuItem = 2;
          break;
        case 3:
          ptrState->currentMenuItem = 0;
          break;
        case 4:
          ptrState->currentMenuItem = 3;
          break;
        default:
          ptrState->currentMenuItem = 0;
      }
      currentMenu = 0;
      start = 0;
      end = menuSize[currentMenu];
      menu();
      break;

    case 10: // Preference sub-menu is 1
      currentMenu = 1;
      start = menuSize[currentMenu-1];
      end = menuSize[currentMenu];
      ptrState->currentMenuItem = start;
      menu();
      break;

    case 11: // Plugins sub-menu is 2
      currentMenu = 2;
      start = menuSize[currentMenu-1];
      end = menuSize[currentMenu];
      ptrState->currentMenuItem = start;
      menu();
      break;

    case 12: // Data set info sub-menu is 3
      currentMenu = 3;
      start = menuSize[currentMenu-1];
      end = menuSize[currentMenu];
      ptrState->currentMenuItem = start;
      menu();
      break;

    case 13: // Max number of particles
      _pref.maxParticles = _pref.maxParticles + 100;
      if (ptrState->smoke == 1)
      {
        //free(_cloud.particle);
        //_cloud.n = _pref.maxParticles;
        //_cloud.particle = (Particle*)malloc(_cloud.n * sizeof(Particle));
        //initParticles();
      }
      menu();
      break;

    case 14: // Smoke particle radius
      _pref.smokeRadius = 1.1*_pref.smokeRadius;
      menu();
      break;

    case 15: // Smoke time step
      _pref.smokeTimeStep = 1.1*_pref.smokeTimeStep;
      menu();
      break;

    case 16: // Smoke emission zone radius
      _pref.emissionRadius = 1.1*_pref.emissionRadius;
      menu();
      break;

    case 17: // Select change
      if (_pref.selectNextZone == NULL)
      {
        _pref.selectNextZone = _plugins.selectNextZone;
        _pref.selectPreviousZone = _plugins.selectPreviousZone;
      }
      else
      {
        _pref.selectNextZone = _pref.selectNextZone->next;
        _pref.selectPreviousZone = _pref.selectPreviousZone->next;
        if (_pref.selectNextZone == NULL)
        {
          _pref.selectNextZone = _plugins.selectNextZone;
          _pref.selectPreviousZone = _plugins.selectPreviousZone;
        }
      }
      menu();
      break;

    case 18: // Colormap change
      if (_pref.colorMap == NULL)
        _pref.colorMap = _plugins.colorMap;
      else
        _pref.colorMap = _pref.colorMap->next;
      if (_pref.colorMap == NULL) _pref.colorMap = _plugins.colorMap;
      menu();
      break;

    case 19: // Screen dump change
      if (_pref.screenDump == NULL)
        _pref.screenDump = _plugins.screenDump;
      else
        _pref.screenDump = _pref.screenDump->next;
      if (_pref.screenDump == NULL) _pref.screenDump = _plugins.screenDump;
      menu();
      break;

    case 20: // Mouse click change
      if (_pref.mouseClick == NULL)
        _pref.mouseClick = _plugins.mouseClick;
      else
      {
        _pref.mouseClick = _pref.mouseClick->next;
        if (_pref.mouseClick == NULL)
          _pref.mouseClick = _plugins.mouseClick;
      }
      menu();
      break;

    case 21: // Texture size
      _pref.sizeTexi = _pref.sizeTexi + 1;
      menu();
      break;

    case 24: // Texture size
      _pref.sizeTexj = _pref.sizeTexj + 1;
      menu();
      break;
  }
}

//=============================================================================
/*
  Change the value of the selected item.
*/
//=============================================================================
void Data::changeMenuItemLeft()
{
  switch (order[ptrState->currentMenuItem])
  {
    case 0: // Keyboard settings
      currentMenu = 4;
      start = menuSize[currentMenu-1];
      end = menuSize[currentMenu];
      ptrState->currentMenuItem = start;
      menu();
      break;

    case 1: // Increase/decrease drawing accuracy
      _pref.speed--;
      if (_pref.speed < 0) _pref.speed = 0;
      menu();
      break;

    case 2: // Change the current blanking function
      if (_pref.blanking == NULL)
      {
        _pref.blanking = _plugins.blanking;
        if (_pref.blanking != NULL)
          findBlankedZones();
      }
      else
      {
        _pref.blanking = _pref.blanking->next;
        findBlankedZones();
      }
      menu();
      break;

    case 3: // Current look-for function
      if (_pref.lookFor == NULL)
        _pref.lookFor = _plugins.lookFor;
      else
        _pref.lookFor = _pref.lookFor->next;
      if (_pref.lookFor == NULL) _pref.lookFor = _plugins.lookFor;
      menu();
      break;

    case 4: // Current add-a-variable function
      if (_pref.addAVariable == NULL)
        _pref.addAVariable = _plugins.addAVariable;
      else
        _pref.addAVariable = _pref.addAVariable->next;
      menu();
      break;

    case 5: // Display the cast

      break;

    case 6: // Exit menu
      ptrState->inmenu = 0;
      glutKeyboardFunc(gkeyboard);
      glutSpecialFunc(garrows);
      glutTimerFunc(ptrState->timeStep, gtimer, 0);
      if (ptrState->smoke == 1)
        glutIdleFunc(gdisplay);
      display();
      break;

    case 9: // Back from sub menu
      switch (currentMenu)
      {
        case 1:
          ptrState->currentMenuItem = 1;
          break;
        case 2:
          ptrState->currentMenuItem = 2;
          break;
        case 3:
          ptrState->currentMenuItem = 0;
          break;
        case 4:
          ptrState->currentMenuItem = 3;
          break;
        default:
          ptrState->currentMenuItem = 0;
      }
      currentMenu = 0;
      start = 0;
      end = menuSize[currentMenu];
      menu();
      break;

    case 10: // Preference sub-menu is 1
      currentMenu = 1;
      start = menuSize[currentMenu-1];
      end = menuSize[currentMenu];
      ptrState->currentMenuItem = start;
      menu();
      break;

    case 11: // Plugins sub-menu is 2
      currentMenu = 2;
      start = menuSize[currentMenu-1];
      end = menuSize[currentMenu];
      ptrState->currentMenuItem = start;
      menu();
      break;

    case 12: // Data set info sub-menu is 3
      currentMenu = 3;
      start = menuSize[currentMenu-1];
      end = menuSize[currentMenu];
      ptrState->currentMenuItem = start;
      menu();
      break;

    case 13: // Max number of particles
      _pref.maxParticles = _pref.maxParticles - 100;
      if (_pref.maxParticles < 100) _pref.maxParticles = 100;
      if (ptrState->smoke == 1)
      {
        //free(_cloud.particle);
        //_cloud.n = _pref.maxParticles;
        //_cloud.particle = (Particle*)malloc(_cloud.n * sizeof(Particle));
        //initParticles();
      }
      menu();
      break;

    case 14: // Smoke particles radius
      _pref.smokeRadius = 0.9*_pref.smokeRadius;
      menu();
      break;

    case 15: // Smoke time step
      _pref.smokeTimeStep = 0.9*_pref.smokeTimeStep;
      menu();
      break;

    case 16: // Smoke emission zone radius
      _pref.emissionRadius = 0.9*_pref.emissionRadius;
      menu();
      break;

    case 17: // Select change
      if (_pref.selectNextZone == NULL)
      {
        _pref.selectNextZone = _plugins.selectNextZone;
        _pref.selectPreviousZone = _plugins.selectPreviousZone;
      }
      else
      {
        _pref.selectNextZone = _pref.selectNextZone->next;
        _pref.selectPreviousZone = _pref.selectPreviousZone->next;
        if (_pref.selectNextZone == NULL)
        {
          _pref.selectNextZone = _plugins.selectNextZone;
          _pref.selectPreviousZone = _plugins.selectPreviousZone;
        }
      }
      menu();
      break;

    case 18: // colorMap change
      if (_pref.colorMap == NULL)
        _pref.colorMap = _plugins.colorMap;
      else
        _pref.colorMap = _pref.colorMap->next;
      if (_pref.colorMap == NULL) _pref.colorMap = _plugins.colorMap;
      menu();
      break;

    case 19: // Screen dump change
      if (_pref.screenDump == NULL)
        _pref.screenDump = _plugins.screenDump;
      else
        _pref.screenDump = _pref.screenDump->next;
      if (_pref.screenDump == NULL) _pref.screenDump = _plugins.screenDump;
      menu();
      break;

    case 20: // Mouse click change
      if (_pref.mouseClick == NULL)
        _pref.mouseClick = _plugins.mouseClick;
      else
      {
        _pref.mouseClick = _pref.mouseClick->next;
        if (_pref.mouseClick == NULL)
          _pref.mouseClick = _plugins.mouseClick;
      }
      menu();
      break;

    case 21: // Texture size
      _pref.sizeTexi = _pref.sizeTexi - 1;
      //if (_pref.sizeTexi == 0) _pref.sizeTexi = 1;
      menu();
      break;

    case 24: // Texture size
      _pref.sizeTexj = _pref.sizeTexj - 1;
      //if (_pref.sizeTexj == 0) _pref.sizeTexj = 1;
      menu();
      break;
  }

}

//=============================================================================
/*
  Change the value of the selected item.
*/
//=============================================================================
void Data::changeMenuItemSpaceBar()
{
  int stateHeader, stateInfo, stateMenu, stateBB;

  switch (order[ptrState->currentMenuItem])
  {
    case 0: // Keyboard settings
      currentMenu = 4;
      start = menuSize[currentMenu-1];
      end = menuSize[currentMenu];
      ptrState->currentMenuItem = start;
      menu();
      break;

    case 1: // Increase
      _pref.speed++;
      if (_pref.speed > 5) _pref.speed = 5;
      menu();
      break;

    case 2: // Current blanking function

      break;

    case 3: // Current look-for function

      break;

    case 4: // Triggers the real creation of variable
      if (_pref.addAVariable != NULL)
        _pref.addAVariable->f(this);
      break;

    case 5: // Display credits

      break;

    case 6: // Exit cplot
      exit(0);
      break;

    case 9: // Back from sub menu
      switch (currentMenu)
      {
        case 1:
          ptrState->currentMenuItem = 1;
          break;
        case 2:
          ptrState->currentMenuItem = 2;
          break;
        case 3:
          ptrState->currentMenuItem = 0;
          break;
        case 4:
          ptrState->currentMenuItem = 3;
          break;
        default:
          ptrState->currentMenuItem = 0;
      }
      currentMenu = 0;
      start = 0;
      end = menuSize[currentMenu];
      menu();
      break;

    case 10: // Preference sub-menu is 1
      currentMenu = 1;
      start = menuSize[currentMenu-1];
      end = menuSize[currentMenu];
      ptrState->currentMenuItem = start;
      menu();
      break;

    case 11: // Plugins sub-menu is 2
      currentMenu = 2;
      start = menuSize[currentMenu-1];
      end = menuSize[currentMenu];
      ptrState->currentMenuItem = start;
      menu();
      break;

    case 12: // Data set info sub-menu is 3
      currentMenu = 3;
      start = menuSize[currentMenu-1];
      end = menuSize[currentMenu];
      ptrState->currentMenuItem = start;
      menu();
      break;

    case 18: // colorMap
      if (ptrState->mode == MESH || ptrState->mode == SOLID)
        if (_zones[0]->nfield >= 1) ptrState->mode = 3;
      menu();
      break;

    case 19: // Screen dump
      if (_pref.screenDump != NULL)
      {
        stateHeader = ptrState->header;
        stateInfo = ptrState->info;
        stateMenu = ptrState->menu;
        stateBB = ptrState->bb;
        ptrState->header = 0;
        ptrState->info = 0;
        ptrState->menu = 0;
        ptrState->bb = 0;
        display();
        dumpWindow();
        ptrState->header = stateHeader;
        ptrState->info = stateInfo;
        ptrState->menu = stateMenu;
        ptrState->bb = stateBB;
        //printTmpMessage("Image dumped to file");
      }
      menu();
      break;
  }
}

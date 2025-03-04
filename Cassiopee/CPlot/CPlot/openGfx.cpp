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

#include "Data.h"

//=============================================================================
void Data::setBgColor()
{
  // Background color
  switch (ptrState->bgColor)
  {
    case 0: // noir
      glClearColor(0.0, 0.0, 0.0, 0.0); break;
 
    case 1: // blanc
      glClearColor(1.0, 1.0, 1.0, 0.0); break;

    case 3: // jaune
      glClearColor(1., 1., 0.75, 0.5); break;

    default: // blanc
      glClearColor(1.0, 1.0, 1.0, 0.0); break; 
  }
}

//=============================================================================
/* Init the GL */
//=============================================================================
void Data::init(void)
{ 
#ifdef __SHADERS__
  glewInit();
#endif
  
  // Cree les textures
  E_Int width, height;
  createNodeTexture();
  createNoise3DTexture();
  createFrameBufferTexture();
  createPngTexture("windtunnel.png", _texEnviron1, width, height);
  createVoxelTexture();

  // Set bg color
  setBgColor();

  // Shade model
  glShadeModel(GL_SMOOTH);
  
  // Enable the depth tests
  glEnable(GL_DEPTH_TEST);

  // Init des shaders
#ifdef __SHADERS__
  _shaders.init();
  _shaders.load();
#endif
}

//=============================================================================
void Data::openGfx()
{
  /* Init */
  glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH | GLUT_MULTISAMPLE);

  // Window size base sur la taille de l'ecran
  int screenWidth = glutGet(GLUT_SCREEN_WIDTH); 
  int screenHeight = glutGet(GLUT_SCREEN_HEIGHT);
  if (_view.w == 0 || _view.h == 0)
  {
    _view.w = screenWidth-320;
    _view.h = screenHeight;
  }
  
  /* Window */
  if (ptrState->offscreen >= 2) glutInitWindowSize(1, 1);
  else glutInitWindowSize(_view.w, _view.h);
  glutInitWindowPosition(320, 0);
  _winId = glutCreateWindow(ptrState->winTitle);
  if (ptrState->offscreen >= 2) glutHideWindow();
  init();

  glutReshapeFunc(reshape);
  glutKeyboardFunc(gkeyboard);
  //glutKeyboardUpFunc(gkeyboardup); // other keys than arrows
  if (ptrState->kkeysActivated == 0) glutSpecialUpFunc(gkeyboardup); // arrows
  glutDisplayFunc(fdisplay);
  glutSpecialFunc(garrows);
  glutMouseFunc(gmouseButton);
  glutMotionFunc(gmouseMotion);
  //glutPassiveMotionFunc(gmousePassiveMotion);

#ifdef __APPLE__
  glutIdleFunc(gidle);
#else
  glutTimerFunc(ptrState->timeStep, gtimer, 0);
#endif

  if (_view.w > 5000 || _view.h > 5000)
  {
    glutFullScreen();
    ptrState->fullScreen = 1;
  }
}

//=============================================================================
void Data::closeGfx()
{
  //glutDestroyWindow(_view.winId);
}

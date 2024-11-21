/*    
    Copyright 2013-2024 Onera.

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
// Nettoie le display.
// Affiche le fond d'ecran
//=============================================================================
void Data::clearDisplay()
{
  setBgColor();
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  switch (ptrState->bgColor)
  {
    case 2: // grey gradient
      glMatrixMode(GL_MODELVIEW); glPushMatrix(); glLoadIdentity(); 
      glMatrixMode(GL_PROJECTION); glPushMatrix(); glLoadIdentity();
      glDisable(GL_DEPTH_TEST);
      glBegin(GL_QUADS);
      glColor3f(255./255., 255./255., 255./255.);
      glVertex3i(-1, -1, -1); 
      glColor3f(233./255., 234./255., 238./255.);
      glVertex3i (1, -1, -1); 
      glColor3f(233./255., 234./255., 238./255.);
      glVertex3i(1, 1, -1); 
      glColor3f(255./255., 255./255., 255./255.);
      glVertex3i(-1, 1, -1); 
      glEnd();
      glEnable(GL_DEPTH_TEST);
      glPopMatrix(); glMatrixMode(GL_MODELVIEW); glPopMatrix();
      break;

    case 4: // bleu gradient
      glMatrixMode(GL_MODELVIEW); glPushMatrix(); glLoadIdentity(); 
      glMatrixMode(GL_PROJECTION); glPushMatrix(); glLoadIdentity();
      glDisable(GL_DEPTH_TEST);
      glBegin(GL_QUADS);
      glColor3f(14./255., 175./255., 250./255.);
      glVertex3i(-1, -1, -1); 
      glVertex3i(1, -1, -1); 
      glColor3f(93./255., 207./255., 236./255.);
      glVertex3i(1, 1, -1); 
      glVertex3i(-1, 1, -1); 
      glEnd();
      glEnable(GL_DEPTH_TEST);
      glPopMatrix(); glMatrixMode(GL_MODELVIEW); glPopMatrix();
      break;

    case 5: // rouge-orange gradient
      glMatrixMode(GL_MODELVIEW); glPushMatrix(); glLoadIdentity(); 
      glMatrixMode(GL_PROJECTION); glPushMatrix(); glLoadIdentity();
      glDisable(GL_DEPTH_TEST);
      glBegin (GL_QUADS);
      glColor3f(243./255., 130./255., 10./255.);
      glVertex3i(-1, -1, -1); 
      glColor3f(114./255., 6./255., 4./255.);
      glVertex3i(1, -1, -1); 
      glColor3f(67./255., 41./255., 42./255.);
      glVertex3i(1, 1, -1); 
      glColor3f(127./255., 109./255., 59./255.);
      glVertex3i(-1, 1, -1); 
      glEnd();
      glEnable(GL_DEPTH_TEST);
      glPopMatrix(); glMatrixMode(GL_MODELVIEW); glPopMatrix();
      break;

    case 6: // textured backgrounds
    case 7:
    case 8:
    case 9:
    case 10:
    case 11:
    case 12:
    case 13: 
      glMatrixMode(GL_TEXTURE); glLoadIdentity(); glMatrixMode(GL_MODELVIEW);
      setOrthographicProjection();
      glPushMatrix(); glLoadIdentity();

      glActiveTexture(GL_TEXTURE0);
      glEnable(GL_TEXTURE_2D);
      glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
      glBindTexture(GL_TEXTURE_2D, _texBackground);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
      glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

      glDisable(GL_DEPTH_TEST);
      glBegin(GL_QUADS);
      glTexCoord2f(0.0, 1.0); glVertex3d(0., 0., 0);
      glTexCoord2f(1.0, 1.0); glVertex3d(_view.w, 0, 0); 
      glTexCoord2f(1.0, 0.0); glVertex3d(_view.w, _view.h, 0); 
      glTexCoord2f(0.0, 0.0); glVertex3d(0, _view.h, 0);
      glEnd(); 
      glEnable(GL_DEPTH_TEST);
      glDisable(GL_TEXTURE_2D);
      
      glPopMatrix();
      resetPerspectiveProjection();
      break;

    default: break;
  }
}

//=============================================================================
/* 
   Fonction principale d'affichage.
   Son comportement est determine par le state.
   - state.alpha 
   Parametre d'alpha blending. Si alpha = 0, render sans alpha blending.
   - state.render 
   Si 1, flush le render a la fin. 0 pas de flush.
   - state.bb
   Si 1, render de la bounding box
   - state.header 
   Si 1, render du header
   - state.info
   Si 1, render de la fenetre info selected zone

   Effectue le render de :
   - bounding box
   - un mode
   - header
   - Menus secondaires
   - Menus principaux
*/
//=============================================================================
void fdisplay()
{
  Data* d = Data::getInstance();
  d->ptrState->render = 1;
  gdisplay();
}
void gdisplay()
{
  Data* d = Data::getInstance();
  if (d->ptrState->farClip == 1) {d->farClipping(); d->ptrState->farClip = 0; }

  d->ptrState->lockDisplay();
  // Free GPU ressources command
  if (d->ptrState->freeGPURes == 1)
  {
    d->ptrState->lockGPURes();
    d->ptrState->freeGPUResources();
    d->ptrState->unlockGPURes();
  }

  if (d->ptrState->render >= 0)
  {
    // Check si display iso solid
    if (d->ptrState->mode == SCALARFIELD) d->createIsoGPURes(d->ptrState->scalarField);
    if (d->ptrState->mode == VECTORFIELD) d->createIsoGPURes(d->ptrState->vectorField1,
                                                             d->ptrState->vectorField2,
                                                             d->ptrState->vectorField3);
    if (d->ptrState->mode == RENDER) d->createIsoGPUResForRender();

    // Update mouse cursor eventuellement
    if (d->ptrState->updateCursor != -1)
    {
      d->setCursor(d->ptrState->updateCursor);
      d->ptrState->updateCursor = -1;
    }

    // Update des textures eventuellement
    if (d->ptrState->updateEnvmap == 1)
    { 
      E_Int width, height;
      d->createImageTexture(d->ptrState->envmapFile, d->_texEnviron1, width, height, true);
      d->ptrState->updateEnvmap = 0; 
    }
    if (d->ptrState->updateBackground == 1)
    { 
      E_Int width, height;
      if (d->ptrState->bgColor == 6)
       d->createPngTexture("paperBackground1.png", d->_texBackground, width, height, false);
      else if (d->ptrState->bgColor == 7)
       d->createPngTexture("paperBackground2.png", d->_texBackground, width, height, false);
      else if (d->ptrState->bgColor == 8)
       d->createPngTexture("paperBackground3.png", d->_texBackground, width, height, false);
      else if (d->ptrState->bgColor == 9)
       d->createPngTexture("paperBackground4.png", d->_texBackground, width, height, false);
      else if (d->ptrState->bgColor == 10)
       d->createPngTexture("paperBackground5.png", d->_texBackground, width, height, false);
      else if (d->ptrState->bgColor == 11)
       d->createPngTexture("paperBackground6.png", d->_texBackground, width, height, false);
      else if (d->ptrState->bgColor == 12)
       d->createPngTexture("paperBackground7.png", d->_texBackground, width, height, false);
      else if (d->ptrState->bgColor == 13)
       d->createImageTexture(d->ptrState->backgroundFile, d->_texBackground, width, height, false);
      else
       d->createPngTexture("paperBackground1.png", d->_texBackground, width, height, false);
      d->ptrState->updateBackground = 0;
    }

    // Display
    if (d->ptrState->offscreen == 0)
    {
      if (d->ptrState->stereo == 0) d->display();
      else d->displayAnaglyph();
    }
    d->ptrState->render = -1;
    // time
    if (d->ptrState->ktimer > 0)
    { d->ptrState->ktime++; 
      d->ptrState->ktimer = (d->ptrState->ktimer>0 ? d->ptrState->ktimer-1 : 0); 
      d->ptrState->render = 1; }
    // Screen shoot eventuel
    if (d->ptrState->shootScreen == 1) d->exportFile();
  }
  else
  {
    // update DL quand on a le temps
    //d->ptrState->lockDisplay();
    d->createGPURes();
    //d->ptrState->unlockDisplay();
  }
  d->ptrState->unlockDisplay();
}

//=============================================================================
void Data::display()
{
  ptrState->lockDisplay();

  // Preprocessing pour le shadow mapping
  if (ptrState->shadow == 1)
  {
#include "shadow.h"
  } // Fin preprocessing

  // Clear the display
  clearDisplay();

  // Set drawing color
  glColor4f(1.0, 1.0, 1.0, ptrState->alpha); // draw color

  // Define camera position
  glLoadIdentity();
  gluLookAt(_view.xcam, _view.ycam, _view.zcam, 
            _view.xeye, _view.yeye, _view.zeye, 
            _view.dirx, _view.diry, _view.dirz);
  computeFrustumPlanes(_view);
  
  glEnable(GL_MULTISAMPLE);
  
  // Display following mode
  // 2D or 3D rendering
  switch (ptrState->mode)
  {
    case MESH:
      displaySMesh(); displayUMesh(); displaySEdges(); displayUEdges();
      break;
      
    case SOLID:
      displaySSolid(); displayUSolid(); displaySEdges(); displayUEdges(); displayNodes();
      break;
      
    case RENDER:
      displaySSolid(); displayUSolid(); displaySEdges(); displayUEdges(); displayNodes();
      displayAllBillBoards(); // must be the last
      break;
      
    case SCALARFIELD:
      displaySIsoSolid(); displayUIsoSolid(); displaySEdges(); displayUEdges(); displayNodes();
      break;
      
    case VECTORFIELD:
      displaySIsoSolid(); displayUIsoSolid(); displaySEdges(); displayUEdges(); displayNodes();
      break;
      
    default: // SOLID
      displaySSolid(); displayUSolid(); displaySEdges(); displayUEdges(); displayNodes();
      break;
  }

  // 1D overlay display
  displayPlots();

#ifdef __SHADERS__
  // Update du frame buffer pour les shaders le necessitant (glass par ex)
  int update = 0;
  for (E_Int i = 0; i < _numberOfZones; i++)
  {
    if (_zones[i]->material == 1) { update = 1; break; }
  }
  if (update == 1)
  {
    // Il faut normalement faire la photo sans les objets glass
    E_Int w = _view.w; E_Int h = _view.h;
    w = w < _frameBufferSize[ptrState->frameBuffer] ? w : (int) _frameBufferSize[ptrState->frameBuffer];
    h = h < _frameBufferSize[ptrState->frameBuffer] ? h : (int) _frameBufferSize[ptrState->frameBuffer];
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, _texFrameBuffer[ptrState->frameBuffer]);
    //glCopyTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, 0, 0, w, h);
    glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 0, 0, w, h, 0);
  }
#endif

  // Post-processing
  E_Int post = 0; double sobelThreshold = -0.5; 
  if (ptrState->DOF == 1) { post = 1; sobelThreshold = ptrState->sobelThreshold; }
  if (ptrState->gamma != 1.) { post = 1; sobelThreshold = ptrState->sobelThreshold; }
  if (ptrState->toneMapping != 0) { post = 1; sobelThreshold = ptrState->sobelThreshold; }
  if (ptrState->mode == SOLID && ptrState->solidStyle == 4) { post = 1; sobelThreshold = 0.5; }
  if (ptrState->panorama == 1) post =1;

  if (post == 1)
  { 
    // Recupere l'image standard dans _texRight
    if (_texRight != 0) glDeleteTextures(1, &_texRight);
    glActiveTexture(GL_TEXTURE0);
    glGenTextures(1, &_texRight);
    glBindTexture(GL_TEXTURE_2D, _texRight);
    glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 0, 0, _view.w, _view.h, 0);
    // Recupere le depth buffer dans _texLeft
    if (_texLeft != 0) glDeleteTextures(1, &_texLeft);
    glActiveTexture(GL_TEXTURE1);
    glGenTextures(1, &_texLeft);
    glBindTexture(GL_TEXTURE_2D, _texLeft);
    glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_DEPTH_COMPONENT, 0, 0, _view.w, _view.h, 0);
    // Z buffer suivi du point actif
    GLint viewport[4];
    GLdouble modelview[16];
    GLdouble projection[16];
    glGetDoublev(GL_MODELVIEW_MATRIX, modelview);
    glGetDoublev(GL_PROJECTION_MATRIX, projection);
    glGetIntegerv(GL_VIEWPORT, viewport);
    GLdouble winX, winY, winZ;
    gluProject(ptrState->activePointX, ptrState->activePointY, ptrState->activePointZ, 
               modelview, projection, viewport, 
               &winX, &winY, &winZ);
    ptrState->activePointZBuf = winZ;

    // Render frame tex
    glClear(GL_COLOR_BUFFER_BIT);
    setOrthographicProjection();
    glPushMatrix();
    glLoadIdentity();
    if (ptrState->panorama == 0)
      displayFrameTex(2, sobelThreshold); // DOF+sobel
    else 
      displayFrameTex(1, sobelThreshold); // panorama
    glPopMatrix();
    resetPerspectiveProjection();
  }

  // Info + legende
  if (ptrState->dim != 1 &&
      ptrState->offscreen != 1 && ptrState->offscreen != 5 &&
      ptrState->offscreen != 6 && ptrState->offscreen != 7)
  {
    // All those functions use glut font and so X
    if (ptrState->header == 1) printHeader();
    if (ptrState->isoLegend > 0) displayIsoLegend(ptrState->isoLegend);
    if (ptrState->info == 1) { displayInfo(); displayAxis(); }
    if (ptrState->info == 1) displayActivePoint();
  }

  // Message overlay
  if (ptrState->message != NULL) printTmpMessage(ptrState->message);

  ptrState->unlockDisplay();

  // Rendering
  if (ptrState->render == 1) glutSwapBuffers();
  else glFlush(); // Force flush (usefull with osmesa)
  
  glDisable(GL_MULTISAMPLE);
}

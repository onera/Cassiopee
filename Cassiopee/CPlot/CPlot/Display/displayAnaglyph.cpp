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

//=============================================================================
void Data::displayAnaglyph()
{
  ptrState->lockDisplay();

  // Preprocessing pour le shadow mapping
  if (ptrState->shadow == 1)
  {
#include "shadow.h"
  } // Fin preprocessing

  //=========
  // Left eye
  //=========

  // Clear the display
  clearDisplay();

  // Set drawing color
  glColor4f(1.0, 1.0, 1.0, ptrState->alpha); // draw color

  // Define camera position
  glLoadIdentity();
  double rx, ry, rz;
  rx = _view.xcam - _view.xeye;
  ry = _view.ycam - _view.yeye;
  rz = _view.zcam - _view.zeye;
  double dx, dy, dz;
  dx = ry*_view.dirz - rz*_view.diry;
  dy = rz*_view.dirx - rx*_view.dirz;
  dz = rx*_view.diry - ry*_view.dirx;
  double dd = dx*dx + dy*dy + dz*dz;

  dd = 20.*epsup*ptrState->stereoDist/sqrt(dd);

  // variation
  E_Float distx = xmax-xmin;
  E_Float disty = ymax-ymin;
  E_Float distz = zmax-zmin;
  E_Float dist = MAX(distx, disty);
  dist = MAX(dist, distz);
  dd = 3*dd*sqrt(rx*rx+ry*ry+rz*rz)/dist;
  // fin variation

  dx = dd*dx; dy = dd*dy; dz = dd*dz;
  //printf("offset: %f %f %f\n", dx, dy, dz);
  gluLookAt(_view.xcam+dx, _view.ycam+dy, _view.zcam+dz, 
            _view.xeye, _view.yeye, _view.zeye, 
            _view.dirx, _view.diry, _view.dirz);
  computeFrustumPlanes(_view);
  
  glEnable(GL_MULTISAMPLE);

  // Display following mode
  if (ptrState->dim != 1)
  {
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
  }
  else
  {
//     // 1D rendering
//     displayPlot();
  }

  // Recuperation dans la texture
  if (_texLeft != 0) glDeleteTextures(1, &_texLeft);
  glActiveTexture(GL_TEXTURE0);
  glGenTextures(1, &_texLeft);
  glBindTexture(GL_TEXTURE_2D, _texLeft);
  glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 0, 0, _view.w, _view.h, 0);
  
  //==========
  // Right eye
  //==========
  
  // Clear the display
  clearDisplay();
  glColor4f(1.0, 1.0, 1.0, ptrState->alpha); // draw color

  // Set drawing color
  glColor4f(1.0, 1.0, 1.0, ptrState->alpha); // draw color

  // Define camera position
  glLoadIdentity();
  gluLookAt(_view.xcam, _view.ycam, _view.zcam, 
            _view.xeye, _view.yeye, _view.zeye, 
            _view.dirx, _view.diry, _view.dirz);
  computeFrustumPlanes(_view);
  
  // Display following mode
  if (ptrState->dim != 1)
  {
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
  }
  else
  {
//     // 1D rendering
//     displayPlot();
  }

#ifdef __SHADERS__
  // Update du frame buffer pour les shaders le necessitant
  // Je pense qu'il faudrait le faire tout le temps
  if (_shaders.currentShader() == _shaders.shader_id(3))
  {
     E_Int w = _view.w; E_Int h = _view.h;
     w = w < _frameBufferSize[ptrState->frameBuffer] ? w : (int) _frameBufferSize[ptrState->frameBuffer];
     h = h < _frameBufferSize[ptrState->frameBuffer] ? h : (int) _frameBufferSize[ptrState->frameBuffer];
    glBindTexture(GL_TEXTURE_2D, _texFrameBuffer[ptrState->frameBuffer]);
    //glCopyTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, 0, 0, w, h);
    glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 0, 0, w, h, 0);
  }
#endif
  
  // Recuperation dans la texture
  if (_texRight != 0) glDeleteTextures(1, &_texRight);
  glActiveTexture(GL_TEXTURE1);
  glGenTextures(1, &_texRight);
  glBindTexture(GL_TEXTURE_2D, _texRight);
  glCopyTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 0, 0, _view.w, _view.h, 0);

  // Render frame tex
  glClear(GL_COLOR_BUFFER_BIT);
  setOrthographicProjection();
  glPushMatrix();
  glLoadIdentity();
  displayFrameTex(0);
  glPopMatrix();
  resetPerspectiveProjection();

  // Info + legende
  if (ptrState->dim != 1 &&
      ptrState->offscreen != 1 && ptrState->offscreen != 5 &&
      ptrState->offscreen != 6 && ptrState->offscreen != 7)
  {
    if (ptrState->header == 1) printHeader();
    if (ptrState->info == 1) 
    { displayInfo(); displayActivePoint(); displayAxis(); }
    if (ptrState->isoLegend > 0) displayIsoLegend(ptrState->isoLegend);
  }

  ptrState->unlockDisplay();
  
  if (ptrState->render == 1) glutSwapBuffers();
  else glFlush(); // Force flush (usefull with osmesa)

  glDisable(GL_MULTISAMPLE);
}

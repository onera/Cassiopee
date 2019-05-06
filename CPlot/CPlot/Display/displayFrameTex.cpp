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
  Display une texture de la taille de l'ecran dans l'ecran.
  mode=0: anaglyph
  mode=2: DOF+GAMMA+SOBEL
*/
//============================================================================
void Data::displayFrameTex(int mode, double sobelThreshold)
{
  glColor3f(1., 0, 0);

#ifdef __SHADERS__
  if (mode == 0) // anaglyph
  {
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, _texLeft);
    glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, _texRight);
    glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

    int shader = _shaders.shader_id(12 + ptrState->stereo);
    if (_shaders.currentShader() != shader) _shaders.activate(shader);
    _shaders[shader]->setUniform("leftEyeTexture", (int)0);
    _shaders[shader]->setUniform("rightEyeTexture", (int)1);
  }
  else if (mode == 2) // DOF
  { 
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, _texRight);
    glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, _texLeft);
    glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
 
    int shader = _shaders.shader_id(20);
    if (_shaders.currentShader() != shader) _shaders.activate(shader);
    _shaders[shader]->setUniform("FrameBuffer", (int)0);
    _shaders[shader]->setUniform("depthMap", (int)1);
    _shaders[shader]->setUniform("focalDepth", (float)ptrState->activePointZBuf);
    _shaders[shader]->setUniform("radius", (float)ptrState->dofPower);
    _shaders[shader]->setUniform("ext", (float)1.);
    _shaders[shader]->setUniform("gamma", (float)ptrState->gamma);
    _shaders[shader]->setUniform("sobelThreshold", (float)sobelThreshold); 
  }
#endif

  glDisable(GL_DEPTH_TEST);
  glBegin(GL_QUADS);
  glTexCoord3f(0, 1, 0);
  glVertex3d(0, 0, 0);
  glTexCoord3f(1, 1, 0);
  glVertex3d(_view.w, 0, 0);
  glTexCoord3f(1, 0, 0);
  glVertex3d(_view.w, _view.h, 0);
  glTexCoord3f(0, 0, 0);
  glVertex3d(0, _view.h, 0);
  glEnd();
  glEnable(GL_DEPTH_TEST);
  glColor3f(1., 1., 1.);
#ifdef __SHADERS__
  _shaders.activate((short unsigned int)0);
#endif

}

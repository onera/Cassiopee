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

//============================================================================
/*
  Display une texture de la taille de l'ecran dans l'ecran.
  mode=0: anaglyph
  mode=2: DOF+GAMMA+SOBEL
*/
//============================================================================
void Data::displayFrameTex(E_Int mode, double sobelThreshold)
{
  glColor3f(1., 0, 0);
  GLuint tex0, tex1, tex2, tex3, tex4, tex5;
    
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

    E_Int shader = _shaders.shader_id(12 + ptrState->stereo);
    if (_shaders.currentShader() != shader) _shaders.activate(shader);
    _shaders[shader]->setUniform("leftEyeTexture", (int)0);
    _shaders[shader]->setUniform("rightEyeTexture", (int)1);
  }
  else if (mode == 1) // panorama
  {
    E_Int w, h;
    glActiveTexture(GL_TEXTURE0);
    createImageTexture("cube_left.png", tex0, w, h, false);
    glBindTexture(GL_TEXTURE_2D, tex0);
    glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glActiveTexture(GL_TEXTURE1);
    createImageTexture("cube_right.png", tex1, w, h, false);
    glBindTexture(GL_TEXTURE_2D, tex1);
    glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glActiveTexture(GL_TEXTURE2);
    createImageTexture("cube_bottom.png", tex2, w, h, false);
    glBindTexture(GL_TEXTURE_2D, tex2);
    glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glActiveTexture(GL_TEXTURE3);
    createImageTexture("cube_top.png", tex3, w, h, false);
    glBindTexture(GL_TEXTURE_2D, tex3);
    glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glActiveTexture(GL_TEXTURE4);
    createImageTexture("cube_back.png", tex4, w, h, false);
    glBindTexture(GL_TEXTURE_2D, tex4);
    glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glActiveTexture(GL_TEXTURE5);
    createImageTexture("cube_front.png", tex5, w, h, false);
    glBindTexture(GL_TEXTURE_2D, tex5);
    glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    
    E_Int shader = _shaders.shader_id(37);
    if (_shaders.currentShader() != shader) _shaders.activate(shader);
    _shaders[shader]->setUniform("cube_left", (int)0);
    _shaders[shader]->setUniform("cube_right", (int)1);
    _shaders[shader]->setUniform("cube_bottom", (int)2);
    _shaders[shader]->setUniform("cube_top", (int)3);
    _shaders[shader]->setUniform("cube_back", (int)4);
    _shaders[shader]->setUniform("cube_front", (int)5);
  }
  else if (mode == 2) // DOF + sobel + gamma
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
 
    E_Int shader = _shaders.shader_id(20);
    if (_shaders.currentShader() != shader) _shaders.activate(shader);
    _shaders[shader]->setUniform("FrameBuffer", (int)0);
    _shaders[shader]->setUniform("depthMap", (int)1);
    _shaders[shader]->setUniform("focalDepth", (float)ptrState->activePointZBuf);
    _shaders[shader]->setUniform("radius", (float)ptrState->dofPower);
    _shaders[shader]->setUniform("ext", (float)1.);
    _shaders[shader]->setUniform("toneMapping", (int)ptrState->toneMapping);
    _shaders[shader]->setUniform("gamma", (float)ptrState->gamma);
    _shaders[shader]->setUniform("sobelThreshold", (float)sobelThreshold);
    _shaders[shader]->setUniform("sharpenCoeff", (float)ptrState->sharpenPower);
    _shaders[shader]->setUniform("ssaoRadius", (float)ptrState->ssaoPower);
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

  if (mode == 1)
  {
    glDeleteTextures(1, &tex0);
    glDeleteTextures(1, &tex1);
    glDeleteTextures(1, &tex2);
    glDeleteTextures(1, &tex3);
    glDeleteTextures(1, &tex4);
    glDeleteTextures(1, &tex5);

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
  }

  glEnable(GL_DEPTH_TEST);
  glColor3f(1., 1., 1.);
#ifdef __SHADERS__
  _shaders.activate((short unsigned int)0);
#endif

}

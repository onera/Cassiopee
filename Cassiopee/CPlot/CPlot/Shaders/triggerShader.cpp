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
#include "shaders_id.h"
#include <iostream>

//=============================================================================
// Trigger le shader approprie en fonction de material, de scale et de color
//=============================================================================
void Data::triggerShader(Zone& z, int material, float scale, float* color)
{
  int shader = 1;
  float dx = K_FUNC::E_max(1.e-6, 1./(z.xmax-z.xmin));
  float dy = K_FUNC::E_max(1.e-6, 1./(z.ymax-z.ymin));
  float dz = K_FUNC::E_max(1.e-6, 1./(z.zmax-z.zmin));
  int t;
  
  switch (material)
  {
    case 1: // Glass
      shader = _shaders.shader_id(shader::glass);
      glActiveTexture(GL_TEXTURE0);
      glBindTexture(GL_TEXTURE_2D, _texEnviron1); // environnement
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
      glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

      glActiveTexture(GL_TEXTURE1);
      glBindTexture(GL_TEXTURE_2D, _texFrameBuffer[ptrState->frameBuffer]); // refraction
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      if (_shaders.currentShader() != shader) _shaders.activate(shader);
      // mix base color et le reste
      _shaders[shader]->setUniform("MixRatio", (float)0.5*z.shaderParam1);
      // mix envmap et refraction
      _shaders[shader]->setUniform("MixRatio2", (float)0.4*z.shaderParam2);
      _shaders[shader]->setUniform("FrameWidth", (float)_frameBufferSize[ptrState->frameBuffer]);
      _shaders[shader]->setUniform("FrameHeight", (float)_frameBufferSize[ptrState->frameBuffer]);
      _shaders[shader]->setUniform("EnvMap", (int)0);
      _shaders[shader]->setUniform("RefractionMap", (int)1);
      break;

    case 2: // Chrome
      shader = _shaders.shader_id(shader::chrome);
      SHADOWTEXTURE;
      glActiveTexture(GL_TEXTURE1);
      glBindTexture(GL_TEXTURE_2D, _texEnviron1); // environnement
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
      glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
      glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
      if (_shaders.currentShader() != shader) _shaders.activate(shader);
      _shaders.activate(shader);
      _shaders[shader]->setUniform("MixRatio", (float)0.6*z.shaderParam1);
      _shaders[shader]->setUniform("intensity", (float)z.shaderParam2);
      _shaders[shader]->setUniform("EnvMap", (int)1);
      _shaders[shader]->setUniform("shadow", (int)ptrState->shadow);
      _shaders[shader]->setUniform("ShadowMap", (int)0);
      break;

    case 3: // Metal
       shader = _shaders.shader_id(shader::metal);
       SHADOWTEXTURE;
       glActiveTexture(GL_TEXTURE1);
       glBindTexture(GL_TEXTURE_3D, _texNoise);
       if (_shaders.currentShader() != shader) _shaders.activate(shader);
       _shaders[shader]->setUniform("Scale", (float)0.3*scale*z.shaderParam2);
       _shaders[shader]->setUniform("bump", (float)(20.*(z.shaderParam2)+10.));
       _shaders[shader]->setUniform("Noise", (int)1);
       _shaders[shader]->setUniform("intensity", (float)(3.*z.shaderParam1+0.5), (float)(3.*z.shaderParam1+0.5), (float)(3.*z.shaderParam1+0.5));
       _shaders[shader]->setUniform("shadow", (int)ptrState->shadow);
       _shaders[shader]->setUniform("ShadowMap", (int)0);
      break;

    case 4: // Wood
      shader = _shaders.shader_id(shader::wood);
      SHADOWTEXTURE;
      glActiveTexture(GL_TEXTURE1);
      glBindTexture(GL_TEXTURE_3D, _texNoise);
      if (_shaders.currentShader() != shader) _shaders.activate(shader);
      _shaders[shader]->setUniform("Scale", (float)0.2*scale*z.shaderParam1);
      _shaders[shader]->setUniform("bump", (float)(20.*z.shaderParam2-20.));
      _shaders[shader]->setUniform(
        "LightWoodColor", (float)color[0], (float)color[1], (float)color[2]);
      float c[3];
      c[0] = color[0]*0.5; c[1] = color[1]*0.5; c[2] = color[2]*0.5;
      if (c[0] + c[1] + c[2] == 0.) {c[0] = 0.2; c[1] = 0.2; c[2] = 0.2;}
      _shaders[shader]->setUniform("DarkWoodColor", 
                                   (float)c[0], (float)c[1], (float)c[2]);
      _shaders[shader]->setUniform("RingFreq", (float)10.);
      _shaders[shader]->setUniform("NoiseScale", 
                                   (float)0.5, (float)0.1, (float)0.1);
      _shaders[shader]->setUniform("Noisiness", (float)4.);
      _shaders[shader]->setUniform("Noise", (int)1);
      _shaders[shader]->setUniform("GrainThreshold", (float)0.2);
      _shaders[shader]->setUniform("LightGrains", (float)1.);
      _shaders[shader]->setUniform("DarkGrains", (float)0.1);
      _shaders[shader]->setUniform("GrainScale", (float)100.);
      if (ptrState->dim == 3) _shaders[shader]->setUniform("lightOn", (int)1);
      else _shaders[shader]->setUniform("lightOn", (int)0);
      _shaders[shader]->setUniform("shadow", (int)ptrState->shadow);
      _shaders[shader]->setUniform("ShadowMap", (int)0);
      break;
      
    case 5: // Marble
      shader = _shaders.shader_id(shader::marble);
      SHADOWTEXTURE;
      glActiveTexture(GL_TEXTURE1);
      glBindTexture(GL_TEXTURE_3D, _texNoise);
      if (_shaders.currentShader() != shader) _shaders.activate(shader);
      _shaders[shader]->setUniform("Scale", (float)0.2*scale*z.shaderParam1);
      _shaders[shader]->setUniform(
        "MarbleColor",  (float)color[0], (float)color[1], (float)color[2]);
      
      _shaders[shader]->setUniform("VeinColor", 
                                   (float)0.1*z.shaderParam2, (float)0.1*z.shaderParam2, (float)0.1*z.shaderParam2);
      _shaders[shader]->setUniform("Noise", (int)1);
      if (ptrState->dim == 3) _shaders[shader]->setUniform("lightOn", (int)1);
      else _shaders[shader]->setUniform("lightOn", (int)0);
      _shaders[shader]->setUniform("shadow", (int)ptrState->shadow);
      _shaders[shader]->setUniform("ShadowMap", (int)0);
      break;
      
    case 6: // Smoke: volumetric shader
      shader = _shaders.shader_id(shader::smoke);
      glEnable(GL_CULL_FACE); glDepthMask(GL_FALSE);
      glActiveTexture(GL_TEXTURE0);
      glBindTexture(GL_TEXTURE_3D, _texVoxelBuffer);
#ifdef __SHADERS__ 
      if (glewIsSupported("GL_EXT_texture3D") != 0)
      {
        glTexImage3D(GL_TEXTURE_3D, 0, GL_LUMINANCE,
                     _voxelBufferSize, _voxelBufferSize, _voxelBufferSize, 
                     0, GL_LUMINANCE, GL_UNSIGNED_BYTE, z._voxelArray);
      }
#endif
      if (_shaders.currentShader() != shader) _shaders.activate(shader);
      //_shaders[shader]->setUniform("BumpDensity", (float)2.);
      //_shaders[shader]->setUniform("BumpSize", (float)60.*scale);
      //_shaders[shader]->setUniform("SpecularFactor", (float)0.6);
      _shaders[shader]->setUniform("DensityTex", (int)0);
      _shaders[shader]->setUniform("Absorption", (float)12.*z.shaderParam2);
      _shaders[shader]->setUniform("xo", (float)z.xmin, (float)z.ymin, (float)z.zmin);
      _shaders[shader]->setUniform("dx", (float)dx, (float)dy, (float)dz);
      _shaders[shader]->setUniform("poscam", (float)_view.xcam, (float)_view.ycam, (float)_view.zcam);
      break;

    case 7: // XRay
      shader = _shaders.shader_id(shader::xray);
      /*glEnable(GL_CULL_FACE);*/ glDepthMask(GL_FALSE);
      if (_shaders.currentShader() != shader) _shaders.activate(shader);
      _shaders[shader]->setUniform("EdgeFalloff", (float)0.9*z.shaderParam1);
      _shaders[shader]->setUniform("intensity", (float)z.shaderParam2);
      break;

    case 8: // Granite
      shader = _shaders.shader_id(shader::granite);
      SHADOWTEXTURE;
      glActiveTexture(GL_TEXTURE1);
      glBindTexture(GL_TEXTURE_3D, _texNoise);
      if (_shaders.currentShader() != shader) _shaders.activate(shader);
      _shaders[shader]->setUniform("Scale", (float)0.3*scale*z.shaderParam1);
      _shaders[shader]->setUniform("bump", (float)(20.*z.shaderParam2+10.));
      _shaders[shader]->setUniform("Noise", (int)1);
      if (ptrState->dim == 3) _shaders[shader]->setUniform("lightOn", (int)1);
      else _shaders[shader]->setUniform("lightOn", (int)0);
      _shaders[shader]->setUniform("shadow", (int)ptrState->shadow);
      _shaders[shader]->setUniform("ShadowMap", (int)0);
      break;

    case 9: // billboarding
      if (z.shaderParam2 < 0.1) // Sphere billboard
      {
        shader = _shaders.shader_id(shader::sphere_billboarding);
        ptrState->billBoardNi = 1;
        ptrState->billBoardNj = 1;
        ptrState->billBoardWidth = 1;
        ptrState->billBoardHeight = 1;
        if (_shaders.currentShader() != shader) _shaders.activate(shader);
        _shaders[shader]->setUniform("specularFactor", (float)1.);
        _shaders[shader]->setUniform("diffuseFactor", (float)1.);
      }
      else if (z.shaderParam2 >= 0.1) // texture billboards
      {
        int type = floor( (z.shaderParam2-0.1)*(_nBillBoards-1)/1.9 +0.5 );
        int t = type;
        t = std::max(t, 0);
        t = std::min(t, (int)_nBillBoards-1);
        if (_billBoardTexs[t] == 0) 
        {
          createImageTexture(_billBoardFiles[t], _billBoardTexs[t], _billBoardWidths[t], _billBoardHeights[t], true);
          //printf("loading %d %s\n", type, _billBoardFiles[t]);
        }
        ptrState->billBoardNi = _billBoardNis[t];
        ptrState->billBoardNj = _billBoardNjs[t];
        ptrState->billBoardWidth = _billBoardWidths[t];
        ptrState->billBoardHeight = _billBoardHeights[t];
        _texBillBoard = _billBoardTexs[t];

        shader = _shaders.shader_id(shader::billboard);
        //glEnable(GL_CULL_FACE);
        glDepthMask(GL_FALSE);
        
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, _texBillBoard);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
        if (_shaders.currentShader() != shader) _shaders.activate(shader);
        _shaders[shader]->setUniform("BTex", (int)0);
        _shaders[shader]->setUniform("Ni", (int)ptrState->billBoardNi);
        _shaders[shader]->setUniform("Nj", (int)ptrState->billBoardNj);
      }
      break;

    case 10: // Bricks
      shader = _shaders.shader_id(shader::brick);
      SHADOWTEXTURE;
      if (_shaders.currentShader() != shader) _shaders.activate(shader);
      _shaders[shader]->setUniform("MortarColor", (float)1., (float)1., (float)1.);
      _shaders[shader]->setUniform("BrickSize", (float)2./scale*z.shaderParam1, (float)1./scale*z.shaderParam1);
      _shaders[shader]->setUniform("BrickPct", (float)0.9*z.shaderParam2, (float)0.85*z.shaderParam2);
      _shaders[shader]->setUniform("shadow", (int)ptrState->shadow);
      _shaders[shader]->setUniform("ShadowMap", (int)0);
      break;

    case 11: // Clouds
      shader = _shaders.shader_id(shader::cloud);
      SHADOWTEXTURE;
      glActiveTexture(GL_TEXTURE1);
      glBindTexture(GL_TEXTURE_3D, _texNoise);
      if (_shaders.currentShader() != shader) _shaders.activate(shader);
      _shaders[shader]->setUniform("Scale", (float)0.05*scale*z.shaderParam1);
      _shaders[shader]->setUniform("Noise", (int)1);
      _shaders[shader]->setUniform("Offset", (float)0., (float)0., (float)0.);
      _shaders[shader]->setUniform("CloudColor", (float)0.8, (float)0.8, (float)0.8);
      _shaders[shader]->setUniform("shadow", (int)ptrState->shadow);
      _shaders[shader]->setUniform("ShadowMap", (int)0);
      break;

    case 12: // gooch
      shader = _shaders.shader_id(shader::gooch);
      SHADOWTEXTURE;
      if (_shaders.currentShader() != shader) _shaders.activate(shader);
      _shaders[shader]->setUniform("specularFactor", (float)z.shaderParam2);
      _shaders[shader]->setUniform("exponent", (float)z.shaderParam1);
      _shaders[shader]->setUniform("shadow", (int)ptrState->shadow);
      _shaders[shader]->setUniform("ShadowMap", (int)0);
      break;

    case 13: // flat shader (no light)
      shader = _shaders.shader_id(shader::flat);
      SHADOWTEXTURE;
      if (_shaders.currentShader() != shader) _shaders.activate(shader);
      _shaders[shader]->setUniform("shadow", (int)ptrState->shadow);
      _shaders[shader]->setUniform("ShadowMap", (int)0);
      break;

    case 14: // textured
      if (_nMaterials > 0)
      {
        shader = _shaders.shader_id(shader::textured_material);
        SHADOWTEXTURE;
        // choix de la texture
        t = (int)round( (z.shaderParam2)*((_nMaterials-1)*0.5) );
        t = std::max(t, 0);
        t = std::min(t, (int)_nMaterials-1);
        glActiveTexture(GL_TEXTURE1);
        if (_materialTexs[t] == 0) 
        {
           createImageTexture(_materialFiles[t], _materialTexs[t], _materialWidths[t], _materialHeights[t], true);
          //printf("loading %d %s\n", 0, _materialFiles[0]);
        }
        glBindTexture(GL_TEXTURE_2D, _materialTexs[t]);
        
        // choix de la bump map (if any)
        bool hasBump = false;
        if (t < _nBumpMaps)
        {
          glActiveTexture(GL_TEXTURE2);
          if (_bumpMapFiles[t] != NULL)
          {
            hasBump = true;
            if (_bumpMapTexs[t] == 0) 
            {
              createImageTexture(_bumpMapFiles[t], _bumpMapTexs[t], _bumpMapWidths[t], _bumpMapHeights[t], true);
            }
          }
          glBindTexture(GL_TEXTURE_2D, _bumpMapTexs[t]);
        }
        if (hasBump)
        {
          glActiveTexture(GL_TEXTURE2);
          glBindTexture(GL_TEXTURE_2D, _bumpMapTexs[t]);
        }

        if (_shaders.currentShader() != shader) _shaders.activate(shader);
        _shaders[shader]->setUniform("specularFactor", (float)z.shaderParam1);
        _shaders[shader]->setUniform("shadow", (int)ptrState->shadow);
        _shaders[shader]->setUniform("hasBump", (int)hasBump);
        float blend = 1.;
        if (ptrState->selectionStyle == 1) // alpha selection
        {
          if (ptrState->selectedZone <= 0) blend = 1;
          else if (z.selected == 0) blend = 0.12;
        }
        if (z.blending != -1.) blend = std::min((float)z.blending, blend);
        _shaders[shader]->setUniform("blend", blend);
        _shaders[shader]->setUniform("ShadowMap", (int)0);
        _shaders[shader]->setUniform("Texmat0", (int)1);
        if (hasBump) _shaders[shader]->setUniform("Texbump0", (int)2);
      }
      else
      {
        _shaders.activate((short unsigned int)0); // no shader
      }
      break;

    default: // phong
      if (ptrState->dim == 3)
      {        
        if (ptrState->solidStyle == 0) // one side
        {
          shader = _shaders.shader_id(shader::unidirectionnel_phong);
          SHADOWTEXTURE;
          if (_shaders.currentShader() != shader)
            _shaders.activate((short unsigned int)shader);
          if (ptrState->mode == RENDER)
          {
            _shaders[shader]->setUniform("specularFactor", (float)z.shaderParam1);
            _shaders[shader]->setUniform("diffuseFactor", (float)z.shaderParam2); 
          }
          else
          {
            _shaders[shader]->setUniform("specularFactor", (float)1.);
            _shaders[shader]->setUniform("diffuseFactor", (float)1.);  
          }
          _shaders[shader]->setUniform("shadow", (int)ptrState->shadow);
          _shaders[shader]->setUniform("ShadowMap", (int)0);
        }
        else if (ptrState->solidStyle == 1 || ptrState->solidStyle == 3 || ptrState->solidStyle == 4) // two sides
        {
          shader = _shaders.shader_id(shader::bidirectionnel_phong);
          SHADOWTEXTURE;
          if (_shaders.currentShader() != shader)
            _shaders.activate((short unsigned int)shader);
          if (ptrState->mode == RENDER)
          {
            _shaders[shader]->setUniform("specularFactor", (float)z.shaderParam1);
            _shaders[shader]->setUniform("diffuseFactor", (float)z.shaderParam2);
          }
          else
          {
            _shaders[shader]->setUniform("specularFactor", (float)1.);
            _shaders[shader]->setUniform("diffuseFactor", (float)1.);
          }
          _shaders[shader]->setUniform("shadow", (int)ptrState->shadow);
          _shaders[shader]->setUniform("ShadowMap", (int)0);
        }
        else _shaders.activate((short unsigned int)shader::None);
      }
      else _shaders.activate((short unsigned int)shader::None);
  }
}

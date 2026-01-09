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
// Shaders pour les isos en mode render
E_Float scale = MAX(zonep->xmax-zonep->xmin, zonep->ymax-zonep->ymin);
scale = MAX(scale, zonep->zmax-zonep->zmin);
scale = 100./(scale+1.e-12);
E_Int shader;
//float alpha = 0.004;
//float dx = (_view.xeye - _view.xcam)*alpha;
//float dy = (_view.yeye - _view.ycam)*alpha;
//float dz = (_view.zeye - _view.zcam)*alpha;
E_Int nofield = -E_Int(zonep->colorR)-2;
E_Int isoLight;
E_Int colormap = (E_Int)(_isoColormap[nofield]);
if (colormap != -1) // isoScales prevails
{
  isoLight = colormap%2;
  colormap = (E_Int)(colormap*0.5);
}
else
{
  colormap = (E_Int)((ptrState->colormap)*0.5);
  isoLight = ptrState->isoLight;
}
//printf("I use colormap %d and isoLight %d\n", colormap, isoLight); fflush(stdout);

switch (zonep->material)
{ 
  case 1: // iso+glass
  {
    shader = _shaders.shader_id(shader::iso_glass);
    glActiveTexture(GL_TEXTURE0);
    if (_texColormap == 0) createColormapTexture();
    fillColormapTexture(colormap);    
    glBindTexture(GL_TEXTURE_1D, _texColormap);
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, _texEnviron1); // environnement
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D, _texFrameBuffer[ptrState->frameBuffer]); // refraction
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

    if (_shaders.currentShader() != shader) 
      _shaders.activate((short unsigned int)shader); 
      
    _shaders[shader]->setUniform("colormap", (int)0);
      
    if (_niso[nofield] == -1)
    {
      _shaders[shader]->setUniform("niso", (float)ptrState->niso);
      _shaders[shader]->setUniform("alpha", (float)1.);
      _shaders[shader]->setUniform("beta", (float)0.);
      _shaders[shader]->setUniform("amin", (float)0.);
      _shaders[shader]->setUniform("amax", (float)1.);
    }
    else 
    { 
      double rmin, rmax, alpha, beta;
      double deltai = MAX(maxf[nofield]-minf[nofield], ISOCUTOFF);
      rmin = (_isoMin[nofield] -minf[nofield])/deltai;
      rmax = (_isoMax[nofield] -minf[nofield])/deltai;
      alpha = 1./MAX(rmax-rmin, ISOCUTOFF); beta = -rmin*alpha;
      _shaders[shader]->setUniform("niso", (float)_niso[nofield]);
      _shaders[shader]->setUniform("alpha", (float)alpha);
      _shaders[shader]->setUniform("beta", (float)beta);
      double amin = (_isoAlphaMin[nofield] - minf[nofield])/deltai;
      double amax = (_isoAlphaMax[nofield] - minf[nofield])/deltai;
      _shaders[shader]->setUniform("amin", (float)amin);
      _shaders[shader]->setUniform("amax", (float)amax); 
    }
    double resf = std::max(_view.w / 1080., 1.);
    _shaders[shader]->setUniform("edgeStyle", (float)(resf*ptrState->isoEdges));
  }
  // mix couleur de base
  _shaders[shader]->setUniform("MixRatio", (float)0.6*zonep->shaderParam1);
  // mix reflection / refraction
  _shaders[shader]->setUniform("MixRatio2", (float)0.4*zonep->shaderParam2);
  _shaders[shader]->setUniform("FrameWidth", (float)_frameBufferSize[ptrState->frameBuffer]);
  _shaders[shader]->setUniform("FrameHeight", (float)_frameBufferSize[ptrState->frameBuffer]);
  _shaders[shader]->setUniform("EnvMap", (int)1);
  _shaders[shader]->setUniform("RefractionMap", (int)2);
  break;

  case 2: // iso+envmap (chrome) 
  {
    shader = _shaders.shader_id(shader::iso_chrome);
    SHADOWTEXTURE;
    glActiveTexture(GL_TEXTURE1);
    if (_texColormap == 0) createColormapTexture();
    fillColormapTexture(colormap);
    glBindTexture(GL_TEXTURE_1D, _texColormap);
    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_2D, _texEnviron1); // environnement
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);

    if (_shaders.currentShader() != shader)
      _shaders.activate((short unsigned int)shader);

    _shaders[shader]->setUniform("colormap", (int)1);
      
    if (_niso[nofield] == -1)
    {
      _shaders[shader]->setUniform("niso", (float)ptrState->niso);
      _shaders[shader]->setUniform("alpha", (float)1.);
      _shaders[shader]->setUniform("beta", (float)0.);
      _shaders[shader]->setUniform("amin", (float)0.);
      _shaders[shader]->setUniform("amax", (float)1.);
    }
    else
    {
      double rmin, rmax, alpha, beta;
      double deltai = MAX(maxf[nofield]-minf[nofield], ISOCUTOFF);
      rmin = (_isoMin[nofield] -minf[nofield])/deltai;
      rmax = (_isoMax[nofield] -minf[nofield])/deltai;
      alpha = 1./MAX(rmax-rmin, ISOCUTOFF); beta = -rmin*alpha;
      _shaders[shader]->setUniform("niso", (float)_niso[nofield]);
      _shaders[shader]->setUniform("alpha", (float)alpha);
      _shaders[shader]->setUniform("beta", (float)beta);
      double amin = (_isoAlphaMin[nofield] - minf[nofield])/deltai;
      double amax = (_isoAlphaMax[nofield] - minf[nofield])/deltai;
      _shaders[shader]->setUniform("amin", (float)amin);
      _shaders[shader]->setUniform("amax", (float)amax);
    }
    double resf = std::max(_view.w / 1080., 1.);
    _shaders[shader]->setUniform("edgeStyle", (float)(resf*ptrState->isoEdges));
  }
  _shaders[shader]->setUniform("MixRatio", (float)0.9*zonep->shaderParam1);
  _shaders[shader]->setUniform("EnvMap", (int)2);
  _shaders[shader]->setUniform("shadow", (int)ptrState->shadow);
  _shaders[shader]->setUniform("ShadowMap", (int)0);
  break;

  case 8: // iso+granite 
  {
    shader = _shaders.shader_id(shader::iso_granite);
    SHADOWTEXTURE;
    glActiveTexture(GL_TEXTURE1);
    if (_texColormap == 0) createColormapTexture();
    fillColormapTexture(colormap);
    glBindTexture(GL_TEXTURE_1D, _texColormap);
    glActiveTexture(GL_TEXTURE2);
    glBindTexture(GL_TEXTURE_3D, _texNoise);
    if (_shaders.currentShader() != shader)
      _shaders.activate((short unsigned int)shader);

    _shaders[shader]->setUniform("Scale", 
                                 (float)(0.3*scale*zonep->shaderParam1));
    _shaders[shader]->setUniform("bump", (float)(20.*zonep->shaderParam2+10.));
    _shaders[shader]->setUniform("Noise", (int)2);
    _shaders[shader]->setUniform("colormap", (int)1);
    
    if (_niso[nofield] == -1)
    {
      _shaders[shader]->setUniform("niso", (float)ptrState->niso);
      _shaders[shader]->setUniform("alpha", (float)1.);
      _shaders[shader]->setUniform("beta", (float)0.);
      _shaders[shader]->setUniform("amin", (float)0.);
      _shaders[shader]->setUniform("amax", (float)1.);
    }
    else 
    { 
      double rmin, rmax, alpha, beta;
      double deltai = MAX(maxf[nofield]-minf[nofield], ISOCUTOFF);
      rmin = (_isoMin[nofield] -minf[nofield])/deltai;
      rmax = (_isoMax[nofield] -minf[nofield])/deltai;
      alpha = 1./MAX(rmax-rmin, ISOCUTOFF); beta = -rmin*alpha;
      _shaders[shader]->setUniform("niso", (float)_niso[nofield]);
      _shaders[shader]->setUniform("alpha", (float)alpha);
      _shaders[shader]->setUniform("beta", (float)beta);
      double amin = (_isoAlphaMin[nofield] - minf[nofield])/deltai;
      double amax = (_isoAlphaMax[nofield] - minf[nofield])/deltai;
      _shaders[shader]->setUniform("amin", (float)amin);
      _shaders[shader]->setUniform("amax", (float)amax);
    }
    double resf = std::max(_view.w / 1080., 1.);
    _shaders[shader]->setUniform("edgeStyle", (float)(resf*ptrState->isoEdges));
    if (isoLight == 1 && ptrState->dim == 3)
      _shaders[shader]->setUniform("lightOn", (int)1);
    else _shaders[shader]->setUniform("lightOn", (int)0);
    _shaders[shader]->setUniform("shadow", (int)ptrState->shadow);
    _shaders[shader]->setUniform("ShadowMap", (int)0);
  }
  break;

 case 10: // iso+brick 
 {
    shader = _shaders.shader_id(shader::iso_brick);
    glActiveTexture(GL_TEXTURE0);
    if (_texColormap == 0) createColormapTexture();
    fillColormapTexture(colormap);
    glBindTexture(GL_TEXTURE_1D, _texColormap);
    if (_shaders.currentShader() != shader)
      _shaders.activate((short unsigned int)shader);
    _shaders[shader]->setUniform("colormap", (int)0);
    
    if (_niso[nofield] == -1)
    {
      _shaders[shader]->setUniform("niso", (float)ptrState->niso);
      _shaders[shader]->setUniform("alpha", (float)1.);
      _shaders[shader]->setUniform("beta", (float)0.);
      _shaders[shader]->setUniform("amin", (float)0.);
      _shaders[shader]->setUniform("amax", (float)1.);
    }
    else 
    { 
      double rmin, rmax, alpha, beta;
      double deltai = MAX(maxf[nofield]-minf[nofield], ISOCUTOFF);
      rmin = (_isoMin[nofield] -minf[nofield])/deltai;
      rmax = (_isoMax[nofield] -minf[nofield])/deltai;
      alpha = 1./MAX(rmax-rmin, ISOCUTOFF); beta = -rmin*alpha;
      _shaders[shader]->setUniform("niso", (float)_niso[nofield]);
      _shaders[shader]->setUniform("alpha", (float)alpha);
      _shaders[shader]->setUniform("beta", (float)beta);
      double amin = (_isoAlphaMin[nofield] - minf[nofield])/deltai;
      double amax = (_isoAlphaMax[nofield] - minf[nofield])/deltai;
      _shaders[shader]->setUniform("amin", (float)amin);
      _shaders[shader]->setUniform("amax", (float)amax);
    }
    _shaders[shader]->setUniform("edgeStyle", (float)zonep->shaderParam1);
    if (isoLight == 1 && ptrState->dim == 3)
      _shaders[shader]->setUniform("lightOn", (int)1);
    else _shaders[shader]->setUniform("lightOn", (int)0);
  }
  break;

 case 12: // iso+gooch 
 {
    shader = _shaders.shader_id(shader::iso_gooch);
    glActiveTexture(GL_TEXTURE0);
    if (_texColormap == 0) createColormapTexture();
    fillColormapTexture(colormap);
    glBindTexture(GL_TEXTURE_1D, _texColormap);
    if (_shaders.currentShader() != shader)
      _shaders.activate((short unsigned int)shader);
    _shaders[shader]->setUniform("colormap", (int)0);
    
    if (_niso[nofield] == -1)
    {
      _shaders[shader]->setUniform("niso", (float)ptrState->niso);
      _shaders[shader]->setUniform("alpha", (float)1.);
      _shaders[shader]->setUniform("beta", (float)0.);
      _shaders[shader]->setUniform("amin", (float)0.);
      _shaders[shader]->setUniform("amax", (float)1.);
    }
    else 
    { 
      float rmin, rmax, alpha, beta;
      float deltai = MAX(maxf[nofield]-minf[nofield], ISOCUTOFF);
      rmin = (_isoMin[nofield] -minf[nofield])/deltai;
      rmax = (_isoMax[nofield] -minf[nofield])/deltai;
      alpha = 1./MAX(rmax-rmin, ISOCUTOFF); beta = -rmin*alpha;
      _shaders[shader]->setUniform("niso", (float)_niso[nofield]);
      _shaders[shader]->setUniform("alpha", (float)alpha);
      _shaders[shader]->setUniform("beta", (float)beta);
      float amin = (_isoAlphaMin[nofield] - minf[nofield])/deltai;
      float amax = (_isoAlphaMax[nofield] - minf[nofield])/deltai;
      _shaders[shader]->setUniform("amin", (float)amin);
      _shaders[shader]->setUniform("amax", (float)amax);
    }
    _shaders[shader]->setUniform("edgeStyle", (float)zonep->shaderParam1);
    //if (isoLight == 1 && ptrState->dim == 3)
    //  _shaders[shader]->setUniform("lightOn", (int)1);
    //else _shaders[shader]->setUniform("lightOn", (int)0);
  }
  break;

  case 13: // iso+flat 
  {
    shader = _shaders.shader_id(shader::iso_flat);
    SHADOWTEXTURE;
    glActiveTexture(GL_TEXTURE1);
    if (_texColormap == 0) createColormapTexture();
    fillColormapTexture(colormap);
    glBindTexture(GL_TEXTURE_1D, _texColormap);
    if (_shaders.currentShader() != shader)
      _shaders.activate((short unsigned int)shader);
    _shaders[shader]->setUniform("colormap", (int)1);
    
    if (_niso[nofield] == -1)
    {
      _shaders[shader]->setUniform("niso", (float)ptrState->niso);
      _shaders[shader]->setUniform("alpha", (float)1.);
      _shaders[shader]->setUniform("beta", (float)0.);
      _shaders[shader]->setUniform("amin", (float)0.);
      _shaders[shader]->setUniform("amax", (float)1.);
    }
    else 
    { 
      double rmin, rmax, alpha, beta;
      double deltai = MAX(maxf[nofield]-minf[nofield], ISOCUTOFF);
      rmin = (_isoMin[nofield] -minf[nofield])/deltai;
      rmax = (_isoMax[nofield] -minf[nofield])/deltai;
      alpha = 1./MAX(rmax-rmin, ISOCUTOFF); beta = -rmin*alpha;
      _shaders[shader]->setUniform("niso", (float)_niso[nofield]);
      _shaders[shader]->setUniform("alpha", (float)alpha);
      _shaders[shader]->setUniform("beta", (float)beta);
      double amin = (_isoAlphaMin[nofield] - minf[nofield])/deltai;
      double amax = (_isoAlphaMax[nofield] - minf[nofield])/deltai;
      _shaders[shader]->setUniform("amin", (float)amin);
      _shaders[shader]->setUniform("amax", (float)amax);
    }
    double resf = std::max(_view.w / 1080., 1.);
    _shaders[shader]->setUniform("edgeStyle", (float)(resf*ptrState->isoEdges));
    _shaders[shader]->setUniform("shadow", (int)ptrState->shadow);
    _shaders[shader]->setUniform("ShadowMap", (int)0);
  }
  break;

  case 7: // iso + Xray
  {
    shader = _shaders.shader_id(shader::iso_xray);
    glActiveTexture(GL_TEXTURE1);
    if (_texColormap == 0) createColormapTexture();
    fillColormapTexture(colormap);
    glBindTexture(GL_TEXTURE_1D, _texColormap);
    if (_shaders.currentShader() != shader)
      _shaders.activate((short unsigned int)shader);
    _shaders[shader]->setUniform("colormap", (int)1);
    
    if (_niso[nofield] == -1)
    {
      _shaders[shader]->setUniform("niso", (float)ptrState->niso);
      _shaders[shader]->setUniform("alpha", (float)1.);
      _shaders[shader]->setUniform("beta", (float)0.);
      _shaders[shader]->setUniform("amin", (float)0.);
      _shaders[shader]->setUniform("amax", (float)1.);
    }
    else 
    { 
      double rmin, rmax, alpha, beta;
      double deltai = MAX(maxf[nofield]-minf[nofield], ISOCUTOFF);
      rmin = (_isoMin[nofield] -minf[nofield])/deltai;
      rmax = (_isoMax[nofield] -minf[nofield])/deltai;
      alpha = 1./MAX(rmax-rmin, ISOCUTOFF); beta = -rmin*alpha;
      _shaders[shader]->setUniform("niso", (float)_niso[nofield]);
      _shaders[shader]->setUniform("alpha", (float)alpha);
      _shaders[shader]->setUniform("beta", (float)beta);
      double amin = (_isoAlphaMin[nofield] - minf[nofield])/deltai;
      double amax = (_isoAlphaMax[nofield] - minf[nofield])/deltai;
      _shaders[shader]->setUniform("amin", (float)amin);
      _shaders[shader]->setUniform("amax", (float)amax);
    }
    double resf = std::max(_view.w / 1080., 1.);
    _shaders[shader]->setUniform("edgeStyle", (float)(resf*ptrState->isoEdges));
    _shaders[shader]->setUniform("EdgeFalloff", (float)0.9*zonep->shaderParam1);
    //_shaders[shader]->setUniform("shadow", (int)ptrState->shadow);
    //_shaders[shader]->setUniform("ShadowMap", (int)0);
    _shaders[shader]->setUniform("blend", (float)1.);
    glDepthMask(GL_FALSE);
  }
  break; 

  case 3: // iso+metal (anisotropic) 
  {
    shader = _shaders.shader_id(shader::iso_metal);
    SHADOWTEXTURE;
    glActiveTexture(GL_TEXTURE1);
    if (_texColormap == 0) createColormapTexture();
    fillColormapTexture(colormap);
    glBindTexture(GL_TEXTURE_1D, _texColormap);
    if (_shaders.currentShader() != shader)
      _shaders.activate((short unsigned int)shader);
    _shaders[shader]->setUniform("intensity", (float)(2.*zonep->shaderParam1+0.5), (float)(2.*zonep->shaderParam1+0.5), (float)(2.*zonep->shaderParam1+0.5));
    _shaders[shader]->setUniform("colormap", (int)1);
    
    if (_niso[nofield] == -1)
    {
      _shaders[shader]->setUniform("niso", (float)ptrState->niso);
      _shaders[shader]->setUniform("alpha", (float)1.);
      _shaders[shader]->setUniform("beta", (float)0.);
      _shaders[shader]->setUniform("amin", (float)0.);
      _shaders[shader]->setUniform("amax", (float)1.);
    }
    else 
    { 
      double rmin, rmax, alpha, beta;
      double deltai = MAX(maxf[nofield]-minf[nofield], ISOCUTOFF);
      rmin = (_isoMin[nofield] -minf[nofield])/deltai;
      rmax = (_isoMax[nofield] -minf[nofield])/deltai;
      alpha = 1./MAX(rmax-rmin, ISOCUTOFF); beta = -rmin*alpha;
      _shaders[shader]->setUniform("niso", (float)_niso[nofield]);
      _shaders[shader]->setUniform("alpha", (float)alpha);
      _shaders[shader]->setUniform("beta", (float)beta);
      double amin = (_isoAlphaMin[nofield] - minf[nofield])/deltai;
      double amax = (_isoAlphaMax[nofield] - minf[nofield])/deltai;
      _shaders[shader]->setUniform("amin", (float)amin);
      _shaders[shader]->setUniform("amax", (float)amax);
    }
    double resf = std::max(_view.w / 1080., 1.);
    _shaders[shader]->setUniform("edgeStyle", (float)(resf*ptrState->isoEdges));
    _shaders[shader]->setUniform("shadow", (int)ptrState->shadow);
    _shaders[shader]->setUniform("ShadowMap", (int)0);
  }
  break;

  default: // iso+solid
  {
    shader = _shaders.shader_id(shader::iso_banded_colormap);
    SHADOWTEXTURE;
    glActiveTexture(GL_TEXTURE1);
    if (_texColormap == 0) createColormapTexture();
    fillColormapTexture(colormap);
    glBindTexture(GL_TEXTURE_1D, _texColormap);
    
    if (_shaders.currentShader() != shader)
      _shaders.activate((short unsigned int)shader);
    _shaders[shader]->setUniform("colormap", (int)1);
    
    if (_niso[nofield] == -1)
    {
      _shaders[shader]->setUniform("niso", (float)ptrState->niso);
      _shaders[shader]->setUniform("alpha", (float)1.);
      _shaders[shader]->setUniform("beta", (float)0.);
      _shaders[shader]->setUniform("amin", (float)0.);
      _shaders[shader]->setUniform("amax", (float)1.); 
    }
    else 
    { 
      double rmin, rmax, alpha, beta;
      double deltai = MAX(maxf[nofield]-minf[nofield], ISOCUTOFF);
      rmin = (_isoMin[nofield] -minf[nofield])/deltai;
      rmax = (_isoMax[nofield] -minf[nofield])/deltai;
      alpha = 1./MAX(rmax-rmin, ISOCUTOFF); beta = -rmin*alpha;
      _shaders[shader]->setUniform("niso", (float)_niso[nofield]);
      _shaders[shader]->setUniform("alpha", (float)alpha);
      _shaders[shader]->setUniform("beta", (float)beta);
      double amin = (_isoAlphaMin[nofield] - minf[nofield])/deltai;
      double amax = (_isoAlphaMax[nofield] - minf[nofield])/deltai;
      _shaders[shader]->setUniform("amin", (float)amin);
      _shaders[shader]->setUniform("amax", (float)amax); 
    }
    double resf = std::max(_view.w / 1080., 1.);
    _shaders[shader]->setUniform("edgeStyle", (float)(resf*ptrState->isoEdges));
    if (isoLight == 1 && ptrState->dim == 3)
      _shaders[shader]->setUniform("lightOn", (int)1);
    else _shaders[shader]->setUniform("lightOn", (int)0);
    _shaders[shader]->setUniform("shadow", (int)ptrState->shadow);
    _shaders[shader]->setUniform("ShadowMap", (int)0);
  }
}

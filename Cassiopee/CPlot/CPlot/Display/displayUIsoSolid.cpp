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
#include "DataDL.h"
#include "ZoneImplDL.h"

//============================================================================
/*
  Display a unstructured mesh solution using fully colored planes.
*/
//============================================================================
void DataDL::displayUIsoSolid()
{
  E_Int zone, zonet;
  if (_numberOfUnstructZones == 0) return;

  // Blending
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_POLYGON_OFFSET_FILL);
  glEnable(GL_POLYGON_OFFSET_LINE);

  // Find colormap and isoLight
  E_Int isoLight;
  E_Int nofield = ptrState->scalarField;
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


#ifdef __SHADERS__
  SHADOWTEXTURE;
  if (ptrState->mode == SCALARFIELD)
  { // shader pour les isos scalaires
    glActiveTexture(GL_TEXTURE1);
    if (_texColormap == 0) createColormapTexture();
    fillColormapTexture(colormap);
    glBindTexture(GL_TEXTURE_1D, _texColormap);
    int s = _shaders.shader_id(shader::iso_banded_colormap);
    if (ptrState->scalarStyle == 2 || ptrState->scalarStyle == 3) s = _shaders.shader_id(shader::iso_colored_lines);
    if (_shaders.currentShader() != s) _shaders.activate((short unsigned int)s);
    _shaders[s]->setUniform("colormap", (int)1);
    if (_niso[nofield] == -1)
    {
      _shaders[s]->setUniform("alpha", (float)1.);
      _shaders[s]->setUniform("beta", (float)0.);
      _shaders[s]->setUniform("niso", (float)ptrState->niso);
      _shaders[s]->setUniform("amin", (float)0.);
      _shaders[s]->setUniform("amax", (float)1.);
    }
    else 
    {
      double rmin, rmax, alpha, beta;
      double deltai = MAX(maxf[nofield]-minf[nofield], ISOCUTOFF);
      rmin = (_isoMin[nofield] -minf[nofield])/deltai;
      rmax = (_isoMax[nofield] -minf[nofield])/deltai;
      alpha = 1./MAX(rmax-rmin, ISOCUTOFF); beta = -rmin*alpha;
      _shaders[s]->setUniform("niso", (float)_niso[nofield]);
      _shaders[s]->setUniform("alpha", (float)alpha);
      _shaders[s]->setUniform("beta", (float)beta);
      double amin = (_isoAlphaMin[nofield] - minf[nofield])/deltai;
      double amax = (_isoAlphaMax[nofield] - minf[nofield])/deltai;
      _shaders[s]->setUniform("amin", (float)amin);
      _shaders[s]->setUniform("amax", (float)amax); 
    }
    double resf = std::max(_view.w / 1080., 1.);
    _shaders[s]->setUniform("edgeStyle", (float)(resf*ptrState->isoEdges));
    _shaders[s]->setUniform("lightOn", (int)0);
    _shaders[s]->setUniform("shadow", (int)ptrState->shadow);
    _shaders[s]->setUniform("ShadowMap", (int)0);
  }
  else
  { // shader pour les isos vectoriels
    int s = _shaders.shader_id(shader::vector_rgb);
    if (ptrState->vectorStyle == 1) s = _shaders.shader_id(shader::vector_arrow);
    else if (ptrState->vectorStyle == 2) s = _shaders.shader_id(shader::vector_line);
    
    if (_shaders.currentShader() != s) _shaders.activate((short unsigned int)s);
    if (s == _shaders.shader_id(shader::vector_rgb))
    {
      _shaders[s]->setUniform("lightOn", (int)0);
      _shaders[s]->setUniform("shadow", (int)ptrState->shadow);
      _shaders[s]->setUniform("ShadowMap", (int)0);
    }
    else
    {
      _shaders[s]->setUniform("lightOn", (int)0);
      _shaders[s]->setUniform("shadow", (int)ptrState->shadow);
      _shaders[s]->setUniform("ShadowMap", (int)0);
      //double diag = 0.01*sqrt((xmax-xmin)*(xmax-xmin)+(ymax-ymin)*(ymax-ymin)+(zmax-zmin)*(zmax-zmin));
      double sc = ptrState->vectorScale/100.;
      double ed = sqrt( (_view.xcam-_view.xeye)*(_view.xcam-_view.xeye)+(_view.ycam-_view.yeye)*(_view.ycam-_view.yeye)+(_view.zcam-_view.zeye)*(_view.zcam-_view.zeye) )*0.1;
      _shaders[s]->setUniform("scale", float(sc*ed));
      _shaders[s]->setUniform("fix_length", (int)ptrState->vectorNormalize);
      _shaders[s]->setUniform("density", float(ptrState->vectorDensity));
      glActiveTexture(GL_TEXTURE1);
      if (_texColormap == 0) createColormapTexture();
      fillColormapTexture((int)_pref.colorMap->varName[0]-48);
      _shaders[s]->setUniform("colormap", (int)1);
      if ( s == _shaders.shader_id(shader::vector_arrow) )
      {
        _shaders[s]->setUniform("show_surface", (int)ptrState->vectorShowSurface);
        _shaders[s]->setUniform("project_vectors", (int)ptrState->vector_projection);
        _shaders[s]->setUniform("style_arrow", (int)ptrState->vectorShape);
      }
    }
  }
#endif 

  // lumiere
  if (isoLight == 1 && ptrState->dim == 3) 
  {
    light(3);
#ifdef __SHADERS__
    if (ptrState->mode == SCALARFIELD)
    {
      int s = _shaders.shader_id(shader::iso_banded_colormap);
      if (ptrState->scalarStyle == 2 || ptrState->scalarStyle == 3) s = _shaders.shader_id(shader::iso_colored_lines);
      _shaders[s]->setUniform("lightOn", (int)1);
    }
    else 
    {
        int s = _shaders.shader_id(shader::vector_rgb);
        if (ptrState->vectorStyle == 1) s = _shaders.shader_id(shader::vector_arrow);
        else if (ptrState->vectorStyle == 2) s = _shaders.shader_id(shader::vector_line);
  	    _shaders[s]->setUniform("lightOn", (int)1);
    }
#endif
  }

  zone = 0;
  while (zone < _numberOfUnstructZones)
  {
    glPolygonOffset(1., zone%10+1.);
    zonet = zone + _numberOfStructZones;
    UnstructZone* zonep = _uzones[zone];
    ZoneImplDL* zoneImpl = static_cast<ZoneImplDL*>(zonep->ptr_impl);

    // if zone is active and in frustum
    if ((zonep->active == 1  ||
         (zonep->active == 0 && ptrState->ghostifyDeactivatedZones == 1)) 
        && isInFrustum(zonep, _view) == 1)
    {
      if (ptrState->simplifyOnDrag == 1 && ptrState->ondrag == 1) displayUBBZone(zonep);
      else
      {
        if (ptrState->mode == SCALARFIELD)
        {
          if (zoneImpl->_DLiso != 0)
            renderUIsoSolidZone(zonep, zonet, ptrState->scalarField);
          else
            displayUIsoSolidZone(zonep, zonet, ptrState->scalarField);
        }
        else // VECTORFIELD
        {
          if (zoneImpl->_DLiso != 0)
            renderUIsoSolidZone(zonep, zone, ptrState->vectorField1, 
                                ptrState->vectorField2, ptrState->vectorField3);
          else
            displayUIsoSolidZone(zonep, zone, ptrState->vectorField1, 
                                ptrState->vectorField2, ptrState->vectorField3);
        }
      }
    }
    zone++;
  }

  noLight();
#ifdef __SHADERS__
  _shaders.activate((short unsigned int)0);
#endif
  glActiveTexture(GL_TEXTURE1);
  glDisable(GL_TEXTURE_2D);
  glActiveTexture(GL_TEXTURE0);
  glDisable(GL_TEXTURE_2D);
  
  // Surimposition maillage eventuellement
  if (ptrState->scalarStyle == 1 || ptrState->scalarStyle == 3)
  {
    glColor4f(0.,0.,0.,1.);
    glPolygonOffset(1., 1.);
    glEnable(GL_LINE_SMOOTH);
    glHint(GL_LINE_SMOOTH_HINT, GL_NICEST);
    glDepthMask(GL_FALSE);
    zone = 0;
    while (zone < _numberOfUnstructZones)
    {
      zonet = zone + _numberOfStructZones;
      UnstructZone* zonep = _uzones[zone];
      // if zone is activated and in frustum
      if ((zonep->active == 1 ||
           (zonep->active == 0 && ptrState->ghostifyDeactivatedZones == 1))
          && isInFrustum(zonep, _view) == 1)
      {
        if (ptrState->simplifyOnDrag == 1 && ptrState->ondrag == 1) displayUBBZone(zonep);
        else displayUMeshZone(zonep, zone, zonet);
      }
      zone++;
    }
    glDepthMask(GL_TRUE);
    glDisable(GL_LINE_SMOOTH);
  } // fin maillage

  glDisable(GL_BLEND);
  glDisable(GL_POLYGON_OFFSET_FILL);
  glDisable(GL_POLYGON_OFFSET_LINE);
  glColor4f(1., 1., 1., 1.);
}

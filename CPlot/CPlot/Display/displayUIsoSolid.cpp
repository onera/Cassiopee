/*    
    Copyright 2013-2018 Onera.

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
#include "../DataDL.h"
#include "../ZoneImplDL.h"

//============================================================================
/*
  Display a unstructured mesh solution using fully colored planes.
*/
//============================================================================
void DataDL::displayUIsoSolid()
{
  int zone, zonet;
  if (_numberOfUnstructZones == 0) return;

  // Blending
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  glEnable(GL_POLYGON_OFFSET_FILL);
  glEnable(GL_POLYGON_OFFSET_LINE);

#ifdef __SHADERS__
  SHADOWTEXTURE;
  if (ptrState->mode == SCALARFIELD)
  { // shader pour les isos scalaires
    glActiveTexture(GL_TEXTURE1);
    if (_texColormap == 0) createColormapTexture();
    fillColormapTexture((int)_pref.colorMap->varName[0]-48);
    glBindTexture(GL_TEXTURE_1D, _texColormap);
    int s = 10;
    if (ptrState->scalarStyle == 2 || ptrState->scalarStyle == 3) s = 29;
    if (_shaders.currentShader() != s) _shaders.activate((short unsigned int)s);
    _shaders[s]->setUniform("colormap", (int)1);
    int nofield = ptrState->scalarField;
    if (_niso[nofield] == -1)
    {
      _shaders[s]->setUniform("alpha", (float)1.);
      _shaders[s]->setUniform("beta", (float)0.);
      _shaders[s]->setUniform("niso", (float)ptrState->niso);
    }
    else 
    {
      float rmin, rmax, alpha, beta;
      float deltai = MAX(maxf[nofield]-minf[nofield], 1.e-6);
      rmin = (_isoMin[nofield] -minf[nofield])/deltai;
      rmax = (_isoMax[nofield] -minf[nofield])/deltai;
      deltai = MAX(rmax-rmin, 1.e-6);
      alpha = 1./deltai; beta = -rmin/deltai;
      _shaders[s]->setUniform("niso", (float)_niso[nofield]);
      _shaders[s]->setUniform("alpha", (float)alpha);
      _shaders[s]->setUniform("beta", (float)beta);
    }
    _shaders[s]->setUniform("edgeStyle", (float)ptrState->isoEdges);
    _shaders[s]->setUniform("lightOn", (int)0);
    _shaders[s]->setUniform("shadow", (int)ptrState->shadow);
    _shaders[s]->setUniform("ShadowMap", (int)0);
  }
  else
  { // shader pour les isos vectoriels
    int s = 27;
    if ( ptrState->vectorStyle == 3) s = 35;
    if ( ptrState->vectorStyle == 2) s = 33;
    if ( ptrState->vectorStyle == 1) s = 34;
    if ( ptrState->vectorStyle == 4) s = 36;
    if (_shaders.currentShader() != s)
      _shaders.activate((short unsigned int)s);
    if ( s == 27 ) {
      _shaders[s]->setUniform("lightOn", (int)0);
      _shaders[s]->setUniform("shadow", (int)ptrState->shadow);
      _shaders[s]->setUniform("ShadowMap", (int)0);
    }
    if ((s == 33)||(s==34)||(s==35)) {
      // ##############################################################
      // Faudra rajouter une option scale pour la taille des vecteurs !
      // ##############################################################
      _shaders[s]->setUniform("lightOn", (int)0);
      _shaders[s]->setUniform("shadow", (int)ptrState->shadow);
      _shaders[s]->setUniform("ShadowMap", (int)0);
      double diag = 0.01*sqrt((xmax-xmin)*(xmax-xmin)+(ymax-ymin)*(ymax-ymin)+(zmax-zmin)*(zmax-zmin));
      double sc = ptrState->vectorScale/100.;
      _shaders[s]->setUniform("scale", float(sc*diag));
      _shaders[s]->setUniform("fix_length", ptrState->vectorNormalize);
      if ((s==34)||(s==35)) _shaders[s]->setUniform("show_surface",ptrState->vectorShowSurface);
    }
    if (s==36) {
      _shaders[s]->setUniform("lightOn", (int)0);
      _shaders[s]->setUniform("shadow", (int)ptrState->shadow);
      _shaders[s]->setUniform("ShadowMap", (int)0);
      double diag = 0.01*sqrt((xmax-xmin)*(xmax-xmin)+(ymax-ymin)*(ymax-ymin)+(zmax-zmin)*(zmax-zmin));
      double sc = ptrState->vectorScale/100.f;
      _shaders[s]->setUniform("scale", float(sc*diag));
      _shaders[s]->setUniform("density", float(ptrState->vectorDensity/diag));
      _shaders[s]->setUniform("fix_length", ptrState->vectorNormalize);
    }

  }
#endif 

  // lumiere
  if (ptrState->isoLight == 1 && ptrState->dim == 3) 
  {
    light(3);
#ifdef __SHADERS__
    if (ptrState->mode == SCALARFIELD)
    {
      int s = 10;
      if (ptrState->scalarStyle == 2 || ptrState->scalarStyle == 3) s = 29;
      _shaders[s]->setUniform("lightOn", (int)1);
    }
    else {
	int s = 27;
	if (ptrState->vectorStyle == 2) s = 33;
	if (ptrState->vectorStyle == 1) s = 34;
	if (ptrState->vectorStyle == 3) s = 35;
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
    zone++;
  }

  noLight();
#ifdef __SHADERS__
  _shaders.activate((short unsigned int)0);
#endif

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
	displayUMeshZone(zonep, zone, zonet);
      }
      zone++;
    }
    glDepthMask(GL_TRUE);
    glDisable(GL_LINE_SMOOTH);
  } // fin maillage

  glDisable(GL_BLEND);
  glDisable(GL_POLYGON_OFFSET_FILL);
  glDisable(GL_POLYGON_OFFSET_LINE);
  glColor3f(1., 1., 1.);
}

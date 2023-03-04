/*    
    Copyright 2013-2023 Onera.

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

# include "kcore.h"
# include "cplot.h"
# include "Data.h"

using namespace K_FLD;
using namespace std;

//=============================================================================
PyObject* K_CPLOT::displayAgainFBO(PyObject* self, PyObject* args)
{
  #include "display1.h"

  // Construction de la chaine de toutes les variables
  E_Int referenceNfield;
  char** referenceVarNames;
  d->getAllVars(structVarString, unstrVarString,
                referenceNfield, referenceVarNames);

  // Init et remplace
  d->initZoneData(structF, structVarString, nit, njt, nkt,
                  unstrF, unstrVarString, cnt, eltType, 
                  zoneNames, renderTags,
                  referenceNfield, referenceVarNames);

  for (size_t i = 0; i < zoneNames.size(); i++) delete [] zoneNames[i];
  for (size_t i = 0; i < renderTags.size(); i++) delete [] renderTags[i];

  for (E_Int i = 0; i < referenceNfield; i++) delete [] referenceVarNames[i];
  delete [] referenceVarNames;
  d->ptrState->clearDeactivatedZones();
  
  // enforce given data
  E_Int mode = getMode(modeObject);
  E_Int scalarField = getScalarField(scalarFieldObject);
  E_Int vectorField1 = getScalarField(vectorFieldObject1);
  E_Int vectorField2 = getScalarField(vectorFieldObject2);
  E_Int vectorField3 = getScalarField(vectorFieldObject3);
  d->enforceGivenData(dim, mode, scalarField, vectorField1, vectorField2,
                      vectorField3, displayBB, displayInfo, displayIsoLegend);
  d->autoPlugins();
  d->enforceGivenData2(xcam, ycam, zcam,
                       xeye, yeye, zeye,
                       dirx, diry, dirz, viewAngle, 
                       meshStyle, solidStyle, scalarStyle, 
                       vectorStyle, vectorScale, vectorDensity, vectorNormalize, vectorShowSurface,
                       vectorShape, vectorProjection,
                       colormap, colormapC1, colormapC2, colormapC3, colormapC,
                       niso, isoEdges, isoScales,
                       bgColor, backgroundFile,
                       -1, -1, -1, shadow, dof,
                       exportFile, exportResolution);

  if (stereo != -1) d->ptrState->stereo = stereo;
  if (stereoDist != -1.) d->ptrState->stereoDist = stereoDist;

  // offscreen rendering?
  if (offscreen > 0) { d->ptrState->offscreen = offscreen; d->ptrState->shootScreen = 1; }
  if (frameBuffer >= 0 && frameBuffer < 10) d->ptrState->frameBuffer = frameBuffer;

  // Free the input arrays
  E_Int structFSize = structF.size();
  for (E_Int i = 0; i < structFSize; i++) RELEASESHAREDS(objs[i], structF[i]);

  E_Int unstrFSize = unstrF.size();
  for (E_Int i = 0; i < unstrFSize; i++) RELEASESHAREDU(obju[i], unstrF[i], cnt[i]);

    d->ptrState->farClip = 1;
    d->ptrState->render = 1;

  // Retourne le hook
  return Py_BuildValue("l", d);
}

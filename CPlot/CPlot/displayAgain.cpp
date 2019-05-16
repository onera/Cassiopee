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

# include "kcore.h"
# include "cplot.h"
# include "Data.h"

using namespace K_FLD;
using namespace std;

//=============================================================================
/* display arrays (again) */
//=============================================================================
PyObject* K_CPLOT::displayAgain(PyObject* self, PyObject* args)
{
  PyObject* arrays;
  int dim;
  PyObject* modeObject;
  PyObject* scalarFieldObject;
  PyObject* vectorFieldObject1, *vectorFieldObject2, *vectorFieldObject3;
  int displayBB, displayInfo, displayIsoLegend;
  int winx, winy;
  int meshStyle, solidStyle, scalarStyle, vectorStyle, colormap, niso;
  E_Float xcam, ycam, zcam, xeye, yeye, zeye, dirx, diry, dirz, isoEdges;
  E_Float viewAngle, stereoDist, vectorScale, vectorDensity;
  int vectorNormalize, vectorShowSurface, vectorShape, vectorProjection;
  char* exportFile; char* exportResolution;
  PyObject* zoneNamesObject;
  PyObject* renderTagsObject;
  PyObject* isoScales;
  int bgColor, shadow, dof, offscreen, stereo;
  if (!PyArg_ParseTuple(args, 
			"OiOOOOOiiiiiiiddiiiiiidO(ii)(ddd)(ddd)(ddd)diiiidssOOi",
                        &arrays, &dim, &modeObject, &scalarFieldObject,
                        &vectorFieldObject1, &vectorFieldObject2, &vectorFieldObject3,
                        &displayBB, &displayInfo, &displayIsoLegend,
                        &meshStyle, &solidStyle, &scalarStyle, 
                        &vectorStyle, &vectorScale, &vectorDensity, &vectorNormalize, 
                        &vectorShowSurface, &vectorShape, &vectorProjection, &colormap,
                        &niso, &isoEdges, &isoScales,
                        &winx, &winy, &xcam, &ycam, &zcam,
                        &xeye, &yeye, &zeye,
                        &dirx, &diry, &dirz, &viewAngle,
                        &bgColor, &shadow, &dof, &stereo, &stereoDist, 
                        &exportFile, &exportResolution,
                        &zoneNamesObject, &renderTagsObject,
                        &offscreen))
  {
    return NULL;
  }

  // Recuperation des noms de zones (eventuellement)
  vector<char*> zoneNames;
  getStringsFromPyObj(zoneNamesObject, zoneNames);

  // Recuperation des tags de render (eventuellement)
  vector<char*> renderTags;
  getStringsFromPyObj(renderTagsObject, renderTags);

  // Recuperation du container de donnees
  Data* d = Data::getInstance();
  
  // Lecture des arrays
  vector<E_Int> res;
  vector<char*> structVarString;
  vector<char*> unstrVarString;
  vector<FldArrayF*> structF;  vector<FldArrayF*> unstrF;
  vector<E_Int> nit; vector<E_Int> njt; vector<E_Int> nkt;
  vector<FldArrayI*> cnt;
  vector<char*> eltType;
  vector<PyObject*> objs, obju;
  E_Boolean skipNoCoord = true;
  E_Boolean skipStructured = false;
  E_Boolean skipUnstructured = false;
  E_Boolean skipDiffVars = false;

  E_Int isOk = K_ARRAY::getFromArrays(
    arrays, res, structVarString, unstrVarString,
    structF, unstrF, nit, njt, nkt, cnt, eltType, objs, obju, 
    skipDiffVars, skipNoCoord, skipStructured, skipUnstructured, true);

  if (isOk == -1)
  {
    PyErr_SetString(PyExc_TypeError,
                    "display: invalid list of arrays.");
    E_Int structFSize = structF.size();
    for (E_Int i = 0; i < structFSize; i++) 
      RELEASESHAREDS(objs[i], structF[i]);

    E_Int unstrFSize = unstrF.size();
    for (E_Int i = 0; i < unstrFSize; i++)
      RELEASESHAREDU(obju[i], unstrF[i], cnt[i]);
    
    return NULL;
  }

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

  for (E_Int i = 0; i < referenceNfield; i++) delete [] referenceVarNames[i];
  delete [] referenceVarNames;

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
                       colormap,
                       niso, isoEdges, isoScales,
                       bgColor, -1, -1, -1, shadow, dof,
                       exportFile, exportResolution);

  if (stereo != -1) d->ptrState->stereo = stereo;
  if (stereoDist != -1.) d->ptrState->stereoDist = stereoDist;

  // offscreen rendering?
  if (offscreen > 0) d->ptrState->offscreen = offscreen;

  if (d->ptrState->offscreen == 1) // MESA offscreen
  {
#ifdef __MESA__
  //printf("Creating OS context...");
  OSMesaContext ctx; 
  ctx = OSMesaCreateContext(OSMESA_RGBA, NULL);
  d->ptrState->offscreenBuffer = (char*)malloc(d->_view.w * d->_view.h * 4 * 
                                               sizeof(GLubyte));
  OSMesaMakeCurrent(ctx, d->ptrState->offscreenBuffer, GL_UNSIGNED_BYTE, 
                    d->_view.w, d->_view.h);

  glShadeModel(GL_SMOOTH);
  glEnable(GL_DEPTH_TEST);
  d->setBgColor();
  d->ptrState->farClip = 1;
  d->ptrState->render = 0;
  d->display();
  d->exportFile();
  //printf("done.\n");
  free(d->ptrState->offscreenBuffer);
  OSMesaDestroyContext(ctx);
#else
  printf("Error: CPlot: MESA offscreen unavailable.\n");
#endif
  }
  else
  { // Direct ou FBO offscreen
    d->ptrState->farClip = 1;
    d->ptrState->render = 1;
  }

  // Free the input arrays
  E_Int structFSize = structF.size();
  for (E_Int i = 0; i < structFSize; i++) 
    RELEASESHAREDS(objs[i], structF[i]);

  E_Int unstrFSize = unstrF.size();
  for (E_Int i = 0; i < unstrFSize; i++)
    RELEASESHAREDU(obju[i], unstrF[i], cnt[i]);

  // Retourne le hook
  return Py_BuildValue("l", d);
}

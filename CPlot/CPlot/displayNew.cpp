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

#ifndef _WIN32
#include <X11/Xlib.h>
#define INITTHREADS XInitThreads()
#else
#define INITTHREADS
#endif

using namespace K_FLD;
using namespace std;

//=============================================================================
// Cree la boucle glut dans une thread
//=============================================================================
static void* threadFunc(void* v)
{
  Data* d = Data::getInstance();
  /* Gfx setup */
  int argc = 0;
  char* com = NULL;
  glutInit(&argc, &com);
  d->openGfx();
  glutMainLoop();
  return NULL;
}

//=============================================================================
/* display arrays */
//=============================================================================
PyObject* K_CPLOT::displayNew(PyObject* self, PyObject* args)
{
  PyObject* arrays;
  int dim;
  PyObject* modeObject;
  PyObject* scalarFieldObject;
  PyObject* vectorFieldObject1, *vectorFieldObject2, *vectorFieldObject3;
  int winx, winy;
  int displayBB, displayInfo, displayIsoLegend;
  int meshStyle, solidStyle, scalarStyle, vectorStyle, colormap, niso;
  E_Float xcam, ycam, zcam, xeye, yeye, zeye, dirx, diry, dirz, isoEdges;
  E_Float stereoDist, viewAngle, vectorScale, vectorDensity;
  int vectorNormalize, vectorShowSurface, vectorShape, vectorProjection;
  char* exportFile; char* exportResolution;
  PyObject* zoneNamesObject;
  PyObject* renderTagsObject;
  PyObject* isoScales;
  int bgColor, shadow, dof, offscreen, stereo;
  if (!PyArg_ParseTuple(args, "OiOOOOOiiiiiiiddiiiiiidO(ii)(ddd)(ddd)(ddd)diiiidssOOi",
                        &arrays, &dim, &modeObject, &scalarFieldObject,
                        &vectorFieldObject1, &vectorFieldObject2, &vectorFieldObject3,
                        &displayBB, &displayInfo, &displayIsoLegend,
                        &meshStyle, &solidStyle, &scalarStyle,
                        &vectorStyle, &vectorScale, &vectorDensity, &vectorNormalize, 
                        &vectorShowSurface, &vectorShape, &vectorProjection, &colormap,
                        &niso, &isoEdges, &isoScales,
                        &winx, &winy, &xcam, &ycam, &zcam,
                        &xeye, &yeye, &zeye,
                        &dirx, &diry, &dirz, &viewAngle, &bgColor,
                        &shadow, &dof, &stereo, &stereoDist,
                        &exportFile, &exportResolution, 
                        &zoneNamesObject, &renderTagsObject, &offscreen))
  {
    return NULL;
  }

  // Recuperation des noms de zones (eventuellement)
  vector<char*> zoneNames;
  getStringsFromPyObj(zoneNamesObject, zoneNames);

  // Recuperation des tags de render (eventuellement)
  vector<char*> renderTags;
  getStringsFromPyObj(renderTagsObject, renderTags);


  // Lecture des arrays
  vector<E_Int> res;
  vector<char*> structVarString; vector<char*> unstrVarString;
  vector<FldArrayF*> structF; vector<FldArrayF*> unstrF;
  vector<E_Int> nit; vector<E_Int> njt; vector<E_Int> nkt;
  vector<FldArrayI*> cnt;
  vector<char*> eltType;
  vector<PyObject*> objs, obju;
  E_Boolean skipNoCoord = true;
  E_Boolean skipStructured = false;
  E_Boolean skipUnstructured = false;
  E_Boolean skipDiffVars = false;
  E_Int isOk = K_ARRAY::getFromArrays(arrays, res, structVarString, unstrVarString,
                                      structF, unstrF, nit, njt, nkt, cnt,
                                      eltType, objs, obju, 
                                      skipDiffVars, skipNoCoord, skipStructured,
                                      skipUnstructured, true);

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

  Data* d = Data::getInstance();

  // Construction de la chaine de toutes les variables
  E_Int referenceNfield;
  char** referenceVarNames;
  d->getAllVars(structVarString, unstrVarString,
                referenceNfield, referenceVarNames);

  d->_CDisplayIsLaunched = 1;

  d->initZoneData(structF, structVarString, nit, njt, nkt,
		  unstrF, unstrVarString, cnt, eltType, 
		  zoneNames, renderTags,
      referenceNfield, referenceVarNames);

  for (size_t i = 0; i < zoneNames.size(); i++) delete [] zoneNames[i];

  for (E_Int i = 0; i < referenceNfield; i++) delete [] referenceVarNames[i];
  delete [] referenceVarNames;

  // Initialisation des Data restantes
  E_Int mode = getMode(modeObject);
  E_Int scalarField = getScalarField(scalarFieldObject);
  E_Int vectorField1 = getScalarField(vectorFieldObject1);
  E_Int vectorField2 = getScalarField(vectorFieldObject2);
  E_Int vectorField3 = getScalarField(vectorFieldObject3);
  d->enforceGivenData(dim, mode, scalarField, vectorField1, vectorField2,
                      vectorField3, displayBB, displayInfo, displayIsoLegend);
  d->initCam();
  d->loadPlugins();
  d->loadPrefs();
  d->autoPlugins();

  // Enforce given data
  d->enforceGivenData2(xcam, ycam, zcam,
                       xeye, yeye, zeye,
                       dirx, diry, dirz, viewAngle,
                       meshStyle, solidStyle, scalarStyle, 
                       vectorStyle, vectorScale, vectorDensity, vectorNormalize, vectorShowSurface,
                       vectorShape, vectorProjection, colormap, 
                       niso, isoEdges, isoScales, 
                       bgColor, -1, -1, -1, shadow, dof,
                       exportFile, exportResolution);

  if (stereo != -1) d->ptrState->stereo = stereo;
  if (stereoDist != -1.) d->ptrState->stereoDist = stereoDist;

  // offscreen rendering?
  if (offscreen > 0) d->ptrState->offscreen = offscreen;

  // Assure la taille de la fenetre
  if (winx != -1) d->_view.w = winx;
  if (winy != -1) d->_view.h = winy;
  d->ptrState->render = 1;

  // Free the input arrays
  E_Int structFSize = structF.size();
  for (E_Int i = 0; i < structFSize; i++) RELEASESHAREDS(objs[i], structF[i]);
  
  E_Int unstrFSize = unstrF.size();
  for (E_Int i = 0; i < unstrFSize; i++) 
    RELEASESHAREDU(obju[i], unstrF[i], cnt[i]);

  if (d->ptrState->offscreen == 1) // MESA offscreen
  {
    int argc = 0;
    char* com = NULL;
    INITTHREADS;
    glutInit(&argc, &com);
#ifdef __MESA__
    /* Init */
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    
    // Window size base sur la taille de l'ecran
    int screenWidth = glutGet(GLUT_SCREEN_WIDTH);
    int screenHeight = glutGet(GLUT_SCREEN_HEIGHT);
    d->_view.w = screenWidth-320;
    d->_view.h = screenHeight;
 
    //printf("Creating OS context..."); fflush(stdout);
    OSMesaContext ctx; 
    ctx = OSMesaCreateContext(OSMESA_RGBA, NULL);
    d->ptrState->offscreenBuffer = (char*)malloc(d->_view.w * d->_view.h * 4 * 
                                                 sizeof(GLubyte));
    OSMesaMakeCurrent(ctx, d->ptrState->offscreenBuffer, GL_UNSIGNED_BYTE, 
                      d->_view.w, d->_view.h);
    d->init();
    d->ptrState->farClip = 1;
    d->ptrState->render = 0;
    d->display();  
    d->exportFile();
    //printf("done.\n");
    free(d->ptrState->offscreenBuffer);
    OSMesaDestroyContext(ctx);
#else
    printf("Error: CPlot: mesa offscreen unavailable.\n");
#endif
  }// if (d->ptrState->offscreen == 1)
  else
  { // direct ou offscreen FBO
    // thread en python
    Py_BEGIN_ALLOW_THREADS;
    Data* d = Data::getInstance();
    d->_save = _save;
    /* Gfx setup */
    int argc = 0;
    char* com = NULL;
    INITTHREADS;
    glutInit(&argc, &com);
    d->openGfx();
    glutMainLoop();
    Py_END_ALLOW_THREADS;
  }

  // Retourne le hook
  return Py_BuildValue("l", d);
}

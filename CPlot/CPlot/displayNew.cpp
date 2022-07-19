/*    
    Copyright 2013-2022 Onera.

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
// Cree la boucle glut dans une thread. Cette fonction n'est plus
// utilise car le threading est fait en python
//=============================================================================
/*
static void* threadFunc(void* v)
{
  Data* d = Data::getInstance();
  int argc = 0;
  char* com = NULL;
  glutInit(&argc, &com);
  d->openGfx();
  glutMainLoop();
  return NULL;
}
*/
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
  char* colormapC1; char* colormapC2; char* colormapC3; PyObject* colormapC;
  E_Float xcam, ycam, zcam, xeye, yeye, zeye, dirx, diry, dirz, isoEdges;
  E_Float stereoDist, viewAngle, vectorScale, vectorDensity;
  int vectorNormalize, vectorShowSurface, vectorShape, vectorProjection;
  char* exportFile; char* exportResolution;
  PyObject* zoneNamesObject;
  PyObject* renderTagsObject;
  PyObject* isoScales;
  char* backgroundFile;
  int bgColor, shadow, dof, offscreen, stereo, frameBuffer;
  if (!PyArg_ParseTuple(args, "OiOOOOOiiiiiiiddiiiiisssOidO(ii)(ddd)(ddd)(ddd)disiiidssOOii",
                        &arrays, &dim, &modeObject, &scalarFieldObject,
                        &vectorFieldObject1, &vectorFieldObject2, &vectorFieldObject3,
                        &displayBB, &displayInfo, &displayIsoLegend,
                        &meshStyle, &solidStyle, &scalarStyle,
                        &vectorStyle, &vectorScale, &vectorDensity, &vectorNormalize, 
                        &vectorShowSurface, &vectorShape, &vectorProjection, 
                        &colormap, &colormapC1, &colormapC2, &colormapC3, &colormapC,
                        &niso, &isoEdges, &isoScales,
                        &winx, &winy, &xcam, &ycam, &zcam,
                        &xeye, &yeye, &zeye,
                        &dirx, &diry, &dirz, &viewAngle, 
                        &bgColor, &backgroundFile,
                        &shadow, &dof, &stereo, &stereoDist,
                        &exportFile, &exportResolution, 
                        &zoneNamesObject, &renderTagsObject, &frameBuffer, &offscreen))
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

  // Read from arrays using Array2/3
  /*
  char* varStringl; char* eltTypel;
  E_Int nil, njl, nkl, resl;
  FldArrayF* fl; FldArrayI* cnl;
  PyObject* o;
  E_Int nz = PyList_Size(arrays);
  for (E_Int i = 0; i < nz; i++)
  {
    o = PyList_GetItem(arrays, i);
    resl = K_ARRAY::getFromArray3(o, varStringl, fl, 
                                  nil, njl, nkl, cnl, eltTypel);
    if (resl == 1)
    {
        structVarString.push_back(varStringl);
        structF.push_back(fl);
        nit.push_back(nil); njt.push_back(njl), nkt.push_back(nkl);
        objs.push_back(o);
    }
    else if (resl == 2)
    {
        unstrVarString.push_back(varStringl);
        unstrF.push_back(fl);
        eltType.push_back(eltTypel); cnt.push_back(cnl);
        obju.push_back(o);
    }
    else printf("Warning: display: array %d is invalid.\n", i);
  } */

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

  // Assure la taille de la fenetre
  if (winx != -1) d->_view.w = winx;
  if (winy != -1) d->_view.h = winy;
  if (offscreen > 0)
  {
    if (d->ptrState->exportWidth == -1) d->ptrState->exportWidth = 1920;
    if (d->ptrState->exportHeight == -1) d->ptrState->exportHeight = 1080;
    d->_view.w = d->ptrState->exportWidth; d->_view.h = d->ptrState->exportHeight;
  }
  d->ptrState->render = 1;

  // Free the input arrays
  E_Int structFSize = structF.size();
  for (E_Int i = 0; i < structFSize; i++) RELEASESHAREDS(objs[i], structF[i]);

  E_Int unstrFSize = unstrF.size();
  for (E_Int i = 0; i < unstrFSize; i++) RELEASESHAREDU(obju[i], unstrF[i], cnt[i]);

  if (d->ptrState->offscreen == 1 ||
      d->ptrState->offscreen == 5 ||
      d->ptrState->offscreen == 6 ||
      d->ptrState->offscreen == 7) // MESA offscreen
  {
    // Dans ce cas, on ne fait pas de glutInit, car il requiert
    // un serveur X
#ifdef __MESA__
    /* Init */

    // Window size base sur l'export
    if (d->ptrState->exportWidth == -1) d->ptrState->exportWidth = 1920;
    if (d->ptrState->exportHeight == -1) d->ptrState->exportHeight = 1080;
    d->_view.w = d->ptrState->exportWidth; d->_view.h = d->ptrState->exportHeight;
    //printf("%d %d\n", d->ptrState->exportWidth, d->ptrState->exportHeight);
    
    //printf("Creating OS context..."); fflush(stdout);
    OSMesaContext* ctx = new OSMesaContext();
    //ctx = OSMesaCreateContext(OSMESA_RGBA, NULL);
    (*ctx) = OSMesaCreateContextExt(OSMESA_RGBA, 32, 0, 0, NULL);
    d->ptrState->ctx = ctx;
    d->ptrState->offscreenBuffer[d->ptrState->frameBuffer] = 
    (char*)malloc(d->_view.w * d->_view.h * 4 * sizeof(GLubyte));
    OSMesaMakeCurrent(*ctx, d->ptrState->offscreenBuffer[d->ptrState->frameBuffer], 
                      GL_UNSIGNED_BYTE, d->_view.w, d->_view.h);
    d->init();
    d->ptrState->farClip = 1;
    d->ptrState->render = 0;
    d->ptrState->shootScreen = 0;
    gdisplay(); // build DL
    d->display();
    d->exportFile();
    //printf("done.\n");
    // use finalizeExport to free OSMesaContext
#else
    printf("Error: CPlot: mesa offscreen unavailable.\n");
#endif
  }
  else
  { // direct ou offscreen FBO
    d->ptrState->farClip = 1;
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

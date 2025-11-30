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
/* display arrays (generic) */
//=============================================================================
PyObject* K_CPLOT::displayNew(PyObject* self, PyObject* args)
{
  #include "display1.h"
    
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
  for (size_t i = 0; i < renderTags.size(); i++) delete [] renderTags[i];

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
                       exportFile, exportResolution, exportAA);

  if (stereo != -1) d->ptrState->stereo = stereo;
  if (stereoDist != -1.) d->ptrState->stereoDist = stereoDist;
  if (panorama != -1) d->ptrState->panorama = panorama;
  if (dofPower != -1.) d->ptrState->dofPower = dofPower;
  if (gamma != -1.) d->ptrState->gamma = gamma;
  if (toneMapping != -1) d->ptrState->toneMapping = toneMapping;
  if (lightOffsetX != -999.) d->ptrState->lightOffsetX = lightOffsetX;
  if (lightOffsetY != -999.) d->ptrState->lightOffsetY = lightOffsetY;

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
    
    // print current renderer
    const char* renderer = (const char*)glGetString(GL_RENDERER);
    printf("Renderer: %s\n", renderer);

    if (posCamList == Py_None)
    {
      //if (d->ptrState->stereo == 0) d->display();
      //else d->displayAnaglyph();
      d->exportFile(); // perfoms display
      // use finalizeExport to free OSMesaContext
    }
    else // list of posCams: ODS
    {
      d->ptrState->odsRun = true;
      E_Int nslits = PyList_Size(posCamList)/9; // xyz+front/top/bot
      d->ptrState->odsNSlits = nslits;
      E_Int height = d->ptrState->exportHeight;
      d->ptrState->odsImage = new char [nslits*height*3];
      d->ptrState->odsFrontImage = new char [2*height*3];
      d->ptrState->odsTopImage = new char [2*height*3];
      
      for (E_Int i = 0; i < 3*nslits; i++)
      {
        //printf("%d / %d\n", i, 3*nslits);
        d->ptrState->odsSlit = i;
        PyObject* v = PyList_GetItem(posCamList, 3*i); 
        d->_view.xcam = PyFloat_AsDouble(v);
        v = PyList_GetItem(posCamList, 3*i+1);
        d->_view.ycam = PyFloat_AsDouble(v);
        v = PyList_GetItem(posCamList, 3*i+2); 
        d->_view.zcam = PyFloat_AsDouble(v);
        v = PyList_GetItem(posEyeList, 3*i); 
        d->_view.xeye = PyFloat_AsDouble(v);
        v = PyList_GetItem(posEyeList, 3*i+1); 
        d->_view.yeye = PyFloat_AsDouble(v);
        v = PyList_GetItem(posEyeList, 3*i+2); 
        d->_view.zeye = PyFloat_AsDouble(v);
        v = PyList_GetItem(dirCamList, 3*i); 
        d->_view.dirx = PyFloat_AsDouble(v);
        v = PyList_GetItem(dirCamList, 3*i+1); 
        d->_view.diry = PyFloat_AsDouble(v);
        v = PyList_GetItem(dirCamList, 3*i+2); 
        d->_view.dirz = PyFloat_AsDouble(v);
        //d->display(); // done in export
        d->exportFile(); // performs display
      }
      delete [] d->ptrState->odsImage;
      delete [] d->ptrState->odsFrontImage;
      delete [] d->ptrState->odsTopImage;
    }
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

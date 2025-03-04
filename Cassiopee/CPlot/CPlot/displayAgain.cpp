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

using namespace K_FLD;
using namespace std;

//=============================================================================
/* display arrays (again) */
//=============================================================================
PyObject* K_CPLOT::displayAgain(PyObject* self, PyObject* args)
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
  if (panorama != -1) d->ptrState->panorama = panorama;

  // offscreen rendering?
  if (offscreen > 0) { d->ptrState->offscreen = offscreen; d->ptrState->shootScreen = 1; }
  if (frameBuffer >= 0 && frameBuffer < 10) d->ptrState->frameBuffer = frameBuffer;

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
#ifdef __MESA__
  // Free context if resolution change
  if ((d->ptrState->exportWidth != -1 && d->_view.w != d->ptrState->exportWidth) || 
      (d->ptrState->exportHeight != -1 && d->_view.h != d->ptrState->exportHeight))  
  {
    d->_view.w = d->ptrState->exportWidth; d->_view.h = d->ptrState->exportHeight;
    free(d->ptrState->offscreenBuffer[d->ptrState->frameBuffer]);
    d->ptrState->offscreenBuffer[d->ptrState->frameBuffer] = NULL;
    OSMesaDestroyContext(*(OSMesaContext*)(d->ptrState->ctx));
    d->ptrState->ctx = NULL;
  }

  if (d->ptrState->ctx == NULL)
  {
    //printf("recreating context\n");
    OSMesaContext* ctx = new OSMesaContext();
    (*ctx) = OSMesaCreateContextExt(OSMESA_RGBA, 32, 0, 0, NULL);
    d->ptrState->ctx = ctx;

    if (d->ptrState->offscreenBuffer[d->ptrState->frameBuffer] == NULL)
      d->ptrState->offscreenBuffer[d->ptrState->frameBuffer] = 
        (char*)malloc(d->_view.w * d->_view.h * 4 * sizeof(GLubyte));
    OSMesaMakeCurrent(*ctx, d->ptrState->offscreenBuffer[d->ptrState->frameBuffer], 
                      GL_UNSIGNED_BYTE, d->_view.w, d->_view.h);
    //d->createNodeTexture();
    //d->createNoise3DTexture();
    //d->createFrameBufferTexture();
    //d->createPngTexture("windtunnel.png", _texEnviron1, width, height);
    //d->createVoxelTexture();
    d->_texColormap = 0; // textures may be lost when destroying context
    d->setBgColor();
    glShadeModel(GL_SMOOTH);
    glEnable(GL_DEPTH_TEST);
    d->_shaders.init(); // shader are attached to context
    d->_shaders.load();
  }
  else
  {
    OSMesaContext& ctx = *((OSMesaContext*)(d->ptrState->ctx));
    if (d->ptrState->offscreenBuffer[d->ptrState->frameBuffer] == NULL)
      d->ptrState->offscreenBuffer[d->ptrState->frameBuffer] = 
      (char*)malloc(d->_view.w * d->_view.h * 4 * sizeof(GLubyte));
    OSMesaMakeCurrent(ctx, d->ptrState->offscreenBuffer[d->ptrState->frameBuffer], 
                      GL_UNSIGNED_BYTE, d->_view.w, d->_view.h);
  }

  d->ptrState->farClip = 1;
  d->ptrState->render = 0; // 1 ou pas?
  d->ptrState->shootScreen = 0;
  gdisplay(); // build DL

  if (posCamList == Py_None)
  {
    //if (d->ptrState->stereo == 0) d->display();
    //else d->displayAnaglyph();
    d->exportFile(); // performs display
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
      //d->display(); // done in export file
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
  { // Direct ou FBO offscreen
    d->ptrState->farClip = 1;
    d->ptrState->render = 1;
  }

  // Retourne le hook
  return Py_BuildValue("l", d);
}

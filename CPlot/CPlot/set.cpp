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
#include "cplot.h"
#include "Data.h"
#include <time.h>
#if defined(_WIN32) || defined(_WIN64)
#  include <winsock.h>
#endif

int findFace(double xp, double yp, double zp, int elt, 
             UnstructZone* zone, double& dist);

//=============================================================================
PyObject* K_CPLOT::setState(PyObject* self, PyObject* args)
{
  int dim;
  PyObject* modeObject;
  PyObject* scalarFieldObject;
  int vectorField1, vectorField2, vectorField3;
  int winx, winy;
  int displayBB, displayInfo, displayIsoLegend;
  int meshStyle, solidStyle, scalarStyle, vectorStyle, colormap, niso;
  E_Float xcam, ycam, zcam, xeye, yeye, zeye, viewAngle, dirx, diry, dirz;
  E_Float isoEdges, vectorScale, vectorDensity;
  int vectorNormalize, vectorShowSurface, vectorShape, vectorProjection;
  int bgColor, shadow, dof;
  int ghostifyDeactivatedZones;
  int edgifyActivatedZones;
  int edgifyDeactivatedZones;
  int stereo; E_Float stereoDist;
  int cursor;
  char* exportFile; char* exportResolution;
  char* envmap; char* message;
  PyObject* isoScales; PyObject* billBoards; 
  PyObject* materials; PyObject* bumpMaps;
  int gridSizeI, gridSizeJ;
  E_Float lightOffsetX, lightOffsetY; 
  E_Float dofPower; E_Float gamma; E_Float sobelThreshold;
  int timer; int selectionStyle;
  int continuousExport; int activateShortCuts;
  E_Float billBoardSize;
  if (!PyArg_ParseTuple(
        args, 
	"iOOiiiiiiiiiiddiiiiiidO(ii)(ddd)(ddd)(ddd)d(dd)iiidddiiississidi(ii)iiiOdOO",
        &dim, &modeObject, &scalarFieldObject,
        &vectorField1, &vectorField2, &vectorField3,
        &displayBB, &displayInfo, &displayIsoLegend,
        &meshStyle, &solidStyle, &scalarStyle, 
        &vectorStyle, &vectorScale, &vectorDensity, &vectorNormalize, 
        &vectorShowSurface, &vectorShape, &vectorProjection, &colormap,
        &niso, &isoEdges, &isoScales,
        &winx, &winy, &xcam, &ycam, &zcam,
        &xeye, &yeye, &zeye, 
        &dirx, &diry, &dirz, &viewAngle,
        &lightOffsetX, &lightOffsetY,
        &bgColor, &shadow, &dof, &dofPower, &gamma, &sobelThreshold,
        &ghostifyDeactivatedZones, &edgifyActivatedZones,
        &edgifyDeactivatedZones,
        &exportFile, &exportResolution, &continuousExport,
        &envmap, &message,
        &stereo, &stereoDist, &cursor,
        &gridSizeI, &gridSizeJ, &timer, &selectionStyle,
        &activateShortCuts, &billBoards, &billBoardSize, 
        &materials, &bumpMaps))
  {
    return NULL;
  }
  Data* d = Data::getInstance();
  E_Int mode = getMode(modeObject);
  E_Int scalarField = getScalarField(scalarFieldObject);
  d->enforceGivenData(dim, mode, scalarField, vectorField1, vectorField2,
                      vectorField3, displayBB, displayInfo, displayIsoLegend);
  d->enforceGivenData2(xcam, ycam, zcam,
                       xeye, yeye, zeye,
                       dirx, diry, dirz, viewAngle,
                       meshStyle, solidStyle, scalarStyle, 
                       vectorStyle, vectorScale, vectorDensity, vectorNormalize,
                       vectorShowSurface, vectorShape, vectorProjection, colormap, 
                       niso, isoEdges, isoScales, bgColor,
                       ghostifyDeactivatedZones, edgifyActivatedZones,
                       edgifyDeactivatedZones,
                       shadow, dof, exportFile, exportResolution);
  if (stereo != -1) d->ptrState->stereo = stereo;
  if (stereoDist != -1.) d->ptrState->stereoDist = stereoDist;
  if (K_STRING::cmp(envmap, "None") != 0)
  {
    strcpy(d->ptrState->envmapFile, envmap);
    d->ptrState->updateEnvmap = 1;
  }
  if (bgColor >= 6) // requires a background texture
  {
    d->ptrState->updateBackground = 1;
  }
  if (K_STRING::cmp(message, "None") == 0) 
  {
    // nothing to do
  }
  else if (K_STRING::cmp(message, "Clear") == 0) 
  {
    delete [] d->ptrState->message;
    d->ptrState->message = NULL;
  }
  else
  {
    delete [] d->ptrState->message;
    d->ptrState->message = new char [strlen(message)+1];
    strcpy(d->ptrState->message, message);
  }

  if (cursor != -1) { int type = int(cursor); d->ptrState->updateCursor = type; }

  if (continuousExport != -1) d->ptrState->continuousExport = continuousExport;
  if (gridSizeI != -1) d->ptrState->gridSizeI = gridSizeI;
  if (gridSizeJ != -1) d->ptrState->gridSizeJ = gridSizeJ;
  if (dofPower != -1.) d->ptrState->dofPower = dofPower;
  if (gamma != -1.) d->ptrState->gamma = gamma;
  if (lightOffsetX != -999.) d->ptrState->lightOffsetX = lightOffsetX;
  if (lightOffsetY != -999.) d->ptrState->lightOffsetY = lightOffsetY;
  //if (viewAngle != -1) d->ptrState->farClip = 1;

  if (timer != -1) d->ptrState->ktimer = timer;
  if (selectionStyle != -1) d->ptrState->selectionStyle = selectionStyle;

  if (activateShortCuts != -1) d->ptrState->kkeysActivated = activateShortCuts;

  // Getting billBoards (files, ni, nj)
  if (billBoards != Py_None) 
  {
    // delete billBoardStorage
    for (int i = 0; i < d->_nBillBoards; i++)
    {
      delete [] d->_billBoardFiles[i];
      //if (d->_billBoardTexs[i] != 0) glDeleteTextures(1, &d->_billBoardTexs[i]);
    }
    delete [] d->_billBoardTexs;
    delete [] d->_billBoardNis; delete [] d->_billBoardNjs;
    delete [] d->_billBoardWidths; delete [] d->_billBoardHeights;
    delete [] d->_billBoardFiles;
    int nb = PyList_Size(billBoards)/3;
    d->_nBillBoards = nb;
    d->_billBoardFiles = new char* [nb];
    d->_billBoardNis = new int [nb];
    d->_billBoardNjs = new int [nb];
    d->_billBoardWidths = new int [nb];
    d->_billBoardHeights = new int [nb];
    d->_billBoardTexs = new GLuint [nb];
    for (int i = 0; i < nb; i++)
    {
      PyObject* o = PyList_GetItem(billBoards, 3*i);
      char* file = NULL;
      if (PyString_Check(o)) file = PyString_AsString(o);
#if PY_VERSION_HEX >= 0x03000000
      else if (PyUnicode_Check(o)) file = PyBytes_AsString(PyUnicode_AsUTF8String(o)); 
#endif
      o = PyList_GetItem(billBoards, 3*i+1);
      int ni = PyLong_AsLong(o);
      o = PyList_GetItem(billBoards, 3*i+2);
      int nj = PyLong_AsLong(o);
      d->_billBoardFiles[i] = new char [128];  
      strcpy(d->_billBoardFiles[i], file);
      //printf("%s %d %d\n", file, ni,nj);
      d->_billBoardTexs[i] = 0;
      d->_billBoardNis[i] = ni; d->_billBoardNjs[i] = nj;
    }
  }
  if (billBoardSize != -1)
  {
    d->ptrState->billBoardSize = billBoardSize;
  }

  if (materials != Py_None)
  {
    for (int i = 0; i < d->_nMaterials; i++)
    {
      delete [] d->_materialFiles[i];
      // DeleteTex must be done by gfx thread
      //if (d->_materialTexs[i] != 0) glDeleteTextures(1, &d->_materialTexs[i]);
    }
    
    delete [] d->_materialTexs;
    delete [] d->_materialFiles;
    delete [] d->_materialWidths;
    delete [] d->_materialHeights;
    int nb = PyList_Size(materials);
    
    d->_nMaterials = nb;
    d->_materialFiles = new char* [nb];
    d->_materialWidths = new int [nb];
    d->_materialHeights = new int [nb];
    d->_materialTexs = new GLuint [nb];
    for (int i = 0; i < nb; i++)
    {
      PyObject* o = PyList_GetItem(materials, i);
      char* file = NULL;
      if (PyString_Check(o)) file = PyString_AsString(o);
#if PY_VERSION_HEX >= 0x03000000
      else if (PyUnicode_Check(o)) file = PyBytes_AsString(PyUnicode_AsUTF8String(o)); 
#endif
      d->_materialFiles[i] = new char [128];  
      strcpy(d->_materialFiles[i], file);
      d->_materialTexs[i] = 0;
    }
  }

  if (bumpMaps != Py_None)
  {
    for (int i = 0; i < d->_nBumpMaps; i++)
    {
      delete [] d->_bumpMapFiles[i];
      //if (d->_bumpMapTexs[i] != 0) glDeleteTextures(1, &d->_bumpMapTexs[i]);
    }
    delete [] d->_bumpMapTexs;
    delete [] d->_bumpMapFiles;
    delete [] d->_bumpMapWidths;
    delete [] d->_bumpMapHeights;
    int nb = PyList_Size(bumpMaps);
    d->_nBumpMaps = nb;
    d->_bumpMapFiles = new char* [nb];
    d->_bumpMapWidths = new int [nb];
    d->_bumpMapHeights = new int [nb];
    d->_bumpMapTexs = new GLuint [nb];
    for (int i = 0; i < nb; i++)
    {
      PyObject* o = PyList_GetItem(bumpMaps, i);
      if (o == Py_None)
      {
        d->_bumpMapFiles[i] = NULL;
      }
      else
      {
        char* file = NULL;
        if (PyString_Check(o)) file = PyString_AsString(o);
#if PY_VERSION_HEX >= 0x03000000
        else if (PyUnicode_Check(o)) file = PyBytes_AsString(PyUnicode_AsUTF8String(o)); 
#endif
        d->_bumpMapFiles[i] = new char [128];  
        strcpy(d->_bumpMapFiles[i], file);
      }
      d->_bumpMapTexs[i] = 0;
    }
  }

  // force render
  d->ptrState->render = 1;
  return Py_BuildValue("i", KSUCCESS);
}

//=============================================================================
/* Set CPlot mode 
   0: mesh
   1: solid
   2: render
   3: scalar
   4: vector
   Retourne 0 (KFAILED), 1 (KSUCCESS)
 */
//=============================================================================
PyObject* K_CPLOT::setMode(PyObject* self, PyObject* args)
{
  PyObject* modeObject;
  if (!PyArg_ParseTuple(args, "O", &modeObject)) return NULL;

  // Get mode from int or string
  E_Int mode = getMode(modeObject);
  if (mode == -1) return Py_BuildValue("i", KFAILED);
  
  Data* d = Data::getInstance();
  if (mode <= 2) d->ptrState->mode = mode;
  else if (d->_numberOfZones > 0) // check variables
  {
    int nv = d->_zones[0]->nfield;
    if (mode == 3 && d->ptrState->scalarField < nv) d->ptrState->mode = mode;
    if (mode == 4 && d->ptrState->vectorField1 < nv && 
        d->ptrState->vectorField2 < nv && d->ptrState->vectorField3 < nv) 
      d->ptrState->mode = mode;
  }
  else d->ptrState->mode = mode;
  
  if (d->_numberOfZones > 0) d->ptrState->render = 1;
  return Py_BuildValue("i", KSUCCESS);
}

//=============================================================================
/* 
   Change CPlot scalar variable 
   Retourne 0 (KFAILED), 1 (KSUCCESS)
 */
//=============================================================================
PyObject* K_CPLOT::changeVariable(PyObject* self, PyObject* args)
{
  Data* d = Data::getInstance();
  if (d->_numberOfZones == 0) return Py_BuildValue("i", KFAILED);

  int scalarField = d->ptrState->scalarField;
  int nv = d->_zones[0]->nfield;
  if (nv == 0) return Py_BuildValue("i", KSUCCESS); // no field
  scalarField++;
  if (scalarField >= nv) scalarField = 0;
  d->ptrState->scalarField = scalarField;
  d->ptrState->render = 1;
  return Py_BuildValue("i", KSUCCESS);
}

//=============================================================================
/* 
   Change CPlot style 
 */
//=============================================================================
PyObject* K_CPLOT::changeStyle(PyObject* self, PyObject* args)
{
  Data* d = Data::getInstance();
  d->changeAppearance();
  d->ptrState->render = 1;
  return Py_BuildValue("i", KSUCCESS);
}

//=============================================================================
/* 
   Change Info display 
 */
//=============================================================================
PyObject* K_CPLOT::changeInfoDisplay(PyObject* self, PyObject* args)
{
  Data* d = Data::getInstance();
  if (d->ptrState->info == 1) d->ptrState->info = 0;
  else d->ptrState->info = 1;
  if (d->ptrState->bb == 1) d->ptrState->bb = 0;
  else d->ptrState->bb = 1;
  if (d->ptrState->header == 1) d->ptrState->header = 0;
  else d->ptrState->header = 1;
  d->ptrState->render = 1;
  return Py_BuildValue("i", KSUCCESS);
}

//=============================================================================
/* 
   Change CPlot autoblanking 
   Retourne 0 (KFAILED), 1 (KSUCCESS)
 */
//=============================================================================
PyObject* K_CPLOT::changeBlanking(PyObject* self, PyObject* args)
{
  Data* d = Data::getInstance();
  if (d->_numberOfZones == 0) return Py_BuildValue("i", KFAILED);

  // Attends que la ressource GPU soit libre pour la modifier :
  d->ptrState->syncGPURes();
  // Libere les ressources utilisees par les zones :
  for (int i = 0; i < d->_numberOfZones; i++)
  {
    Zone* z = d->_zones[i];
    z->freeGPURessources( true, false );
  }

  int v = d->ptrState->autoblank;
  if (v == 0) d->ptrState->autoblank = 1;
  else d->ptrState->autoblank = 0;
  d->autoPlugins();
  d->ptrState->render = 1;
  return Py_BuildValue("i", KSUCCESS);
}

//=============================================================================
/* 
   Set dim (1D/2D/3D)
 */
//=============================================================================
PyObject* K_CPLOT::setDim(PyObject* self, PyObject* args)
{
  int dim;
  if (!PyArg_ParseTuple(args, "i", &dim)) return NULL;
  
  if (dim < 1 || dim > 3)
  {
    PyErr_SetString(PyExc_TypeError,
                    "setDim: dim must be 1, 2 or 3.");
    return NULL;
  }

  Data* d = Data::getInstance();
  if (d->ptrState->dim == 1 && dim == 2) 
  {
    d->ptrState->dim = dim;
    d->ptrState->syncGPURes();
    for (int i = 0; i < d->_numberOfZones; i++)
    {
      Zone* z = d->_zones[i];
      z->freeGPURessources(true, false);
    }
    d->_view.xcam1D = d->_view.xcam;
    d->_view.ycam1D = d->_view.ycam;
    d->_view.zcam1D = d->_view.zcam;
    d->_view.xeye1D = d->_view.xeye;
    d->_view.yeye1D = d->_view.yeye;
    d->_view.zeye1D = d->_view.zeye;
    d->_view.dirx1D = d->_view.dirx;
    d->_view.diry1D = d->_view.diry;
    d->_view.dirz1D = d->_view.dirz;
    d->_view.xcam = d->_view.xcam2D;
    d->_view.ycam = d->_view.ycam2D;
    d->_view.zcam = d->_view.zcam2D;
    d->_view.xeye = d->_view.xeye2D;
    d->_view.yeye = d->_view.yeye2D;
    d->_view.zeye = d->_view.zeye2D;
    d->_view.dirx = d->_view.dirx2D;
    d->_view.diry = d->_view.diry2D;
    d->_view.dirz = d->_view.dirz2D;
    //d->roll1Dto2D();
  }
  else if (d->ptrState->dim == 1 && dim == 3) 
  {
    d->ptrState->dim = dim;
    d->ptrState->syncGPURes();
    for (int i = 0; i < d->_numberOfZones; i++)
    {
      Zone* z = d->_zones[i];
      z->freeGPURessources(true, false);
    }
    d->_view.xcam1D = d->_view.xcam;
    d->_view.ycam1D = d->_view.ycam;
    d->_view.zcam1D = d->_view.zcam;
    d->_view.xeye1D = d->_view.xeye;
    d->_view.yeye1D = d->_view.yeye;
    d->_view.zeye1D = d->_view.zeye;
    d->_view.dirx1D = d->_view.dirx;
    d->_view.diry1D = d->_view.diry;
    d->_view.dirz1D = d->_view.dirz;
    d->_view.xcam = d->_view.xcam3D;
    d->_view.ycam = d->_view.ycam3D;
    d->_view.zcam = d->_view.zcam3D;
    d->_view.xeye = d->_view.xeye3D;
    d->_view.yeye = d->_view.yeye3D;
    d->_view.zeye = d->_view.zeye3D;
    d->_view.dirx = d->_view.dirx3D;
    d->_view.diry = d->_view.diry3D;
    d->_view.dirz = d->_view.dirz3D;
    //d->roll1Dto3D();
  }
  else if (d->ptrState->dim == 2 && dim == 3) 
  {
    d->ptrState->dim = dim;
    d->ptrState->syncGPURes();
    for (int i = 0; i < d->_numberOfZones; i++)
    {
      Zone* z = d->_zones[i];
      z->freeGPURessources(true, false);
    }
    d->_view.xcam2D = d->_view.xcam;
    d->_view.ycam2D = d->_view.ycam;
    d->_view.zcam2D = d->_view.zcam;
    d->_view.xeye2D = d->_view.xeye;
    d->_view.yeye2D = d->_view.yeye;
    d->_view.zeye2D = d->_view.zeye;
    d->_view.dirx2D = d->_view.dirx;
    d->_view.diry2D = d->_view.diry;
    d->_view.dirz2D = d->_view.dirz;
    d->_view.xcam = d->_view.xcam3D;
    d->_view.ycam = d->_view.ycam3D;
    d->_view.zcam = d->_view.zcam3D;
    d->_view.xeye = d->_view.xeye3D;
    d->_view.yeye = d->_view.yeye3D;
    d->_view.zeye = d->_view.zeye3D;
    d->_view.dirx = d->_view.dirx3D;
    d->_view.diry = d->_view.diry3D;
    d->_view.dirz = d->_view.dirz3D;
    //d->roll2Dto3D();
  }
  else if (d->ptrState->dim == 2 && dim == 1) 
  {
    d->ptrState->dim = dim;
    d->ptrState->syncGPURes();
    for (int i = 0; i < d->_numberOfZones; i++)
    {
      Zone* z = d->_zones[i];
      z->freeGPURessources(true, false);
    }
    d->_view.xcam2D = d->_view.xcam;
    d->_view.ycam2D = d->_view.ycam;
    d->_view.zcam2D = d->_view.zcam;
    d->_view.xeye2D = d->_view.xeye;
    d->_view.yeye2D = d->_view.yeye;
    d->_view.zeye2D = d->_view.zeye;
    d->_view.dirx2D = d->_view.dirx;
    d->_view.diry2D = d->_view.diry;
    d->_view.dirz2D = d->_view.dirz;
    d->_view.xcam = d->_view.xcam1D;
    d->_view.ycam = d->_view.ycam1D;
    d->_view.zcam = d->_view.zcam1D;
    d->_view.xeye = d->_view.xeye1D;
    d->_view.yeye = d->_view.yeye1D;
    d->_view.zeye = d->_view.zeye1D;
    d->_view.dirx = d->_view.dirx1D;
    d->_view.diry = d->_view.diry1D;
    d->_view.dirz = d->_view.dirz1D;
    //d->roll2Dto1D();
  }
  else if (d->ptrState->dim == 3 && dim == 1) 
  {
    d->ptrState->dim = dim;
    d->ptrState->syncGPURes();
    for (int i = 0; i < d->_numberOfZones; i++)
    {
      Zone* z = d->_zones[i];
      z->freeGPURessources(true);
      //d->ptrState->syncGPURes();
    }
    d->_view.xcam3D = d->_view.xcam;
    d->_view.ycam3D = d->_view.ycam;
    d->_view.zcam3D = d->_view.zcam;
    d->_view.xeye3D = d->_view.xeye;
    d->_view.yeye3D = d->_view.yeye;
    d->_view.zeye3D = d->_view.zeye;
    d->_view.dirx3D = d->_view.dirx;
    d->_view.diry3D = d->_view.diry;
    d->_view.dirz3D = d->_view.dirz;
    d->_view.xcam = d->_view.xcam1D;
    d->_view.ycam = d->_view.ycam1D;
    d->_view.zcam = d->_view.zcam1D;
    d->_view.xeye = d->_view.xeye1D;
    d->_view.yeye = d->_view.yeye1D;
    d->_view.zeye = d->_view.zeye1D;
    d->_view.dirx = d->_view.dirx1D;
    d->_view.diry = d->_view.diry1D;
    d->_view.dirz = d->_view.dirz1D;
    //d->roll3Dto1D();
  }
  else if (d->ptrState->dim == 3 && dim == 2) 
  {
    d->ptrState->dim = dim;
    d->ptrState->syncGPURes();
    for (int i = 0; i < d->_numberOfZones; i++)
    {
      Zone* z = d->_zones[i];
      z->freeGPURessources(true, false);
    }
    d->_view.xcam3D = d->_view.xcam;
    d->_view.ycam3D = d->_view.ycam;
    d->_view.zcam3D = d->_view.zcam;
    d->_view.xeye3D = d->_view.xeye;
    d->_view.yeye3D = d->_view.yeye;
    d->_view.zeye3D = d->_view.zeye;
    d->_view.dirx3D = d->_view.dirx;
    d->_view.diry3D = d->_view.diry;
    d->_view.dirz3D = d->_view.dirz;
    d->_view.xcam = d->_view.xcam2D;
    d->_view.ycam = d->_view.ycam2D;
    d->_view.zcam = d->_view.zcam2D;
    d->_view.xeye = d->_view.xeye2D;
    d->_view.yeye = d->_view.yeye2D;
    d->_view.zeye = d->_view.zeye2D;
    d->_view.dirx = d->_view.dirx2D;
    d->_view.diry = d->_view.diry2D;
    d->_view.dirz = d->_view.dirz2D;
    //d->roll3Dto2D();
  }

  // Determine la variable secondaire (en 2D, c'est le repere)
  if (dim == 2)
  {
    if (d->_view.diry == 0. && d->_view.dirz == 0.) d->ptrState->var2D = 0; // xy
    else
    {
      double vx = d->_view.xcam - d->_view.xeye;
      double vy = d->_view.ycam - d->_view.yeye;
      double vz = d->_view.zcam - d->_view.zeye;
      if (vx == 0. && vz == 0.) d->ptrState->var2D = 1; // xz
      if (vy == 0. && vz == 0.) d->ptrState->var2D = 2; // yz
    }
  }
  
  d->ptrState->render = 1;
  return Py_BuildValue("i", KSUCCESS);
}

//=============================================================================
/* 
   setSelectedZones
 */
//=============================================================================
PyObject* K_CPLOT::setSelectedZones(PyObject* self, PyObject* args)
{
  PyObject* o;
  if (!PyArg_ParseTuple(args, "O", &o)) return NULL;
  
  if (PyList_Check(o) == 0)
  {
    PyErr_SetString(
      PyExc_TypeError, 
      "setSelectedZones: arg must be a list.");
    return NULL;
  }
  Data* d = Data::getInstance();
  PyObject* tpl;
  int n = PyList_Size(o);
  for (int i = 0; i < n; i++)
  {
    tpl = PyList_GetItem(o, i);
    // tpl must be a tuple of two ints (no, 1)
    if (PyTuple_Check(tpl) == 0)
    {
      PyErr_SetString(
        PyExc_TypeError, 
        "setSelectedZones: arg must be a list of tuples (noz, 1).");
      return NULL;
    }
    if (PyTuple_Size(tpl) != 2) 
    {
      PyErr_SetString(
        PyExc_TypeError, 
        "setSelectedZones: arg must be a list of tuples (noz, 1).");
      return NULL;
    }
    int noz = PyLong_AsLong(PyTuple_GetItem(tpl, 0));
    int status = PyLong_AsLong(PyTuple_GetItem(tpl, 1));
    if (noz < 0 || noz > d->_numberOfZones-1)
    {
      printf("Warning: setSelectedZones: number of zone is invalid.\n");
    }
    else if (status != 1 && status != 0)
    {
      printf("Warning: setSelectedZones: status of zone is invalid.\n");
    }
    else
    {
      if (status == 1)
      {
        d->ptrState->selectedZone = noz+1;
        Zone* z = d->_zones[noz];
        if (z->selected == 0) // not previously selected
        {
          d->ptrState->activePointX = z->x[0];
          d->ptrState->activePointY = z->y[0];
          d->ptrState->activePointZ = z->z[0];
        }
      }
      d->_zones[noz]->selected = status;
    }
  }
  d->ptrState->render = 1;
  return Py_BuildValue("i", KSUCCESS);
}

//=============================================================================
/*
  Unselect all zones
*/
//=============================================================================
PyObject* K_CPLOT::unselectAllZones(PyObject* self, PyObject* args)
{
  Data* d = Data::getInstance();
  for (int i = 0; i < d->_numberOfZones; i++)
    d->_zones[i]->selected = 0;

  d->ptrState->selectedZone = 0;
  d->ptrState->render = 1;
  return Py_BuildValue("i", KSUCCESS);
}

//=============================================================================
/* 
   setActiveZones
 */
//=============================================================================
PyObject* K_CPLOT::setActiveZones(PyObject* self, PyObject* args)
{
  PyObject* o;
  if (!PyArg_ParseTuple(args, "O", &o)) return NULL;
  
  if (PyList_Check(o) == 0)
  {
    PyErr_SetString(
      PyExc_TypeError, 
      "setActiveZones: arg must be a list.");
    return NULL;
  }
  Data* d = Data::getInstance();
  PyObject* tpl;
  int n = PyList_Size(o);
  for (int i = 0; i < n; i++)
  {
    tpl = PyList_GetItem(o, i);
    // tpl must be a tuple of two ints (no, 1)
    if (PyTuple_Check(tpl) == 0)
    {
      PyErr_SetString(
        PyExc_TypeError, 
        "setActiveZones: arg must be a list of tuples (noz, 1).");
      return NULL;
    }
    if (PyTuple_Size(tpl) != 2)
    {
      PyErr_SetString(
        PyExc_TypeError, 
        "setActiveZones: arg must be a list of tuples (noz, 1).");
      return NULL;
    }
    int noz = PyLong_AsLong(PyTuple_GetItem(tpl, 0));
    int status = PyLong_AsLong(PyTuple_GetItem(tpl, 1));
    if (noz < 0 || noz > d->_numberOfZones-1)
    {
      // Je supprime ce message, car il arrive que l'interface
      // essaie sans succes de reactiver certains blocs et que ce n'est pas
      // bloquant
      //printf("Warning: setActiveZones: number of zone is invalid (%d).\n", 
      //       noz);
    }
    else if (status != 1 && status != 0)
    {
      printf("Warning: setActiveZones: status of zone is invalid (%d).\n",
             status);
    } 
    else d->_zones[noz]->active = status;
  }
  d->ptrState->render = 1;
  return Py_BuildValue("i", KSUCCESS);
}

//=============================================================================
/* 
   setZoneNames
 */
//=============================================================================
PyObject* K_CPLOT::setZoneNames(PyObject* self, PyObject* args)
{
  PyObject* o;
  if (!PyArg_ParseTuple(args, "O", &o)) return NULL;
  
  if (PyList_Check(o) == 0)
  {
    PyErr_SetString(
      PyExc_TypeError, 
      "setZoneNames: arg must be a list.");
    return NULL;
  }
  Data* d = Data::getInstance();
  PyObject* tpl;
  int n = PyList_Size(o);
  for (int i = 0; i < n; i++)
  {
    tpl = PyList_GetItem(o, i);
    // tpl must be a tuple of int, string of type (no, 'newName')
    if (PyTuple_Check(tpl) == 0)
    {
      PyErr_SetString(
        PyExc_TypeError, 
        "setZoneNames: arg must be a list of tuples (noz, 'name').");
      return NULL;
    }
    if (PyTuple_Size(tpl) != 2)
    {
      PyErr_SetString(
        PyExc_TypeError, 
        "setZoneNames: arg must be a list of tuples (noz, 'name').");
      return NULL;
    }
    int noz = PyLong_AsLong(PyTuple_GetItem(tpl, 0));
    
    char* name = NULL;
    PyObject* l = PyTuple_GetItem(tpl, 1);
    if (PyString_Check(l)) name = PyString_AsString(l);
#if PY_VERSION_HEX >= 0x03000000
    else if (PyUnicode_Check(l)) name = PyBytes_AsString(PyUnicode_AsUTF8String(l)); 
#endif
    
    if (noz < 0 || noz > d->_numberOfZones-1)
    {
      printf("Warning: setZoneNames: number of zone is invalid.\n");
    }
    else strcpy(d->_zones[noz]->zoneName, name);
  }
  d->ptrState->render = 1;
  return Py_BuildValue("i", KSUCCESS);
}

//==============================================================================
/*
  Set active point
*/
//==============================================================================
PyObject* K_CPLOT::setActivePoint(PyObject* self, PyObject* args)
{
  E_Float px, py, pz;
#if defined E_DOUBLEREAL
  if (!PyArg_ParseTuple(args, "ddd", &px, &py, &pz)) return NULL;
#else
  if (!PyArg_ParseTuple(args, "fff", &px, &py, &pz)) return NULL;
#endif
  Data* d = Data::getInstance();
  
  int ret; int zone, ind, indE; double dist;
  ret = d->findBlockContaining(px, py, pz, zone, ind, indE, dist);
  Zone* z = d->_zones[zone];
  if (ret == 1)
  {
    d->ptrState->activePointX = px;
    d->ptrState->activePointY = py;
    d->ptrState->activePointZ = pz;
    d->ptrState->activePointZBuf = 1.;

    if (zone < d->_numberOfStructZones)
    {
      StructZone* zz = (StructZone*)z;
      int ni = zz->ni; 
      int nj = zz->nj;
      int k = ind / (ni*nj);
      int j = (ind - k*ni*nj)/ni;
      int i = ind - k*ni*nj - j*ni;
      d->ptrState->activePointI = i+1;
      d->ptrState->activePointJ = j+1;
      d->ptrState->activePointK = k+1;
    }
    else
    {
      d->ptrState->activePointI = ind; // indice du noeud le plus proche
      d->ptrState->activePointJ = indE; // indice de l'element contenant P
      UnstructZone* zz = (UnstructZone*)z;
      if (zz->eltType != 10) // autre que NGON
      {
        int* c = zz->connect;
        int size = zz->eltSize;
        int ne = zz->ne;
        int v = 0;
        for (v = 0; v < size; v++)
        {
          if (c[indE+v*ne] == ind+1) break;
        }
        d->ptrState->activePointK = -v-1;
      }
      else d->ptrState->activePointK = findFace(px, py, pz, indE, zz, dist);
    }
    for (int n = 0; n < z->nfield; n++)
    {
      double* f = z->f[n];
      d->ptrState->activePointF[n] = f[ind];
    }
    return Py_BuildValue("i", KSUCCESS);
  }
  return Py_BuildValue("i", KFAILED);
}

//=============================================================================
/* 
   Look for selected zone
*/
//=============================================================================
PyObject* K_CPLOT::lookFor(PyObject* self, PyObject* args)
{
  Data* d = Data::getInstance();
  if (d->_pref.lookFor != NULL) d->_pref.lookFor->f(d);
  d->ptrState->render = 1;
  return Py_BuildValue("i", KSUCCESS);
}

//=============================================================================
// Remet la keyboard string (keys) a 0
//=============================================================================
PyObject* K_CPLOT::resetKeyboard(PyObject* self, PyObject* args)
{
  Data* d = Data::getInstance();
  d->ptrState->kcursor = 0;
  return Py_BuildValue("i", KSUCCESS);
}

//=============================================================================
// Set shaderPath (chemin des shaders)
//=============================================================================
PyObject* K_CPLOT::setShaderPath(PyObject* self, PyObject* args)
{
  char* path;
  if (!PyArg_ParseTuple(args, "s", &path)) return NULL;
  Data* d = Data::getInstance();
  strcpy(d->ptrState->shaderPath, path);
  return Py_BuildValue("i", KSUCCESS);
}

//=============================================================================
// Set window title (file name + file path)
//=============================================================================
PyObject* K_CPLOT::setWindowTitle(PyObject* self, PyObject* args)
{
  char* file; char* path;
  if (!PyArg_ParseTuple(args, "ss", &file, &path)) return NULL;
  Data* d = Data::getInstance();
  strcpy(d->ptrState->winTitle, "CPlot - ");
  strcat(d->ptrState->winTitle, file);
  strcpy(d->ptrState->file, file);
  strcpy(d->ptrState->filePath, path);
  if (d->_winId != 0)
  {
    glutSetWindowTitle(d->ptrState->winTitle);
    glutSetIconTitle(d->ptrState->winTitle);
  }
  return Py_BuildValue("i", KSUCCESS);
}

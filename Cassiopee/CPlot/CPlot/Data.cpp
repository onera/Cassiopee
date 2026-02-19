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

#include "cplot.h"
#include "Data.h"
#if defined(_WIN32) ||defined(_WIN64)
# include <direct.h>
#endif

Data::RenderID Data::_renderID = Data::DL;
Data* Data::_instance = NULL;

//=============================================================================
/* Ce constructeur est incomplet.
   La classe Data est complete quand :
   xmin, ymin, zmin, xmax, ymax, zmax ont ete initialises
   et tous les init appeles.
*/
//=============================================================================
Data::Data(CPlotState* ptState)
{
  _instance = NULL;
  _numberOfZones = 0;
  _numberOfStructZones = 0;
  _numberOfUnstructZones = 0;
  _zones = NULL;
  _szones = NULL;
  _uzones = NULL;
  _niter = 0;
  _winId = 0;
  _texEnviron1 = 0;
  _texLeft = 0;
  _texRight = 0;
  _texSupp = 0;
  _shadowMap = 0;
  _texColormap = 0;
  _texColormapType = -1;
  _texColormapMinMax = -1;
  for (E_Int i = 0; i < 16; i++) _bias[i] = 0.;
  _bias[0] = 0.5; _bias[5]= 0.5; _bias[10] = 0.5; _bias[12] = 0.5;
  _bias[13] = 0.5; _bias[14] = 0.5; _bias[15] = 1.;

  // Initial window size
  _view.w = 0;
  _view.h = 0;

  // Camera view angle
  _view.angle = 50.;

  // frameBufferSize doit etre > taille du viewport
  for (E_Int i = 0; i < 10; i++) _frameBufferSize[i] = 2048;
  // Nbre de voxels dans les voxels buffer
  _voxelBufferSize = 96;
  // Init variables
  ptrState = ptState;
  initState();
  // Initial global bounding box
  xmax = 1; xmin = 0; ymax = 1; ymin = 0; zmax = 1; zmin = 0;
  // Initial iso values
  _nfield = 0;
  minf = NULL;
  maxf = NULL;
  _isoMin = NULL;
  _isoMax = NULL;
  _isoAlphaMin = NULL;
  _isoAlphaMax = NULL;
  _niso = NULL;
  _isoColormap = NULL;

  // Colormaps (init on demand)
  _colormapSizeViridis = 0;
  _colormapRViridis = NULL;
  _colormapGViridis = NULL; 
  _colormapBViridis = NULL; 
  _colormapSizeInferno = 0;
  _colormapRInferno = NULL;
  _colormapGInferno = NULL; 
  _colormapBInferno = NULL; 
  _colormapSizeMagma = 0;
  _colormapRMagma = NULL;
  _colormapGMagma = NULL; 
  _colormapBMagma = NULL; 
  _colormapSizePlasma = 0;
  _colormapRPlasma = NULL;
  _colormapGPlasma = NULL; 
  _colormapBPlasma = NULL; 
  _colormapSizeJet = 0;
  _colormapRJet = NULL;
  _colormapGJet = NULL; 
  _colormapBJet = NULL; 
  _colormapSizeGreys = 0;
  _colormapRGreys = NULL;
  _colormapGGreys = NULL; 
  _colormapBGreys = NULL; 
  _colormapSizeNiceBlue = 0;
  _colormapRNiceBlue = NULL;
  _colormapGNiceBlue = NULL; 
  _colormapBNiceBlue = NULL; 
  _colormapSizeGreens = 0;
  _colormapRGreens = NULL;
  _colormapGGreens = NULL; 
  _colormapBGreens = NULL; 

  // billBoards settings
  E_Int nb = 5; E_Int c = 0;
  _nBillBoards = nb;
  _billBoardFiles = new char* [nb];
  _billBoardNis = new E_Int [nb];
  _billBoardNjs = new E_Int [nb];
  _billBoardWidths = new E_Int [nb]; // init from files
  _billBoardHeights = new E_Int [nb];
  _billBoardTexs = new GLuint [nb];
  _billBoardFiles[c] = new char [128];
  strcpy(_billBoardFiles[c], "smoke1.png");
  _billBoardTexs[c] = 0;
  _billBoardNis[c] = 4; _billBoardNjs[c] = 4; c++;
  _billBoardFiles[c] = new char [128]; 
  strcpy(_billBoardFiles[c], "smoke2.png");
  _billBoardTexs[c] = 0;
  _billBoardNis[c] = 4; _billBoardNjs[c] = 4; c++;
  _billBoardFiles[c] = new char [128]; 
  strcpy(_billBoardFiles[c], "smoke3.png");
  _billBoardTexs[c] = 0;
  _billBoardNis[c] = 4; _billBoardNjs[c] = 4; c++;
  _billBoardFiles[c] = new char [128];
  strcpy(_billBoardFiles[c], "fire1.png");
  _billBoardTexs[c] = 0;
  _billBoardNis[c] = 2; _billBoardNjs[c] = 2; c++;
  _billBoardFiles[c] = new char [128];
  strcpy(_billBoardFiles[c], "rock1.png");
  _billBoardTexs[c] = 0;
  _billBoardNis[c] = 4; _billBoardNjs[c] = 4; c++;
  
  // material settings
  _nMaterials = 0;
  _materialFiles = NULL;
  _materialWidths = NULL;
  _materialHeights = NULL;
  _materialTexs = NULL;
  
  _nBumpMaps = 0; // must be equal to nMaterials
  _bumpMapFiles = NULL;
  _bumpMapWidths = NULL;
  _bumpMapHeights = NULL;
  _bumpMapTexs = NULL;
  
  // Fonts
  _font1Size = 12;
  _font2Size = 18;
  _font3Size = 12;
  _oglText1 = NULL;
  _oglText2 = NULL;
  _oglText3 = NULL;
  
  // Init cam
  initCam();
  _CDisplayIsLaunched = 0;

}

//=============================================================================
Data::~Data()
{
  if (minf != NULL) delete [] minf;
  if (maxf != NULL) delete [] maxf;
  if (_isoMin != NULL) delete [] _isoMin;
  if (_isoMax != NULL) delete [] _isoMax;
  if (_isoAlphaMin != NULL) delete [] _isoAlphaMin;
  if (_isoAlphaMax != NULL) delete [] _isoAlphaMax;
  if (_isoColormap != NULL) delete [] _isoColormap;
  if (_colormapRViridis != NULL) delete [] _colormapRViridis;
  if (_colormapGViridis != NULL) delete [] _colormapGViridis;
  if (_colormapBViridis != NULL) delete [] _colormapBViridis;

  size_t ns = _slots1D.size();
  for (size_t i = 0; i < ns; i++) delete _slots1D[i];
  pthread_mutex_destroy(&ptrState->lock_mutex);
  pthread_mutex_destroy(&ptrState->gpures_mutex);
  // delete all allocated data of ptrState
  delete [] ptrState->message;
  delete [] ptrState->colormapR;
  delete [] ptrState->colormapG;
  delete [] ptrState->colormapB;
  delete ptrState;
  // delete billBoardStorage
  for (E_Int i = 0; i < _nBillBoards; i++)
  {
    delete [] _billBoardFiles[i];
    if (_billBoardTexs[i] != 0) glDeleteTextures(1, &_billBoardTexs[i]);
  }
  delete [] _billBoardTexs;
  delete [] _billBoardNis; delete [] _billBoardNjs;
  delete [] _billBoardFiles;

  // delete material storage
  for (E_Int i = 0; i < _nMaterials; i++)
  {
    delete [] _materialFiles[i];
    if (_materialTexs[i] != 0) glDeleteTextures(1, &_materialTexs[i]);
  }
  delete [] _materialTexs;
  delete [] _materialFiles;

  for (E_Int i = 0; i < _nBumpMaps; i++)
  {
    delete [] _bumpMapFiles[i];
    if (_bumpMapTexs[i] != 0) glDeleteTextures(1, &_bumpMapTexs[i]);
  }
  delete [] _bumpMapTexs;
  delete [] _bumpMapFiles;
} 

//=============================================================================
/*
   Initialise l'etat.
   Pas de parametres.
*/
//=============================================================================
void Data::initState()
{
  _numberOfZones = 0;
  _numberOfStructZones = 0;
  _numberOfUnstructZones = 0;

  // State variables
  ptrState->mode = MESH;
  ptrState->dim = 3;

  ptrState->scalarField = 0;
  ptrState->vectorField1 = 0; ptrState->vectorField2 = 0; ptrState->vectorField3 = 0;

  ptrState->ijk1D = 0;
  ptrState->var1D = 0;
  ptrState->var2D = 0;

  ptrState->ghostifyDeactivatedZones = 0;
  ptrState->edgifyDeactivatedZones = 0;
  ptrState->edgifyActivatedZones = 0;
  ptrState->simplifyOnDrag = 0;

  ptrState->alpha = 1.;
  ptrState->currentMenuItem = 0;
  ptrState->smoke = 0;
  ptrState->inmenu = 0;

  // Selected zone
  ptrState->selectedZone = 0;
  ptrState->deactivatedZones = NULL;

  // Active point (=selected point)
  ptrState->activePointX = 0.;
  ptrState->activePointY = 0.;
  ptrState->activePointZ = 0.;
  ptrState->activePointZBuf = 0.999;
  ptrState->activePointI = 0;
  ptrState->activePointJ = 0;
  ptrState->activePointK = 0;
  ptrState->activePointL = 0;
  ptrState->activePointF = NULL;

  // Render
  ptrState->farClip = 0;
  ptrState->render = 1;
  ptrState->bb = 0;
  ptrState->header = 1;
  ptrState->info = 1;
  ptrState->menu = 1;
  ptrState->texture = 0;
  ptrState->isoLegend = 0;
  ptrState->offscreen = 0;
  ptrState->frameBuffer = 0; 
  for (E_Int i = 0; i < 10; i++) ptrState->offscreenBuffer[i] = NULL;
  for (E_Int i = 0; i < 10; i++) ptrState->offscreenDepthBuffer[i] = NULL;

  // overlay message
  ptrState->message = NULL;

  // Stereo
  ptrState->stereo = 0;
  ptrState->stereoDist = 1./30.;

  // Post-processing
  ptrState->lightOffsetX = 0.; // entre -1 et 1
  ptrState->lightOffsetY = 0.; // entre 0 et 5
  ptrState->DOF = 0;
  ptrState->shadow = 0;
  ptrState->dofPower = 0.; // inactif par defaut
  ptrState->gamma = 1.; // inactif par defaut
  ptrState->toneMapping = 0; // rien par defaut
  ptrState->sobelThreshold = -0.5; // inactif par defaut
  ptrState->sharpenPower = -0.5; // inactif par defaut
  ptrState->ssaoPower = -0.5; // inactif par defaut
  ptrState->panorama = 0; // trigger panorama rendering

  strcpy(ptrState->winTitle, "CPlot - array/pyTree display");
  strcpy(ptrState->file, "tmp");
  strcpy(ptrState->filePath, ".");
  strcpy(ptrState->localPathName, ".");
  
  // Local directory path name
#if defined(_WIN32) ||defined(_WIN64)
  _getcwd(ptrState->localPathName, 120);
#else
   getcwd(ptrState->localPathName, 120);
#endif

  // Chemin des shaders
  ptrState->shaderPath[0] = '\0';

  // Mouse
  ptrState->activeMouseButton = GLUT_MIDDLE_BUTTON;
  ptrState->activeMouseX = 0;
  ptrState->activeMouseY = 0;
  ptrState->currentMouseButton = 5; // unclicked
  ptrState->ondrag = 0; // not drag moving

  // Keyboard
  ptrState->kcursor = 0; // position in keys
  ptrState->kkeysActivated = 1; // activate short cuts in CPlot

  // Styles
  ptrState->meshStyle = 2;
  ptrState->solidStyle = 0;
  ptrState->scalarStyle = 0;
  ptrState->vectorStyle = 0;
  ptrState->vectorScale = 100.f;
  ptrState->vectorDensity = 0.f;
  ptrState->vectorNormalize = 0;
  ptrState->vectorShowSurface = 1;
  ptrState->vectorShape = 0;
  ptrState->vector_projection = 0;
  ptrState->selectionStyle = 0;
  ptrState->colormap = 0;
  ptrState->colormapR1 = 0.;
  ptrState->colormapG1 = 0.;
  ptrState->colormapB1 = 0.; 
  ptrState->colormapR2 = 1.;
  ptrState->colormapG2 = 1.; 
  ptrState->colormapB2 = 1.;
  ptrState->colormapR3 = 0.5; 
  ptrState->colormapG3 = 0.5;
  ptrState->colormapB3 = 0.5;
  ptrState->colormapR = NULL;
  ptrState->colormapG = NULL;
  ptrState->colormapB = NULL;
  ptrState->isoLight = 1;
  ptrState->niso = 25;
  ptrState->isoEdges = -0.5;

  // Lock + timestep pour le display
  pthread_mutex_init(&ptrState->lock_mutex, NULL);
  //ptrState->lock_mutex = PTHREAD_MUTEX_INITIALIZER;
  pthread_cond_init(&ptrState->unlocked_display, NULL);
  //ptrState->unlocked_display = PTHREAD_COND_INITIALIZER;
  ptrState->lock = 0;

  ptrState->timeStep = 50; // en ms
  ptrState->mouseFirstClicked = 0;

  // Export
  strcpy(ptrState->exportFile, "CPlot");
  ptrState->exportWidth = -1;
  ptrState->exportHeight = -1;
  ptrState->exportAA = 0;
  ptrState->continuousExport = 0;
  ptrState->ptrFile = NULL;
  ptrState->context = NULL;
  ptrState->exportNumber = 0;
  ptrState->shootScreen = 0;
  pthread_mutex_init(&ptrState->export_mutex, NULL);
  pthread_cond_init(&ptrState->unlocked_export, NULL);
  ptrState->_mustExport = 0;

  ptrState->odsRun = false;
  ptrState->odsNSlits = 0;
  ptrState->odsSlit = 0;
  ptrState->odsImage = NULL;
  ptrState->odsFrontImage = NULL;
  ptrState->odsTopImage = NULL;

  // Others
  ptrState->fullScreen = 0;
  ptrState->bgColor = 0;
  strcpy(ptrState->backgroundFile, "paperBackground1.png");
  ptrState->autoblank = 1;

  // Textures
  strcpy(ptrState->envmapFile, "windtunnel.png");
  ptrState->updateEnvmap = 0;
  ptrState->updateBackground = 0;

  // Cursor
  ptrState->updateCursor = -1;
  ptrState->cursorType = -1;

  // BillBoard
  ptrState->billBoardType = -1; // aucun loade
  ptrState->billBoardNi = 1;
  ptrState->billBoardNj = 1;
  ptrState->billBoardD = -1; // taille constante
  ptrState->billBoardT = -1; // pas de champ ref pour la turbulence
  ptrState->billBoardSize = -0.5; // size calculee par la distance

  // 1D overlay display
  ptrState->gridSizeI = 1;
  ptrState->gridSizeJ = 1;

  // time
  ptrState->ktime = 0;
  ptrState->ktimer = 0;

  // High order :
  ptrState->outer_tesselation = 2;
  ptrState->inner_tesselation = 1;

  // locks
  ptrState->freeGPURes = 0;

  pthread_mutex_init(&ptrState->gpures_mutex, NULL);
  //ptrState->gpures_mutex = PTHREAD_MUTEX_INITIALIZER;
  pthread_cond_init(&ptrState->unlocked_gpures, NULL);
  //ptrState->unlocked_gpures = PTHREAD_COND_INITIALIZER;
  ptrState->_lockGPURes = 0;
}
//=============================================================================
void Data::freeGPUResources(int mode, int start, int end, int permanent)
{
  CPlotState& state = (*ptrState);
  state.freeGPUData[0] = mode; state.freeGPUData[1] = start;
  state.freeGPUData[2] = end; state.freeGPUData[3] = permanent;
  state.freeGPURes = 1;
  freeGPURes(state.freeGPUData[0], state.freeGPUData[1], state.freeGPUData[2], state.freeGPUData[3]);
}
//=============================================================================
void Data::updateGPUResources(int mode, int size, int permanent, void* updatedPointer)
{
  CPlotState& state = (*ptrState);
  state.freeGPUData[0] = mode; state.freeGPUData[1] = size;
  state.freeGPUData[2] = -1; state.freeGPUData[3] = permanent;
  state.freeGPUPtr = static_cast<int*>(updatedPointer);
  ptrState->freeGPURes = 1;
  freeGPURes(state.freeGPUData[0], state.freeGPUData[1],
	         state.freeGPUPtr, state.freeGPUData[3]);
}
//=============================================================================
/*
   InitCam
   Init the Camera postion, the eye position and
   the camera direction.
   Called only once at the beginning.

   IN: xmin, ymin, zmin, xmax, ymax, zmax (globals)
   OUT: xeye, yeye, zeye, xcam, ycam, zcam, dirx, diry, dirz
         id. for saved 3D and 2D modes.
*/
//=============================================================================
void Data::initCam()
{
  double distx, disty, distz, dist;

  // Clipping planes to far
  _view.clipping = 0;

  // Eye position is initialiazed at the center of model
  _view.xeye = 0.5*(xmax+xmin);
  _view.yeye = 0.5*(ymax+ymin);
  _view.zeye = 0.5*(zmax+zmin);

  // Camera is initialized at 1.2 distance of model in x
  distx = xmax-xmin;
  disty = ymax-ymin;
  distz = zmax-zmin;

  dist = MAX(distx, disty);
  dist = MAX(dist, distz);
  dist = dist*1.2;

  if (ptrState->dim == 3)
  {
    _view.xcam = _view.xeye + dist;
    _view.ycam = _view.yeye;
    _view.zcam = _view.zeye;

    // Camera direction is z axis in 3D mode
    _view.dirx = 0.;
    _view.diry = 0.;
    _view.dirz = 1.;

    // Auto set best view (added)
    if (fabs(xmin+1.e-4) < 1.e-11 && fabs(xmax-1.e-4) < 1.e-11) //YZ
    {
     _view.xcam = _view.xeye + dist;
     _view.ycam = _view.yeye;
     _view.zcam = _view.zeye;
     _view.dirx = 0.;
     _view.diry = 0.;
     _view.dirz = 1.;
    }
    else if (fabs(ymin+1.e-4) < 1.e-11 && fabs(ymax-1.e-4) < 1.e-11) // XZ
    {
      _view.xcam = _view.xeye;
      _view.ycam = _view.yeye + dist;
      _view.zcam = _view.zeye;
      _view.dirx = 0.;
      _view.diry = 0.;
      _view.dirz = 1.;
    }
    else if (distz < 1.e-11) // XY
    {
      _view.xcam = _view.xeye;
      _view.ycam = _view.yeye;
      _view.zcam = _view.zeye + dist;
      _view.dirx = 0.;
      _view.diry = 1.;
      _view.dirz = 0.;
    }
  }
  else
  {
    _view.xcam = _view.xeye;
    _view.ycam = _view.yeye;
    _view.zcam = _view.zeye + dist;

    // Camera direction is z axis in 3D mode
    _view.dirx = 0.;
    _view.diry = 1.;
    _view.dirz = 0.;
  }

  // Set the 3D saved positions
  _view.xcam3D = _view.xeye + dist;
  _view.ycam3D = _view.yeye;
  _view.zcam3D = _view.zeye;
  _view.xeye3D = _view.xeye;
  _view.yeye3D = _view.yeye;
  _view.zeye3D = _view.zeye;
  _view.dirx3D = 0.;
  _view.diry3D = 0.;
  _view.dirz3D = 1.;

  // Set the 2D saved positions
  _view.xcam2D = _view.xeye;
  _view.ycam2D = _view.yeye;
  _view.zcam2D = _view.zeye + dist;
  _view.xeye2D = _view.xeye;
  _view.yeye2D = _view.yeye;
  _view.zeye2D = _view.zeye;
  _view.dirx2D = 0.;
  _view.diry2D = 1.;
  _view.dirz2D = 0.;

  // Set the 1D saved positions
  _view.xcam1D = _view.xeye;
  _view.ycam1D = _view.yeye;
  _view.zcam1D = _view.zeye + dist;
  _view.xeye1D = _view.xeye;
  _view.yeye1D = _view.yeye;
  _view.zeye1D = _view.zeye;
  _view.dirx1D = 0.;
  _view.diry1D = 1.;
  _view.dirz1D = 0.;
}

//============================================================================
/* Init 2D camera for the different secondary variables var2D */
//============================================================================
void Data::init2DCam()
{
  double distx, disty, distz, dist;

  // Camera is initialiazed at four distance of model in x
  distx = xmax-xmin;
  disty = ymax-ymin;
  distz = zmax-zmin;

  dist = MAX(distx, disty);
  dist = MAX(dist, distz);

  switch (ptrState->var2D)
  {
    case 0: // (x,y)
      _view.xcam2D = _view.xeye;
      _view.ycam2D = _view.yeye;
      _view.zcam2D = _view.zeye + dist;
      _view.xeye2D = _view.xeye;
      _view.yeye2D = _view.yeye;
      _view.zeye2D = _view.zeye;
      _view.dirx2D = 0.;
      _view.diry2D = 1.;
      _view.dirz2D = 0.;
      break;

    case 1:
      _view.xcam2D = _view.xeye;
      _view.ycam2D = _view.yeye - dist;
      _view.zcam2D = _view.zeye;
      _view.xeye2D = _view.xeye;
      _view.yeye2D = _view.yeye;
      _view.zeye2D = _view.zeye;
      _view.dirx2D = 0.;
      _view.diry2D = 0.;
      _view.dirz2D = 1.;
      break;

    case 2:
      _view.xcam2D = _view.xeye - dist;
      _view.ycam2D = _view.yeye;
      _view.zcam2D = _view.zeye;
      _view.xeye2D = _view.xeye;
      _view.yeye2D = _view.yeye;
      _view.zeye2D = _view.zeye;
      _view.dirx2D = 0.;
      _view.diry2D = 0.;
      _view.dirz2D = 1.;
      break;
  }
}

//============================================================================
/* Init 1D camera for the different secondary variables var1D */
//============================================================================
void Data::init1DCam()
{
  double distx, disty, distz, dist;

  distx = xmax-xmin;
  disty = ymax-ymin;
  distz = zmax-zmin;

  dist = MAX(distx, disty);
  dist = MAX(dist, distz);

  switch (ptrState->var1D)
  {
    case 0: // (x,y)
      _view.xcam1D = _view.xeye;
      _view.ycam1D = _view.yeye;
      _view.zcam1D = _view.zeye + dist;
      _view.xeye1D = _view.xeye;
      _view.yeye1D = _view.yeye;
      _view.zeye1D = _view.zeye;
      _view.dirx1D = 0.;
      _view.diry1D = 1.;
      _view.dirz1D = 0.;
      break;

    case 1: // (y,z)
      _view.xcam1D = _view.xeye + dist;
      _view.ycam1D = _view.yeye;
      _view.zcam1D = _view.zeye;
      _view.xeye1D = _view.xeye;
      _view.yeye1D = _view.yeye;
      _view.zeye1D = _view.zeye;
      _view.dirx1D = 0.;
      _view.diry1D = 0.;
      _view.dirz1D = 1.;
      break;

    case 2: // (z,x)
      _view.xcam1D = _view.xeye;
      _view.ycam1D = _view.yeye - dist;
      _view.zcam1D = _view.zeye;
      _view.xeye1D = _view.xeye;
      _view.yeye1D = _view.yeye;
      _view.zeye1D = _view.zeye;
      _view.dirx1D = 1.;
      _view.diry1D = 0.;
      _view.dirz1D = 0.;
      break;

    case 3: // s
      _view.xcam1D = _view.xeye;
      _view.ycam1D = _view.yeye;
      _view.zcam1D = _view.zeye + dist;
      _view.xeye1D = _view.xeye;
      _view.yeye1D = _view.yeye;
      _view.zeye1D = _view.zeye;
      _view.dirx1D = 0.;
      _view.diry1D = 1.;
      _view.dirz1D = 0.;
      break;
  }
}

//=============================================================================
/* Set the default preference settings */
//=============================================================================
void Data::loadPrefs()
{
  // -- Performance --
  // Speed of gfx cards
#ifdef _CARD0_
  _pref.speed = 1;
  _pref.nroll = 10;

  _pref.maxParticles = 100;
  _pref.smokeRadius = 0.1;
  _pref.emissionRadius = 1.;
  _pref.smokeTimeStep = 1.;
#endif

  // Mac
#ifdef _CARD1_
  _pref.speed = 1;
  _pref.nroll = 20;

  _pref.maxParticles = 100;
  _pref.smokeRadius = 0.5;
  _pref.emissionRadius = 1.;
  _pref.smokeTimeStep = 1.;
#endif

  // NVidia
#ifdef _CARD2_
  _pref.speed = 1;
  _pref.nroll = 20;

  _pref.maxParticles = 100;
  _pref.smokeRadius = 0.5;
  _pref.emissionRadius = 1.;
  _pref.smokeTimeStep = 1.;
#endif

  // -- Textures --
  _pref.sizeTexi = 3;
  _pref.sizeTexj = 3;

  // -- Plugin functions --
  _pref.colorMap = NULL;
  _pref.lookFor = NULL;
  _pref.addAVariable = NULL;
  _pref.blanking = NULL;
  _pref.selectNextZone = NULL;
  _pref.selectPreviousZone = NULL;
  _pref.mouseClick = NULL;
  _pref.mouseMultipleClick = NULL;
  _pref.mouseAccurateClick = NULL;
  _pref.mouseRightClick = NULL;
  _pref.screenDump = NULL;

  // Default connected plugins
  _pref.colorMap = _plugins.colorMap;
  _pref.screenDump = _plugins.screenDump;
  _pref.selectNextZone = _plugins.selectNextZone;
  _pref.selectPreviousZone = _plugins.selectPreviousZone;
  _pref.mouseClick = _plugins.mouseClick;
  _pref.mouseMultipleClick = _plugins.mouseMultipleClick;
  _pref.mouseAccurateClick = _plugins.mouseAccurateClick;
  _pref.mouseRightClick = _plugins.mouseRightClick;
  _pref.lookFor = _plugins.lookFor;
}

//=============================================================================
void Data::enforceGivenData(
  E_Int dim, E_Int mode, E_Int scalarField,
  E_Int vectorField1, E_Int vectorField2, E_Int vectorField3,
  E_Int displayBB, E_Int displayInfo, E_Int displayIsoLegend)
{
  if (dim != -1) ptrState->dim = dim;
  if (mode != -1)
  {
    if (mode == SCALARFIELD && _numberOfZones > 0)  // check valid field
    {
      E_Int nf = _zones[0]->nfield;
      if (scalarField < nf && scalarField >= 0)
      { ptrState->mode = mode; ptrState->scalarField = scalarField; }
      else if (scalarField == -1 && ptrState->scalarField < nf)
      { ptrState->mode = mode; }
    }
    else if (mode == VECTORFIELD && _numberOfZones > 0) // check valid fields
    {
      E_Int nf = _zones[0]->nfield;
      if (vectorField1 < nf && vectorField1 >= 0)
        ptrState->vectorField1 = vectorField1;
      if (vectorField2 < nf && vectorField2 >= 0)
        ptrState->vectorField2 = vectorField2;
      if (vectorField3 < nf && vectorField3 >= 0)
        ptrState->vectorField3 = vectorField3;

      if (ptrState->vectorField1 < nf && ptrState->vectorField2 < nf &&
          ptrState->vectorField3 < nf)
      { ptrState->mode = mode; }
    }
    else if (mode == RENDER)
    {
      if (_numberOfZones > 0)
      {
        E_Int nf = _zones[0]->nfield;
        if (scalarField < nf && scalarField >= 0) ptrState->scalarField = scalarField;
        if (vectorField1 < nf && vectorField1 >= 0) ptrState->vectorField1 = vectorField1;
        if (vectorField2 < nf && vectorField2 >= 0) ptrState->vectorField2 = vectorField2;
        if (vectorField3 < nf && vectorField3 >= 0) ptrState->vectorField3 = vectorField3;
      }
      ptrState->mode = mode;
    }
    else ptrState->mode = mode;
  }
  if (displayBB != -1) ptrState->bb = displayBB;
  if (displayInfo != -1)
  {
    ptrState->info = displayInfo; ptrState->header = displayInfo;
  }
  if (displayIsoLegend != -1) ptrState->isoLegend = displayIsoLegend;
}

//=============================================================================
void Data::enforceGivenData2(float xcam, float ycam, float zcam,
                             float xeye, float yeye, float zeye,
                             float dirx, float diry, float dirz,
                             float viewAngle,
                             E_Int meshStyle, E_Int solidStyle, E_Int scalarStyle,
                             E_Int vectorStyle, float vectorScale, float vectorDensity, E_Int vectorNormalize, 
                             E_Int vectorShowSurface, E_Int vectorShape, E_Int vector_projection, 
                             E_Int colormap, 
                             char* colormapC1, char* colormapC2, char* colormapC3, PyObject* colormapC,
                             E_Int niso, float isoEdges, PyObject* isoScales,
                             E_Int bgColor, char* backgroundFile, E_Int ghostifyDeactivatedZones,
                             E_Int edgifyActivatedZones,
                             E_Int edgifyDeactivatedZones,
                             E_Int shadow, E_Int dof,
                             char* exportFile, char* exportResolution, E_Int exportAA)
{
  if (xcam != -999) _view.xcam = xcam;
  if (ycam != -999) _view.ycam = ycam;
  if (zcam != -999) _view.zcam = zcam;
  if (xeye != -999) _view.xeye = xeye;
  if (yeye != -999) _view.yeye = yeye;
  if (zeye != -999) _view.zeye = zeye;
  if (dirx != -999) _view.dirx = dirx;
  if (diry != -999) _view.diry = diry;
  if (dirz != -999) _view.dirz = dirz;
  if (viewAngle != -1) { _view.angle = viewAngle; }
  if (meshStyle != -1) ptrState->meshStyle = meshStyle;
  if (solidStyle != -1) ptrState->solidStyle = solidStyle;
  if (colormap != -1)
  {
    ptrState->colormap = colormap;
    ptrState->isoLight = 0;
    _pref.colorMap = _plugins.colorMap;
    for (E_Int i = 0; i < colormap/2; i++)
    {
      if (_pref.colorMap->next == NULL)
        _pref.colorMap = _plugins.colorMap;
      else
        _pref.colorMap = _pref.colorMap->next;
    }
    if (2*(colormap/2)-colormap != 0) ptrState->isoLight = 1;
  }
  if (strlen(colormapC1) > 1)
  {
    colorString2RGB(colormapC1, ptrState->colormapR1, ptrState->colormapG1, ptrState->colormapB1);
  }
  if (strlen(colormapC2) > 1)
  {
    colorString2RGB(colormapC2, ptrState->colormapR2, ptrState->colormapG2, ptrState->colormapB2);
  }
  if (strlen(colormapC3) > 1)
  {
    colorString2RGB(colormapC3, ptrState->colormapR3, ptrState->colormapG3, ptrState->colormapB3);
  }
  if (colormapC != Py_None)
  {
    // colormapC must be a regular colormap of rgb values between 0 and 1.
    FldArrayF colors;
    if (PyArray_Check(colormapC)) // numpy
    {
      if (PyArray_ISFLOAT((PyArrayObject*)colormapC)) 
        K_ARRAY::getFromList(colormapC, colors);
      else
      {
        FldArrayI colorsI;
        K_ARRAY::getFromList(colormapC, colorsI);
        E_Int s = colorsI.getSize();
        colors.malloc(s);
        for (E_Int i = 0; i < s; i++) colors[i] = colorsI[i];
      }
    }
    else // list
    {
      K_ARRAY::getFromList(colormapC, colors);
    }

    E_Int size = colors.getSize()/3;

    ptrState->colormapSize = size;

    if (ptrState->colormapR != NULL) delete [] ptrState->colormapR;
    if (ptrState->colormapG != NULL) delete [] ptrState->colormapG;
    if (ptrState->colormapB != NULL) delete [] ptrState->colormapB;
    
    ptrState->colormapR = new float [size];
    ptrState->colormapG = new float [size];
    ptrState->colormapB = new float [size];

    float* r = ptrState->colormapR;
    float* g = ptrState->colormapG;
    float* b = ptrState->colormapB;

    float cmax = -1.;
    for (E_Int i = 0; i < size; i++) 
    {
      r[i] = colors[3*i];
      g[i] = colors[3*i+1];
      b[i] = colors[3*i+2];
      cmax = MAX(cmax, r[i]);
      cmax = MAX(cmax, g[i]);
      cmax = MAX(cmax, b[i]);
    }
    if (cmax > 1.5)
    {
      for (E_Int i = 0; i <size; i++)
      {
        r[i] = r[i]/255.;
        g[i] = g[i]/255.;
        b[i] = b[i]/255.;
      }
    }
  }

  if (scalarStyle != -1) ptrState->scalarStyle = scalarStyle;
  if (vectorStyle != -1) ptrState->vectorStyle = vectorStyle;
  if (vectorScale > 0.) ptrState->vectorScale = vectorScale;
  if (vectorDensity > -0.5) ptrState->vectorDensity = vectorDensity;
  if (vectorNormalize != -1) ptrState->vectorNormalize = vectorNormalize;
  if (vectorShowSurface != -1) ptrState->vectorShowSurface = vectorShowSurface;
  if ((vector_projection > -1) and (vector_projection < 2)) ptrState->vector_projection = vector_projection;
  if ((vectorShape > -1) and (vectorShape < 3)) ptrState->vectorShape = vectorShape;
  if (niso != -1) ptrState->niso = niso;
  if (isoEdges != -1) ptrState->isoEdges = isoEdges;

  if (isoScales != NULL)
  {
    if (PyList_Check(isoScales) == true)
    {
      E_Int size = PyList_Size(isoScales);
      
      // double liste (nouvelle methode)
      if (size > 0 && PyList_Check(PyList_GetItem(isoScales, 0)) == true) 
      {
        for (E_Int i = 0; i < size; i++)
        {
          PyObject* l = PyList_GetItem(isoScales, i);
          E_Int nelts = PyList_Size(l);
          for (E_Int j = 0; j < nelts; j++)
          {
            PyObject* f = PyList_GetItem(l, 0); // nfield
            PyObject* n = PyList_GetItem(l, 1); // nombre d'isos
            PyObject* min = PyList_GetItem(l, 2); // min for isos
            PyObject* max = PyList_GetItem(l, 3); // max for isos
            E_Int nfield = getScalarField(f);
            if (nfield == -1) nfield = 0; // not found
            E_Int niso = 10;
            if (PyLong_Check(n) == true) niso = PyLong_AsLong(n);
            else niso = E_Int(PyFloat_AsDouble(n));
            if (_nfield > nfield)
            {
              _niso[nfield] = niso;
              _isoMin[nfield] = PyFloat_AsDouble(min);
              _isoMax[nfield] = PyFloat_AsDouble(max);
              _isoAlphaMin[nfield] = -1.e38;
              _isoAlphaMax[nfield] = 1.e38;

              if (nelts == 5)
              {
                PyObject* cmap = PyList_GetItem(l, 4); // colormap for isos
                _isoColormap[nfield] = (int)(PyLong_AsLong(cmap));
              }
              else if (nelts == 6)
              {
                PyObject* amin = PyList_GetItem(l, 4); // alpha min for isos
                PyObject* amax = PyList_GetItem(l, 5); // alpha max for isos
                _isoAlphaMin[nfield] = PyFloat_AsDouble(amin);
                _isoAlphaMax[nfield] = PyFloat_AsDouble(amax);
              }
              else if (nelts == 7)
              {
                PyObject* amin = PyList_GetItem(l, 4); // alpha min for isos
                PyObject* amax = PyList_GetItem(l, 5); // alpha max for isos
                _isoAlphaMin[nfield] = PyFloat_AsDouble(amin);
                _isoAlphaMax[nfield] = PyFloat_AsDouble(amax);
                PyObject* cmap = PyList_GetItem(l, 6); // colormap for isos
                if (PyLong_Check(cmap) == true) _isoColormap[nfield] = (int)(PyLong_AsLong(cmap));
                else _isoColormap[nfield] = (int)(PyFloat_AsDouble(cmap));
              }
            }
          }
        } 
      }
      else
      {
        // liste a plat (ancienne methode)
        E_Int cpt = 0;
        while (cpt < size)
        {
          PyObject* f = PyList_GetItem(isoScales, cpt); // nfield
          PyObject* n = PyList_GetItem(isoScales, cpt+1); // nombre d'isos
          PyObject* min = PyList_GetItem(isoScales, cpt+2); // min for isos
          PyObject* max = PyList_GetItem(isoScales, cpt+3); // max for isos
          E_Int nfield = getScalarField(f);
          if (nfield == -1) nfield = 0; // not found
          E_Int niso = 10;
          if (PyLong_Check(n) == true) niso = PyLong_AsLong(n);
          else niso = E_Int(PyFloat_AsDouble(n));
          if (_nfield > nfield)
          {
            _niso[nfield] = niso;
            _isoMin[nfield] = PyFloat_AsDouble(min);
            _isoMax[nfield] = PyFloat_AsDouble(max);
            _isoAlphaMin[nfield] = -1.e38;
            _isoAlphaMax[nfield] = 1.e38;
          }
          cpt += 4;
        }
      }
    }
  }

  if (bgColor != -1) ptrState->bgColor = bgColor;
  if (strcmp(backgroundFile, "None") != 0) strcpy(ptrState->backgroundFile, backgroundFile);
  if (bgColor >= 6) // requires a background texture
    ptrState->updateBackground = 1;

  if (ghostifyDeactivatedZones != -1)
    ptrState->ghostifyDeactivatedZones = ghostifyDeactivatedZones;

  if (edgifyActivatedZones != -1)
    ptrState->edgifyActivatedZones = edgifyActivatedZones;
  if (edgifyDeactivatedZones != -1)
    ptrState->edgifyDeactivatedZones = edgifyDeactivatedZones;

  if (shadow != -1) ptrState->shadow = shadow;

  if (dof != -1) ptrState->DOF = dof;

  if (strcmp(exportResolution, "None") != 0)
  {
    // look for x
    E_Int s = strlen(exportResolution);
    E_Int i = 0;
    while (i < s)
    {
      if (exportResolution[i] == 'x') break;
      i++;
    }
    if (i != s)
    {
      char number[256];
      for (E_Int j = 0; j < i; j++) number[j] = exportResolution[j];
      number[i] = '\0';
      E_Int w = atoi(number);
      for (E_Int j = i+1; j < s; j++) number[j-i-1] = exportResolution[j];
      number[s-i-1] = '\0';
      E_Int h = atoi(number);
      ptrState->exportWidth = w;
      ptrState->exportHeight = h;
    }
  }

  if (exportAA != -1) ptrState->exportAA = exportAA;

  if (strcmp(exportFile, "None") != 0)
  {
    // Choisi le plugins aproprie
    strcpy(ptrState->exportFile, exportFile);
    char ext[5];
    getFileExt(exportFile, ext);

    struct chain_function_void2* p = _plugins.screenDump;
    while (p != NULL)
    {
      if (strcmp(p->extension, ext) == 0)
      {
        _pref.screenDump = p; break;
      }
      p = p->next;
    }
    ptrState->shootScreen = 1;
  }
}

//=============================================================================
void gidle() // utilise sous mac, car timer ne semble pas marcher
{
  gdisplay();

  // pour la detection du double click
  Data* d = Data::getInstance();
  if (d->ptrState->mouseFirstClicked > 0) d->ptrState->mouseFirstClicked++;
  if (d->ptrState->mouseFirstClicked > 5) d->ptrState->mouseFirstClicked = 0;
}
//=============================================================================
void gtimer(int val)
{
  Data* d = Data::getInstance();
  if (d->ptrState->inmenu == 0)
  {
    // lance le render
    gdisplay();

    // pour la detection du double click
    if (d->ptrState->mouseFirstClicked > 0) d->ptrState->mouseFirstClicked++;
    if (d->ptrState->mouseFirstClicked > 5) d->ptrState->mouseFirstClicked = 0;

    // relance le timer
    glutTimerFunc(d->ptrState->timeStep, gtimer, 0);
  }
}

//=============================================================================
/* Realloue si necessaire les vecteurs globaux de Data et du _state
   dependant de nfield */
//=============================================================================
void Data::reallocNFieldArrays(E_Int nfield)
{
  //printf("allocating field %d %d\n", nfield, _nfield);
  if (nfield > 0 && nfield > _nfield)
  {
    double* n = new double [nfield];
    if (minf != NULL) delete [] minf;
    minf = n;

    n = new double [nfield];
    if (maxf != NULL) delete [] maxf;
    maxf = n;

    n = new double [nfield];
    if (_isoMin != NULL) delete [] _isoMin;
    _isoMin = n;

    n = new double [nfield];
    if (_isoMax != NULL) delete [] _isoMax;
    _isoMax = n;
    
    n = new double [nfield];
    if (_isoAlphaMin != NULL) delete [] _isoAlphaMin;
    _isoAlphaMin = n;
    
    n = new double [nfield];
    if (_isoAlphaMax != NULL) delete [] _isoAlphaMax;
    _isoAlphaMax = n;

    E_Int* ni = new E_Int [nfield];
    for (E_Int i = 0; i < nfield; i++) ni[i] = -1;
    if (_isoColormap != NULL) delete [] _isoColormap;
    _isoColormap = ni;

    n = new double [nfield];
    for (E_Int i = 0; i < nfield; i++) n[i] = 0.;
    if (ptrState->activePointF != NULL) delete [] ptrState->activePointF;
    ptrState->activePointF = n;

    E_Int* m = new E_Int [nfield];
    if (_niso != NULL) delete [] _niso;
    _niso = m;
    for (E_Int i = 0; i < nfield; i++) _niso[i] = -1;

    _nfield = nfield;
  }
}

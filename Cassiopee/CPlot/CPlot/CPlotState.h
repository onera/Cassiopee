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

//=============================================================================
// Keep information on the transcient software state
//=============================================================================
#ifndef _CPLOT_STATE_H_
#define _CPLOT_STATE_H_
#include <pthread.h>
#include <sys/time.h>
#include <cstdio>
#include <vector>
#include "Def/DefTypes.h"

#if defined( _WIN32 ) || defined( _WIN64 )
#include <winsock.h>
#endif

/* Define a chained integer list */
struct chain_int 
{
  E_Int      value;
  chain_int* next;
};

/* Define informations on the state of software */
struct CPlotState 
{
  E_Int mode;  // mode=MESH, SOLID, RENDER, SCALARFIELD, VECTORFIELD
  E_Int dim;   // 1D, 2D, 3D display

  E_Int scalarField;                               // variable number for scalar field
  E_Int vectorField1, vectorField2, vectorField3;  // variables for vector field

  E_Int ijk1D;  // i,j,k mode (1D): i(0), j(1), k(2)
  E_Int var1D;  // 1D secondary variable (inutilise)
  E_Int var2D;  // 2D secondary variable (toggle entre xy, yz, xz, pour dim=2)

  E_Int ghostifyDeactivatedZones;  // display deactivated zones as blended
  E_Int edgifyDeactivatedZones;    // display deactivated zones as edges (NI)
  E_Int edgifyActivatedZones;      // display activated zones as edges (NI)
  E_Int simplifyOnDrag;            // simplify display on drag

  float alpha;            // parameter for alpha blending (1: no alpha blending)
  int currentMenuItem;    // current displayed menu item
  int inmenu;             // we are in menu mode or not?
  int smoke;              // if particles smoke is to be displayed

  E_Int farClip;            // if 1 perform farClipping before display
  E_Int render;             // -1: dont push data, dont flush buffer
                            // 0: push but dont flush
                            // 1, push and flush at end of display
  E_Int bb;                 // if 1, display BB in global display
  E_Int header;             // if 1, display header in global display
  E_Int info;               // if 1, display the info on selected zone
  E_Int texture;            // if 1, display textures
  E_Int menu;               // if 1, display bottom menu
  E_Int isoLegend;          // if 1, display iso legend
  E_Int offscreen;          // 1: offscreen rendering with mesa, 2: with FBO, 
                            // 3: with FBO composite incomplete, 4: with FBO composite finalized
  E_Int frameBuffer;        // numero du frame buffer (0-9)
  char* offscreenBuffer[10];  // buffer for rendering (offscreen)
  float* offscreenDepthBuffer[10]; // Depth buffer storage for offscreen rendering 3 et 4
  void* ctx;
  
  // overlay message
  char* message;

  // Stereo
  E_Int  stereo;        // stereo mode (anaglyph)
  double stereoDist;    // stereo distance

  // Shadow mapping
  E_Int shadow;

  // Light position (for light and shadow)
  double lightOffsetX;
  double lightOffsetY;

  // Depth of field
  E_Int DOF;
  double dofPower;
    
  // Gamma correction
  double gamma;
    
  // Tone Mapping type
  E_Int toneMapping;

  // Sobel threshold
  double sobelThreshold;

  // Sharpen power
  double sharpenPower;

  // ssao power
  double ssaoPower;

  // panorama mode
  E_Int panorama;

  // Last selected zone
  E_Int selectedZone;
  std::vector<E_Int> ambSelections; // keep track of eventual overlapping zones when selecting
  E_Int ambSelSet; // keep track of chosen selected zone in ambSelections;

  // List of deactivated zones
  chain_int* deactivatedZones;
  // Clicked point
  double activePointX, activePointY, activePointZ, activePointZBuf;
  E_Int activePointI, activePointJ, activePointK, activePointL;
  double* activePointF;
  // Pressed keys
  E_Int kkeysActivated;  // if 1, short cut keys are activated, 0: deactivated
  E_Int kcursor;
  unsigned char keys[128];  // store pressed keys

  // Textures
  char shaderPath[1028];   // chemin pour les fichiers + shaders
  char envmapFile[1028];   // fichier image pour les envmaps
  E_Int updateEnvmap;      // dit a display d'updater l'envmap
  E_Int updateBackground;  // dit a display d'updater la texture de background

  // Billboard (material=9)
  E_Int billBoardType;  // type du billboard en memoire (texture)
  E_Int billBoardNi;    // nbre de samples en i
  E_Int billBoardNj;    // nbre de samples en j
  E_Int billBoardD;     // no du champ a utiliser pour la taille des billboards
  E_Int billBoardT;     // no du champ a utiliser pour le choix du billboard
  E_Int billBoardWidth; // nbre de pixel de l'image
  E_Int billBoardHeight; // idem
  double billBoardSize; // Reference size for billboards. If = 0., automatic from cam

  // Autres
  char winTitle[120];        // titre de la fenetre
  char file[120];            // fichier en cours de display (si il y en a)
  char filePath[1028];       // local directory of current file
  char localPathName[1028];  // local directory name (where cassiopee is launched)

  // cursor
  E_Int cursorType;    // current cursor type
  E_Int updateCursor;  // dit a display d'updater le curseur de la souris

  // Position de la souris (no E_Int here)
  int activeMouseButton;  // bouton de la souris presse (dernier click)
  int activeMouseX;       // position de la souris (au dernier click)
  int activeMouseY;
  int currentMouseButton;  // bouton courant
  float currentMousePosX;    // position courante de la souris
  float currentMousePosY;    // quand on est en drag mode (shift + click)
  float currentMousePosZ;
  int ondrag;    // 1 if we are dragging/moving view
  int modifier;  // modifieur quand la souris est pressee

  // Styles
  E_Int meshStyle;       // style pour le mesh mode
  E_Int solidStyle;      // style pour le solid mode
  E_Int scalarStyle;     // style pour le scalar mode
  E_Int vectorStyle;     // style pour le vector mode
  float vectorScale;     // scale pour le vector mode
  float vectorDensity;   // Density of vectors in vector mode
  E_Int vectorNormalize; // Normalize all vectors before display them ( 1 : yes, 0 : no )
  E_Int vectorShowSurface;// Show the triangle emmiting the vector field.
  E_Int vectorShape;     // Shape of the arrow
  E_Int vector_projection; // Project ( 1 ) or not ( 0 ) the vector on the surface of the obstacle
  E_Int selectionStyle;  // style pour la selection (0: bleue, 1: alpha)
  E_Int colormap;        // colormap type
  E_Int colormapSize;    // number of colors in colormap
  float* colormapR; // colormap red
  float* colormapG; // colormap green
  float* colormapB; // colormap blue

  float colormapR1, colormapG1, colormapB1; // starting color for bi/tricolor colormaps
  float colormapR2, colormapG2, colormapB2; // ending color for bi/tricolor colormaps  
  float colormapR3, colormapG3, colormapB3; // mid color for tricolor colormaps

  E_Int isoLight;        // light ou pas light -> isoLight
  E_Int niso;            // nbre d'isos (global)
  float isoEdges;        // edge entre les isos

  // Options pour le high order
  int inner_tesselation; // Degre de raffinement pour les triangles internes
  int outer_tesselation; // Degre de raffinement pour les triangles au bord.

  void clearDeactivatedZones();
  void insertDeactivatedZones(E_Int i);
  void removeDeactivatedZones(E_Int i);
  void printDeactivatedZones();

  CPlotState() : lock(0), _lockGPURes(0), _mustExport(0), _isExporting(0) {}

  virtual ~CPlotState( ) {
      for (E_Int i = 0; i < 10; i++) {
      if (offscreenBuffer[i] != NULL) free(offscreenBuffer[i]); 
      if (offscreenDepthBuffer[i] != NULL) free(offscreenDepthBuffer[i]); } 
      clearDeactivatedZones(); }
    
  // lock=1 pendant le display, les donnees ne doivent alors
  // pas etre modifiees!
  volatile int lock;
  pthread_mutex_t lock_mutex;
  pthread_cond_t unlocked_display;
  void lockDisplay() {
      pthread_mutex_lock(&lock_mutex);
      lock = 1;
      pthread_mutex_unlock(&lock_mutex); }
  void unlockDisplay() {
      pthread_mutex_lock(&lock_mutex);
      lock = 0;
      pthread_cond_signal(&unlocked_display);
      pthread_mutex_unlock(&lock_mutex); }
  void syncDisplay(int time=5) {
      pthread_mutex_lock(&lock_mutex);
      if (lock == 1) pthread_cond_wait(&unlocked_display, &lock_mutex);
      pthread_mutex_unlock(&lock_mutex); }

  int    freeGPURes;      // dit a display de liberer les ressources du GPU
  E_LONG freeGPUData[4];  // Etats de liberation des ressources GPUs, independant de la methode de rendering
  // 0: mode, 1: start ou size (ptr), 2: end ou -1 (ptr),
  // 3: permanent
  int* freeGPUPtr;  // ptr
  virtual void freeGPUResources() = 0;
  pthread_mutex_t gpures_mutex;
  pthread_cond_t unlocked_gpures;
  volatile int _lockGPURes;
  // lockGPURes=display1: On bloque les ressources GPUs qui ne doivent pas etre
  //                      modifiees pendant ce temps
  void lockGPURes() {
      pthread_mutex_lock(&gpures_mutex);
      _lockGPURes = 1;
      pthread_mutex_unlock(&gpures_mutex); }
  void unlockGPURes() {
      pthread_mutex_lock(&gpures_mutex);
      _lockGPURes = 0;
      pthread_cond_signal(&unlocked_gpures);
      pthread_mutex_unlock(&gpures_mutex); }
  virtual void syncGPURes(int dt = 5) {
      pthread_mutex_lock(&gpures_mutex);
      if (_lockGPURes == 1) pthread_cond_wait(&unlocked_gpures, &gpures_mutex);
      pthread_mutex_unlock(&gpures_mutex); }

  E_Int timeStep;           // temps en ms entre 2 appels de display
  E_Int mouseFirstClicked;  // true if mouse has been clicked already before
                            // redisplay

  // Export
  char exportFile[1028];            // nom du fichier pour l'export
  E_Int exportWidth, exportHeight;  // resolution de l'export
  E_Int exportAA;                   // antialiasing at export (1 or 0)
  E_Int continuousExport;           // =0, pas d'export continue, =1, export continue
  FILE* ptrFile;                    // ptr fichier pour les exports continus
  void* context;                    // contexte pour les exports continus
  E_Int exportNumber;               // no ajoute au nom du fichier
  E_Int shootScreen;                // dit a display de shooter
  pthread_mutex_t export_mutex;
  volatile int _mustExport, _isExporting;
  pthread_cond_t  unlocked_export;

  bool odsRun; // true if ODS run
  E_Int odsNSlits; // the number of slits
  E_Int odsSlit; // the number of current slit for ODS
  char* odsImage; // pointer to ODS final image
  char* odsFrontImage; // pointer to ods front slit image
  char* odsTopImage; // pointer to ods top slit image

  // Others
  E_Int fullScreen;  // 1: full screen, 0: window
  E_Int bgColor;     // background color
  char backgroundFile[1028]; // image file name for background (bgColor=11)
  E_Int autoblank;   // if 1, cplot try to set automatically the
                     // blanking plugin

  // 1D overlay display size
  E_Int gridSizeI;  // nombre de slots 1D affichables sur le display
  E_Int gridSizeJ;

  // animated shaders
  E_Int ktime;   // the current time, incremented at each display
  E_Int ktimer;  // if ktimer!=0, redisplay is triggered, timer
                 // is decremented at each display
};

// Init software state
void initState( );
#endif

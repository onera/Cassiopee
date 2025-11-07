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

// All data for one shot display
#ifndef _CPLOT_DATA_H_
#define _CPLOT_DATA_H_

#define _CARD0_
#define MAXSTRINGLENGTH 128  /* max length of strings */

// Les differents mode
#define MESH 0
#define SOLID 1
#define RENDER 2
#define SCALARFIELD 3
#define VECTORFIELD 4

// cutoff pour les isos
#define ISOCUTOFF 1.e-30

#ifdef __SHADERS__
#include "GL/glew.h"
#include "ShaderManager.h"
#include "TesselationShaderManager.hpp"
#endif

#ifndef __SHADERS__
// La fonction glActiveTexture est une extension ARB
#define glActiveTexture(x)
#endif 

#ifdef __MESA__
#define GLAPI extern
#ifdef _WIN32
  #define APIENTRY __stdcall
#else
  #define APIENTRY
#endif
#include <GL/osmesa.h>
#endif

#ifdef _MPI
#include <mpi.h>
#endif

#ifdef _WIN32
#define GLAPI extern
#endif

#include "GL/glut.h"
#include "Zone.h"
#include "StructZone.h"
#include "UnstructZone.h"
#include "Preferences.h"
#include "Plugins/Plugins.h"
#include "ViewInfo.h"
#include "CPlotState.h"
#include "Particles.h"
#include "Functions.h"
#include "Slot1D.h"
#include "Zone1D.h"
#include "Fonts/OpenGLText.h"
class Data
{   
  public:
    Data(CPlotState* pt_state);
    virtual ~Data();

  // le getInstance ici n'est valable que si deja appele dans les
  // classes filles :
  static Data* getInstance(); // <--- put in DataDL and DataVBO class.
  enum RenderID { Direct=0, VBO=1, DL=2, END_GUARD };
  static RenderID _renderID;
  static Data* _instance;

  Preferences _pref;
  Plugins _plugins;
  ViewInfo _view;
  CPlotState* ptrState;
  volatile int _CDisplayIsLaunched;
    
  E_Int _numberOfZones; // all zones
  Zone** _zones;
    
  E_Int _numberOfStructZones; // struct zones
  StructZone** _szones;

  E_Int _numberOfUnstructZones; // unstruct zones
  UnstructZone** _uzones;

  // Global data
  double xmin, ymin, zmin, xmax, ymax, zmax; // global BB
  double dmoy; // taille moyenne de la BB

  double epsup;     // minimum up
  double epsstrafe; // minimum strafe

  E_Int _niter; // iteration (si calcul) ou transmis
  E_Int _winId; // main window id
  PyThreadState* _save;

  // Field data
  E_Int _nfield; // dim of variables field arrays
  double* minf; // global min-max
  double* maxf;
  double* _isoMin; // min pour chaque champ (si specifie)
  double* _isoMax; // max pour chaque champ
  double* _isoAlphaMin; // min pour alpha iso
  double* _isoAlphaMax; // min pour alpha iso
  E_Int* _niso; // nbre d'iso pour chaque champ (si specifie)
  E_Int* _isoColormap; // colormap pour chaque champ (si specifie)

  // Colormap data
  E_Int _colormapSizeViridis; // number of colors in colormap
  float* _colormapRViridis; // colormap red
  float* _colormapGViridis; // colormap green
  float* _colormapBViridis; // colormap blue
  E_Int _colormapSizeInferno; // number of colors in colormap
  float* _colormapRInferno; // colormap red
  float* _colormapGInferno; // colormap green
  float* _colormapBInferno; // colormap blue
  E_Int _colormapSizeMagma; // number of colors in colormap
  float* _colormapRMagma; // colormap red
  float* _colormapGMagma; // colormap green
  float* _colormapBMagma; // colormap blue
  E_Int _colormapSizePlasma; // number of colors in colormap
  float* _colormapRPlasma; // colormap red
  float* _colormapGPlasma; // colormap green
  float* _colormapBPlasma; // colormap blue
  E_Int _colormapSizeJet; // number of colors in colormap
  float* _colormapRJet; // colormap red
  float* _colormapGJet; // colormap green
  float* _colormapBJet; // colormap blue
  E_Int _colormapSizeGreys; // number of colors in colormap
  float* _colormapRGreys; // colormap red
  float* _colormapGGreys; // colormap green
  float* _colormapBGreys; // colormap blue
  E_Int _colormapSizeNiceBlue; // number of colors in colormap
  float* _colormapRNiceBlue; // colormap red
  float* _colormapGNiceBlue; // colormap green
  float* _colormapBNiceBlue; // colormap blue
  E_Int _colormapSizeGreens; // number of colors in colormap
  float* _colormapRGreens; // colormap red
  float* _colormapGGreens; // colormap green
  float* _colormapBGreens; // colormap blue
  
  // Texture data
  GLuint _texNodes; // texture pour les nodes
  GLuint _texNoise; // texture pour le noise
  GLuint _texBillBoard; // texture for the billboard
  GLuint _texBackground; // texture for background
  GLuint _texColormap; // texture for colormap
  E_Int _texColormapType; // type stored in colormap
  double _texColormapMinMax; // Min max of colormap stored 
  E_Int _frameBufferSize[10]; // size of frame buffer
  GLuint _texFrameBuffer[10]; // texture for frame buffer
  GLuint _texEnviron1; // texture environnement 1
  E_Int _voxelBufferSize; // size of voxel buffer
  GLuint _texVoxelBuffer; // texture pour le voxel buffer
  GLuint _texLeft; // texture for left eye rendered frame (anaglyph)
  GLuint _texRight; // texture for right eye rendered frame (anaglyph)
  GLuint _texSupp; // texture supp pour le post-processing
  GLuint _shadowMap; // texture pour le shadow mapping
  double _bias[16];
  double _lightModelView[16];
  double _lightProjection[16];

  // billboards image files and texture storage
  E_Int _nBillBoards; 
  char** _billBoardFiles;
  E_Int* _billBoardNis;
  E_Int* _billBoardNjs;
  E_Int* _billBoardWidths;
  E_Int* _billBoardHeights;
  GLuint* _billBoardTexs;

  // Material image files and texture storage
  E_Int _nMaterials; 
  char** _materialFiles;
  E_Int* _materialWidths;
  E_Int* _materialHeights;
  GLuint* _materialTexs;

  // BumpMaps image files and texture storage
  E_Int _nBumpMaps; // must be equal to nMaterials, but accepts NULL
  char** _bumpMapFiles;
  E_Int* _bumpMapWidths;
  E_Int* _bumpMapHeights;
  GLuint* _bumpMapTexs;

#ifdef __SHADERS__
  CPlot::ShaderManager _shaders; // Shaders
  void triggerShader(Zone& z, int material, float scale, float* color);
#endif
  
  // Data for overlay plot display
  std::vector<Slot1D*> _slots1D;

public:
  // Init openGL
  void init();
  void setBgColor();
  
  // Init zone data
  virtual E_Int initZoneData(std::vector<K_FLD::FldArrayF*>& structF,
                std::vector<char*>& structVarString,
                std::vector<E_Int>& nit,
                std::vector<E_Int>& njt, 
                std::vector<E_Int>& nkt,
                std::vector<K_FLD::FldArrayF*>& unstrF,
                std::vector<char*>& unstrVarString,
                std::vector<K_FLD::FldArrayI*>& cnt,
                std::vector<char*>& eltType,
                std::vector<char*>& zoneNames,
                std::vector<char*>& zoneTags,
                E_Int referenceNfield=-1,
                char** referenceVarNames=NULL);
  StructZone* createStructZone(K_FLD::FldArrayF* structF, char* varString,
             E_Int posx, E_Int posy, E_Int posz,
             E_Int ni, E_Int nj, E_Int nk,
             char* zoneName, char* zoneTags,
             E_Int referenceNfield=-1, char** referenceVarNames=NULL,
             E_Int mustComplete=0);
  UnstructZone* createUnstrZone(K_FLD::FldArrayF* unstrF, char* varString,
                E_Int posx, E_Int posy, E_Int posz,
                K_FLD::FldArrayI* cn, char* eltType,
                char* zoneName, char* zoneTags, 
                E_Int referenceNfield=-1, char** referenceVarNames=NULL,
                E_Int mustComplete=0);
  void reallocNFieldArrays(E_Int nfield);
  // Init _state
  virtual void initState();
  // Init camera
  void initCam();
  void init2DCam();
  void init1DCam();
  // Load preferences
  void loadPrefs();
  // enforce data
  void enforceGivenData(E_Int dim, E_Int mode, E_Int scalarField,
            E_Int vectorField1, E_Int vectorField2, E_Int vectorField3,
            E_Int displayBB, E_Int displayInfo,
            E_Int displayIsoLegend);
  void enforceGivenData2(float xcam, float ycam, float zcam,
             float xeye, float yeye, float zeye,
             float dirx, float diry, float dirz,
             float viewAngle,
             E_Int meshStyle, E_Int solidStyle,
             E_Int scalarStyle, E_Int vectorStyle, 
             float vectorScale, float vectorDensity, E_Int vectorNormalize,
             E_Int vectorShowSurface, E_Int vectorShape, E_Int vector_projection,
             E_Int colormap, 
             char* colormapC1, char* colormapC2, char* colormapC3, PyObject* colormapC,
             E_Int niso, float isoEdges, PyObject* isoScales,
             E_Int bgColor, char* backgroundFile,
             E_Int ghostifyDeactivatedZones,
             E_Int edgifyActivatedZones, 
             E_Int edgifyDeactivatedZones,
             E_Int shadow, E_Int dof,
             char* exportFile, char* exportResolution);
  void rgb2hsv(float r, float g, float b, float& h, float& s, float& v);
  void hsv2rgb(float h, float s, float v, float& r, float& g, float& b);
  void colorString2RGB(char* color, float& colorR, float& colorG, float& colorB);
  void initViridis();
  void initInferno();
  void initMagma();
  void initPlasma();
  void initJet();
  void initGreys();
  void initGreens();
  void codeFromRenderTag(Zone& z, char* tag, 
             float& colorR, float& colorG, float& colorB,
             E_Int& material, double& blending, E_Int& meshOverlay,
             float& meshColorR, float& meshColorG, float& meshColorB, float& meshWidth,
             float& shaderParam1, float& shaderParam2);
  void getAllVars(std::vector<char*>& structVarString,
                  std::vector<char*>& unstrVarString,
                  E_Int& referenceNfield,
                  char**& referenceVarNames);
  void replaceVolumetricZones();
  
  // openGfx
  void openGfx();
  void closeGfx();

  // Misc
  E_Int findBlockContaining(double x, double y, double z,
                            E_Int& zone, E_Int& ind, E_Int& indE, double& dist, E_Int& ncon);

  // Create textures
  E_Int createNodeTexture();
  E_Int createNoise3DTexture();
  E_Int createColormapTexture();
  void fillColormapTexture(E_Int type);
  E_Int createFrameBufferTexture();
  E_Int createImageTexture(const char* filename, GLuint &tex,
                           E_Int& width, E_Int& height, 
                           bool mipmap=true);
  E_Int createPngTexture(const char* filename, GLuint &tex,
                         E_Int& width, E_Int& height, 
                         bool mipmap=true);
  E_Int createJpgTexture(const char* filename, GLuint &tex,
                         E_Int& width, E_Int& height, 
                         bool mipmap=true);
  E_Int createVoxelTexture();
  void voxelize(UnstructZone& zn, UnstructZone& z);
  void voxelize(StructZone& zn, StructZone& z);

  // Keys
  void keyboard(unsigned char key, E_Int x, E_Int y);
  void arrows(unsigned char key, E_Int x, E_Int y);
  void moveDown(double alpha, double dx, double dy, double dz, double d,
                double dirx, double diry, double dirz);
  void strafeDown(double alpha, double dx, double dy, double dz, double d,
                  double dirx, double diry, double dirz);
  void moveUp(double alpha, double dx, double dy, double dz, double d,
              double dirx, double diry, double dirz);
  void strafeUp(double alpha, double dx, double dy, double dz, double d,
                double dirx, double diry, double dirz);
  void moveRight(double alpha, double dx, double dy, double dz, double d,
                 double dirx, double diry, double dirz);
  void strafeRight(double alpha, double dx, double dy, double dz, double d,
                   double dirx, double diry, double dirz);
  void moveLeft(double alpha, double dx, double dy, double dz, double d,
                double dirx, double diry, double dirz);
  void strafeLeft(double alpha, double dx, double dy, double dz, double d,
                  double dirx, double diry, double dirz);
  void tiltRight(double alpha, double dx, double dy, double dz, double d,
                 double dirx, double diry, double dirz);
  void tiltLeft(double alpha, double dx, double dy, double dz, double d,
                double dirx, double diry, double dirz);
  void rotateHeadRight(double alpha, double dx, double dy, double dz, 
                       double d, double dirx, double diry, double dirz);
  void rotateHeadLeft(double alpha, double dx, double dy, double dz, 
                      double d,
                      double dirx, double diry, double dirz);
  void rotateHeadUp(double alpha, double dx, double dy, double dz, double d,
                    double dirx, double diry, double dirz);
  void rotateHeadDown(double alpha, double dx, double dy, double dz, 
                      double d,
                      double dirx, double diry, double dirz);
  void changeIPlanePlus();
  void changeJPlanePlus();
  void changeKPlanePlus();
  void changeIPlaneMinus();
  void changeJPlaneMinus();
  void changeKPlaneMinus();
  void changeSecondaryVariablePlus();
  void changeSecondaryVariableMinus();
  void changeBlankingFunction();
  void changeAppearance();
  void changeAmbSelection();
  void mouseButton(E_Int button, E_Int etat, E_Int x, E_Int y);
  void mouseMotion(E_Int x, E_Int y);
  void mousePassiveMotion(E_Int x, E_Int y);
  void mouseDrag(E_Int x, E_Int y);

  // Local display
  void fog();
  void light(E_Int type);
  void noLight();
  void setCursor(E_Int type);
  double dist2BB(double x, double y, double z,
    double xmin, double ymin, double zmin,
    double xmax, double ymax, double zmax);
  void computeSteps0(StructZone* zonep, 
                     E_Int& stepi, E_Int& stepj, E_Int& stepk);
  void computeSteps1(StructZone* zonep, 
                     E_Int& stepi, E_Int& stepj, E_Int& stepk);
  void computeSteps(StructZone* zonep, 
                    E_Int& stepi, E_Int& stepj, E_Int& stepk);
  void activateZone();
  void deactivateZone();
  void clearDisplay();
  virtual void createGPURes() = 0;
  virtual void createIsoGPURes(E_Int nofield) = 0; // scalaire
  virtual void createIsoGPURes(E_Int nofield1, E_Int nofield2, E_Int nofield3) = 0; // vecteur
  virtual void createIsoGPUResForRender() = 0; // scalaire for render mode
  virtual void freeGPURes(int mode, int start, int end, int permanent) = 0;
  virtual void freeGPURes(int mode, int size, int* ptr, int permanent) = 0;
  void display();
  void displayBB();
  void displayBB2();
  void displayFrameTex(E_Int mode, double sobelThreshold=-0.5);
  void displayAnaglyph();
  void displayActivePoint();
  void displaySEdges();
  void displayUEdges();
  void displaySBBZone(StructZone* z);
  void displayUBBZone(UnstructZone* z);
  virtual void displaySMesh() = 0;
  void displaySMeshZone(StructZone* zonep, E_Int zone);
  virtual void displayUMesh() = 0;
  void displayUMeshZone(UnstructZone* zonep, E_Int zone, E_Int zonet);
  void displayUMeshZone_ho(UnstructZone* zonep, E_Int zone, E_Int zonet);
  virtual void displaySSolid() = 0;
  void displaySSolidZone(StructZone* zonep, E_Int zone);
  virtual void displayUSolid() = 0;
  void displayUSolidZone(UnstructZone* zonep, E_Int zone, E_Int zonet);
  void displayUSolidHOZone(UnstructZone* zonep, E_Int zone, E_Int zonet);
  virtual void displaySIsoSolid() = 0;
  void displaySIsoSolidZone(StructZone* zonep, E_Int zone, E_Int nofield);
  void displaySIsoSolidZone(StructZone* zonep, E_Int zone, E_Int nofield1,
                            E_Int nofield2, E_Int nofield3);
  virtual void displayUIsoSolid() = 0;
  void displayUIsoSolidZone(UnstructZone* zonep, E_Int zonet, E_Int nofield);
  void displayUIsoSolidZone(UnstructZone* zonep, E_Int zonet, 
                            E_Int nofield1, E_Int nofield2, E_Int nofield3);
  void displayNodes();
  void displayBillBoards(Zone* zonep, E_Int zone);
  void displayAllBillBoards();
  void createGPUSMeshZone(StructZone* zonep, E_Int zone);
  virtual void renderGPUSMeshZone(StructZone* zonep, E_Int zone) = 0;
  void createGPUUMeshZone(UnstructZone* zonep, E_Int zone, E_Int zonet);
  virtual void renderGPUUMeshZone(UnstructZone* zonep, E_Int zone, E_Int zonet) = 0;
  void createGPUSSolidZone(StructZone* zonep, E_Int zone);
  virtual void renderGPUSSolidZone(StructZone* zonep, E_Int zone) = 0;
  virtual void createGPUUSolidZone(UnstructZone* zonep, E_Int zone, E_Int zonet) = 0;
  virtual void renderGPUUSolidZone(UnstructZone* zonep, E_Int zone, E_Int zonet) = 0;
  virtual void createGPUUSolidHOZone(UnstructZone* zonep, E_Int zone, E_Int zonet) = 0;
  virtual void renderGPUUSolidHOZone(UnstructZone* zonep, E_Int zone, E_Int zonet) = 0;
  virtual void createGPUSIsoSolidZone(StructZone* zonep, E_Int zone, E_Int nofield) = 0;
  virtual void createGPUSIsoSolidZone(StructZone* zonep, E_Int zone, E_Int nofield1,
                                      E_Int nofield2, E_Int nofield3) = 0;
  virtual void createGPUUIsoSolidZone(UnstructZone* zonep, E_Int zone, E_Int zonet, E_Int nofield) = 0;
  virtual void createGPUUIsoSolidZone(UnstructZone* zonep, E_Int zone, E_Int zonet, 
                                      E_Int nofield1, E_Int nofield2, E_Int nofield3) = 0;
  virtual void renderSIsoSolidZone(StructZone* zonep, E_Int zone, E_Int nofield) = 0;
  virtual void renderSIsoSolidZone(StructZone* zonep, E_Int zone, E_Int nofield1,
                                   E_Int nofield2, E_Int nofield3) = 0;
  virtual void renderUIsoSolidZone(UnstructZone* zonep, E_Int zonet, E_Int nofield) = 0;
  virtual void renderUIsoSolidZone(UnstructZone* zonep, E_Int zonet, E_Int nofield1,
                                   E_Int nofield2, E_Int nofield3) = 0;

  // roll
  void roll3Dto2D();
  void roll2Dto1D();
  void roll1Dto3D();
  void roll1Dto2D();
  void roll2Dto3D();
  void roll3Dto1D();
  void rollto2Dws();
  void rollto1Dws();

  // Text
  E_Int textWidth(E_Int fontSize, char* string);
  E_Int textWidth1(E_Int fontSize, char* string);
  E_Int textWidth2(E_Int fontSize, char* string);
  E_Int textHeight(E_Int fontSize);
  void renderBitmapString(float x, float y, float z,
                          E_Int fontSize, char *string,
                          float colorR=1., float colorG=1., float colorB=1., float colorA=1.,
                          float nx=0., float ny=1., float nz=0.,
                          float r=1.);
  void* getGlutFont(E_Int fontSize);
  void renderBitmapString1(float x, float y, float z,
                           E_Int fontSize, char *string,
                           float colorR=1., float colorG=1., float colorB=1., float colorA=1.,
                           float nx=0., float ny=1., float nz=0.,
                           float r=1.);
  OpenGLText* getOpenGLText(E_Int fontSize);
  void renderBitmapString2(float x, float y, float z,
                           E_Int fontSize, char *string,
                           float colorR=1., float colorG=1., float colorB=1., float colorA=1.,
                           float nx=0., float ny=1., float nz=0.,
                           float r=1.);
                            
  void renderStringWithShadow(float x, float y, float z,
      E_Int fontSize, char *myString,
      float fgColorR, float fgColorG,
      float fgColorB, float fgColorA,
      float shColorR, float shColorG, 
      float shColorB, float shColorA,
      double offtx=1., double offty=0., double offtz=0.,
      double offnx=0., double offny=1., double offnz=0., 
      double r=1.);
  void resetPerspectiveProjection();
  void setOrthographicProjection();
  void displayText(char* text);
  void displayBigText(E_Int posx, E_Int posy, char* text);
  void displaySmallText(E_Int posx, E_Int posy, char* text);
  void printHeader();
  void printTmpMessage(const char* text);
  void displayInfoWindow(char* text, E_Int l);
  void displayInfo();

  // Legend
  void displayIsoLegend(E_Int dir);

  // 3D Axis
  void displayAxis();

  // 1D
  void displayPlots();
  void displayPlot(Slot1D* s);
  void plotZone(Slot1D* s, Zone1D* z, E_Float posx, E_Float posy,
                E_Float dx, E_Float dy, E_Int var1, E_Int var2);
  E_Float getTick(E_Float rmin, E_Float rmax);
  void plot1DAxis(Slot1D* s, E_Float posx, E_Float posy,
                  E_Float dx, E_Float dy, E_Float blend);
  void getCharsFromVarName(char* varName, char& c1, char& c2);
  E_Int getActivePointIndex(Zone1D* z, E_Int var1, E_Int var2,
                            E_Int& e1, E_Int& e2, double& alpha);
  E_Int display1DActivePoint(Slot1D* s, Zone1D* z, E_Float posx, E_Float posy,
                             E_Float dx, E_Float dy, E_Int var1, E_Int var2);
  E_Int link2View(Zone1D* z, E_Int var1, E_Int var2,
                  E_Float& r1min, E_Float& r1max, 
                  E_Float& r2min, E_Float& r2max);

  // Menu
  void displayMenu();
  void displayDimensionMenu(E_Int* x);
  void displayVariableMenu(E_Int* x);
  void displayAxisMenu(E_Int* x);
  void displayZoneMenu(E_Int* x);
  void menu();
  void menuKeyboard(unsigned char key, int x, int y);
  void menuArrows(int key, int x, int y);
  void displayMenuString(E_Int no, char* msg, E_Int* l, E_Int* sizeMax,
                         E_Int type, void* value);
  void changeMenuItemRight();
  void changeMenuItemLeft();
  void changeMenuItemSpaceBar();

  // Plugins
  void loadPlugins();
  void loadPluginsPath();
  void autoPlugins();
  E_Int checkVariable(E_Int zone, const char* varName);
  void findBlankedZones();
  void dumpWindow();
  char* export2Image(E_Int exportWidth, E_Int exportHeight);
  void superSample(E_Int w, E_Int h, char* im1, char* im2, E_Int factor);
  void gaussianBlur(E_Int w, E_Int h, char* im1, char* im2, E_Int r, double eps);
  void mixImages(E_Int w, E_Int h, char* im1, char* im2, 
                 double alpha, double beta);
  void sharpenImage(E_Int w, E_Int h, char* im1, char* im2, double amount,
                    E_Int radius, E_Int threshold);
  void localBlur(E_Int w, E_Int h, char* im1, char* im2);
  void specPostProcess(char* in, E_Int ni, E_Int nj, float* depth, char* out);
  FILE* fopenw(const char* path, const char* mode);                                                              
  void exportFile();
  void finalizeExport();
  void dataMouseClickSelect(E_Int button, E_Int etat, E_Int x, E_Int y, 
                            E_Int multiple, E_Int accurate);
  void dataMouseRightClickSelect(E_Int button, E_Int etat, E_Int x, E_Int y);

  // Clipping
  void farClipping();
  void closeClipping();
  void veryCloseClipping();
  void veryVeryCloseClipping();
  void adaptiveClipping(double d);

  // Tools
  inline E_Int isBlanked(Zone& zone, E_Int zonet, E_Int n1, E_Int n2, E_Int n3) 
  {
    E_Int ret1 = _pref.blanking->f(this, n1, zone.blank, zonet);
    E_Int ret2 = _pref.blanking->f(this, n2, zone.blank, zonet);
    E_Int ret3 = _pref.blanking->f(this, n3, zone.blank, zonet);
    return (ret1*ret2*ret3 != 0 ? 1 : 0);
  }

  inline E_Int isBlanked(Zone& zone, E_Int zonet, E_Int n1, E_Int n2, E_Int n3, E_Int n4) 
  {
    E_Int ret1 = _pref.blanking->f(this, n1, zone.blank, zonet);
    E_Int ret2 = _pref.blanking->f(this, n2, zone.blank, zonet);
    E_Int ret3 = _pref.blanking->f(this, n3, zone.blank, zonet);
    E_Int ret4 = _pref.blanking->f(this, n3, zone.blank, zonet);
    return (ret1*ret2*ret3*ret4 != 0 ? 1 : 0);
  }

  E_Int numberOfNonBlankedCells(UnstructZone& zone, E_Int zonet);

  E_Int _font1Size; // size of font1
  E_Int _font2Size; // size of font2
  E_Int _font3Size; // size of font3
  OpenGLText* _oglText1; // texture font ptr
  OpenGLText* _oglText2; // texture font ptr
  OpenGLText* _oglText3; // texture font ptr
    
 protected:
  virtual void freeGPUResources(int mode, int start, int end, int permanent);
  virtual void updateGPUResources(int mode, int size, int permanent, void* updatedPointer);
  virtual ZoneImpl* createZoneImpl() = 0;
};

// global functions (no E_Int here to be compatible with glut)
void gdisplay();
void fdisplay();
void reshape(int w, int h);
void gkeyboard(unsigned char key, int x, int y);
void gkeyboardup(int key, int x, int y);
void garrows(int key, int x, int y);
void gidle();
void gtimer(int val);
void gmouseButton(int button, int etat, int x, int y);
void gmouseMotion(int x, int y);
void gmousePassiveMotion(int x, int y);
void gmenuKeyboard(unsigned char key, int x, int y);
void gmenuArrows(int key, int x, int y);

// frustum
void computeFrustumPlanes(ViewInfo& view);
//void computeFrustumPlanes2();
int isInFrustum(Zone* z, ViewInfo& view);

// files
void getFileExt(char* file, char* ext);
#endif

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

// Les differentes fonts
#define FONT1 GLUT_BITMAP_HELVETICA_12
#define FONTSIZE1 12
#define FONT2 GLUT_BITMAP_HELVETICA_18
#define FONTSIZE2 18
#define FONT3 GLUT_BITMAP_HELVETICA_10
#define FONTSIZE3 10

#ifdef __SHADERS__
#include "GL/glew.h"
#include "Shaders/ShaderManager.h"
#include "Shaders/TesselationShaderManager.hpp"
#endif

#ifndef __SHADERS__
// La fonction glActiveTexture est une extension ARB
#define glActiveTexture(x)
#endif 

#ifdef __MESA__
#define GLAPI extern
#include <GL/osmesa.h>
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

class Data
{   
  public:
    Data(CPlotState* pt_state);
    virtual ~Data();

  // le getInstance ici n'est valable que si deja appele dans les
  // classes filles :
  static Data* getInstance(); // <--- put in DataDL and DataVBO class.
  enum RenderID {
    Direct = 0, VBO = 1, DL = 2, END_GUARD
  };
  static RenderID _renderID;
  static Data* _instance;

  Preferences _pref;
  Plugins _plugins;
  ViewInfo _view;
  CPlotState* ptrState;
  volatile int _CDisplayIsLaunched;
    
  int _numberOfZones; // all zones
  Zone** _zones;
    
  int _numberOfStructZones; // struct zones
  StructZone** _szones;

  int _numberOfUnstructZones; // unstruct zones
  UnstructZone** _uzones;

  double xmin, ymin, zmin, xmax, ymax, zmax; // global BB
  double dmoy; // taille moyenne de la BB

  int _nfield; // dim of variables field arrays
  double* minf; // global min-max
  double* maxf;
  double* _isoMin; // min pour chaque champ (si specifie)
  double* _isoMax; // max pour chaque champ
  int* _niso; // nbre d'iso pour chaque champ (si specifie)

  double epsup;     // minimum up
  double epsstrafe; // minimum strafe
  
  int _niter; // iteration (si calcul) ou transmis
  int _winId; // main window id
  PyThreadState* _save;
  
  GLuint _texNodes; // texture pour les nodes
  GLuint _texNoise; // texture pour le noise
  GLuint _texBillBoard; // texture for the billboard
  GLuint _texBackground; // texture for background
  GLuint _texColormap; // texture for colormap
  int _texColormapType; // type stored in colormap
  int _frameBufferSize; // size of frame buffer
  GLuint _texFrameBuffer; // texture for frame buffer
  GLuint _texEnviron1; // texture environnement 1
  int _voxelBufferSize; // size of voxel buffer
  GLuint _texVoxelBuffer; // texture pour le voxel buffer
  GLuint _texLeft; // texture for left eye rendered frame (anaglyph)
  GLuint _texRight; // texture for right eye rendered frame (anaglyph)
  GLuint _texSupp; // texture supp pour le post-processing
  GLuint _shadowMap; // texture pour le shadow mapping
  double _bias[16];
  double _lightModelView[16];
  double _lightProjection[16];
    
  // billboards image files and texture storage
  int _nBillBoards; 
  char** _billBoardFiles;
  int* _billBoardNis;
  int* _billBoardNjs;
  int* _billBoardWidths;
  int* _billBoardHeights;
  GLuint* _billBoardTexs;

  // Material image files and texture storage
  int _nMaterials; 
  char** _materialFiles;
  int* _materialWidths;
  int* _materialHeights;
  GLuint* _materialTexs;

  // BumpMaps image files and texture storage
  int _nBumpMaps; // must be equal to nMaterials, but accepts NULL
  char** _bumpMapFiles;
  int* _bumpMapWidths;
  int* _bumpMapHeights;
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
  virtual int initZoneData( std::vector<K_FLD::FldArrayF*>& structF,
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
  void reallocNFieldArrays(int nfield);
  // Init _state
  virtual void initState();
  // Init camera
  void initCam();
  void init2DCam();
  void init1DCam();
  // Load preferences
  void loadPrefs();
  // enforce data
  void enforceGivenData(int dim, int mode, int scalarField,
			int vectorField1, int vectorField2, int vectorField3,
			int displayBB, int displayInfo, 
			int displayIsoLegend);
  void enforceGivenData2(float xcam, float ycam, float zcam,
			 float xeye, float yeye, float zeye,
			 float dirx, float diry, float dirz,
			 float viewAngle,
			 int meshStyle, int solidStyle, 
			 int scalarStyle, int vectorStyle, float vectorScale, float vectorDensity, int vectorNormalize,
       int vectorShowSurface, int vectorShape, int vector_projection,
       int colormap, int niso, float isoEdges, PyObject* isoScales,
			 int bgColor, int ghostifyDeactivatedZones,
			 int edgifyActivatedZones, 
			 int edgifyDeactivatedZones,
			 int shadow, int dof,
			 char* exportFile, char* exportResolution);
  void codeFromRenderTag(Zone& z, char* tag, 
			 double& colorR, double& colorG, double& colorB,
			 int& material, double& blending, int& meshOverlay,
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
    int findBlockContaining(double x, double y, double z,
                            int& zone, int& ind, int& indE, double& dist);

    // Create textures
    int createNodeTexture();
    int createNoise3DTexture();
    int createColormapTexture();
    void fillColormapTexture(int type);
    int createFrameBufferTexture();
    int createPngTexture(const char* filename, GLuint &tex,
                         int& width, int& height, 
                         bool mipmap=true);
    int createVoxelTexture();
    void voxelize(UnstructZone& zn, UnstructZone& z);
    void voxelize(StructZone& zn, StructZone& z);

    // Keys
    void keyboard(unsigned char key, int x, int y);
    void arrows(int key, int x, int y);
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
                         double d,
                         double dirx, double diry, double dirz);
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
    void mouseButton(int button, int etat, int x, int y);
    void mouseMotion(int x, int y);
    void mousePassiveMotion(int x, int y);
    void mouseDrag(int x, int y);

    // Local display
    void fog();
    void light(int type);
    void noLight();
    void setCursor(int type);
    double dist2BB(double x, double y, double z,
    double xmin, double ymin, double zmin,
    double xmax, double ymax, double zmax);
    void computeSteps0(StructZone* zonep, 
                       int& stepi, int& stepj, int& stepk);
    void computeSteps1(StructZone* zonep, 
                       int& stepi, int& stepj, int& stepk);
    void computeSteps(StructZone* zonep, 
                      int& stepi, int& stepj, int& stepk);
    void activateZone();
    void deactivateZone();
    void clearDisplay();
    virtual void createGPURes() = 0;
    virtual void createIsoGPURes(int nofield) = 0; // scalaire
    virtual void createIsoGPURes(int nofield1, int nofield2, int nofield3) = 0; // vecteur
    virtual void createIsoGPUResForRender() = 0; // scalaire for render mode
    virtual void freeGPURes(int mode, int start, int end, int permanent) = 0;
    virtual void freeGPURes(int mode, int size, int* ptr, int permanent) = 0;
    void display();
    void displayBB();
    void displayBB2();
    void displayFrameTex(int mode, double sobelThreshold=-0.5);
    void displayAnaglyph();
    void displayActivePoint();
    void displaySEdges();
    void displayUEdges();
    virtual void displaySMesh() = 0;
    void displaySMeshZone(StructZone* zonep, int zone);
    virtual void displayUMesh() = 0;
    void displayUMeshZone(UnstructZone* zonep, int zone, int zonet);
    void displayUMeshZone_ho(UnstructZone* zonep, int zone, int zonet);
    virtual void displaySSolid() = 0;
    void displaySSolidZone(StructZone* zonep, int zone);
    virtual void displayUSolid() = 0;
    void displayUSolidZone(UnstructZone* zonep, int zone, int zonet);
    void displayUSolidHOZone(UnstructZone* zonep, int zone, int zonet);
    virtual void displaySIsoSolid() = 0;
    void displaySIsoSolidZone(StructZone* zonep, int zone, int nofield);
    void displaySIsoSolidZone(StructZone* zonep, int zone, int nofield1,
                              int nofield2, int nofield3);
    virtual void displayUIsoSolid() = 0;
    void displayUIsoSolidZone(UnstructZone* zonep, int zonet, 
                              int nofield);
    void displayUIsoSolidZone(UnstructZone* zonep, int zonet, 
                              int nofield1, int nofield2, int nofield3);
    void displayNodes();
    void displayBillBoards(Zone* zonep, int zone);
    void displayAllBillBoards();
    void createGPUSMeshZone(StructZone* zonep, int zone);
    virtual void renderGPUSMeshZone(StructZone* zonep, int zone) = 0;
    void createGPUUMeshZone(UnstructZone* zonep, int zone, int zonet);
    virtual void renderGPUUMeshZone(UnstructZone* zonep, int zone, int zonet) = 0;
    void createGPUSSolidZone(StructZone* zonep, int zone);
    virtual void renderGPUSSolidZone(StructZone* zonep, int zone) = 0;
    virtual void createGPUUSolidZone(UnstructZone* zonep, int zone, int zonet) = 0;
    virtual void renderGPUUSolidZone(UnstructZone* zonep, int zone, int zonet) = 0;
    virtual void createGPUUSolidHOZone(UnstructZone* zonep, int zone, int zonet) = 0;
    virtual void renderGPUUSolidHOZone(UnstructZone* zonep, int zone, int zonet) = 0;
    void createGPUSIsoSolidZone(StructZone* zonep, int zone, int nofield);
    void createGPUSIsoSolidZone(StructZone* zonep, int zone, int nofield1,
                               int nofield2, int nofield3);
    void createGPUUIsoSolidZone(UnstructZone* zonep, int zone, int zonet, 
                               int nofield);
    void createGPUUIsoSolidZone(UnstructZone* zonep, int zone, int zonet, 
                               int nofield1, int nofield2, int nofield3);
    virtual void renderSIsoSolidZone(StructZone* zonep, int zone, int nofield) = 0;
    virtual void renderSIsoSolidZone(StructZone* zonep, int zone, int nofield1,
				     int nofield2, int nofield3) = 0;
    virtual void renderUIsoSolidZone(UnstructZone* zonep, int zonet, int nofield) = 0;
    virtual void renderUIsoSolidZone(UnstructZone* zonep, int zonet, int nofield1,
                             int nofield2, int nofield3) = 0;

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
    int textWidth(void* font, char* string);
    int textHeight(void* font);
    void renderBitmapString(float x, float y, float z,
                            void *font, char *string,
                            float nx=0., float ny=1., float nz=0.,
                            float r=1.);
    void renderStringWithShadow(
      float x, float y, float z,
      void *font, char *myString,
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
    void displayBigText(int posx, int posy, char* text);
    void displaySmallText(int posx, int posy, char* text);
    void printHeader();
    void printTmpMessage(const char* text);
    void displayInfoWindow(char* text, int l);
    void displayInfo();

    // Legend
    void displayIsoLegend(int dir);

    // 3D Axis
    void displayAxis();

    // 1D
    void displayPlots();
    void displayPlot(Slot1D* s);
    void plotZone(Slot1D* s, Zone1D* z, E_Float posx, E_Float posy,
                  E_Float dx, E_Float dy, int var1, int var2);
    E_Float getTick(E_Float rmin, E_Float rmax);
    void plot1DAxis(Slot1D* s, E_Float posx, E_Float posy,
                    E_Float dx, E_Float dy, E_Float blend);
    void getCharsFromVarName(char* varName, char& c1, char& c2);
    int getActivePointIndex(Zone1D* z, int var1, int var2,
                            int& e1, int& e2, double& alpha);
    int display1DActivePoint(Slot1D* s, Zone1D* z, E_Float posx, E_Float posy,
                             E_Float dx, E_Float dy, int var1, int var2);
    int link2View(Zone1D* z, int var1, int var2,
                  E_Float& r1min, E_Float& r1max, 
                  E_Float& r2min, E_Float& r2max);

    // Menu
    void displayMenu();
    void displayDimensionMenu(int* x);
    void displayVariableMenu(int* x);
    void displayAxisMenu(int* x);
    void displayZoneMenu(int* x);
    void menu();
    void menuKeyboard(unsigned char key, int x, int y);
    void menuArrows(int key, int x, int y);
    void displayMenuString(int no, char* msg, int* l, int* sizeMax,
                           int type, void* value);
    void changeMenuItemRight();
    void changeMenuItemLeft();
    void changeMenuItemSpaceBar();

    // Plugins
    void loadPlugins();
    void loadPluginsPath();
    void autoPlugins();
    int checkVariable(int zone, const char* varName);
    void findBlankedZones();
    void dumpWindow();
    char* export2Image(int exportWidth, int exportHeight);
    void superSample(int w, int h, char* im1, char* im2, int factor);
    void gaussianBlur(int w, int h, char* im1, char* im2, int r, double eps);
    void mixImages(int w, int h, char* im1, char* im2, 
                   double alpha, double beta);
    void sharpenImage(int w, int h, char* im1, char* im2, double amount,
                      int radius, int threshold);
    void exportFile();
    void finalizeExport();
    void dataMouseClickSelect(int button, int etat, int x, int y, 
                              int multiple, int accurate);
    void dataMouseRightClickSelect(int button, int etat, int x, int y);

    // Clipping
    void farClipping();
    void closeClipping();
    void veryCloseClipping();
    void veryVeryCloseClipping();
    void adaptiveClipping(double d);

    // Tools
    inline int isBlanked( Zone& zone, int zonet, int n1, int n2, int n3 ) {
      int ret1 = _pref.blanking->f(this, n1, zone.blank, zonet);
      int ret2 = _pref.blanking->f(this, n2, zone.blank, zonet);
      int ret3 = _pref.blanking->f(this, n3, zone.blank, zonet);
      return (ret1*ret2*ret3 != 0 ? 1 : 0);
    }

    inline int isBlanked( Zone& zone, int zonet, int n1, int n2, int n3, int n4 ) {
      int ret1 = _pref.blanking->f(this, n1, zone.blank, zonet);
      int ret2 = _pref.blanking->f(this, n2, zone.blank, zonet);
      int ret3 = _pref.blanking->f(this, n3, zone.blank, zonet);
      int ret4 = _pref.blanking->f(this, n3, zone.blank, zonet);
      return (ret1*ret2*ret3*ret4 != 0 ? 1 : 0);
    }

    unsigned long numberOfNonBlankedCells( UnstructZone& zone, int zonet );
 protected:
    virtual void freeGPUResources( int mode, int start, int end, int permanent );
    virtual void updateGPUResources( int mode, int size, int permanent, void* updatedPointer );
    virtual ZoneImpl* createZoneImpl( ) = 0;
};

// global functions
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
int isInFrustum(Zone* z, ViewInfo& view);

// files
void getFileExt(char* file, char* ext);
#endif

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

# ifndef _GENERATOR_GENERATOR_H_
# define _GENERATOR_GENERATOR_H_

# include "kcore.h"

# define FldArrayF K_FLD::FldArrayF
# define FldArrayI K_FLD::FldArrayI
# define FldArrayIS K_FLD::FldArrayIS

namespace K_GENERATOR
{
  PyObject* cartStruct(PyObject* self, PyObject* args);
  PyObject* cartHexa(PyObject* self, PyObject* args);
  PyObject* cartTetra(PyObject* self, PyObject* args);
  PyObject* cartPenta(PyObject* self, PyObject* args);
  PyObject* cartPyra(PyObject* self, PyObject* args);
  PyObject* cartNGon(PyObject* self, PyObject* args);
  PyObject* cylinderMesh(PyObject* self, PyObject* args);
  PyObject* cylinderMesh2(PyObject* self, PyObject* args);
  PyObject* cylinderMesh3(PyObject* self, PyObject* args);
  PyObject* delaunay(PyObject* self, PyObject* args);
  PyObject* TFIMesh(PyObject* self, PyObject* args);
  PyObject* TTMMesh(PyObject* self, PyObject* args);
  PyObject* mapMesh(PyObject* self, PyObject* args);
  PyObject* checkMesh(PyObject* self, PyObject* args);
  PyObject* checkPointInCEBBOfMesh(PyObject* self, PyObject* args);
  PyObject* getBBOfCells(PyObject* self, PyObject* args);
  PyObject* obbox(PyObject* self, PyObject* args);
  PyObject* bboxIntersection(PyObject* self, PyObject* args);
  PyObject* _bboxIntersectionZ(PyObject* self, PyObject* args);   
  PyObject* obboxIntersection(PyObject* self, PyObject* args);
  PyObject* _obboxIntersectionZ(PyObject* self, PyObject* args);
  PyObject* crossIntersection(PyObject* self, PyObject* args); 
  PyObject* _crossIntersectionZ(PyObject* self, PyObject* args);
  PyObject* barycenter(PyObject* self, PyObject* args);
  PyObject* getCEBBIntersectionOfArrays(PyObject* self, PyObject* args);
  PyObject* getVolumeMapOfMesh(PyObject* self, PyObject* args);
  PyObject* getOrthogonalityMap(PyObject* self, PyObject* args);
  PyObject* getRegularityMap(PyObject* self, PyObject* args);
  PyObject* getNormalMapOfMesh(PyObject* self, PyObject* args);
  PyObject* getCircumCircleMap(PyObject* self, PyObject* args);
  PyObject* getInCircleMap(PyObject* self, PyObject* args);
  PyObject* enforceLineMesh(PyObject* self, PyObject* args);
  PyObject* enforceXMesh(PyObject* self, PyObject* args);
  PyObject* enforceMoinsXMesh(PyObject* self, PyObject* args);
  PyObject* enforcePlusXMesh(PyObject* self, PyObject* args);
  PyObject* enforceYMesh(PyObject* self, PyObject* args);
  PyObject* enforceMoinsYMesh(PyObject* self, PyObject* args);
  PyObject* enforcePlusYMesh(PyObject* self, PyObject* args);
  PyObject* enforcePoint(PyObject* self, PyObject* args);
  PyObject* enforceCurvature(PyObject* self, PyObject* args);
  PyObject* addPointInDistribution(PyObject* self, PyObject* args);
  PyObject* hyper2DMesh(PyObject* self, PyObject* args);
  PyObject* hyper2D2Mesh(PyObject* self, PyObject* args);
  PyObject* hyper2D3Mesh(PyObject* self, PyObject* args);
  PyObject* hyper2D4Mesh(PyObject* self, PyObject* args);
  PyObject* closeMesh(PyObject* self, PyObject* args);
  PyObject* closeAllMeshes(PyObject* self, PyObject* args);
  PyObject* pointedHat(PyObject* self, PyObject* args);
  PyObject* stitchedHat(PyObject* self, PyObject* args);
  PyObject* growMesh(PyObject* self, PyObject* args);
  PyObject* stackMesh(PyObject* self, PyObject* args);
  PyObject* computeCellPlanarity(PyObject* self, PyObject* args);
  PyObject* enforceMesh(PyObject* self, PyObject* args);
  PyObject* checkDelaunay(PyObject* self, PyObject* args);
  PyObject* selectInsideElts(PyObject* self, PyObject* args);
  PyObject* densifyMesh(PyObject* self, PyObject* args);
  PyObject* T3mesher2D(PyObject* self, PyObject* args);
  PyObject* fittingPlaster(PyObject* self, PyObject* args);
  PyObject* gapfixer(PyObject* self, PyObject* args);
  PyObject* gapsmanager(PyObject* self, PyObject* args);
  PyObject* front2Hexa(PyObject* self, PyObject* args);
  PyObject* front2Struct(PyObject* self, PyObject* args);
  PyObject* snapFront(PyObject* self, PyObject* args);
  PyObject* snapSharpEdges(PyObject* self, PyObject* args);
  PyObject* fillWithStruct(PyObject* self, PyObject* args);
  PyObject* octree(PyObject* self, PyObject* args);
  PyObject* balanceOctree(PyObject* self, PyObject* args);
  PyObject* octree3(PyObject* self, PyObject* args);
  PyObject* octree2AMR(PyObject* self, PyObject* args);
  PyObject* octree2Struct(PyObject* self, PyObject* args);
  PyObject* extendCartGrids(PyObject* self, PyObject* args);
  PyObject* adaptOctree(PyObject* self, PyObject* args);
  PyObject* adaptOctree3(PyObject* self, PyObject* args);
  PyObject* conformOctree3(PyObject* self, PyObject* args);
  PyObject* modifyIndicToExpandLayer(PyObject* self, PyObject* args);
  PyObject* straightenVector(PyObject* self, PyObject* args);
  PyObject* computeEta(PyObject* self, PyObject* args);
  PyObject* getLocalStepFactor(PyObject* self, PyObject* args);
  PyObject* getEdgeRatio(PyObject* self, PyObject* args);
  PyObject* getMaxLength(PyObject* self, PyObject* args);
  PyObject* getTriQualityMap(PyObject* self, PyObject* args);
  PyObject* netgen1(PyObject* self, PyObject* args);
  PyObject* netgen2(PyObject* self, PyObject* args);
  PyObject* tetgen(PyObject* self, PyObject* args);
  PyObject* mmgs(PyObject* self, PyObject* args);
  PyObject* quad2Pyra(PyObject* self, PyObject* args);

  void computeEta(E_Int nic, E_Float* xc, E_Float* yc, E_Float* zc, 
                  E_Float* nxc, E_Float* nyc, E_Float* nzc, 
                  E_Int loop, E_Int niter,
                  E_Float* etax, E_Float* etay, E_Float* etaz);
  void normalize(E_Int npts, E_Float* vx, E_Float* vy, E_Float* vz);
  E_Int getClosestIndex(E_Int ni, E_Float x, E_Float y, E_Float z, 
                        E_Float* xc, E_Float* yc, E_Float* zc);
  void getConstrainedEta(
    E_Float etax0, E_Float etay0, E_Float etaz0, E_Int indA, 
    E_Int ni, E_Float* xc, E_Float* yc, E_Float* zc, 
    E_Float* etam, E_Float toldist);
  void applyConstraintsOnVector(
    E_Int ni, E_Float* xt, E_Float* yt, E_Float* zt,
    E_Int loop, std::vector<E_Int>& constrainedPts, 
    std::vector<FldArrayF*>& constraints,
    E_Float* vx, E_Float* vy, E_Float* vz, E_Float toldist);
  void relaxNearConstraints(
    E_Int ni, E_Float* xt, E_Float* yt, E_Float* zt,
    E_Int loop, std::vector<E_Int>& constrainedPts, 
    E_Float* alp,  E_Float* vx, E_Float* vy, E_Float* vz, E_Float dalphamax);
  void straightenVector(
    E_Int ni, E_Float* xt, E_Float* yt, E_Float* zt,
    E_Int loop, E_Int niter, std::vector<E_Int>& constrainedPts, 
    std::vector<FldArrayF*>& constraints,
    FldArrayF& vect, E_Float toldist);
/* Moteur enforce */
  E_Int enforceCommon(const char* name, char* varString, 
		      E_Int ni, E_Int nj, E_Int nk,
		      FldArrayF& f, E_Float P0, E_Float eh,
		      E_Int supp, E_Int add, FldArrayF& out,
		      E_Int& niout, E_Int& njout, E_Int& nkout,
		      E_Int autoAdd=1);

/* Curvilinear absciss calculation */
  void fa(E_Float x, FldArrayF& s, E_Int* i);

/* Cassiopee Kernel dependent functions 
   Attention : si les options de gencart sont basees sur un DesNumerics
   le DesNum n'est pas detruit a la fin, ni le DesAutomesh */
  PyObject* gencart(PyObject* self, PyObject* args);

  void readOptions(PyObject* options, E_Float& dfar, 
                   E_Float& dfxm, E_Float& dfxp, E_Float& dfym, E_Float& dfyp, 
                   E_Float& dfzm, E_Float& dfzp,
                   E_Int& automesh_type, E_Int& vmin, E_Int& snearType, 
                   E_Float& snearMul, E_Int& maxLevels, E_Int& startLevel, 
                   E_Float& stepFactor, E_Int& patchedGrids);

// void getListOfBlocks(PyObject* O, list<BlkBaseBlock*>& listOfBlks);
// void getListOfMeshes(PyObject* O, list<BlkMesh*>& listOfMeshes);

  void cellPlanarityStructured(E_Int ni, E_Int nj, E_Int nk,
                               FldArrayF& f, E_Int posx, E_Int posy, E_Int posz,
                               FldArrayF& dist, 
                               E_Int& ni1, E_Int& nj1, E_Int& nk1);
  void cellPlanarityUnstructured(FldArrayF& f, FldArrayI& cn,
                                 E_Int posx, E_Int posy, E_Int posz,
                                 FldArrayF& dist);
/* close a structured mesh */
  void closeStructuredMesh(E_Float* xt, E_Float* yt, E_Float* zt,
                           E_Int ni, E_Int nj, E_Int nk, E_Float eps);    
/* close an unstructured mesh */
  void closeUnstructuredMesh(E_Int posx, E_Int posy,E_Int posz,E_Float eps,
                             char* eltType, FldArrayF& f, FldArrayI& cn);
/* close a BAR mesh */
  void closeBARMesh(E_Int posx, E_Int posy, E_Int posz, 
                    FldArrayF& f, FldArrayI& cn);
  E_Int closeBARMeshElt(E_Int posx, E_Int posy, E_Int posz,
                        FldArrayF& f, FldArrayI& cn, E_Int i);

  PyObject* TFI2D(PyObject* arrays);

  PyObject* TFI3D(PyObject* arrays);

  PyObject* TFITRI(PyObject* arrays);

  PyObject* TFITETRA(PyObject* arrays);

  PyObject* TFIPENTA(PyObject* arrays);

/* TFI 2D structure. Retourne 0 si echec.*/
  short TFIstruct2D(E_Int ni, E_Int nj, E_Int nfld, 
                    E_Int imin, E_Int imax, E_Int jmin, E_Int jmax, 
                    std::vector<FldArrayF*>& fields,
                    FldArrayF& coords);
  short TFIstruct2D2(E_Int ni, E_Int nj, E_Int nfld,
                     E_Int posx, E_Int posy, E_Int posz,
                     E_Int imin, E_Int imax, E_Int jmin, E_Int jmax,
                     std::vector<FldArrayF*>& fields,
                     FldArrayF& coords);
/* TFI 3D structure. Retourne 0 si echec.*/
  short TFIstruct3D(E_Int ni, E_Int nj, E_Int nk, E_Int nfld, 
                    E_Int imin, E_Int imax, E_Int jmin, E_Int jmax, E_Int kmin, E_Int kmax, 
                    std::vector<FldArrayF*>& fields,
                    FldArrayF& coords);
  short TFIstruct3D2(
                    E_Int ni, E_Int nj, E_Int nk, E_Int nfld,
                    E_Int posx, E_Int posy, E_Int posz,
                    E_Int imin, E_Int imax, E_Int jmin, E_Int jmax, E_Int kmin, E_Int kmax,
                    std::vector<FldArrayF*>& fields, FldArrayF& coords);
  E_Int reorderTFI2D(E_Int posx, E_Int posy, E_Int posz,
                     std::vector<E_Int>& nit, std::vector<FldArrayF*>& fields, 
                     FldArrayIS& newOrder);

  E_Int reorderTFI3D(E_Int posx, E_Int posy, E_Int posz,
                     std::vector<E_Int>& nit, std::vector<E_Int>& njt,
                     std::vector<FldArrayF*>& fields, 
                     FldArrayIS& newOrder, E_Float eps=1.e-10);

  E_Int reorderTFITRI(E_Int posx, E_Int posy, E_Int posz,
                      std::vector<E_Int>& nit, 
                      std::vector<FldArrayF*>& fields, 
                      FldArrayIS& newOrder);
/* Verifie la continuite des frontieres imin, imax, ...
   aux coins. Vaut 0 si pas continu, 1 sinon  */
  short checkContinuousBnds3D(E_Int ni, E_Int nj, E_Int nk,
                              E_Float* ximin, E_Float* yimin, E_Float* zimin,
                              E_Float* ximax, E_Float* yimax, E_Float* zimax,
                              E_Float* xjmin, E_Float* yjmin, E_Float* zjmin,
                              E_Float* xjmax, E_Float* yjmax, E_Float* zjmax,
                              E_Float* xkmin, E_Float* ykmin, E_Float* zkmin,
                              E_Float* xkmax, E_Float* ykmax, E_Float* zkmax);

  short checkContinuousBndsTRI(E_Int ni,
                               E_Float* xt1, E_Float* yt1, E_Float* zt1,
                               E_Float* xt2, E_Float* yt2, E_Float* zt2,
                               E_Float* xt3, E_Float* yt3, E_Float* zt3);

  short checkContinuousBndTETRA(E_Int n,
                                E_Float* xt1, E_Float* yt1, E_Float* zt1,
                                E_Float* xt2, E_Float* yt2, E_Float* zt2,
                                E_Float* xt3, E_Float* yt3, E_Float* zt3,
                                E_Float* xt4, E_Float* yt4, E_Float* zt4);

  short checkContinuousBndsPENTA(
    E_Int n, E_Int p,
    E_Float* xtmin, E_Float* ytmin, E_Float* ztmin,
    E_Float* xtmax, E_Float* ytmax, E_Float* ztmax,
    E_Float* ximin, E_Float* yimin, E_Float* zimin,
    E_Float* xjmin, E_Float* yjmin, E_Float* zjmin,
    E_Float* xdiag, E_Float* ydiag, E_Float* zdiag);

/* Calcul des coefficients de l 'eq. du plan ax + by + cz + d = 0*/
  E_Int compPlaneCoefficients(E_Int posx, E_Int posy, E_Int posz, 
                              FldArrayF& f, 
                              E_Float& a, E_Float& b, E_Float& c, E_Float& d);

/* Calcul des normales a une courbe */
  void normals(E_Int n, E_Float* x, E_Float* y, E_Float* z,
               E_Float* nx, E_Float *ny, E_Float* nz);

/* Calcul de la distance a une courbe sauf le point ind */
  E_Float dist2Curve(E_Int n, E_Float* x, E_Float* y, E_Float* z, E_Int* found,
                     E_Float px, E_Float py, E_Float pz, E_Int ind,
                     E_Int& nearest);

  E_Int getNeighbourElts(E_Int npts, E_Float* xt, E_Float* yt, E_Float* zt, 
                         FldArrayI& cEV, 
                         std::vector< std::vector<E_Int> >& nghbrs,
                         E_Int type=0, E_Float mindh=1.e-9);
  /* Retourne les voisins par elt ayant une facette en commun */
  E_Int getNeighbourQuads(E_Int npts, E_Float* xt, E_Float* yt, E_Float* zt, 
                          FldArrayI& cEV, 
                          std::vector< std::vector<E_Int> >& nghbrs );
  /* Retourne les voisins par elt ayant une facette en commun */
  E_Int getNeighbourHexas(E_Int npts, E_Float* xt, E_Float* yt, E_Float* zt, 
                          FldArrayI& cEV, 
                          std::vector< std::vector<E_Int> >& nghbrs );
  /* Retourne les voisins par elt ayant une facette ou un point en commun */
  E_Int getNeighbourHexas2(E_Int npts, E_Float* xt, E_Float* yt, E_Float* zt, 
                           FldArrayI& cEV, 
                           std::vector< std::vector<E_Int> >& nghbrs,
                           E_Float mindh=1.e-9);
  E_Int getNeighbourQuads2(E_Int npts, E_Float* xt, E_Float* yt, E_Float* zt, 
                           FldArrayI& cEV, 
                           std::vector< std::vector<E_Int> >& nghbrs,
                           E_Float mindh=1.e-9);
  void splitElement(E_Int et, E_Int npts, E_Int nelts, E_Float indic, 
                    E_Int* indir,
                    FldArrayI& cn, E_Float* xt, E_Float* yt, E_Float* zt, 
                    FldArrayI& cno, FldArrayF& fo, FldArrayF& indicout, 
                    E_Int& no, E_Int& eto);
  void splitElement(E_Int et, E_Int npts, E_Int nelts, E_Float indic, 
                    FldArrayI& cn, E_Float* xt, E_Float* yt, E_Float* zt, 
                    FldArrayI& cno, FldArrayF& fo, FldArrayF& indicout, 
                    E_Int& no, E_Int& eto);
  void splitElement27(E_Int et, E_Int npts, E_Float indic,
                      FldArrayI& cn, E_Float* xt, E_Float* yt, E_Float* zt, 
                      FldArrayI& cno, FldArrayF& fo, FldArrayF& indicout, 
                      E_Int& no, E_Int& eto);  
  /* fonctions de adaptOctree */
  E_Int mergeOctreeElement(E_Int et, E_Int npts, E_Float indic, FldArrayI& cn,
                           E_Float xs, E_Float ys, E_Float zs,
                           E_Float xe, E_Float ye, E_Float ze, 
                           E_Float* xt, E_Float* yt, E_Float* zt, 
                           E_Float* dht, E_Float* indict,
                           FldArrayI& cno, FldArrayF& fo, FldArrayF& indicout,
                           E_Int& no, E_Int& eto, FldArrayIS& dejaVu);
  E_Int mergeOctreeElement27(E_Int et, E_Int npts, E_Float indic, FldArrayI& cn,
                             E_Float xs, E_Float ys, E_Float zs,
                             E_Float xe, E_Float ye, E_Float ze, 
                             E_Float* xt, E_Float* yt, E_Float* zt, 
                             E_Float* dht, E_Float* indict,
                             FldArrayI& cno, FldArrayF& fo, FldArrayF& indicout,
                             E_Int& no, E_Int& eto, FldArrayIS& dejaVu);

  void modifyIndicator(E_Int nelts, FldArrayF& dht, std::vector< std::vector<E_Int> >& cEEN,
                       E_Int posi, E_Float* indict);

  /* Insert patterns to conformize an octree3
     IN : et : element to conformize
     IN : tag : tag (0 or 1) defined for each vertex of the element to determine the pattern
     IN : f : field containing coordinates of the initial octree3
     IN : cn : connectivity of the initial octree
     IN/OUT : fo : fields of the conformized octree after insertion of conformizing elements
     IN/OUT : cno : connectivity of the conformized octree after insertion of conformizing elements
     IN/OUT : eto : number of elements already built in the conformized octree
     IN/OUT : no : nb of pts already built in the conformized octree */
  void insert2DPattern1(E_Int et, E_Int let, E_Int* tag,
                        FldArrayF& f, FldArrayI& cn,
                        FldArrayF& fo, FldArrayI& cno, FldArrayI& levelso, E_Int& eto, E_Int& no);
  void insert2DPattern2(E_Int et, E_Int let, E_Int* tag,
                        FldArrayF& f, FldArrayI& cn,
                        FldArrayF& fo, FldArrayI& cno, FldArrayI& levelso, E_Int& eto, E_Int& no);
  void insert2DPattern3(E_Int et, E_Int let, E_Int* tag,
                        FldArrayF& f, FldArrayI& cn,
                        FldArrayF& fo, FldArrayI& cno, FldArrayI& levelso, E_Int& eto, E_Int& no);
  void insert2DPattern4(E_Int et, E_Int let, E_Int* tag,
                        FldArrayF& f, FldArrayI& cn,
                        FldArrayF& fo, FldArrayI& cno, FldArrayI& levelso, E_Int& eto, E_Int& no);
  E_Int find2Template(E_Int* tags);
  E_Int find3Template(E_Int* tags);
  E_Int find4Template(E_Int* tags);
  void insert3DPattern1(E_Int et, E_Int let, E_Int* tag,
                        FldArrayF& f, FldArrayI& cn,
                        FldArrayF& fo, FldArrayI& cno, FldArrayI& levelso, E_Int& eto, E_Int& no);
  void insert3DPattern3(E_Int et, E_Int let, E_Int* tag,
                        FldArrayF& f, FldArrayI& cn,
                        FldArrayF& fo, FldArrayI& cno, FldArrayI& levelso, E_Int& eto, E_Int& no); 
  void insert3DPattern4(E_Int et, E_Int let, E_Int* tag,
                        FldArrayF& f, FldArrayI& cn,
                        FldArrayF& fo, FldArrayI& cno, FldArrayI& levelso, E_Int& eto, E_Int& no);
  
  E_Int conformizeElements(E_Float mindh, E_Float maxdh, E_Int lmax, E_Int npts,
                           E_Int posx, E_Int posy, E_Int posz,
                           FldArrayF& f, FldArrayI& cn, FldArrayI& levels,
                           FldArrayF& fo, FldArrayI& cno, FldArrayI& levelso, 
                           E_Int& no, E_Int& eto);
  /* Balancing of an octree or an octree3 */

  void checkBalancing2(FldArrayI& cn, FldArrayF& coords);
  void checkBalancing3(FldArrayI& cn, FldArrayF& coords);
  void balanceOctree2(FldArrayF& f, FldArrayI& cn, E_Int corners);

  void getValidNgbrsForMerge(E_Int et, E_Float* indict, E_Float* dht, 
                             E_Float xs, E_Float ys, E_Float zs,
                             E_Float xe, E_Float ye, E_Float ze,
                             E_Float* xt, E_Float* yt, E_Float* zt,
                             FldArrayIS& dejaVu, FldArrayI& cn,
                             std::vector<E_Int>& candidats);
  /* close */
  void closeOneWindow(
    FldArrayF& f1, E_Int ni1, E_Int nj1, E_Int nk1,
    E_Int posx1, E_Int posy1, E_Int posz1,
    E_Int im1, E_Int im2, E_Int jm1, E_Int jm2, E_Int km1, E_Int km2,
    FldArrayF& f2, E_Int ni2, E_Int nj2, E_Int nk2,
    E_Int posx2, E_Int posy2, E_Int posz2,
    E_Int is1, E_Int is2, E_Int js1, E_Int js2, E_Int ks1, E_Int ks2,
    E_Float eps, E_Boolean check);

  void closeAllStructuredMeshes(
    std::vector<FldArrayF*>& structF,
    std::vector<E_Int>& nit, std::vector<E_Int>& njt, std::vector<E_Int>& nkt,
    std::vector<E_Int>& posx, std::vector<E_Int>& posy, std::vector<E_Int>& posz,
    E_Float eps);

/* determination des blocs intersectants noz1 ds bboxes
   IN : noz1 : numero de la zone dans bboxes
   IN : minB : xmin(z1), ymin(z1), zmin(z1) 
   IN : maxB: xmin(z1), ymin(z1), zmin(z1)
   IN : bboxes : liste des bboxes de ttes les zones
   IN : tol tolerance d intersection
   OUT : listOfZones : indices des zones dans bboxes intersectant la zone noz1 */
  void getBlocksIntersecting(E_Int noz1, E_Float* minB, E_Float* maxB,
                             FldArrayF& bboxes, E_Float tol, 
                             std::vector<E_Int>& listOfZones); 

  E_Int detectMatchCart(E_Int nowin1, E_Int ni1, E_Int nj1, E_Int nk1,
                        E_Float* xt1, E_Float* yt1, E_Float* zt1,
                        E_Int ni2, E_Int nj2, E_Int nk2,
                        E_Float* xt2, E_Float* yt2, E_Float* zt2,
                        E_Int& im1, E_Int& ip1, E_Int& jm1, E_Int& jp1, E_Int& km1, E_Int& kp1,
                        E_Int& im2, E_Int& ip2, E_Int& jm2, E_Int& jp2, E_Int& km2, E_Int& kp2);

  /* Calcul de la normale pour l element et par relaxation avec prise en compte des elts voisins */
  void smoothNormalsForTriangle(E_Int et, std::vector<E_Int>& eltsVoisins,
                                E_Float* sx, E_Float* sy, E_Float* sz, 
                                E_Float& sxet, E_Float& syet, E_Float& szet);
  /* Calcul d'un tableau d indirection servant a construire un front de maillage structure selon un critre sur le cellN */
  void computeStructFrontIndi(E_Float* cellN, E_Int ni, E_Int nj, E_Int nk, E_Int var1, E_Int var2, E_Int corners, E_Int* indi, E_Int& np);
  /* Calcul d'un tableau d indirection servant a construire un front de maillage non-structure selon un critre sur le cellN */
  void computeUnstrFrontIndi(E_Float* cellN, E_Int npts, E_Int var1, E_Int var2, E_Int corners, FldArrayI& cn, E_Int* indi, E_Int& np);
  /* Snap un maillage par rapport a une surface */
  void snapMesh(FldArrayF& surface, std::vector<E_Int> sizet,E_Int posx,E_Int posy,E_Int posz,
                std::vector<E_Float*>& coordx, std::vector<E_Float*>& coordy, 
                std::vector<E_Float*>& coordz, std::vector<E_Float*>& indic,
                std::vector<FldArrayI*>& connect);
}
# undef FldArrayF
# undef FldArrayI
# undef FldArrayIS

#endif

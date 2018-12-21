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

#ifndef _TRANSFORM_TRANSFORM_H_
#define _TRANSFORM_TRANSFORM_H_

# include <locale>
# include <cctype>
# include "kcore.h"
# include "Fld/FldArray.h"

# define FldArrayF K_FLD::FldArrayF
# define FldArrayI K_FLD::FldArrayI

namespace K_TRANSFORM
{

/* joint deux arrays 1 et 2 structures
   IN: f1: champs de l array 1 : contient coord
   IN: im1, jm1, km1: dimensions de array1
   IN: posx1, posy1, posz1: position des coordonnees ds f1
   IN: pos1: position des champs communs a f2 ds f1
   IN: f2: champs de l array 2 : contient coord
   IN: im2, jm2, km2: dimensions de array2
   IN: posx2, posy2, posz2: position des coordonnees ds f2
   IN: pos2: position des champs communs a f1 ds f2
   OUT: field: array resultant des champs obtenus par join des 2 arrays
   OUT: im, jm, km: dimensions de l array resultant
   retourne 0 si pas de join possible, 1 si ok 
   ATTENTION: on fait une copie car joinStructured modifie f1 et f2
*/
  E_Int joinStructured(FldArrayF f1, E_Int im1, E_Int jm1, E_Int km1,
                       E_Int posx1, E_Int posy1, E_Int posz1,
                       FldArrayF f2, E_Int im2, E_Int jm2, E_Int km2,
                       E_Int posx2, E_Int posy2, E_Int posz2,
                       std::vector<E_Int>& pos1, std::vector<E_Int>& pos2,
                       FldArrayF& field,
                       E_Int& im, E_Int& jm, E_Int& km, E_Float tol);

  E_Int joinBothStructured(FldArrayF f1, E_Int im1, E_Int jm1, E_Int km1,
                           E_Int posx1, E_Int posy1, E_Int posz1,
                           FldArrayF f2, E_Int im2, E_Int jm2, E_Int km2,
                           E_Int posx2, E_Int posy2, E_Int posz2,
                           FldArrayF fc1, E_Int imc1, E_Int jmc1, E_Int kmc1,
                           FldArrayF fc2, E_Int imc2, E_Int jmc2, E_Int kmc2,
                           std::vector<E_Int>& pos1, std::vector<E_Int>& pos2,
                           std::vector<E_Int>& posc1, std::vector<E_Int>& posc2,
                           FldArrayF& field, E_Int& im, E_Int& jm, E_Int& km, 
                           FldArrayF& fieldc, E_Int& imc, E_Int& jmc, E_Int& kmc, E_Float tol);

  // join structure 3d: memes arguments que joinStructured 
  E_Int joinstructured3d(FldArrayF& f1, E_Int im1, E_Int jm1, E_Int km1,
                         E_Int posx1, E_Int posy1, E_Int posz1,
                         FldArrayF& f2, E_Int im2, E_Int jm2, E_Int km2,
                         E_Int posx2, E_Int posy2, E_Int posz2,
                         std::vector<E_Int>& pos1, std::vector<E_Int>& pos2,
                         FldArrayF& field, E_Int& im, E_Int& jm, E_Int& km, E_Float tol);

  E_Int joinbothstructured3d(FldArrayF& f1, E_Int im1, E_Int jm1, E_Int km1,
                             E_Int posx1, E_Int posy1, E_Int posz1,
                             FldArrayF& f2, E_Int im2, E_Int jm2, E_Int km2,
                             E_Int posx2, E_Int posy2, E_Int posz2,
                             FldArrayF& fc1, E_Int imc1, E_Int jmc1, E_Int kmc1,
                             FldArrayF& fc2, E_Int imc2, E_Int jmc2, E_Int kmc2,
                             std::vector<E_Int>& pos1, std::vector<E_Int>& pos2,
                             std::vector<E_Int>& posc1, std::vector<E_Int>& posc2,
                             FldArrayF& field, E_Int& im, E_Int& jm, E_Int& km, 
                             FldArrayF& fieldc, E_Int& imc, E_Int& jmc, E_Int& kmc, E_Float tol);

  E_Int nextCornerMatchingIndices(E_Int i1, E_Int j1, E_Int k1, 
                                  E_Int im1, E_Int jm1, E_Int km1,  
                                  E_Int im2, E_Int jm2, E_Int km2, 
                                  E_Int posx1, E_Int posy1, E_Int posz1, 
                                  E_Int posx2, E_Int posy2, E_Int posz2,
                                  FldArrayF& f1, FldArrayF& f2, 
                                  E_Int& i2, E_Int& j2, E_Int& k2,
                                  E_Float eps=1.e-12);

  // join structure 2d: memes arguments que joinStructured
  E_Int joinstructured2d(FldArrayF& f1, E_Int im1, E_Int jm1, E_Int km1,
                         E_Int posx1, E_Int posy1, E_Int posz1,
                         FldArrayF& f2, E_Int im2, E_Int jm2, E_Int km2,
                         E_Int posx2, E_Int posy2, E_Int posz2,
                         std::vector<E_Int>& pos1, std::vector<E_Int>& pos2,
                         FldArrayF& field,
                         E_Int& im, E_Int& jm, E_Int& km, E_Float tol);

  E_Int joinbothstructured2d(FldArrayF& f1, E_Int im1, E_Int jm1, E_Int km1,
                             E_Int posx1, E_Int posy1, E_Int posz1,
                             FldArrayF& f2, E_Int im2, E_Int jm2, E_Int km2,
                             E_Int posx2, E_Int posy2, E_Int posz2,
                             FldArrayF& fc1, E_Int imc1, E_Int jmc1, E_Int kmc1,
                             FldArrayF& fc2, E_Int imc2, E_Int jmc2, E_Int kmc2,
                             std::vector<E_Int>& pos1, std::vector<E_Int>& pos2,
                             std::vector<E_Int>& posc1, std::vector<E_Int>& posc2,
                             FldArrayF& field, E_Int& im, E_Int& jm, E_Int& km, 
                             FldArrayF& fieldc, E_Int& imc, E_Int& jmc, E_Int& kmc, E_Float tol);
  // join structure 1d: memes arguments que joinStructured
  E_Int joinstructured1d(FldArrayF& f1, E_Int im1, E_Int jm1, E_Int km1,
                         E_Int posx1, E_Int posy1, E_Int posz1,
                         FldArrayF& f2, E_Int im2, E_Int jm2, E_Int km2,
                         E_Int posx2, E_Int posy2, E_Int posz2,
                         std::vector<E_Int>& pos1, std::vector<E_Int>& pos2,
                         FldArrayF& field,
                         E_Int& im, E_Int& jm, E_Int& km, E_Float tol);
  E_Int joinbothstructured1d(FldArrayF& f1, E_Int im1, E_Int jm1, E_Int km1,
                             E_Int posx1, E_Int posy1, E_Int posz1,
                             FldArrayF& f2, E_Int im2, E_Int jm2, E_Int km2,
                             E_Int posx2, E_Int posy2, E_Int posz2,
                             FldArrayF& fc1, E_Int imc1, E_Int jmc1, E_Int kmc1,
                             FldArrayF& fc2, E_Int imc2, E_Int jmc2, E_Int kmc2,
                             std::vector<E_Int>& pos1, std::vector<E_Int>& pos2,
                             std::vector<E_Int>& posc1, std::vector<E_Int>& posc2,
                             FldArrayF& field, E_Int& im, E_Int& jm, E_Int& km, 
                             FldArrayF& fieldc, E_Int& imc, E_Int& jmc, E_Int& kmc, E_Float tol);

  /* Join 2 arrays non structures */
  E_Int joinUnstructured(FldArrayF& f1, FldArrayI& cn1,
                         E_Int posx1, E_Int posy1, E_Int posz1,
                         FldArrayF& f2, FldArrayI& cn2,
                         E_Int posx2, E_Int posy2, E_Int posz2,
                         E_Int nfld, char* eltType,
                         FldArrayF& field, FldArrayI& cn, E_Float tol);
  E_Int joinBothUnstructured(FldArrayF& f1, FldArrayI& cn1,
                             FldArrayF& fc1, FldArrayI& cnc1,
                             E_Int posx1, E_Int posy1, E_Int posz1,
                             FldArrayF& f2, FldArrayI& cn2,
                             FldArrayF& fc2, FldArrayI& cnc2,
                             E_Int posx2, E_Int posy2, E_Int posz2,
                             std::vector<E_Int>& pos1, std::vector<E_Int>& pos2,
                             std::vector<E_Int>& posc1, std::vector<E_Int>& posc2,
                             char* eltType,
                             FldArrayF& field, FldArrayI& cn, 
                             FldArrayF& fieldc, FldArrayI& cnc,
                             E_Float tol);

  /* Join 2 arrays NGON */
  E_Int joinNGON(FldArrayF& f1, FldArrayI& cn1,
                 E_Int posx1, E_Int posy1, E_Int posz1,
                 FldArrayF& f2, FldArrayI& cn2,
                 E_Int posx2, E_Int posy2, E_Int posz2,
                 E_Int nfld, FldArrayF& field, FldArrayI& cn, E_Float tol);
  E_Int joinBothNGON(FldArrayF& f1, FldArrayI&  cn1,
                     FldArrayF& fc1, FldArrayI&  cnc1,
                     E_Int posx1, E_Int posy1, E_Int posz1,
                     FldArrayF& f2, FldArrayI& cn2,
                     FldArrayF& fc2, FldArrayI& cnc2,
                     E_Int posx2, E_Int posy2, E_Int posz2,
                     std::vector<E_Int>& pos1,  std::vector<E_Int>& pos2,
                     std::vector<E_Int>& posc1,  std::vector<E_Int>& posc2,
                     FldArrayF& field, FldArrayI& cn,
                     FldArrayF& fieldc, FldArrayI& cnc,
                     E_Float tol);

  /* Split d'une courbe structuree  definie par f */
  void splitSplineStruct(E_Float dmax, E_Float cvmax,
                         E_Int posx, E_Int posy, E_Int posz, FldArrayF* f,
                         std::vector<FldArrayF*>& fsplit);


  void collapseMinVertexInTriangle(FldArrayI& cn, 
                                   E_Float* xt, E_Float* yt, E_Float* zt);

  /* Fonctions necessaires a projectSmoothDir
     si oriented != 0 : le projete doit etre dans la direction de (nx,ny,nz) */
  void projectSmoothDirWithoutPrecond(E_Float nx, E_Float ny, E_Float nz,
                                      E_Int im1, E_Int jm1, E_Int km1,
                                      E_Int nelts2, FldArrayI& cn2,
                                      E_Float* fx2, E_Float* fy2, E_Float* fz2,
                                      E_Float* fx, E_Float* fy, E_Float* fz, E_Int oriented=0);
  void projectSmoothDirWithPrecond(E_Float nx, E_Float ny, E_Float nz,
                                   E_Int im1, E_Int jm1, E_Int km1,
                                   E_Int nelts2, FldArrayI& cn2,
                                   E_Float* fx2, E_Float* fy2, E_Float* fz2,
                                   E_Float* fx, E_Float* fy, E_Float* fz, 
                                   E_Int oriented=0);
  void smoothUnprojectedPts(E_Int im1, E_Int jm1, E_Int km1, E_Int npts,
                            E_Float* fx, E_Float* fy, E_Float* fz, short* tagp);

  /* Fonctions necessaires a reorderAll */
  E_Int pivotingBlks(E_Int bg, E_Int bd, FldArrayI& sortedBlks, 
                     FldArrayF& dist);

  void quickSortBlks(E_Int bg, E_Int bd, FldArrayI& sortedBlks, 
                     FldArrayF& dist);

  void sortBlocks(FldArrayF& distMat, 
                  std::vector<FldArrayI*>& listOfSortedBlks);

  // Algorithme recursif de parcours de graphe: Depth First Search   
  //IN: nb_father:  no du bloc pere
  //IN: nb_cur: no du bloc courant
  //IN: rel: relation entre les blocs : distance < eps
  //IN/OUT: dejaVu 
  //IN/OUT: tagOpp : sens de la normale par rapport a la reference 
  void graphPathSearch( E_Int nb_father, E_Int nb_cur, 
                        std::vector<FldArrayF*>& listOfFields, 
                        std::vector<E_Int>& nit, std::vector<E_Int>& njt,
                        std::vector<E_Int>& posxt, std::vector<E_Int>& posyt,
                        std::vector<E_Int>& poszt,
                        std::vector<FldArrayI*>& listOfSortedBlks,
                        FldArrayI& rel, FldArrayF& distMat,
                        FldArrayI& dejaVu,
                        FldArrayI& tagOpp);

/* Reorder a 2D grid */
  void reorder(const E_Int im, const E_Int jm, 
               FldArrayF& field);

/* Return 1 if bbox of blocks intersect, 0 elsewhere */
  E_Int testBBIntersection(E_Int noblk1, E_Int noblk2,
                           FldArrayF& bbox);

/* Return 1 if one negative volume cell exists in the mesh */
  E_Int checkNegativeVolumeCells(E_Int dim, E_Int im, E_Int jm, E_Int km, 
                                 FldArrayF& coords);
  void increaseGradeForRotation(
    E_Int indA1, E_Int indB1, E_Int indC1, E_Int indD1,
    E_Float* xt1, E_Float* yt1, E_Float* zt1,
    E_Int indA2, E_Int indB2, E_Int indC2, E_Int indD2,
    E_Float* xt2, E_Float* yt2, E_Float* zt2,
    E_Int& grade);
  /* Computes the max deviation angle (from 180 deg) between zones 1 and 
     2 sharing an edge of direction dir1 in z1 and dir2 in z2 */
  E_Float compAngleBetweenZones(E_Int dir1, E_Int ni1, E_Int nj1,
                                E_Float* xt1, E_Float* yt1, E_Float* zt1,
                                E_Int dir2, E_Int ni2, E_Int nj2,
                                E_Float* xt2, E_Float* yt2, E_Float* zt2);

  void breakNGonElements(FldArrayF& field, FldArrayI& cFNEF, 
                         std::vector<FldArrayI*>& cEV, 
                         std::vector<FldArrayF*>& fields, 
                         std::vector<E_Int>& eltType);
  
  /* Cree le dual d un maillage NGON 2D ou 3D.
   IN: f: champs contenant les coordonnees localises aux noeuds du 
   maillage primal
   IN: cn: connectivite NGON  du primal
   OUT: fd: champs aux noeuds du maillage dual
   OUT: cNGD: connectivite NGON du dual 
   fd et cNGD sont alloues dans ces fonctions */
  void dualNGON2D(FldArrayF& f, FldArrayI& cn, E_Int extraPoints,
                  FldArrayF& fd, FldArrayI& cNGD);
  void dualNGON3D(FldArrayF& f, FldArrayI& cn, 
                  FldArrayF& fd, FldArrayI& cNGD);
  //E_Int createDegeneratedPrimalMesh3D(E_Int nptsp, FldArrayI& cNG, 
  //                                    FldArrayI& cNGD);
  E_Int createDegeneratedPrimalMesh3D(FldArrayF& fNG, FldArrayI& cNG, 
                                      FldArrayF& fNGD, FldArrayI& cNGD);
  void  umbrella(FldArrayF& coord, FldArrayF& coordo,
                 FldArrayF& move, FldArrayF& proj,
                 FldArrayF* f3, FldArrayI* cn3,
                 E_Int posx1, E_Int posy1, E_Int posz1,
                 E_Int posx3, E_Int posy3, E_Int posz3,
                 E_Int projConstraintOn, E_Float delta, E_Int type,
                 E_Float eps, E_Int niter,
                 std::vector< std::vector<E_Int> > &cVN,
                 E_Float xR, E_Float yR, E_Float zR, E_Float radius);
  
  E_Int extractVectorComponents(char* varString, PyObject* listOfFieldVectors, 
                                std::vector<E_Int>& posvx, std::vector<E_Int>& posvy, std::vector<E_Int>& posvz);

  // Flip Edges dans un maillage TRI
  void flipEdges(FldArrayI& ct, E_Int np, E_Float* x, E_Float* y, E_Float* z, E_Float* indic=NULL, E_Int mode=1);

  // Contract Edges dans un maillage TRI
  void contractEdges(FldArrayI& ct, E_Int np, E_Float* x, E_Float* y, E_Float* z, E_Int mode=1);

  // Check TRI
  void checkTriMesh(FldArrayI& ct, E_Int np, E_Float* x, E_Float* y, E_Float* z,
                    E_Int& ne, E_Int& ni);

 /*------------------------------------------------------------------*/
  PyObject* _cart2CylZ(PyObject* self, PyObject* args);
  PyObject* _cart2CylA(PyObject* self, PyObject* args);
  PyObject* _cyl2CartZ(PyObject* self, PyObject* args);
  PyObject* _cyl2CartA(PyObject* self, PyObject* args);
  PyObject* translate(PyObject* self, PyObject* args);
  PyObject* rotateA1(PyObject* self, PyObject* args);
  PyObject* rotateA2(PyObject* self, PyObject* args);
  PyObject* rotateA3(PyObject* self, PyObject* args);
  PyObject* homothety(PyObject* self, PyObject* args);
  PyObject* contract(PyObject* self, PyObject* args);
  PyObject* symetrize(PyObject* self, PyObject* args);
  PyObject* perturbate(PyObject* self, PyObject* args);
  PyObject* smooth(PyObject* self, PyObject* args);
  PyObject* deform(PyObject* self, PyObject* args);
  PyObject* deform2(PyObject* self, PyObject* args);
  PyObject* deformPoint(PyObject* self, PyObject* args);
  PyObject* deformMeshStruct(PyObject* self, PyObject* args);
  PyObject* computeDeformationVector(PyObject* self, PyObject* args);

  PyObject* join(PyObject* self, PyObject* args);
  PyObject* joinBoth(PyObject* self, PyObject* args);
  PyObject* joinAll(PyObject* self, PyObject* args);
  PyObject* joinAllBoth(PyObject* self, PyObject* args);
  PyObject* patch(PyObject* self, PyObject* args);
  PyObject* patch2(PyObject* self, PyObject* args);

  PyObject* subzoneStruct(PyObject* self, PyObject* args);
  PyObject* subzoneStructInt(PyObject* self, PyObject* args);
  PyObject* subzoneStructIntBoth(PyObject* self, PyObject* args);
  PyObject* subzoneUnstruct(PyObject* self, PyObject* args);
  PyObject* subzoneUnstructBoth(PyObject* self, PyObject* args);
  PyObject* subzoneElements(PyObject* self, PyObject* args);
  PyObject* subzoneElementsBoth(PyObject* self, PyObject* args);
  PyObject* subzoneFaces(PyObject* self, PyObject* args);
  PyObject* subzoneFacesBoth(PyObject* self, PyObject* args);

  PyObject* oneovern(PyObject* self, PyObject* args);
  PyObject* reorder(PyObject* self, PyObject* args);
  PyObject* reorderAll(PyObject* self, PyObject* args);
  PyObject* reorderAllUnstr(PyObject* self, PyObject* args);
  PyObject* addkplane(PyObject* self, PyObject* args);
  PyObject* addkplaneCenters(PyObject* self, PyObject* args);

  PyObject* projectAllDirs(PyObject* self, PyObject* args);
  PyObject* projectDir(PyObject* self, PyObject* args);
  PyObject* projectOrtho(PyObject* self, PyObject* args);
  PyObject* projectOrthoSmooth(PyObject* self, PyObject* args);
  PyObject* projectRay(PyObject* self, PyObject* args);
  PyObject* projectSmoothDir(PyObject* self, PyObject* args);

  PyObject* splitCurvatureAngle(PyObject* self, PyObject* args);
  PyObject* splitCurvatureRadius(PyObject* self, PyObject* args);
  PyObject* splitConnexity(PyObject* self, PyObject* args);
  PyObject* splitConnexityBasics(FldArrayF* f, FldArrayI* cn, 
                                 char* eltType, char* varString,
                                 E_Int posx, E_Int posy, E_Int posz);
  PyObject* splitConnexityNGon(FldArrayF* f, FldArrayI* cn, char* varString,
                               E_Int posx, E_Int posy, E_Int posz);
  PyObject* splitConnexityNODE(FldArrayF* f, FldArrayI* cn,
                               char* eltType, char* varString,
                               E_Int posx, E_Int posy, E_Int posz);

  PyObject* splitSharpEdges(PyObject* self, PyObject* args);
  PyObject* splitSharpEdgesBasics(FldArrayF* f, FldArrayI* cn, 
                                  char* eltType, char* varString,
                                  E_Int posx, E_Int posy, E_Int posz,
                                  E_Int type, E_Float alphaRef, 
                                  E_Float* dirVect);
  PyObject* splitSharpEdgesNGon(FldArrayF* f, FldArrayI* cn, char* varString,
                                E_Int posx, E_Int posy, E_Int posz,
                                E_Float alphaRef, E_Float* dirVect);
  PyObject* splitSharpEdgesList(PyObject* self, PyObject* args);
  PyObject* splitBAR(PyObject* self, PyObject* args);
  PyObject* splitTBranches(PyObject* self, PyObject* args);
  PyObject* splitTRI(PyObject* self, PyObject* args);
  PyObject* splitManifold(PyObject* self, PyObject* args);

  PyObject* collapse(PyObject* self, PyObject* args);
  PyObject* mergeCartGrids(PyObject* self, PyObject* args);
  PyObject* mergeStructGrids(PyObject* self, PyObject* args);
  PyObject* breakElements(PyObject* self, PyObject* args);
  PyObject* splitNGon(PyObject* self, PyObject* args);
  PyObject* splitElement(PyObject* self, PyObject* args);
  PyObject* dualNGon(PyObject* self, PyObject* args);
  PyObject* flipEdges(PyObject* self, PyObject* args);
  PyObject* contractEdges(PyObject* self, PyObject* args);
  PyObject* checkTriMesh(PyObject* self, PyObject* args);
}
# undef FldArrayF
# undef FldArrayI

#endif

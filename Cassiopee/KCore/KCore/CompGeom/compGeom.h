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

#ifndef _KCORE_COMPGEOM_H
#define _KCORE_COMPGEOM_H
# include <vector>
# include <list>
# include "Def/DefTypes.h"
# include "Fld/FldArray.h"
# include "Def/DefFunction.h"
# include "Sort/sort.h"
  
namespace K_COMPGEOM
{
typedef struct {
    E_Int s1, s2;
} Edge;

# include "Triangle.h"

  //======================================================================
  // Barycenters
  //======================================================================
  /* Calcul du barycentre d'un nuage de points
     IN: n: nombre de pts dans le nuage
     IN: x, y, z: coord. du nuage de pts
     OUT: xb, yb, zb: coord. du barycentre
  */
  void barycenter(E_Int n, 
                  E_Float* x, E_Float* y, E_Float* z,
                  E_Float& xb, E_Float& yb, E_Float& zb);

  /* Calcul le barycentre pondere d'un ensemble de points.
     IN: n: nombre de pts dans le nuage
     IN: x, y, z: coord. du nuage de pts
     IN: w: ponderation pour chaque point
     OUT: xb, yb, zb: coord. du barycentre
  */
  void weightedBarycenter(E_Int n,
                          E_Float* x, E_Float* y, E_Float* z, E_Float* w,
                          E_Float& xb, E_Float& yb, E_Float& zb);

  /* Computes the barycenter of a NGON face 
     IN : posface: position of the face in NGon connectivity
     IN : ptrNG: pointer on the connectivity NGon
     IN : xt, yt, zt: coordinates of vertices
     OUT: coordinates of the barycenter of the face */
  void getNGONFaceBarycenter(E_Int posface, E_Int* ptrNG, 
                             E_Float* xt, E_Float* yt, E_Float* zt,
                             E_Float& xbf, E_Float& ybf, E_Float& zbf);
  /* Computes the barycenter of a NGON element
     IN: noet: number of the element
     IN: posEltsp: position of elts in connectivity NFace
     IN: posFacesp: position of faces in connectivity NGon
     IN: posface: position of the face in NGon connectivity
     IN: ptrNF: connectivity Elt/Faces (NFace connectivity)
     IN: xt, yt, zt: coordinates of vertices
     OUT: coordinates of the barycenter of the face */
  void getNGONEltBarycenter(E_Int noet, E_Int* posEltsp, E_Int* posFacesp, 
                            E_Int* ptrNF, E_Int* ptrNG,
                            E_Float* xt, E_Float* yt, E_Float* zt,
                            E_Float& xbg, E_Float& ybg, E_Float& zbg);

  //======================================================================
  // Beziers, splines, nurbs... 
  //======================================================================
  /* Evalue une courbe de Bezier a partir de points de controle
     IN: n: nb de pts de controle
     IN: N: nb de pts totaux dans la courbe de Bezier
     IN: xt, yt, zt: pts de controle de la courbe
     OUT: coord: tableau des coordonnees des pts de la courbe de Bezier 
          resultante (doit etre deja alloue)
  */ 
  void bezier(E_Int n, E_Int N, 
              E_Float* xt, E_Float* yt, E_Float* zt,
              K_FLD::FldArrayF& coord);
  /* Evalue une courbe de Bezier reguliere a partir de pts de controle 
     IN: n: nb de pts de controle
     IN: density: nb de pts par unite de longueur dans la courbe de Bezier
     IN: N: nbre de pts sur la courbe finale (si N>0, density est ignore)
     IN: xt, yt, zt: pts de controle de la courbe
     OUT: coord: tableau des coordonnees de la bezier reguliere
  */
  void regularBezier(E_Int n, E_Int N, E_Float density, 
                     E_Float* xt, E_Float* yt, E_Float* zt,
                     K_FLD::FldArrayF& coord);
  /* Evalue une surface de Bezier a partir du maillage imxjm des points de 
     controle
     IN: imxjm: nb de pts de controle
     IN: NxM: nombre de points sur la surface de Bezier
     IN: xt, yt, zt: pts de controle de la surface
     OUT: coord: tableau des coordonnees des pts de la courbe de Bezier 
          resultante (doit etre deja alloue) 
  */
  void bezier2D(E_Int im, E_Int jm, E_Int N, E_Int M, 
                E_Float* xt, E_Float* yt, E_Float* zt, 
                K_FLD::FldArrayF& coord);
  /* Evalue une surface de Bezier reguliere a partir du maillage imxjm 
     des points de controle
     IN: imxjm: nb de pts de controle
     IN: NxM: nombre de points sur la surface de Bezier
     IN: densite: si > 0, remplace N, M
     IN: xt, yt, zt: pts de controle de la surface
     OUT: coord: tableau des coordonnees des pts de la courbe de Bezier 
          resultante (doit etre deja alloue)
     OUT: niout, njout: nbre de pts effectivement cree dans la surface
  */
  void regularBezier2D(E_Int n, E_Int m, E_Int N, E_Int M, 
                       E_Float density, 
                       E_Float* xt, E_Float* yt, E_Float* zt,
                       K_FLD::FldArrayF& coord, E_Int& niout, E_Int& njout);

  /* interne */
  E_Int factorielle(E_Int i);
  E_Int combinations(E_Int n, E_Int i);
  E_Float Bernstein(E_Int i, E_Int n, E_Float t);

  /* Evalue une courbe B-spline a partir des points de controle
     IN: n: nb de pts de controle
     IN: ordern: Ordre de la spline 
     IN: N: nombre de points totaux sur la spline
     IN: xt,yt,zt: points de controle
     OUT: coord: points de la spline
  */
  void spline(E_Int n, E_Int ordern, E_Int N, 
              E_Float* xt, E_Float* yt, E_Float* zt, 
              K_FLD::FldArrayF& coord);

  void regularSpline(E_Int n, E_Int ordern, E_Int N, E_Float density,
                     E_Float* xt, E_Float* yt, E_Float* zt,
                     K_FLD::FldArrayF& coord);

  /* Evalue une surface B-spline a partir des points de controle
     IN: imxjm: nb de pts de controle
     IN: NxM: nombre de points totaux sur la spline dans chaque direction
     IN: ordern,m: Ordre de la spline dans chaque direction
     IN: xt,yt,zt: points de controle
     OUT: coord: points de la spline
  */
  void spline2D(E_Int im, E_Int jm, E_Int ordern, E_Int N, 
                E_Int orderm, E_Int M, 
                E_Float* xt, E_Float* yt, E_Float* zt, 
                K_FLD::FldArrayF& coord);
  /* Evalue une surface B-spline reguliere a partir des points de controle
     IN: imxjm: nb de pts de controle
     IN: NxM: nombre de points totaux sur la spline dans chaque direction
     IN: ordern,m: Ordre de la spline dans chaque direction
     IN: xt,yt,zt: points de controle
     OUT: coord: points de la spline
  */
  void regularSpline2D(E_Int im, E_Int jm, E_Int ordern, E_Int N, 
                       E_Int orderm, E_Int M, E_Float density, 
                       E_Float* xt, E_Float* yt, E_Float* zt, 
                       K_FLD::FldArrayF& coord, E_Int& niout, E_Int& njout);
  /* interne */
  void evalSpline(E_Float t, K_FLD::FldArrayF& x, E_Int n, E_Int c, 
                  E_Float* xt, E_Float* yt, E_Float* zt,
                  E_Float& xo, E_Float& yo, E_Float& zo);

  /* Evalue une nurbs a partir des points de controle
     IN: im: nb de pts de controle
     IN: N,M: nombre de points totaux sur la nurbs dans chaque direction
     IN: ordern,m: Ordre de la nurbs dans chaque direction
     IN: xt,yt,zt: points de controle
     IN: W: poids des points de controle
     OUT: coord: points de la la nurbs 
  */
  void nurbs(E_Int im, E_Int ordern, E_Int N, 
             E_Float* xt, E_Float* yt, E_Float* zt, E_Float* W, 
             K_FLD::FldArrayF& coord);

  void regularNurbs(E_Int im, E_Int ordern, E_Int N, E_Float density, 
                    E_Float* xt, E_Float* yt, E_Float* zt, E_Float* W, 
                    K_FLD::FldArrayF& coord);
   /* Evalue une surface nurbs a partir des points de controle
     IN: im,jm: nb de pts de controle
     IN: N,M: nombre de points totaux sur la nurbs dans chaque direction
     IN: ordern,m: Ordre de la nurbs dans chaque direction
     IN: xt,yt,zt: points de controle
     IN: W: poids des points de controle
     OUT: coord: points de la la nurbs 
  */
  void nurbs2D(E_Int im, E_Int jm, E_Int ordern, E_Int N, E_Int orderm, 
               E_Int M, E_Float* xt, E_Float* yt, E_Float* zt, E_Float* W,
               K_FLD::FldArrayF& coord);
  /* interne */
  void evalNurbs(E_Float t, K_FLD::FldArrayF& x, E_Int n, E_Int c, 
                 K_FLD::FldArrayF& N);

  void regularNurbs2D(E_Int im, E_Int jm, E_Int ordern, E_Int N, 
                      E_Int orderm, E_Int M, E_Float density, 
                      E_Float* xt, E_Float* yt, E_Float* zt, E_Float* W, 
                      K_FLD::FldArrayF& coord, E_Int& niout, E_Int& njout);

  //======================================================================
  // delaunay triangulation 
  //======================================================================
  /* Triangulation de delaunay dans un plan
     IN: coefa, coefbm coefc, coefd: defini le plan 
     (coefa x + coefb y + coefc z + coefd = 0)
     IN: coord: coordonnees des pts a trianguler
     OUT: connect: connectivite calculee
     IN: keepBB: 1: on garde la BB, 0 on la supprime
   */
  void delaunay(E_Float coefa, E_Float coefb, E_Float coefc, 
                E_Float coefd,
                K_FLD::FldArrayF& coord, K_FLD::FldArrayI& connect, 
                E_Int keepBB);

  /* Interne : tile the bounding box of the cloud of points in triangles */
  void compAndTileBoundingBox(E_Float coefa, E_Float coefb, 
                              E_Float coefc, E_Float coefd,
                              K_FLD::FldArrayF& field, 
                              std::list<K_COMPGEOM::Triangle*>& triangles);

  /* Interne : insert the triangle edges AB, BC and CA in the list of edges */
  void insertTriangleEdgesInList(K_COMPGEOM::Triangle* oneTriangle,
                                 std::list<Edge*>& edges);
  /* Interne : create new Triangles with list of triangles and current point
     + Reinit edges list */
  void insertNewTriangles(E_Float coefa, E_Float coefb, 
                          E_Float coefc, E_Float coefd,
                          E_Int ind, K_FLD::FldArrayF& field, 
                          std::list<Edge*>& edges,
                          std::list<K_COMPGEOM::Triangle*>& triangles);

  /* Interne: remove triangles using bbox vertices */
  void removeTrianglesWithBBoxVertices(
    E_Int sizeIni, std::list<K_COMPGEOM::Triangle*>& triangles);

  /* Interne: build the connectivity */
  void buildConnectivity(std::list<K_COMPGEOM::Triangle*>& triangles,
                         K_FLD::FldArrayI& connect);
  
  //===========================================================================
  // Triangles / quads
  //===========================================================================
  /* Calcul l'angle entre deux triangles : si les triangles
     partagent une arete commune : retourne l'angle en degres dans [0,360].
     Sinon retourne -1000. */
  E_Float getAlphaAngleBetweenTriangles( 
    E_Float* ptA1, E_Float* ptB1, E_Float* ptC1,
    E_Float* ptA2, E_Float* ptB2, E_Float* ptC2);
  /* Calcul l'angle entre deux quads : si les quads
     partagent une arete commune : retourne l'angle en degres dans [0,360].
     Sinon retourne -1000. */
  E_Float getAlphaAngleBetweenQuads(
    E_Float* ptA1, E_Float* ptB1, E_Float* ptC1, E_Float* ptD1, 
    E_Float* ptA2, E_Float* ptB2, E_Float* ptC2, E_Float* ptD2);
  /*  Calcul l'angle entre 2 polygones p1 et p2. Ils doivent etre
      numerotes en tournant. */
  E_Float getAlphaAngleBetweenPolygons(std::vector<E_Float*>& p1, 
                                       std::vector<E_Float*>& p2);

  /* Verifie si un quad est convexe ou non.
     Retourne 0 : quad convexe
     retourne i : quad non convexe au sommet i
     retourne -1 : quad degenere */
  E_Int checkQuadConvexity(E_Float* pt1, E_Float* pt2, 
                           E_Float* pt3, E_Float* pt4);

  /* Calcul l'angle entre deux segments : si les segments
     partagent une arete commune : retourne l'angle en degres dans [0,360].
     dirVect est le vecteur orthogonal (approximatif) au plan contenant les 2 segments
     Si pas de pt commun entre les 2 segments, retourne -1000. */
  E_Float getAlphaAngleBetweenBars(E_Float* ptA1, E_Float* ptB1, 
                                   E_Float* ptA2, E_Float* ptB2,   E_Float* dirVect);

  /* Intersection de triangles (algo de Sam : cf Nuga/include/Triangle.h) 
     Retourne -1 si intersection entre 2 triangles non coplanaires
     Retourne  0 si pas d'intersection 
  */
  E_Int crossIntersectionOfTriangles(
    E_Float* ptA1, E_Float* ptB1, E_Float* ptC1,
    E_Float* ptA2, E_Float* ptB2, E_Float* ptC2,
    E_Float eps);
  /* Calcule si 2 triangles s'intersectent : algo de Moller
     IN : ptA, ptB, ptC sommets des 2 triangles 
     OUT :  0: pas d intersection
            1: intersection sur un segment interne aux triangles
           -1: intersection en un point commun
           -2: intersection sur une arete commune
           -3: coplanaires et intersectant sur des edges
           -4: coplanaires et un triangle est interne a l'autre */
  E_Int trianglesIntersection(E_Float* ptA1, E_Float* ptB1, E_Float* ptC1,
                              E_Float* ptA2, E_Float* ptB2, E_Float* ptC2,
                              E_Float eps=1.e-10);
  /* Determine l'intersection entre deux segments situes sur le meme 
     plan z=cste */
  E_Int getTypeOfSegmentIntersection(E_Float* P1T1, E_Float* P2T1,
                                     E_Float* P1T2, E_Float* P2T2,
                                     E_Float eps);
  /* Teste si 2 triangles coplanaires s'intersectent 
     Retourne  0 si pas d intersection 
     Retourne -1 si intersection en au moins un point et pas sommet 
     Retourne  1 si intersection en un sommet 
     Retourne  2 si intersection en une arete du triangle
     Retourne 10 si un triangle est contenu dans l'autre */
  E_Int testCoplanarTrianglesIntersection(
    E_Float* ptA1, E_Float* ptB1, E_Float* ptC1,
    E_Float* ptA2, E_Float* ptB2, E_Float* ptC2,
    E_Float* N1, E_Float* N2, E_Float eps = 1.e-10);
  /* Interne */
  E_Int testCommonEdgeIntersection(
    E_Float* ptA1, E_Float* ptB1, E_Float* ptC1,
    E_Float* ptA2, E_Float* ptB2, E_Float* ptC2,
    E_Float resA1, E_Float resB1, E_Float resC1, 
    E_Float resA2, E_Float resB2, E_Float resC2, 
    E_Float* N1, E_Float eps);
  /* Interne */
  E_Int testCommonVertexIntersection(
    E_Float* ptA1, E_Float* ptB1, E_Float* ptC1,
    E_Float* ptA2, E_Float* ptB2, E_Float* ptC2,
    E_Float resA1, E_Float resB1, E_Float resC1, 
    E_Float resA2, E_Float resB2, E_Float resC2, 
    E_Float eps);
  /* Calcul si un point est dans un triangle plan (x,y).
     IN: p1, p2, p3: triangle coordinates
     IN: p: tested point
     IN: Ball: null or circum-circle of triangle (xR,yR,R)
     IN: BB: null or bounding box of triangle (xmin,ymin,xmax,ymax)
     
     Retourne: 0: out, 1: in
  */
  E_Int pointInTriangle2D(E_Float* p0, E_Float* p1, E_Float* p2,
                          E_Float* p,
                          E_Float* Ball, E_Float* BB);

  /* Calcul la distance a un triangle.
     IN: p0, p1, p2: triangle coordinates
     IN: p: tested point
     IN: treatment: traitement effectue si p ne se projete pas dans 
     le triangle
     treatment=0, retourne le projete sur le plan du triangle
     treatment=1, retourne le sommet du traiangle le plus proche du projete
     treatment=2, retourne le point de l'edge le plus proche
     OUT: dist2: square distance to projected point
     OUT: in: true, projection falls in the triangle.
                false, projection falls out of triangle.
     OUT: xp, yp, zp: coordinates of projected point.
     return: 0: OK, -1: Failed.
  */
  E_Int distanceToTriangle(E_Float* p0, E_Float* p1, E_Float* p2,
                           E_Float* p, E_Int treatment, 
                           E_Float& dist2, E_Bool& in, 
                           E_Float& xp, E_Float& yp, E_Float& zp,
                           E_Float& sigma0, E_Float& sigma1);
  
  /* Calcul de la distance a un segment (AB)
     IN: p: tested point
     IN: pA, pB: BAR (segment)
     IN: treatement: traitement effectue si p ne se projete pas dans la BAR
     treatment=0, retourne le projete sur la BAR prolongee a l'infini
     treatment=1, retourne le sommet de la BAR le plus proche du projete
     OUT: xp, yp, zp: coord du point projete
     OUT: in: true, le pt projete est dans la BAR,
              false, le pt est dehors
     OUT: dist2: distance au carre du pt p au pt projete.
     Retourne 0: OK
     -1: Failed.
   */
  E_Int distanceToBar(E_Float* pA, E_Float* pB,
                      E_Float* p, E_Int treatment, 
                      E_Float& xp, E_Float& yp, E_Float& zp,
                      E_Bool& in, E_Float& dist2);

  /* Calcule la distance minimale entre deux blocs structurés 
     et retourne les indices correspondants */
  void compMeanDist(const E_Int ni1, const E_Int nj1,
    const E_Float* x1, const E_Float* y1, const E_Float* z1,
    const E_Int ni2, const E_Int nj2,
    const E_Float* x2, const E_Float* y2, const E_Float* z2,
    E_Int& ind1s, E_Int& ind2s, E_Float& dmin);

  /* Analyse l'orientation de deux blocs structurés. 
     Retourne -1 si les normales sont inversées, 1 sinon */
  void rectifyNormals(const E_Int ni1, const E_Int nj1, const E_Int ind1,
    const E_Float* x1, const E_Float* y1, const E_Float* z1,
    const E_Int ni2, const E_Int nj2, const E_Int ind2,
    const E_Float* x2, const E_Float* y2, const E_Float* z2,
    const E_Float distmin,
    E_Int& isopp);

  /* Calcul de l'aire d'un triangle ABC a partir des longueurs a, b et c
     de ses trois cotes (formule de Heron). 
     Retourne 0. si triangle degenere */
  E_Float compTriangleArea(const E_Float a, const E_Float b, const E_Float c);

  /* Calcul le cercle circonscrit a un triangle du plan (x,y)
     IN: p1, p2, p3: triangles coordinates (in plane x,y)
     OUT: pc, R: center and radius of circum-center
     Retourne:
     0: OK
     -1: Failed
  */
  E_Int circumCircle(E_Float* p1, E_Float* p2, E_Float* p3,
                     E_Float* pc, E_Float& R);

  /* Calcul du rayon du cercle circonscrit au triangle P1P2P3. 
     Retourne 0 si le triangle est degenere */
  E_Float circumCircleRadius(E_Float& p1x, E_Float& p1y, E_Float& p1z,
			     E_Float& p2x, E_Float& p2y, E_Float& p2z,
			     E_Float& p3x, E_Float& p3y, E_Float& p3z);
  
  /* Calcul du rayon du cercle inscrit au triangle P1P2P3. 
     Retourne 0 si triangle degenere */
  E_Float inscribedCircleRadius(E_Float* p1, E_Float* p2, E_Float* p3);
  
  /* Calcul l'intersection d'un rayon et d'un triangle.
     IN: p1, p2, p3: triangle coordinates
     IN: pr0, pr1: deux points du rayon.
     Retourne :
     pi: coord. du point d'intersection.
     -1: si le triangle est degenere
     0: pas d'intersection
     1: intersection en un point unique I1.
     2: le rayon et le triangle sont dans le meme plan.
  */
  E_Int intersectRayTriangle(E_Float* p1, E_Float* p2, E_Float* p3,
                             E_Float* pR0, E_Float* pR1,
                             E_Float* pi);

  /* Calcul les coordonnees parametriques d'un point P dans un triangle.
     IN: p1, p2, p3: triangle coordinates
     IN: p: le point
     OUT: s, t: coord. parametriques
     Retourne:
     1: le point est dedans.
     0: le point est dehors
  */
  E_Int computeParamCoord(E_Float* p1, E_Float* p2, E_Float* p3,
                          E_Float* P,
                          E_Float& s, E_Float& t);

  /* Projette selon la direction dirx,diry,dirz un point sur la surface 
     triangulaire definie par ses coordonnees fx2, fy2, fz2 et sa 
     connectivite cn2
     oriented != 0 : the projected point must be in the same direction as the dir vector
     Retourne xo,yo,zo : point de projection
     Retourne le no du triangle si projete trouve, -1 sinon. */
  E_Int projectOneDirWithoutPrecond(E_Float x, E_Float y, E_Float z,
                                    E_Float dirx, E_Float diry, E_Float dirz,
                                    E_Float* fx2, E_Float* fy2, E_Float* fz2,
                                    K_FLD::FldArrayI& cn2, 
                                    E_Float& xo, E_Float& yo, E_Float& zo, E_Int oriented=0);
  /* Idem mais utilise le preconditionement par BBox Tree */
  E_Int projectOneDirWithPrecond(E_Float x, E_Float y, E_Float z,
                                 E_Float dirx, E_Float diry, E_Float dirz,
                                 E_Float* fx2, E_Float* fy2, E_Float* fz2,
                                 K_FLD::FldArrayI& cn2, 
                                 E_Float& xo, E_Float& yo, E_Float& zo, E_Int oriented=0);
  
  /* Projette un point (x,y,z) suivant une direction dirx, diry, dirz
     sur une surface triangulaire definie par ses coordonnees fx2, fy2, fz2, 
     sa connectivite cn2.
     oriented : si vaut 1, le projete doit etre dans le meme sens que (dirx,diry,dirz)
     Retourne xo,yo,zo: point de projection
     Retourne le no du triangle si projete trouve, -1 sinon. */
  E_Int projectDir(E_Float x, E_Float y, E_Float z,
                   E_Float dirx, E_Float diry, E_Float dirz,
                   E_Float* fx2, E_Float* fy2, E_Float* fz2,
                   K_FLD::FldArrayI& cn2, 
                   E_Float& xo, E_Float& yo, E_Float& zo, E_Int oriented=0);
  /* Meme chose mais avec une liste de triangles a parcourir uniqt */
  E_Int projectDir(E_Float x, E_Float y, E_Float z,
                   E_Float dirx, E_Float diry, E_Float dirz,
                   E_Float* fx2, E_Float* fy2, E_Float* fz2,
                   std::vector<E_Int>& indices, K_FLD::FldArrayI& cn2, 
                   E_Float& xo, E_Float& yo, E_Float& zo, E_Int oriented=0);

  /* Projette les pts de coordonnees fx, fy, fz suivant une direction 
     dirx, diry, dirz sur une surface triangulaire definie par ses 
     coordonnees fx2, fy2, fz2 et sa connectivite cn2. 
     oriented : si vaut 1, le projete doit etre dans le meme sens que (dirx,diry,dirz)
     Modifie fx, fy, fz si le projete est trouve. */
  void projectDirWithoutPrecond(
    E_Float nx, E_Float ny, E_Float nz,
    E_Int npts, E_Int nelts2,  K_FLD::FldArrayI& cn2,
    E_Float* fx2, E_Float* fy2, E_Float* fz2,
    E_Float* fx, E_Float* fy, E_Float* fz, E_Int oriented=0);

  /* Meme chose mais avec preconditionnement par bboxTree */
  void projectDirWithPrecond(
    E_Float nx, E_Float ny, E_Float nz,
    E_Int npts, E_Int nelts2, K_FLD::FldArrayI& cn2,
    E_Float* fx2, E_Float* fy2, E_Float* fz2,
    E_Float* fx, E_Float* fy, E_Float* fz, E_Int oriented=0);

  /* Meme chose mais avec une liste de vecteur de pts a projeter.
     sizet est le nombre de pts de chaque vecteur fxt,fyt,fzt */
  void projectDirWithPrecond(
    E_Float nx, E_Float ny, E_Float nz,
    K_FLD::FldArrayI& cn2, E_Float* fx2, E_Float* fy2, E_Float* fz2,
    std::vector<E_Int>& sizet,
    std::vector<E_Float*>& fxt, std::vector<E_Float*>& fyt, std::vector<E_Float*>& fzt,
    E_Int oriented=0);

  /* Meme chose mais avec des directions variables pour chaque point a
   projeter. */
  void projectDirWithPrecond(
    std::vector<E_Float*>& nxt, std::vector<E_Float*>& nyt, std::vector<E_Float*>& nzt,
    K_FLD::FldArrayI& cn2, E_Float* fx2, E_Float* fy2, E_Float* fz2,
    std::vector<E_Int>& sizet,
    std::vector<E_Float*>& fxt, std::vector<E_Float*>& fyt, std::vector<E_Float*>& fzt,
    E_Int oriented=0);

  /* Projette orthogonalement un point (x,y,z) sur une surface 
     triangulaire definie par ses coordonnees fx2, fy2, fz2, sa 
     connectivite cn2.
     Retourne xo,yo,zo: point de projection
     Retourne le no du triangle si projete trouve, -1 sinon. */
  E_Int projectOrtho(E_Float x, E_Float y, E_Float z,
                     E_Float* fx2, E_Float* fy2, E_Float* fz2,
                     K_FLD::FldArrayI& cn2, 
                     E_Float& xo, E_Float& yo, E_Float& zo,
                     E_Float* p0, E_Float* p1, E_Float* p2, E_Float* p);
  /* Meme chose mais indices est la liste des elements de cn2 a 
     parcourir uniquement */
  E_Int projectOrthoPrecond(E_Float x, E_Float y, E_Float z,
                            E_Float* fx2, E_Float* fy2, E_Float* fz2,
                            std::vector<E_Int>& indices, K_FLD::FldArrayI& cn2, 
                            E_Float& xo, E_Float& yo, E_Float& zo,
                            E_Float* p0, E_Float* p1, E_Float* p2, E_Float* p);

  /* Projection orthogonale des pts de coordonnees fx,fy,fz a projeter 
     orthogonalement sur une surface de cn cn2, coords fx2, fy2, fz2. 
     Optimisation par kdtree et bboxtree */
  void projectOrthoWithoutPrecond(E_Int npts, K_FLD::FldArrayI& cn2, 
                                  E_Float* fx2, E_Float* fy2, E_Float* fz2,
                                  E_Float* fx, E_Float* fy, E_Float* fz);
  
  void projectOrthoWithPrecond(E_Int npts, K_FLD::FldArrayI& cn2, 
                               E_Int posx2, E_Int posy2, E_Int posz2,
                               K_FLD::FldArrayF& f2, 
                               E_Float* fx, E_Float* fy, E_Float* fz);

  /* Meme chose mais sur une liste de surfaces 
     IN: posx2, posy2, posz2: position de x,y,z de la surface de projection
     IN: f2, cn2 coordonnees de la surface de projection et connectivite 
     IN: sizet: taille de chaque tableau de fxt, fyt, fzt
     IN: fxt, fyt, fzt: liste des coordonnees pour toutes les zones a projeter */
  void projectOrthoWithPrecond(
    E_Int posx2, E_Int posy2, E_Int posz2, 
    K_FLD::FldArrayI& cn2, K_FLD::FldArrayF& f2, std::vector<E_Int>& sizet,
    std::vector<E_Float*>& fxt, std::vector<E_Float*>& fyt, std::vector<E_Float*>& fzt);

  void projectRay(E_Int npts, E_Int nelts2,
                  E_Float Px, E_Float Py, E_Float Pz,
                  E_Float* fx2, E_Float* fy2, E_Float* fz2,
                  K_FLD::FldArrayI& cn2, 
                  E_Float* fx, E_Float* fy, E_Float* fz);
  void projectRay(E_Int nelts2,
                  E_Float Px, E_Float Py, E_Float Pz,
                  E_Float* fx2, E_Float* fy2, E_Float* fz2,
                  K_FLD::FldArrayI& cn2, 
                  std::vector<E_Int>& sizet, std::vector<E_Float*>& fxt, 
                  std::vector<E_Float*>& fyt, std::vector<E_Float*>& fzt);

  //===========================================================================
  // Lines
  //===========================================================================
  /* Intersecte 2 segments 2D (dans le plan x,y)
     IN: ps1a, ps1b: coordonnees des points du segment 1
     IN: ps2a, ps2b: coordonnees des points du segment 2
     OUT: pi0, pi1: coordonnees des points d'intersections (si il y en a)
     Retourne :
     0: segments disjoints
     1: intersection en un point unique (pi0)
     2: les deux segments se recouvrent de pi0 a pi1.
     Remarque: tous les points doivent etre alloues a E_Float[3].
  */
  E_Int intersect2Segments(E_Float* ps1a, E_Float* ps1b,
                           E_Float* ps2a, E_Float* ps2b,
                           E_Float* pi0, E_Float* pi1);
  /* Interne */
  E_Int inSegment(E_Float* p, E_Float* psa, E_Float* psb);

  /* Point in segment
     IN: p0, p1: segment
     IN: p: pt a tester
     IN: eps: tolerance
     OUT: 1: p est le pt p0
     OUT: 2: p est le pt p1
     OUT: 3: p est dans p0-p1
     OUT: 0: sinon
  */
  E_Int pointInSegment(E_Float* p0, E_Float* p1, E_Float* p, 
                       E_Float eps=1.e-6);

  //===========================================================================
  // Polygones
  //===========================================================================
  /* Calcul si un point est dans un polygone plan (x,y).
     IN: coord: coordonnees du polygone
     IN: connectivite du polygone
     IN: p: tested point
     IN: Ball: null or circum-circle of triangle (xR,yR,R)
     IN: BB: null or bounding box of triangle (xmin,ymin,xmax,ymax)
     
     return: 0: out, 1: in
  */
  E_Int pointInPolygon2D(E_Float* xt, E_Float* yt, E_Float* zt,
                         K_FLD::FldArrayI& connect,
                         E_Float* p,
                         E_Float* Ball, E_Float* BB);
  /* Calcul si un point est dans un polygone plan (x,y) a epsilon pres.
     Si le point est dans le poylgone a epsilon, pres cette routine
     retourne 1 (in).
     IN: coord: coordonnees du polygone
     IN: connectivite du polygone
     IN: p: tested point
     IN: epsilon: tolerance
     IN: Ball: null or circum-circle of triangle (xR,yR,R)
     IN: BB: null or bounding box of triangle (xmin,ymin,xmax,ymax)
     
     return: 0 (out), 1 (in)
  */
  E_Int pointInPolygon2D(E_Float* xt, E_Float* yt, E_Float* zt,
                         K_FLD::FldArrayI& connect,
                         E_Float* p, E_Float epsilon,
                         E_Float* Ball, E_Float* BB);

  //===========================================================================
  // Bounding boxes et Cartesian Elements Bounding boxes (CEBB)
  //===========================================================================
  /* Bounding box d'une grille structuree */
  void boundingBoxStruct(
   E_Int im, E_Int jm, E_Int km, 
   E_Float* x, E_Float* y, E_Float* z,
   E_Float& xmin, E_Float& ymin, E_Float& zmin,
   E_Float& xmax, E_Float& ymax, E_Float& zmax);

  /* Bounding box d'une grille non structuree */ 
  void boundingBoxUnstruct(
   const E_Int npts, const E_Float* xt, const E_Float* yt, const E_Float* zt,
   E_Float& xmin, E_Float& ymin, E_Float& zmin,
   E_Float& xmax, E_Float& ymax, E_Float& zmax);

  /* Bounding box globale d'une liste de grilles */
  void globalBoundingBox(
    std::vector<E_Int>& posxt, std::vector<E_Int>& posyt, 
    std::vector<E_Int>& poszt, std::vector<K_FLD::FldArrayF*>& listOfFields,
    E_Float& xmin, E_Float& ymin, E_Float& zmin,
    E_Float& xmax, E_Float& ymax, E_Float& zmax);

  /* Bounding box de toutes les cellules d'une grille structuree
     IN: im, jm, km: dimensions de l'array definissant la grille
     IN: coord: coordonnees de la grille
     OUT: bbox(ncells,6): xmin, ymin, zmin, xmax,ymax,zmax
     bbox est alloue ici.
  */
  void boundingBoxOfStructCells(E_Int im, E_Int jm, E_Int km,
                                E_Float* x, E_Float* y, E_Float* z,
                                K_FLD::FldArrayF& bbox);
  /* Bounding box de toutes les cellules d'une grille non structuree
     IN: connect: connectivite de la grille
     IN: coord: coordonnees de la grille
     OUT: bbox(nelts,6): xmin, ymin, zmin, xmax,ymax,zmax
     bbox est alloue ici.
  */
  void boundingBoxOfUnstrCells(K_FLD::FldArrayI& connect,
                               E_Float* xt, E_Float* yt, E_Float* zt,
                               K_FLD::FldArrayF& bbox);

   /* Bounding box de toutes les cellules d'une grille NGon
     IN: connect: connectivite de la grille
     IN: coord: coordonnees de la grille
     OUT: bbox(nelts,6): xmin, ymin, zmin, xmax,ymax,zmax
     bbox est alloue ici.
  */
  void boundingBoxOfNGonCells(K_FLD::FldArrayI& connect,
    E_Float* xt, E_Float* yt, E_Float* zt,
    K_FLD::FldArrayF& bbox);

  /* Bounding box d'une cellule issue d'une grille structuree. */
  void boundingBoxOfStructCell(E_Int ind, E_Int im, E_Int jm, E_Int km,
    E_Float* x, E_Float* y, E_Float* z,
    E_Float& xmin, E_Float& ymin, E_Float& zmin,
    E_Float& xmax, E_Float& ymax, E_Float& zmax, E_Int loc);

  /* Bounding box d'une cellule issue d'une grille non structuree. */
  void boundingBoxOfUnstrCell(
    E_Int noet, K_FLD::FldArrayI& connect, 
    E_Float* xt, E_Float* yt, E_Float* zt, 
    E_Float& xmin, E_Float& ymin, E_Float& zmin,
    E_Float& xmax, E_Float& ymax, E_Float& zmax);

  /* Bounding box d'une cellule issue d'une grille NGon. */
  void boundingBoxOfNGonCell(
    E_Int noet, K_FLD::FldArrayI& connect, 
    E_Float* xt, E_Float* yt, E_Float* zt,
    E_Float& xmin, E_Float& ymin, E_Float& zmin,
    E_Float& xmax, E_Float& ymax, E_Float& zmax);

  /* Intersection des bounding boxes de 2 grilles structurees
     IN: ni1, nj1, nk1: dimension de la grille1
     IN: posx1, posy1, posz1: positions de x,y,z ds f1
     IN: f1: champ defini sur la grille 1 contenant les coordonnees
     IN: ni2, nj2, nk2: dimension de la grille2 
     IN: posx2, posy2, posz2: positions de x,y,z ds f2
     IN: f2: champ defini sur la grille 2 contenant les coordonnees
     OUT: xmin1, xmax1, ymin1, ymax1, zmin1, zmax1 : bbox de grille 1
     OUT: xmin2, xmax2, ymin2, ymax2, zmin2, zmax2 : bbox de grille 2
     RETOURNE:
     1 si bbox s'intersectent, 
     0 sinon   */
  E_Int compBoundingBoxIntersection(
    E_Int ni1, E_Int nj1, E_Int nk1, 
    E_Int posx1, E_Int posy1, E_Int posz1, K_FLD::FldArrayF& f1,
    E_Int ni2, E_Int nj2, E_Int nk2, 
    E_Int posx2, E_Int posy2, E_Int posz2, K_FLD::FldArrayF& f2, 
    E_Float& xmin1, E_Float& xmax1, E_Float& ymin1, E_Float& ymax1, 
    E_Float& zmin1, E_Float& zmax1, E_Float& xmin2, E_Float& xmax2, 
    E_Float& ymin2, E_Float& ymax2, E_Float& zmin2, E_Float& zmax2,
    E_Float tol=1.e-10);
 
  /* Recherche si les CEBB  de deux grilles structurees s'intersectent
     IN: ni1, nj1, nk1: dimension de la grille1
     IN: posx1, posy1, posz1: positions de x,y,z ds f1
     IN: f1: champ defini sur la grille 1 contenant les coordonnees
     IN: ni2, nj2, nk2: dimension de la grille2 
     IN: posx2, posy2, posz2: positions de x,y,z ds f2
     IN: f2: champ defini sur la grille 2 contenant les coordonnees
     Retourne 1 si intersection, 0 sinon */
  E_Int compCEBBIntersection(E_Int ni1, E_Int nj1, E_Int nk1, 
                             E_Int posx1, E_Int posy1, E_Int posz1,
                             K_FLD::FldArrayF& f1,
                             E_Int ni2, E_Int nj2, E_Int nk2,
                             E_Int posx2, E_Int posy2, E_Int posz2,
                             K_FLD::FldArrayF& f2,
                             E_Float tol=1.e-10);

  /* Intersection de 2 bounding boxes
     IN: xmin1, xmax1, ymin1, ymax1, zmin1, zmax1: bbox 1
     IN: xmin2, xmax2, ymin2, ymax2, zmin2, zmax2: bbox 2
     RETOURNE:
     1 si bbox s'intersectent, 
     0 sinon   */
  E_Int compBoundingBoxIntersection(E_Float xmin1, E_Float xmax1, 
                                    E_Float ymin1, E_Float ymax1, 
                                    E_Float zmin1, E_Float zmax1, 
                                    E_Float xmin2, E_Float xmax2, 
                                    E_Float ymin2, E_Float ymax2, 
                                    E_Float zmin2, E_Float zmax2,
                                    E_Float tol=1.e-10);

  /* Calcul du tableau des elements cartesiens dans la direction dir
     Comme il s'agit de cellules cartesiennes seules sont stockees 
     [xmin1,ymin1,zmin1,xmax1,ymax1,zmax1]
     IN: dir: direction de la hauteur
     IN: ni, nj, nk: dimensions de la grille
     IN: xmin,...: bbox de la grille
     IN: coord: coordonnees de la grille
     OUT: cartEltArray: tableau contenant les min et max de chq elt cart
     retourne 1 si elts cartesiens possibles ds la direction dir, 0 sinon
  */
  short 
  compCartEltsArray(E_Int dir, E_Int ni, E_Int nj, E_Int nk,
                    E_Float xmin, E_Float ymin, E_Float zmin,
                    E_Float xmax, E_Float ymax, E_Float zmax,
                    K_FLD::FldArrayF& coord, K_FLD::FldArrayF& cartEltArray);
  
  /* Calcul des dimensions dim1 et dim2 des elements cartesiens d'un bloc
     selon la direction dir */
  void compDimOfCartElts(E_Int dir, E_Int ni, E_Int nj, E_Int nk,
                         E_Float xmin, E_Float ymin, E_Float zmin,
                         E_Float xmax, E_Float ymax, E_Float zmax,
                         K_FLD::FldArrayF& coord,
                         E_Int& nbElts1, E_Int& nbElts2);

  /* Calcul de la hauteur des elements cartesiens dans 1 direction donnee
     IN: dir: direction de la hauteur des elts : 1-x , 2-y, 3-z
     IN: nbElts1, nbElts2: nb d elets ds les autres directions
     IN: ni, nj, nk: dimensions de la grille 
     IN: xmin, xmax, ymin, ymax, zmin, zmax: bounding box de la grille
     IN: coord: coordonnees de la grille
     OUT: cartmin, cartmax: elements cartesiens ds la direction 
     retourne 0 si impossible ds la direction consideree
  */
  short compCEBBHeightInDir(E_Int dir, E_Int nbElts1, E_Int nbElts2,
                            E_Int ni, E_Int nj, E_Int nk, 
                            E_Float xmin, E_Float ymin, E_Float zmin,
                            E_Float xmax, E_Float ymax, E_Float zmax,
                            K_FLD::FldArrayF& coord,
                            K_FLD::FldArrayF& cartMin, 
                            K_FLD::FldArrayF& cartMax);

  /* Calcul si un point est dans une bounding box
     IN: xmin, xmax,...: bounding box
     IN: p tested point
     retourne: 0 (out), 1 (in)
  */
  E_Int pointInBB(E_Float xmin, E_Float ymin, E_Float zmin,
                  E_Float xmax, E_Float ymax, E_Float zmax,
                  E_Float* p, E_Float tol=1.e-10);

  //===========================================================================
  // Curvature
  //===========================================================================
  /* Calcul de la courbure pour un i-array:
     IN: npts: nb de pts du i-array
     IN: xt, yt, zt: coordonnees de la courbe
     OUT: courbure aux pts de la courbe. */
  void compCurvature(E_Int npts,
                     E_Float* xt, E_Float* yt, E_Float* zt,
                     K_FLD::FldArrayF& curv);

  /* 
     Calcul de l'angle entre les segments pour un "i-array"
     IN: npts: nbre de pts du i-array
     IN: xt, yt, zt: coord. des pts du i-array
     IN: dirVect est le vecteur orthogonal (approximativement) au plan moyen contenant l'i-array
     OUT: angle: angle entre 0 et 360 degres.
  */
  void compCurvatureAngle(E_Int npts, 
                          E_Float* xt, E_Float* yt, E_Float* zt, 
                          E_Float* dirVect,
                          K_FLD::FldArrayF& angle);

  /* Calcul de l'angle entre les segments pour une "BAR"
     IN: npts: nbre de pts de la BAR
     IN: xt, yt, zt: coord. des pts de la BAR
     IN: cn : connectivite vertex-elts
     IN: dirVect est le vecteur orthogonal (approximativement) au plan moyen contenant la BAR array
     OUT: angle entre 0 et 360 degres.
  */
  void compCurvatureAngleForBar(E_Int npts, 
                                E_Float* xt, E_Float* yt, E_Float* zt,
                                K_FLD::FldArrayI& cn, 
                                E_Float* dirVect,
                                K_FLD::FldArrayF& angle);
  
  /* Calcul de l'angle pour un "TRI"
     IN: npts: nbre de vertex dans le maillage TRI
     IN: xt, yt, zt: coord. des pts du TRI
     IN: cn: connectivite du TRI
     OUT: angle pour chaque element et pour chaque arete. 
     Retourne 1: calcul angle correct
     Retourne 0: echec calcul
  */
  E_Int compCurvatureAngleForTri(E_Int npts,
                                 E_Float* xt, E_Float* yt, E_Float* zt,
                                 K_FLD::FldArrayI& cn, 
                                 K_FLD::FldArrayF& angle);

  /* calcul de la hauteur li�e � la courbure pour des i-arrays */
  void compStructCurvatureHeight1D(E_Int im, E_Float* xt, E_Float* yt, E_Float* zt, 
                                   E_Float* hmaxt);
  /* calcul de la hauteur li�e � la courbure pour des (i,j)-arrays */
  void compStructCurvatureHeight2D(E_Int im, E_Int jm, 
                                   E_Float* xt, E_Float* yt, E_Float* zt, 
                                   E_Float* hmaxt);
  /* calcul de la hauteur li�e � la courbure pour des BAR*/
  void compCurvatureHeightForBAR(E_Int npts, E_Float* xt, E_Float* yt, E_Float* zt, 
                                 K_FLD::FldArrayI& cn, E_Float* hmaxt);
  /* calcul de la hauteur li�e � la courbure pour des surfaces TRI et QUAD */
  void compCurvatureHeightForTRIQUAD(E_Int npts, E_Float* xt, E_Float* yt, E_Float* zt, 
                                     K_FLD::FldArrayI& cn,
                                     E_Float* hmaxt);

  /* calcul du max/min/ratio des edges d'une cellule */
  E_Int getEdgeLength(E_Float* x, E_Float* y, E_Float* z,
                      E_Int im, E_Int jm, E_Int km, 
                      K_FLD::FldArrayI* cn, char* eltType,
                      E_Int dim, E_Int type, E_Float* out);

  // Make a rotation of a mesh
  // IN: dim        : mesh size
  // IN: teta       : angle
  // IN: center     : center of rotation
  // IN: axis       : rotation vector
  // OUT: xo, yo, zo : rotated mesh point
  void rotateMesh(
    const E_Int dim, const E_Float teta,
    const E_Float* center, const E_Float* axis,
    E_Float* xo, E_Float* yo, E_Float* zo);
  // Make a rotation of a mesh - same function as rotateMesh, but different interface and omp present
  // IN: npts       : mesh size
  // IN: teta       : angle
  // IN: xc, yc, zc : center of rotation
  // IN: nx, ny, nz : rotation vector
  // IN: x, y, z    : mesh coordinate
  // OUT: xo, yo, zo : rotated mesh point
  void rotateMesh2(
    const E_Int npts, const E_Float teta,
    const E_Float xc, const E_Float yc, const E_Float zc,
    const E_Float nx, const E_Float ny, const E_Float nz,
    const E_Float* x,const E_Float* y, const E_Float* z,
    E_Float* xo, E_Float* yo, E_Float* zo);

  // Map a 1D distribution over a profile
  // IN: ni         : number of pnts in input line
  // IN: x, y, z    : input line
  // IN: no         : number of pnts in output line
  // IN: d          : distribution
  // OUT: xo, yo, zo : output line
  void onedmap(
    const E_Int ni,
    const E_Float* x, const E_Float* y, const E_Float* z,
    const E_Int no, const E_Float* d,
    E_Float* xo, E_Float* yo, E_Float* zo,
    E_Float* s, E_Float* dx, E_Float* dy, E_Float* dz);

  // Map a 1D distribution over a profile with bar connectivity
  // IN: npts       : number of pnts in input line
  // IN: x, y, z    : input line
  // IN: no         : number of pnts in output line
  // IN: d          : distribution
  // IN: net        : number of elements in input line
  // IN: cn1        : 1st vertex connectivity of input line
  // IN: cn2        : 2nd vertex connectivity of input line
  // IN: neto       : number of element in output line
  // OUT: cn1o      : 1st vertex connectivity of output line
  // OUT: cn2o      : 2nd vertex connectivity of output line
  // OUT: xo, yo, zo : output line
  void onedmapbar(
    const E_Int npts,
    const E_Float* x, const E_Float* y, const E_Float* z,
    const E_Int no, const E_Float* d,
    const E_Int net, const E_Int* cn1, const E_Int* cn2,
    const E_Int neto, E_Int* cn1o, E_Int* cn2o,
    E_Float* xo, E_Float* yo, E_Float* zo,
    E_Float* s, E_Float* dx, E_Float* dy, E_Float* dz);

  // Compute the slope for a 1D line
  // Called internally - onedmap
  // OUT: dx, dy, dz
  void slope(
    const E_Int m,
    const E_Float* x0, const E_Float* y0, const E_Float* z0,
    E_Float* dx, E_Float* dy, E_Float* dz);

  // Compute the slope for a 1D bar with connectivity
  // Called internally - onedmapbar
  // OUT: dx, dy, dz
  void slopebar(
    const E_Int npts, const E_Int net,
    const E_Int* cn1, const E_Int* cn2,
    const E_Float* x0, const E_Float* y0, const E_Float* z0,
    E_Float* dx, E_Float* dy, E_Float* dz);

  // Compute the parametrization for a line
  // Called internally - onedmap
  // OUT: stota, s0
  void paramFunc(
    const E_Int m,
    const E_Float* x0, const E_Float* y0, const E_Float* z0,
    const E_Float* dx, const E_Float* dy, const E_Float* dz,
    E_Float& stota, E_Float* s0);

  // Compute the parametrization for a bar (with connectivity)
  // Called internally - onedmapbar
  // OUT: stota, s0
  void paramFuncBar(
    const E_Int npts, const E_Int net,
    const E_Int* cn1, const E_Int* cn2,
    const E_Float* x0, const E_Float* y0, const E_Float* z0,
    const E_Float* dx, const E_Float* dy, const E_Float* dz,
    E_Float& stota, E_Float* s0);

  // Hermite cubic interpolation helper
  inline E_Float valcub(E_Float a, E_Float b, E_Float c, E_Float d, E_Float t);

  // Interpolate a 1D distribution over a profile
  // Called internally - onedmap
  // OUT: tabx, taby, tabz
  void interp(
    const E_Int im0, const E_Int im, const E_Float stota,
    const E_Float* s0, const E_Float* s,
    const E_Float* tabx0, const E_Float* taby0, const E_Float* tabz0,
    const E_Float* dx0, const E_Float* dy0, const E_Float* dz0,
    E_Float* tabx, E_Float* taby, E_Float* tabz);

  // Interpolate a 1D distribution over a profile with bar connectivity
  // Called internally - onedmapbar
  // OUT: tabx, taby, tabz, cn1, cn2
  void interpbar(
    const E_Int im0, const E_Int im, const E_Float stota,
    const E_Float* s0, const E_Float* s,
    const E_Float* tabx0, const E_Float* taby0, const E_Float* tabz0,
    const E_Float* dx0, const E_Float* dy0, const E_Float* dz0,
    E_Float* tabx, E_Float* taby, E_Float* tabz, 
    const E_Int net, E_Int* cn1, E_Int* cn2);
}

#endif

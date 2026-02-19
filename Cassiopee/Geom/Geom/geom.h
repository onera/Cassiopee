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

#ifndef _GEOM_GEOM_H_
#define _GEOM_GEOM_H_

# include "kcore.h"

namespace K_GEOM
{
  PyObject* naca(PyObject* self, PyObject* args);
  PyObject* line(PyObject* self, PyObject* args);
  PyObject* circle(PyObject* self, PyObject* args);
  PyObject* sphere(PyObject* self, PyObject* args);
  PyObject* cone(PyObject* self, PyObject* args);
  PyObject* torus(PyObject* self, PyObject* args);
  PyObject* triangle(PyObject* self, PyObject* args);
  PyObject* quadrangle(PyObject* self, PyObject* args);
  PyObject* bezier(PyObject* self, PyObject* args);
  PyObject* lineGenerate(PyObject* self, PyObject* args);
  PyObject* lineGenerate2(PyObject* self, PyObject* args);
  PyObject* addSeparationLine(PyObject* self, PyObject* args);
  PyObject* axisym(PyObject* self, PyObject* args);
  PyObject* volumeFromCrossSections(PyObject* self, PyObject* args);
  PyObject* getLength(PyObject* self, PyObject* args);
  PyObject* getDistantIndex(PyObject* self, PyObject* args);
  PyObject* getCurvatureAngle(PyObject* self, PyObject* args);
  PyObject* getCurvatureRadius(PyObject* self, PyObject* args);
  PyObject* getCurvatureHeight(PyObject* self, PyObject* args);
  PyObject* getCurvilinearAbscissa(PyObject* self, PyObject* args);
  PyObject* polyline(PyObject* self, PyObject* args);
  PyObject* spline(PyObject* self, PyObject* args);
  PyObject* nurbs(PyObject* self, PyObject* args);
  PyObject* getSharpestAngleForVertices(PyObject* self, PyObject* args);
  PyObject* getNearestPointIndex(PyObject* self, PyObject* args);
  PyObject* getUV(PyObject* self, PyObject* args);

  /* Calcul des centres des cercles circonscrits de tous les triangles
     IN: coord: coordonnees des vertices
     IN: cn: connectivite element->vertex
     OUT: coordCC: coordonnees des centres des cercles circonscrits aux
     triangles
     coordCC est alloue par la routine */
  void compCCCenters(E_Float* xt, E_Float* yt, E_Float* zt, K_FLD::FldArrayI& cn,
                     K_FLD::FldArrayF& coordCC);

  /* Creation du tableau des aretes du diagramme de Voronoi d une triangulation
     IN: nvertex: nb de vertices dans la triangulation
     IN: cn: connectivite element -> vertex
     OUT: vedges: retourne les 2 elements voisins a une arete du diagramme de Voronoi */
  void compVoronoiEdges(E_Int nv, K_FLD::FldArrayI& cn,
                        K_FLD::FldArrayI& vedges);

  /* Construction des tetraedres formes par une base tri et un vertex de
     type Ti
     type: tableau informant sur le type de tetraedre forme: type0
     IN: type0: type de la triangulation (1 ou 2)
     IN: coordpoly1: coordonnees des sommets du contour 1
     IN: cnpoly1: connectivite BAR du contour1
     IN: cn1: connectivite triangle element->vertex de la triangulation T1
     IN: coordCC1: coordonnees des centres des cercles circonscrits aux triangles de T1
     IN: coord1: coordonnees des points de T1
     IN: coord2: coordonnees des points de T2
     IN/OUT: ne: nb d elements courant
     IN/OUT: nv: nb de vertices courant
     IN/OUT: type: tagge les triangles selon qu ils sont base d un tetraedre de type 1/2 ou non
     IN/OUT: coord: coordonnees des sommets du maillage tetraedrique mis a jour
     IN/OUT: cn: connectivite du maillage tetraedrique
  */
  void compTetraType1(E_Int type0, E_Int posxc1, E_Int posyc1, E_Int poszc1,
                      K_FLD::FldArrayF& fc1, K_FLD::FldArrayI& cnpoly1,
                      K_FLD::FldArrayI& cn1, K_FLD::FldArrayF& coordCC1,
                      K_FLD::FldArrayF& coord1, K_FLD::FldArrayF& coord2,
                      E_Int& nv, E_Int& ne, K_FLD::FldArrayI& type1,
                      K_FLD::FldArrayF& coord, K_FLD::FldArrayI& cn,
                      K_FLD::FldArrayI& type);

  /* Calcul les tetraedres de type T12 : issus de 2 aretes intersectantes
     des diagrammes de Voronoi de T1 et T2
     IN: coordpoly1: coordonnees des sommets du contour 1
     IN: cnpoly1: connectivite BAR du contour1
     IN: cn1: connectivite triangle element->vertex de la triangulation T1
     IN: coordCC1: coordonnees des centres des cercles circonscrits aux triangles de T1
     IN: coord1: coordonnees des points de T1
     IN: coordpoly2: coordonnees des sommets du contour 2
     IN: cnpoly2: connectivite BAR du contour2
     IN: cn2: connectivite triangle element->vertex de la triangulation T2
     IN: coordCC2: coordonnees des centres des cercles circonscrits aux triangles de T2
     IN: coord2: coordonnees des points de T2
     IN: type1: triangles de T1 sont base d un tetraedre de type 1 ou non
     IN: type2: triangles de T2 sont base d un tetraedre de type 2 ou non
     IN/OUT: ne: nb d elements courant
     IN/OUT: nv: nb de vertices courant
     IN/OUT: coord: coordonnees des sommets du maillage tetraedrique mis a jour
     IN/OUT: cn: connectivite du maillage tetraedrique
     IN/OUT: type: type des tetraedres (t1, t2, t12)
  */
  void compTetraType12(E_Int posxc1, E_Int posyc1, E_Int poszc1,
                       K_FLD::FldArrayF& fc1, K_FLD::FldArrayI& cnpoly1,
                       K_FLD::FldArrayI& cn1,
                       K_FLD::FldArrayF& coordCC1, K_FLD::FldArrayF& coord1,
                       E_Int posxc2, E_Int posyc2, E_Int poszc2,
                       K_FLD::FldArrayF& fc2, K_FLD::FldArrayI& cnpoly2,
                       K_FLD::FldArrayI& cn2, K_FLD::FldArrayF& coordCC2,
                       K_FLD::FldArrayF& coord2,
                       K_FLD::FldArrayI& type1, K_FLD::FldArrayI& type2,
                       E_Int& nv, E_Int& ne,
                       K_FLD::FldArrayF& coord, K_FLD::FldArrayI& cn,
                       K_FLD::FldArrayI& type);

  /* Determine les indices des vertices dans la triangulation de Delaunay
     communs aux 2 elements et1 et et2. Retourne les indices des sommets
     ind demarre a 1
     IN: et1, et2: numero des elements
     IN: cEV: connectivite element->vertex
     OUT: ind: tableau des indices des 2 sommets communs aux triangles et1 et et2
  */
  void getDelaunayVertices(E_Int et1, E_Int et2, K_FLD::FldArrayI& cEV,
                           K_FLD::FldArrayI& ind);

  /* Perturbe la triangulation T1 si les coordonnees (x,y) de 2 points,
     l'un etant dans T1, l'autre dans T2, sont identiques.
     IN: posx1, posy1, posz1: position des coordonnees ds f1
     IN: posx2, posy2, posz2: position des coordonnees ds f2
     IN: f1: coordonnees des pts de la triangulation T1
     IN: f2: coordonnees des pts de la triangulation T2
     OUT: coord1: coordonnees de T1 eventuellement perturbees
     OUT: coord2: coordonnees de T2
  */
  void perturbate2D(E_Int posx1, E_Int posy1, E_Int posz1,
                    K_FLD::FldArrayF& f1,
                    E_Int posx2, E_Int posy2, E_Int posz2,
                    K_FLD::FldArrayF& f2,
                    K_FLD::FldArrayF& coord1, K_FLD::FldArrayF& coord2);

  /* Verifie les tetraedres : pour tous les tetraedres de type t1 ou t2,
     si le tetraedre n a pas un voisin t12, le supprimer.
     IN: type : tableau des types des tetraedres
     IN/OUT: cn : connectivite TETRA/VERTEX
     IN/OUT: coord : points dans le maillage tetraedrique
  */
  void checkTetrahedra(K_FLD::FldArrayI& type, K_FLD::FldArrayI& cn,
                       K_FLD::FldArrayF& coord);

  /* Test si les trois aretes d'un triangle sont interieures
     au contour ou sur le contour.
     IN: ind1, ind2, ind3: indices globaux des sommets du triangle
     IN: coord: coordonnees des sommets
     IN: coordc: coordonnees des points du contour
     IN: cnc: connectivite BAR du contour
     Retourne 0 si l'une des aretes est exterieure
     Retourne 1 si les 3 aretes interieures
  */
  E_Int testIfEdgesCentersInContour(
    E_Int ind1, E_Int ind2,E_Int ind3,
    K_FLD::FldArrayF& coord,
    E_Int posxc1, E_Int posyc1, E_Int poszc1,
    K_FLD::FldArrayF& coordpoly1, K_FLD::FldArrayI& cnc);

  /* Calcule l'angle de courbure pour un "i-array" */
  void computationForiArray(E_Int im, E_Float* xt, E_Float* yt,
                            K_FLD::FldArrayF& angle);

  /* Calcule l'angle de courbure pour un "Bar" */
  void computationForBar(E_Int im, E_Int sizef, E_Float* xt, E_Float* yt,
                         K_FLD::FldArrayI* cn, K_FLD::FldArrayF& angle);

  /* Calcule l'angle de courbure pour un "Tri" */
  void computationForTri(E_Int im, E_Int sizef, E_Float* xt, E_Float* yt,
                         E_Float* zt, K_FLD::FldArrayI* cn,
                         K_FLD::FldArrayF& angle);
  /* Calcul des B-spline avec l'algorithme de De Boor */
  void evalSpline(E_Float t, K_FLD::FldArrayF& x, E_Int n, E_Int c,
                  E_Float* xt, E_Float* yt, E_Float* zt,
                  E_Float& xo, E_Float& yo, E_Float& zo);


  /* Calcul des nurbs avec l'algorithme de De Boor */
  void evalNurbs (E_Float t, K_FLD::FldArrayF& x, E_Int n, E_Int c,
                  K_FLD::FldArrayF& N);

  /* Calcule l'angle de courbure pour un i-array */
  void getCurvatureAngleForiArray(
    E_Int im, E_Float* xt, E_Float* yt, E_Float* zt,
    K_FLD::FldArrayF& angle);
  /* Calcule l'angle de courbure pour un BAR-array */
  void getCurvatureAngleForBar(
    E_Int sizef, E_Float* xt, E_Float* yt, E_Float* zt,
    K_FLD::FldArrayI* cn, K_FLD::FldArrayF& angle);
  /* Calcule l'angle de courbure pour un TRI-array */
  E_Int getCurvatureAngleForTri(
    E_Float* xt, E_Float* yt, E_Float* zt,
    K_FLD:: FldArrayI* cn,
    K_FLD::FldArrayF& angle);

  /* Calcul le factoriel */
  E_Int factorielle(E_Int i);

  /* Evalue le ieme polynome de Bernstein de degre n en t */
  E_Float Bernstein(E_Int i,E_Int j,E_Float t);
}
#endif

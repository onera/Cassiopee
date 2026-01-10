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
# ifndef _POST_ZIPPER_INTERSECTION_H_
# define _POST_ZIPPER_INTERSECTION_H_

# include "StructBlock.h"
# include "TriangleZ.h"
# include<list>
# define FldArrayF K_FLD::FldArrayF

/* teste l'intersection des cellules */
E_Bool testIntersectionOfCells(StructBlock* blk1, 
                                  StructBlock* blk2,
                                  E_Int ind1, E_Int ind2,
                                  E_Float eps);

/* calcul des coordonnees des sommets d'une cellule
 IN: blk : bloc contenant la cellule
           ind : indice du premier sommet de la cellule 
 OUT : coordonnees des 4 sommets de la cellule*/
void compCoordinatesOfCellVertices(StructBlock* blk, E_Int ind, 
                                   FldArrayF& ptA, FldArrayF& ptB, 
                                   FldArrayF& ptC, FldArrayF& ptD);

/* Creation de la liste des 4 triangles de ABCD*/
void createListOfTriangles(StructBlock* blk, E_Int ind, 
                           std::list<TriangleZ*>& triangles);

/* Test si le triangle t est degenere ou non */
E_Bool isDegenerated(TriangleZ* t, E_Float eps);

/* Compute epsilon tolerance from real frame to new frame */
void compEpsilonInNewFrame(E_Float eps, 
                           FldArrayF& mat01, FldArrayF& eps01, 
                           FldArrayF& mat10, FldArrayF& eps10);

/*Calcul des matrices de passage du repere R0(x,y,z) vers R1(X,Y,Z), et 
  de R1 vers R0. Retourne false si le determinant est nul.
  R1 est le repere tel que le triangle t est unitaire dans le plan Z=0*/
E_Bool compMatrix(TriangleZ* t, FldArrayF& mat10, FldArrayF& mat01);

/* teste si le triangle est plan dans la direction Z=0 
   Retourne 1 si oui, -1 si plan parallele,0 si secant */
E_Int isPlanar(E_Float epsZ, 
               FldArrayF& ptA, FldArrayF& ptB, FldArrayF& ptC);

/* Teste si le triangle unite dans le plan (u,v) et le triangle ABC 
   s'intersectent dans le cas ou ils sont coplanaires */
E_Bool 
testIntersectionInPlane(FldArrayF& eps12, FldArrayF& ptA2, 
                        FldArrayF& ptB2, FldArrayF& ptC2,
                        FldArrayF& eps21, FldArrayF& ptA1, 
                        FldArrayF& ptB1, FldArrayF& ptC1);

/* Teste si le segment [AB] de l'espace et les aretes du triangle unite
   de R2 s'intersectent*/
E_Bool 
testIntersectionInSpace(FldArrayF& ptA1, FldArrayF& ptB1, FldArrayF& ptC1,
                        FldArrayF& eps01, FldArrayF& intersect1);

/* Teste si le point ptA est dans le triangle unite */
E_Bool testIfPtInTriangle(FldArrayF& eps, FldArrayF& ptA);

/* Teste si le triangle ABC est dans le triangle unite */
E_Bool testIfTriangleInterior(FldArrayF& ptA, FldArrayF& ptB, 
                                 FldArrayF& ptC, FldArrayF& eps);


/* Calcule les coordonnees des sommets du triangle t1 dans le repere (X,Y,Z)
   dont la matrice de passage est mat.
   IN : t1 : triangle unite dans le plan Z=0
        t2 : triangle pour lequel le changement de repere est a effectuer
        mat : matrice de passage du repere (x,y,z) vers (X,Y,Z)
   OUT :pt1, pt2, pt3 : coordonnees des sommets de t1 dans le nouveau repere*/ 
void compCoordInNewFrame(TriangleZ* t1, TriangleZ* t2, FldArrayF& mat,
                         FldArrayF& pt1, FldArrayF& pt2, FldArrayF& pt3);

/* Changement de repere
   IN : coord : coordonnees de pts dans R0
   IN : mat01 : matrice de changement de repere de R0 vers R1
   IN : x0_1  : deplacement, defini dans R1 
   OUT : newCoord : coordonnees dans le repere R1 */
void compPtsInNewFrame(FldArrayF& mat10, FldArrayF& x0_1, 
                       FldArrayF& coord, FldArrayF& newCoord);


/* Etant donnees 2 combinaisons de points sur la meme droite de l espace
   (2 segts, 1segt+1point, 2 points), teste s'ils s'intersectent */
E_Bool 
testIntersectionOfSegments(FldArrayF& mat10, TriangleZ* T1,
                           FldArrayF& mat20, TriangleZ* T2,
                           FldArrayF& eps10, FldArrayF& eps20,
                           FldArrayF& intersectPts1, 
                           FldArrayF& intersectPts2);


void writeTriangles(char* fileName,
                    FldArrayF& ptA1, FldArrayF& ptB1, FldArrayF& C1,
                    FldArrayF& ptA2, FldArrayF& ptB2, FldArrayF& C2);
# undef FldArrayF

#endif

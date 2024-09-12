/*    
    Copyright 2013-2024 Onera.

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
#include "Intersection.h"
#include <stdlib.h>

using namespace K_FUNC;
using namespace std;
using namespace K_FLD;
using namespace K_CONST;

//=============================================================================
/* Teste l'intersection des cellules */
//=============================================================================
E_Boolean testIntersectionOfCells(StructBlock* blk1, StructBlock* blk2,
                                  E_Int ind1, E_Int ind2,
                                  E_Float eps)
{
  //R0 : repere initial/physique (x,y,z)
  //R1 : repere dans lequel T1 est unitaire dans le plan Z1=0 (X1,Y1,Z1)
  //R2 : repere dans lequel T2 est unitaire dans le plan Z2=0 (X2,Y2,Z2)
  E_Int pl;
  E_Boolean test1, test2;
  FldArrayF mat10(3,3); // matrice de passage du repere 1 au repere 0 
  FldArrayF mat20(3,3); // matrice de passage du repere 2 au repere 0
  FldArrayF mat01(3,3); // matrice de passage du repere 0 au repere 1 
  FldArrayF mat02(3,3); // matrice de passage du repere 0 au repere 2
  FldArrayF ptA1(3);
  FldArrayF ptA2(3);
  FldArrayF ptB1(3);
  FldArrayF ptB2(3);
  FldArrayF ptC1(3);
  FldArrayF ptC2(3);
  FldArrayF eps01(3);
  FldArrayF eps02(3);
  FldArrayF eps10(3);
  FldArrayF eps20(3);
  FldArrayF intersect1;
  FldArrayF intersect2;

  // Creation des triangles de  ABCD et EFGH
  list<TriangleZ*> triangles1;
  list<TriangleZ*>::iterator itrT1;
  list<TriangleZ*>::iterator itrT2;
  list<TriangleZ*> triangles2;
  
  createListOfTriangles(blk1, ind1, triangles1);
  createListOfTriangles(blk2, ind2, triangles2);

  itrT1 = triangles1.begin();
  itrT2 = triangles2.begin();
  
  while( itrT1 != triangles1.end())
  {
    if (compMatrix(*itrT1, mat10, mat01))// matrice de passage
    {
      while (itrT2 != triangles2.end())
      {
        if (compMatrix(*itrT2, mat20, mat02))// matrice de passage
        {
          // change les coordonnees de T2 dans le repere (X1,Y1,Z1)
          compCoordInNewFrame(*itrT1, *itrT2, mat01, ptA1, ptB1, ptC1);
          // change les coordonnees de T1 dans le repere (X2,Y2,Z2)
          compCoordInNewFrame(*itrT2, *itrT1, mat02, ptA2, ptB2, ptC2);
          
          compEpsilonInNewFrame(eps, mat01, eps01, mat10, eps10);
          compEpsilonInNewFrame(eps, mat02, eps02, mat20, eps20);

          // teste si t123 est dans le plan Z=0        
          pl = isPlanar(eps01[2], ptA1, ptB1, ptC1);
          
          switch (pl)
          {
            case 1://plan Z=0
              if(testIntersectionInPlane(eps01, ptA1, ptB1, ptC1, eps02, ptA2, ptB2, ptC2))
              {
                for (auto itrT1 = triangles1.begin(); itrT1 != triangles1.end(); itrT1++)
                  delete *itrT1;
                for (auto itrT2 = triangles2.begin(); itrT2 != triangles2.end(); itrT2++)
                  delete *itrT2;
                triangles1.clear(); triangles2.clear();
                return true;
              }
              break;
            case 0://plan secant a Z=0
              // droite D : intersection des 2 plans contenant T1 et T2
              //calcul des intersections de la droite Davec chaque triangle
              test1 = testIntersectionInSpace(ptA1, ptB1, ptC1, eps01, intersect1);
              test2 = testIntersectionInSpace(ptA2, ptB2, ptC2, eps02, intersect2);

              if (test1 && test2)
              {
                //intersection ?
                if (testIntersectionOfSegments(mat10, *itrT1, mat20, *itrT2, 
                                               eps10, eps20, intersect1, intersect2))
                {
                  for (auto itrT1 = triangles1.begin(); itrT1 != triangles1.end(); itrT1++)
                    delete *itrT1;
                  for (auto itrT2 = triangles2.begin(); itrT2 != triangles2.end(); itrT2++)
                    delete *itrT2;
                  triangles1.clear(); triangles2.clear();
                  return true;
                }
              }
              break;
            default:;// plan parallele a Z=0
          }
        }
        itrT2++;
      }
    }
    itrT1++;
  }
  
  for (auto itrT1 = triangles1.begin(); itrT1 != triangles1.end(); itrT1++)
    delete *itrT1;
  for (auto itrT2 = triangles2.begin(); itrT2 != triangles2.end(); itrT2++)
    delete *itrT2;
  triangles1.clear(); triangles2.clear();
  return false;
}

//=============================================================================
/* calcul des coordonnees des sommets d'une cellule
 IN: blk : bloc contenant la cellule
           ind : indice du premier sommet de la cellule 
 OUT : coordonnees des 4 sommets de la cellule */
//=============================================================================
void compCoordinatesOfCellVertices(StructBlock* blk, E_Int ind, 
                                   FldArrayF& ptA, FldArrayF& ptB, 
                                   FldArrayF& ptC, FldArrayF& ptD)
{
  E_Int im = blk->getIm();
  E_Int jm = blk->getJm();  
  FldArrayF& coord = blk->getCoord();
  E_Int inci = 1;
  E_Int incj = im;
  E_Int indA = ind;
  
  E_Int j = indA/im + 1;
  E_Int i = indA - (j-1)*im + 1;
  if ( i == im)
    inci = -1;
  if ( j == jm)
    incj = -im;

  E_Int indB = indA + inci;
  E_Int indC = indB + incj;
  E_Int indD = indA + incj;
  for (E_Int eq = 1; eq <= 3; eq++)
  {
    ptA[eq-1] = coord(indA,eq);
    ptB[eq-1] = coord(indB,eq);
    ptC[eq-1] = coord(indC,eq);
    ptD[eq-1] = coord(indD,eq);
  }
}

//=============================================================================
/* Creation de la liste des 4 triangles de ABCD*/
//=============================================================================
void createListOfTriangles(StructBlock* blk, E_Int ind,
                           list<TriangleZ*>& triangles)
{
  FldArrayF ptA(3);
  FldArrayF ptB(3);
  FldArrayF ptC(3);
  FldArrayF ptD(3);
  compCoordinatesOfCellVertices(blk, ind, ptA, ptB, ptC, ptD);
  
  // insertion des triangles ABC, ABD, BCD, ACD
  TriangleZ* t1 = new TriangleZ(ptA, ptB, ptC);
  triangles.push_back(t1);

  t1 = new TriangleZ(ptA, ptB, ptD);
  triangles.push_back(t1);
  
  t1 = new TriangleZ(ptB, ptC, ptD);
  triangles.push_back(t1);
  
  t1 = new TriangleZ(ptA, ptC, ptD);
  triangles.push_back(t1);
}

//=============================================================================
/* Test si le triangle t est degenere ou non */
//=============================================================================
E_Boolean isDegenerated(TriangleZ* t, E_Float eps)
{
  FldArrayF& ptA = t->getField1();
  FldArrayF& ptB = t->getField2();
  FldArrayF& ptC = t->getField3();

  E_Float xAB = ptB[0]-ptA[0];
  E_Float yAB = ptB[1]-ptA[1];
  E_Float zAB = ptB[2]-ptA[2];
  E_Float xAC = ptC[0]-ptA[0];
  E_Float yAC = ptC[1]-ptA[1];
  E_Float zAC = ptC[2]-ptA[2];

  // Attention : ce n'est pas du produit scalaire qu il faut mettre !!!
  E_Float normAB = sqrt(xAB*xAB + yAB*yAB + zAB*zAB);
  E_Float normAC = sqrt(xAC*xAC + yAC*yAC + zAC*zAC);
  if (normAB < eps || normAC < eps)
    return false;

  normAB = 1./normAB;
  normAC = 1./normAC;
  xAB = normAB * xAB;
  xAC = normAC * xAC;
  yAB = normAB * yAB;
  yAC = normAC * yAC;
  zAB = normAB * zAB;
  zAC = normAC * zAC;
  E_Float ps = xAB*xAC + yAB*yAC +zAB*zAC;

  if ( E_abs(ps) <= eps )
    return true;
  else return false;
}

//=============================================================================
/* Calcul des matrices de passage du repere R0(x,y,z) vers R1(X,Y,Z), et 
   de R1 vers R0.R1 est le repere tel que le triangle t est unitaire dans 
   le plan dir=0
   retourne False si pas de determinant non nul
   dir = 1 : plan X = 0 
   dir = 2 : plan Y = 0
   dir = 3 : plan Z = 0
   mat10 : matrice de R1 vers R0
   mat01 : matrice de R0 vers R1
   Le point pris comme reference (xA,yA,xA) est le 1er pt du triangle
*/
//=============================================================================
E_Boolean compMatrix(TriangleZ* t, FldArrayF& mat10, FldArrayF& mat01)
{
  E_Int dir;
  E_Float invdet;
  E_Float eps = 1.e-12;
  mat01.setAllValuesAtNull();
  mat10.setAllValuesAtNull();
  FldArrayF& pt1 = t->getField1();
  FldArrayF& pt2 = t->getField2();
  FldArrayF& pt3 = t->getField3();

  E_Float dx12, dy12, dz12, dx13, dy13, dz13;
  E_Float dx120 = pt2[0]-pt1[0];
  E_Float dy120 = pt2[1]-pt1[1];
  E_Float dz120 = pt2[2]-pt1[2];
  E_Float dx130 = pt3[0]-pt1[0];
  E_Float dy130 = pt3[1]-pt1[1];
  E_Float dz130 = pt3[2]-pt1[2];
  E_Float n12 = sqrt(dx120*dx120 + dy120*dy120 + dz120*dz120);
  E_Float n13 = sqrt(dx130*dx130 + dy130*dy130 + dz130*dz130);

  if (n12 < eps || n13 < eps)
  {
    printf("Warning: degenerated triangle not treated.\n");
    return false;
  }
  n12 = 1./n12;
  n13 = 1./n13;

  dx12 = dx120*n12;
  dy12 = dy120*n12;
  dz12 = dz120*n12;
  dx13 = dx130*n13;
  dy13 = dy130*n13;
  dz13 = dz130*n13;
 
  //determinant normalise pour bon preconditionnement
  E_Float det1 = dy12*dz13-dy13*dz12;
  E_Float det2 = dx13*dz12-dx12*dz13; 
  E_Float det3 = dx12*dy13-dx13*dy12;
  
  if ( E_abs(det3) >= eps ) dir = 3;
  else if ( E_abs(det1) >= eps ) dir = 1;
  else if ( E_abs(det2) >= eps ) dir = 2;
  else return false;
  
  E_Float lambda1, lambda2, lambda3;
 
  switch (dir)
  {
    case 1:
      invdet = 1./det1;
      lambda1 = n12 * invdet;
      lambda2 = n13 * invdet;
      lambda3 =  1. * invdet;
      //matrice de passage de R0 vers R1
      mat01(0,2) =  dz13 * lambda1;
      mat01(0,3) = -dy13 * lambda1;
      mat01(1,2) = -dz12 * lambda2;
      mat01(1,3) =  dy12 * lambda2;
      mat01(2,1) =  1.;
      mat01(2,2) =  (dx13*dz12 - dx12*dz13) * lambda3;
      mat01(2,3) =  (dx12*dy13 - dx13*dy12) * lambda3;       
      
      //matrice de passage de R1 vers R0
      mat10(0,1) = dx120;
      mat10(0,2) = dx130;
      mat10(0,3) = 1.;
      mat10(1,1) = dy120;
      mat10(1,2) = dy130;
      mat10(1,3) = 0.;
      mat10(2,1) = dz120;
      mat10(2,2) = dz130;
      mat10(2,3) = 0.;  
      break;
    case 2:
      invdet = 1./det2;
      lambda1 = n12 * invdet;
      lambda2 = n13 * invdet;
      lambda3 =   1.* invdet;
      
      //matrice de passage de R0 vers R1      
      mat01(0,1) = -dz13 * lambda1;
      mat01(0,3) =  dx13 * lambda1;
      mat01(1,1) =  dz12 * lambda2;
      mat01(1,3) = -dx12 * lambda2;
      mat01(2,1) =  (dy12*dz13 - dy13*dz12) * lambda3;
      mat01(2,2) =  1.;
      mat01(2,3) =  (dx12*dy13 - dx13*dy12) * lambda3;  
  
      //matrice de passage de R1 vers R0
      mat10(0,1) = dx120;
      mat10(0,2) = dx130;
      mat10(0,3) = 0.;
      mat10(1,1) = dy120;
      mat10(1,2) = dy130;
      mat10(1,3) = 1.;
      mat10(2,1) = dz120;
      mat10(2,2) = dz130;
      mat10(2,3) = 0.; 
      break;

    case 3:
      invdet = 1./det3;
      lambda1 = n12 * invdet;
      lambda2 = n13 * invdet;
      lambda3 =  1. * invdet;

      //matrice de passage de R0 vers R1      
      mat01(0,1) =  dy13 * lambda1;
      mat01(0,2) = -dx13 * lambda1;
      mat01(1,1) = -dy12 * lambda2;
      mat01(1,2) =  dx12 * lambda2;
      mat01(2,1) = (dy12*dz13-dy13*dz12) * lambda3;
      mat01(2,2) = (dx13*dz12-dx12*dz13) * lambda3;
      mat01(2,3) = 1.;

      //matrice de passage de R1 vers R0
      mat10(0,1) = dx120;
      mat10(0,2) = dx130;
      mat10(0,3) = 0.;
      mat10(1,1) = dy120;
      mat10(1,2) = dy130;
      mat10(1,3) = 0.;
      mat10(2,1) = dz120;
      mat10(2,2) = dz130;
      mat10(2,3) = 1.;
      break;
    default:
      printf("Error: Not a valid value for dir in compMatrix.\n");
      exit(0);
  }
  return true;
}

//=============================================================================
/* Calcule les coordonnees des sommets de t2 dans le repere dont la matrice
   de passage est mat. */
//=============================================================================
void compCoordInNewFrame(TriangleZ* t1, TriangleZ* t2, FldArrayF& mat,
                         FldArrayF& ptA, FldArrayF& ptB, FldArrayF& ptC)
{
  FldArrayF& pt1 = t1->getField1();
  FldArrayF& pt4 = t2->getField1();
  FldArrayF& pt5 = t2->getField2();
  FldArrayF& pt6 = t2->getField3();

  for (E_Int eq = 0; eq < 3; eq++)
  {
    ptA[eq] = 
      mat(eq,1) * (pt4[0]-pt1[0]) + 
      mat(eq,2) * (pt4[1]-pt1[1]) +
      mat(eq,3) * (pt4[2]-pt1[2]);
    ptB[eq] =    
      mat(eq,1) * (pt5[0]-pt1[0]) + 
      mat(eq,2) * (pt5[1]-pt1[1]) +
      mat(eq,3) * (pt5[2]-pt1[2]);
    ptC[eq] =    
      mat(eq,1) * (pt6[0]-pt1[0]) + 
      mat(eq,2) * (pt6[1]-pt1[1]) +
      mat(eq,3) * (pt6[2]-pt1[2]);
  }
}

//=============================================================================
/* Teste si le triangle unite A1B1C1 dans le plan (X,Y) et le triangle 
   A2B2C2 s'intersectent dans le cas ou ils sont coplanaires 
   MODIFIER LES NOMS DE VARIABLES->PAS TRES CLAIRES
*/
//=============================================================================
E_Boolean 
testIntersectionInPlane(FldArrayF& eps12, FldArrayF& ptA2, 
                        FldArrayF& ptB2, FldArrayF& ptC2,
                        FldArrayF& eps21, FldArrayF& ptA1, 
                        FldArrayF& ptB1, FldArrayF& ptC1)
{
  E_Float xmin, xmax, ymin, ymax;
  E_Float k;
  E_Float epsX = eps12[0];
  E_Float epsY = eps12[1];
  E_Float epsXY  = epsX + epsY;

  E_Float xA = ptA2[0];
  E_Float yA = ptA2[1];
  E_Float xB = ptB2[0];
  E_Float yB = ptB2[1];
  E_Float xC = ptC2[0];
  E_Float yC = ptC2[1];
 
//   cout << xA << " "<< yA << " " << ptA2[2] << endl;
//   cout << xB << " "<< yB << " " << ptB2[2] << endl;
//   cout << xC << " "<< yC << " " << ptC2[2] << endl;
  E_Float xint = xA;
  E_Float yint = yA;

  /* 0- test si A2,B2 ou C2 est dans le triangle unite */
  if ( testIfTriangleInterior(ptA2, ptB2, ptC2, eps12) == true)
    return true;

  /* 1- test si un des sommets du triangle unite est interieur a ABC*/
  if ( testIfTriangleInterior(ptA1, ptB1, ptC1, eps21) == true)
    return true;
  
  /* 2-  test intersection des droites avec X=0*/
  //droite (AB)
  xmin = E_min(xB,xA);
  xmax = E_max(xB,xA); 
  ymin = E_min(yB,yA);
  ymax = E_max(yB,yA); 
  if (E_abs(xB-xA) <= epsX) // droite X = xA
  {
    if (E_abs(xB) <= epsX && ymin <=1.+epsY && ymax >= -epsY) // X=0?
      return true;
  }
  else // droite y = a x + b
  {
    yint = yA - (yB-yA)/(xB-xA) * xA;
    if ( yint >=-epsY && yint <= 1.+epsY && 
         xmin <= epsX && xmax >=-epsX )
      return true;
  }
  
  //droite (AC)
  xmin = E_min(xC,xA);
  xmax = E_max(xC,xA);
  ymin = E_min(yC,yA);
  ymax = E_max(yC,yA); 
  if (E_abs(xC-xA) <= epsX ) // droite X = xA
  {
    if (E_abs(xA)<=epsX && ymin <=1.+epsY && ymax >= -epsY) 
      return true;
  }
  else 
  {
    yint= yA - (yC-yA)/(xC-xA) * xA;
    if ( yint >=-epsY && yint <= 1.+epsY && 
         xmin <=epsX && xmax >=-epsX)
      return true;
  }
  
  //droite(BC)
  xmin = E_min(xB,xC);
  xmax = E_max(xB,xC);
  ymin = E_min(yC,yB);
  ymax = E_max(yC,yB);
  if (E_abs(xC-xB) <= epsX) // droite x = xB
  {
    if (E_abs(xB) <= epsX && ymin <= 1.+epsY && ymax >= -epsY) 
      return true;
  }
  else 
  {
    yint= yB - (yC-yB)/(xC-xB) * xB;
    if ( yint >=-epsY && yint <= 1.+epsY && 
         xmin <= epsX && xmax >= -epsX)
      return true;
  }

  /* 3-  test intersection des droites avec Y=0*/
  //droite (AB)
  xmin = E_min(xB,xA);
  xmax = E_max(xB,xA); 
  ymin = E_min(yB,yA);
  ymax = E_max(yB,yA);
  if (E_abs(yB-yA) <= epsY ) // droite y = yA
  {
    if (E_abs(yB) <= epsY && xmin <=1.+epsX && xmax >= -epsX)//Y=0? 
      return true;
  }
  else// droite x = a y + b
  {
    xint= xA - (xB-xA)/(yB-yA) * yA;
    if ( xint >=-epsX && xint <= 1.+epsX && 
         ymin <= epsY && ymax >= -epsY)
     return true;
  }
   
  //droite (AC)
  xmin = E_min(xC,xA);
  xmax = E_max(xC,xA);
  ymin = E_min(yC,yA);
  ymax = E_max(yC,yA);
  if (E_abs(yC-yA) <= epsY  && xmin <=1.+epsX && xmax >=-epsX) // droite y = yA
  {
    if (E_abs(yA) <= epsY) 
      return true;
  }
  else 
  {
    xint= xA - (xC-xA)/(yC-yA) * yA;
    if ( xint >=-epsX && xint <= 1.+epsX && 
         ymin <= epsY && ymax >= -epsY)
      return true;
  }
  
  //droite(BC)
  xmin = E_min(xB,xC);
  xmax = E_max(xB,xC);
  ymin = E_min(yC,yB);
  ymax = E_max(yC,yB);
  if (E_abs(yC-yB) <= epsY && xmin <=1.+epsX && xmax >= -epsX) // droite y = yB
  {
    if (E_abs(yB) <= epsY) 
      return true;
  }
  else 
  {
    xint= xB - (xC-xB)/(yC-yB) * yB;
    if ( xint >=-epsX && xint <= 1.+epsX && 
         ymin <= epsY && ymax >= -epsY)
      return true;
  }

  /* 4-  test intersection des droites avec Y=1-X*/
  //droite (AB)
  ymin = E_min(yA,yB);
  ymax = E_max(yA,yB);
  if (E_abs(xB-xA) < epsX)// X = XB
  {
    yint = 1-xB;
    if (yint >= ymin-epsY && yint <= ymax+epsY &&
        xB >= -epsX && xB <= 1.+epsX)
      return true;
  }
  else//droite ax+b
  {
    k = (yB-yA)/(xB-xA);
    if ( E_abs(k+1.) <= epsY ) // k= -1 
    {
      if(E_abs(xB+yB-1) <= epsXY)// droites confondues ?
      {
        if (ymin <= 1.+epsY && ymax >= -epsY )//segments s intersectent
          return true;
      }
    }
    else 
    {
      xint = (1.+k*xB-yB)/(1.+k);
      if ( yint >= -epsY && yint <= 1.+epsY &&
           yint >= ymin-epsY && yint <= ymax+epsY)
        return true;
    }
  }
 //droite (BC)
  ymin = E_min(yB,yC);
  ymax = E_max(yB,yC);
  if (E_abs(xC-xB) <= epsX)// X = XC
  {
    yint = 1-xC;
    if (yint >= ymin-epsY && yint <= ymax+epsY &&
        xC >= -epsX && xC <= 1.+epsX)
      return true;
  }
  else//droite ax+b
  {
    k = (yC-yB)/(xC-xB);
    if ( E_abs(k+1.) <= epsY ) // k= -1 
    {
      if(E_abs(xC+yC-1) <= epsXY)// droites confondues ?
      {
        if (ymin <= 1.+epsY && ymax >= -epsY )//segments s intersectent
          return true;
      }
    }
    else 
    {
      xint = (1.+k*xC-yC)/(1.+k);
      if ( yint >= -epsY && yint <= 1.+epsY &&
           yint >= ymin-epsY && yint <= ymax+epsY)
        return true;
    }
  }
  //droite (AC)
  ymin = E_min(yA,yC);
  ymax = E_max(yA,yC);
  if (E_abs(xC-xA) <= epsX)// X = XC
  {
    yint = 1-xC;
    if (yint >= ymin-epsY && yint <= ymax+epsY &&
        xC >= -epsX && xC <= 1.+epsX)
      return true;
  }
  else//droite ax+b
  {
    k = (yC-yA)/(xC-xA);
    if ( E_abs(k+1.) <= epsY ) // k= -1 
    {
      if(E_abs(xC+yC-1.) <= epsXY)// droites confondues ?
      {
        if (ymin <= 1.+epsY && ymax >= -epsY )//segments s intersectent
          return true;
      }
    }
    else 
    {
      xint = (1.+k*xC-yC)/(1.+k);
      if ( yint >= -epsY && yint <= 1.+epsY &&
           yint >= ymin-epsY && yint <= ymax+epsY)
        return true;
    }
  }
  return false;
}

//=============================================================================
/* Teste l'intersection du triangle T1 dans son repere R1 
   avec la droite d'intersection des plans contenant T1 et T2:
 IN : coordonnees des pts du triangle T2 dans le repere R1: ptA1,ptB1,ptC1
 IN : eps01 : tolerance geometrique dans R1
OUT : intersect1 : coordonnees des pts d intersection
*/
//=============================================================================
E_Boolean 
testIntersectionInSpace(FldArrayF& ptA1, FldArrayF& ptB1, FldArrayF& ptC1,
                        FldArrayF& eps01, FldArrayF& intersect1)
{
  E_Float a1,b1,d1;
  E_Float xint, yint;
  E_Float xAB = ptB1[0]-ptA1[0];
  E_Float yAB = ptB1[1]-ptA1[1];
  E_Float zAB = ptB1[2]-ptA1[2];
  E_Float xAC = ptC1[0]-ptA1[0];
  E_Float yAC = ptC1[1]-ptA1[1];
  E_Float zAC = ptC1[2]-ptA1[2];

  E_Float epsa1 = E_max(eps01[1], eps01[2]);
  E_Float epsb1 = E_max(eps01[0], eps01[2]);
  E_Float epsd1 = E_max(epsb1, eps01[1]);
  E_Float epsa1b1 = E_max(epsa1,epsb1);
  E_Float epsa1d1 = E_max(epsa1,epsd1);

  FldArrayF tmpPts(9,3);
  E_Int cnt = 0;

  /*-------------*/
  /* Repere R1   */
  /*-------------*/
  // 1-equation de la droite d'intersection D des deux plans dans R1
  // a1 X + b1 Y + c1 Z + d1 = 0, Z = 0
  E_Float det = xAB*yAC - xAC*yAB;
  a1 = zAB*yAC - zAC*yAB;
  b1 = zAC*xAB - zAB*xAC;
  d1 = -ptA1[0] * a1 - ptA1[1] * b1 + det * ptA1[2];

//   cout << "a1 "<< a1 <<" "<< b1 << " " << d1 << endl;
//   cout << "eps :"<< eps01[0] << " "<< eps01[1]<< " " << eps01[2]<<endl;
//   cout << "epsa:"<< epsa1 << " "<< epsb1 << " "<< epsd1 << " "<< epsa1b1<< " "<< epsa1d1<<endl;

  if (E_abs(a1)+E_abs(b1) <= epsa1 + epsb1) 
  {intersect1.malloc(0); return false;}

  // 2-recherche des points d'intersection de D avec cote y=[0,1],x=0
  if (E_abs(b1) <= epsb1)
  {
    yint = 0.;
    if (E_abs(d1) <= epsd1)
    {
      tmpPts(cnt,1) = 0.;
      tmpPts(cnt,2) = 0.;
      tmpPts(cnt,3) = 0.;
      
      tmpPts(cnt+1,1) = 0.;
      tmpPts(cnt+1,2) = 1.;
      tmpPts(cnt+1,3) = 0.;
      cnt = cnt+2;
    }
  }
  else
  {
    yint = - d1/b1;

    if (yint >= -eps01[1] && yint <= 1.+eps01[1])
    {
      tmpPts(cnt,1) = 0.;
      tmpPts(cnt,2) = yint;
      tmpPts(cnt,3) = 0.;
      cnt++;
    }
  }

  // 3-recherche des points d'intersection de D avec cote x=[0,1],y=0
  if (E_abs(a1) <= epsa1)
  {
    if (E_abs(d1) <= epsd1)
    {
      tmpPts(cnt,1) = 0.;
      tmpPts(cnt,2) = 0.;
      tmpPts(cnt,3) = 0.;
      cnt++;
      tmpPts(cnt,1) = 0.;
      tmpPts(cnt,2) = yint;
      tmpPts(cnt,3) = 0.;
      cnt++;
    }
  }
  else 
  {
    xint = - d1/a1;
//    cout << "intersect2 "<< xint << " "<< 0. <<endl;

    if (xint >= -eps01[0] && xint <= 1.+eps01[0])
    {
      tmpPts(cnt,1) = xint;
      tmpPts(cnt,2) = 0.;
      tmpPts(cnt,3) = 0.;
      cnt++;
    }
  }

  // 4-recherche des points d'intersection de D avec cote x+y=1
  if (E_abs(a1-b1) <= epsa1b1)
  {
    if (E_abs(a1+d1) <= epsa1d1)
    {
      tmpPts(cnt,1) = 0.;
      tmpPts(cnt,2) = 1.;
      tmpPts(cnt,3) = 0.;
      cnt++;
      tmpPts(cnt,1) = 1.;
      tmpPts(cnt,2) = 0.;
      tmpPts(cnt,3) = 0.;
      cnt++;
    }
  }
  else 
  {
    xint = (b1+d1)/(b1-a1);
    yint = 1.-xint;

    if (xint>=-eps01[0] && yint>=-eps01[1])
    {
      tmpPts(cnt,1) = xint;
      tmpPts(cnt,2) = yint;
      tmpPts(cnt,3) = 0.;
      cnt++;
    }
  }
 
  // tri de la liste des pts d'intersection : il doit en rester 2
  if (cnt == 0 ) return false;

  FldArrayF work(cnt,3);
  E_Int c = 0;
  E_Float xi, yi, zi;
  for (E_Int i = 0; i < cnt; i++)
  {
    E_Boolean found = false;
    xi = tmpPts(i,1);
    yi = tmpPts(i,2);
    zi = tmpPts(i,3);
    for (E_Int j = 0; j < i; j++)
    {
      E_Float dxij = E_abs(tmpPts(j,1)-xi);
      E_Float dyij = E_abs(tmpPts(j,2)-yi);
      E_Float dzij = E_abs(tmpPts(j,3)-zi);

      if (dxij <=eps01[0] && dyij <=eps01[1] && dzij<=eps01[2])
      {
        found = true;
        break;
      }
    }
    if (found == false)
    {
      work(c,1)=xi; work(c,2)=yi, work(c,3)=zi;
      c++;
    }
  }
  if (c > 0)
  {
    work.reAllocMat(c,3);
    intersect1 = work; // copy
    return true;
  }
  else 
  {
    intersect1.malloc(0);
    return false;
  }
}

//=========================================================================
/* teste si le triangle est plan dans la direction Z=0 
   Retourne 1 si oui, -1 si plan parallele,0 si secant */
//=========================================================================
E_Int isPlanar(E_Float epsZ, 
               FldArrayF& ptA, FldArrayF& ptB, FldArrayF& ptC)
{                 
  E_Float dAB = E_abs(ptB[2]-ptA[2]);
  E_Float dAC = E_abs(ptC[2]-ptA[2]);
  E_Float dBC = E_abs(ptC[2]-ptB[2]);
 
  if ( dAB <= epsZ && dAC <= epsZ && dBC <= epsZ)
  {
    if(E_abs(ptA[2]) <= epsZ)
    {
      ptA[2] = 0.;
      ptB[2] = 0.;
      ptC[2] = 0.;
      return 1;// plan Z=0
    }
    else 
      return -1;//plan parallele a Z=0
  }
  return 0;
}

//=============================================================================
/* Compute epsilon tolerance eps01 from real frame R0 to frame 
   R1 and eps10 from R1 to R0 */
//=============================================================================
void compEpsilonInNewFrame(E_Float eps, 
                           FldArrayF& mat01, FldArrayF& eps01, 
                           FldArrayF& mat10, FldArrayF& eps10) 
{
  E_Float epsX, epsY, epsZ;

  // compute epsilon tolerance from R0 to R1 
  epsX = E_abs(mat01(0,1)) * eps;
  epsX = E_max(epsX, E_abs(mat01(0,2))*eps);
  epsX = E_max(epsX, E_abs(mat01(0,3))*eps);

  epsY = E_abs(mat01(1,1)) * eps;
  epsY = E_max(epsY, E_abs(mat01(1,2))*eps);
  epsY = E_max(epsY, E_abs(mat01(1,3))*eps);

  epsZ = E_abs(mat01(2,1)) * eps;
  epsZ = E_max(epsZ, E_abs(mat01(2,2))*eps);
  epsZ = E_max(epsZ, E_abs(mat01(2,3))*eps);

  eps01[0] = epsX;
  eps01[1] = epsY;
  eps01[2] = epsZ;

  // compute eps to go back to R0 from R1
  eps10[0] = E_abs(mat10(0,1))*epsX + E_abs(mat10(0,2))*epsY + 
    E_abs(mat10(0,3))*epsZ;
  eps10[0] = E_max(eps10[0],eps);

  eps10[1] = E_abs(mat10(1,1))*epsX + E_abs(mat10(1,2))*epsY + 
    E_abs(mat10(2,3))*epsZ; 
  eps10[1] = E_max(eps10[1],eps);
  
  eps10[2] = E_abs(mat10(2,1))*epsX + E_abs(mat10(2,2))*epsY + 
    E_abs(mat10(2,3))*epsZ; 
  eps10[2] = E_max(eps10[2],eps);
  return;
}

//=============================================================================
/* Teste si le point ptA est dans le triangle unite */
//=============================================================================
E_Boolean testIfPtInTriangle(FldArrayF& eps, FldArrayF& ptA)
{
  E_Float epsX = eps[0];
  E_Float epsY = eps[1];
  E_Float epsXY = epsX + epsY;
  E_Float xA = ptA[0];
  E_Float yA = ptA[1];
 
  if ( xA >= -epsX && yA >=-epsY &&  xA+yA <= 1.+epsXY)
    return true;

  return false;
}

//=============================================================================
/* Teste si le triangle ABC est dans le triangle unite */
//=============================================================================
E_Boolean testIfTriangleInterior(FldArrayF& ptA, FldArrayF& ptB, 
                                 FldArrayF& ptC, FldArrayF& eps)
{
  if ( testIfPtInTriangle(eps, ptA) == true)
    return true;

  if ( testIfPtInTriangle(eps, ptB) == true)
    return true;

  if ( testIfPtInTriangle(eps, ptC) == true)
    return true;

  return false;
}

//=============================================================================
/* Changement de repere
   IN : coord : coordonnees de pts dans R0
   IN : mat01 : matrice de changement de repere de R0 vers R1
   IN : x0_1  : deplacement, defini dans R1 
   OUT : newCoord : coordonnees dans le repere R1 */
//=============================================================================
void compPtsInNewFrame(FldArrayF& mat10, FldArrayF& x0_1, 
                       FldArrayF& coord, FldArrayF& newCoord)
{
  E_Int size = coord.getSize();
  newCoord.malloc(size,3);

  for (E_Int ind = 0; ind < size; ind++)
  {
    for(E_Int eq = 0; eq < 3; eq++)
    {
      newCoord(ind,eq+1) = x0_1[eq] + mat10(eq,1) * coord(ind,1) +
        mat10(eq,2) * coord(ind,2) + mat10(eq,3) * coord(ind,3);  
    }
  }
}

//=============================================================================
/* Etant donnees 2 combinaisons de points sur la meme droite de l espace
   (2 segts, 1segt+1point, 2 points), teste s'ils s'intersectent */
//=============================================================================
E_Boolean testIntersectionOfSegments(FldArrayF& mat10, TriangleZ* T1,
                                     FldArrayF& mat20, TriangleZ* T2,
                                     FldArrayF& eps10, 
                                     FldArrayF& eps20,
                                     FldArrayF& intersectPts1, 
                                     FldArrayF& intersectPts2)
{
  FldArrayF coordI1;
  FldArrayF coordI2;

  // Coordonnees des points d intersection dans repere R0
  compPtsInNewFrame(mat10, T1->getField1(), intersectPts1, coordI1);
  
  // Coordonnees des points d intersection dans repere R0
  compPtsInNewFrame(mat20, T2->getField1(), intersectPts2, coordI2);

  E_Int n1 = intersectPts1.getSize();
  E_Int n2 = intersectPts2.getSize();
  
  E_Float xA1, xB1, xA2, xB2;
  E_Float xmin1, xmax1, xmin2, xmax2;

  E_Float epsX = eps10[0] + eps20[0];
//   cout << " n1 :"<< n1 << " " << n2<<endl;

  if (n1 == 2 && n2 == 2)
  {
    xA1 = coordI1(0,1);
    xB1 = coordI1(1,1);
    xA2 = coordI2(0,1);
    xB2 = coordI2(1,1);
    xmin1 = E_min(xA1, xB1);
    xmax1 = E_max(xA1, xB1);
    xmin2 = E_min(xA2, xB2);
    xmax2 = E_max(xA2, xB2);
   
    if (xmax1 >= xmin2-epsX && xmin1 <= xmax2 + epsX)
      return true;
    else return false;
  }
  else if (n1 == 1 && n2 == 2)
  {
   //  cout << "x1 :"<< coordI1(0,1) << " "<<coordI1(0,2)<< " "<<coordI1(0,3)<<endl;
//     cout << "x3:"<< coordI2(0,1) << " "<<coordI2(0,2)<< " "<<coordI2(0,3) <<endl;
//     cout << "x4 :"<< coordI2(1,1) << " "<<coordI2(1,2)<< " "<<coordI2(1,3) <<endl;
    xA1 = coordI1(0,1);
    xA2 = coordI2(0,1);
    xB2 = coordI2(1,1);
    xmin2 = E_min(xA2, xB2);
    xmax2 = E_max(xA2, xB2);
    if (xA1 > xmin2-epsX && xA1 < xmax2+epsX)
      return true;
    else return false;
  }
  else if (n1 == 2 && n2 == 1)
  {
//     cout << "x1 :"<< coordI1(0,1) << " "<<coordI1(0,2)<< " "<<coordI1(0,3)<<endl;
//     cout << "x2:"<< coordI1(1,1) << " "<<coordI1(1,2) <<" "<<coordI1(1,3) <<endl;
//     cout << "x3:"<< coordI2(0,1) << " "<<coordI2(0,2)<< " "<<coordI2(0,3) <<endl;
    
    xA2 = coordI2(0,1);
    xA1 = coordI1(0,1);
    xB1 = coordI1(1,1);
    xmin1 = E_min(xA1, xB1);
    xmax1 = E_max(xA1, xB1);    
    if (xA2 > xmin1-epsX && xA2 < xmax1+epsX)
      return true;
    else return false;
  }
  else if (n1 == 1 && n2 == 1)
  {
  //   cout << "x1 :"<< coordI1(0,1) << " "<<coordI1(0,2)<< " "<<coordI1(0,3)<<endl;
//     cout << "x3:"<< coordI2(0,1) << " "<<coordI2(0,2)<< " "<<coordI2(0,3) <<endl;
    
    xA2 = coordI2(0,1);
    xA1 = coordI1(0,1);
    if (E_abs(xA1-xA2) < epsX) return true;
    else return false;
  }
  else
  {
 //    for (E_Int i = 0; i < n1; i++)
//       cout << "x1 :"<< coordI1(i,1) << " "<<coordI1(i,2)<< " "<<coordI1(i,3)<<endl;
//     for (E_Int i = 0; i < n2; i++)
//       cout << "x3:"<< coordI2(i,1) << " "<<coordI2(i,2)<< " "<<coordI2(i,3) <<endl;
    return false;
  }
  return false;
}

//=============================================================================
/* Ecriture des triangles dans un fichier */
//=============================================================================
void writeTriangles(char* fileName,
                    FldArrayF& ptA, FldArrayF& ptB, FldArrayF& ptC,
                    FldArrayF& ptD, FldArrayF& ptE, FldArrayF& ptF)
{
  FldArrayF f(6,3);
  f(0,1) = ptA[0]; f(0,2) = ptA[1]; f(0,3) = ptA[2];
  f(1,1) = ptB[0]; f(1,2) = ptB[1]; f(1,3) = ptB[2];
  f(2,1) = ptC[0]; f(2,2) = ptC[1]; f(2,3) = ptC[2];
  f(3,1) = ptD[0]; f(3,2) = ptD[1]; f(3,3) = ptD[2];
  f(4,1) = ptE[0]; f(4,2) = ptE[1]; f(4,3) = ptE[2];
  f(5,1) = ptF[0]; f(5,2) = ptF[1]; f(5,3) = ptF[2];
  FldArrayI connect(2,3);
  connect(0,1) = 1; connect(0,2) = 2; connect(0,3) = 3;
  connect(1,1) = 4; connect(1,2) = 5; connect(1,3) = 6;

  //K_IO::GenIO::getInstance()->tpwriteTriangles(fileName, f, connect, 
  //                                             false);
}

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

# include "compGeom.h"
# include <stdio.h>
# include "Nuga/include/Triangle.h"

using namespace std;

//=============================================================================
// Test intersection de 2 triangles par Sam L.
//=============================================================================
E_Int K_COMPGEOM::crossIntersectionOfTriangles(
  E_Float* ptA1, E_Float* ptB1, E_Float* ptC1,
  E_Float* ptA2, E_Float* ptB2, E_Float* ptC2,
  E_Float eps)
{
  E_Bool tol_is_absolute = false;
  E_Float u00, u01; E_Int tx;
  E_Bool overlap;
  E_Float tolu = 1.e-6; // tolerance sur u00 correspondant a AP = u00.A2B2 
  //test triangle A1B1C1 avec edge A2B2
  bool intersect = K_MESH::Triangle::intersect<3>(ptA1, ptB1, ptC1, ptA2, ptB2,
                                                  eps, tol_is_absolute, 
                                                  u00, u01, tx, overlap);
  if (intersect == true) 
  {
    
    if (overlap == false && u00 > tolu && u00 < 1-tolu) 
    {
      return -1;
    }
  }
  //test triangle A1B1C1 avec edge B2C2
  intersect = K_MESH::Triangle::intersect<3>(ptA1, ptB1, ptC1, ptB2, ptC2, 
                                             eps, tol_is_absolute, 
                                             u00, u01, tx, overlap);
  if (intersect == true) 
  {
    if (overlap == false && u00 > tolu && u00 < 1-tolu) 
    {
      return -1;
    }
  }
  //test triangle A1B1C1 avec edge C2A2
  intersect = K_MESH::Triangle::intersect<3>(ptA1, ptB1, ptC1, ptC2, ptA2, 
                                             eps, tol_is_absolute, 
                                             u00, u01, tx, overlap);
  if (intersect == true)
  {
    if (overlap == false && u00 > tolu && u00 < 1-tolu) 
    {
      return -1;
    }
  }

  //test triangle A2B2C2 avec edge A1B1
  intersect = K_MESH::Triangle::intersect<3>(ptA2, ptB2, ptC2, ptA1, ptB1, 
                                             eps, tol_is_absolute, 
                                             u00, u01, tx, overlap);
  if (intersect == true)
  {
    if (overlap == false && u00 > tolu && u00 < 1-tolu) return -1;
  }
  //test triangle A2B2C2 avec edge B1C1
  intersect = K_MESH::Triangle::intersect<3>(ptA2, ptB2, ptC2, ptB1, ptC1, 
                                             eps, tol_is_absolute, 
                                             u00, u01, tx, overlap);
  if (intersect == true)
  {
    if (overlap == false && u00 > tolu && u00 < 1-tolu) return -1;
  }
  //test triangle A2B2C2 avec edge C1A1
  intersect = K_MESH::Triangle::intersect<3>(ptA2, ptB2, ptC2, ptC1, ptA1, 
                                             eps, tol_is_absolute, 
                                             u00, u01, tx, overlap);
  if (intersect == true) 
  {
    if ( overlap == false && u00 > tolu && u00 < 1-tolu) return -1;
  }
  return 0;
}

//=============================================================================
/* Compute triangle-triangle intersection par l'algorithme de Moller
   IN: ptA, ptB, ptC sommets des 2 triangles 
   OUT:  0: pas d'intersection
          1: intersection sur une ligne: les 2 intervalles s'intersectent
         -1: intersection en un point commun
         -2: intersection sur une arete commune
         -3: coplanaires et intersectants */
//=============================================================================
E_Int K_COMPGEOM::trianglesIntersection(
  E_Float* ptA1, E_Float* ptB1, E_Float* ptC1,
  E_Float* ptA2, E_Float* ptB2, E_Float* ptC2,
  E_Float eps)
{
  //E_Float eps2 = eps*eps;
  E_Float lref = K_CONST::E_MAX_FLOAT; //distance de reference
  E_Float dx, dy, dz;
  dx = K_FUNC::E_abs(ptB1[0]-ptA1[0]);
  dx = K_FUNC::E_min(dx, K_FUNC::E_abs(ptC1[0]-ptA1[0]));
  dx = K_FUNC::E_min(dx, K_FUNC::E_abs(ptC1[0]-ptB1[0]));
  dx = K_FUNC::E_min(dx, K_FUNC::E_abs(ptB2[0]-ptA2[0]));
  dx = K_FUNC::E_min(dx, K_FUNC::E_abs(ptC2[0]-ptA2[0]));
  dx = K_FUNC::E_min(dx, K_FUNC::E_abs(ptB2[0]-ptC2[0]));

  dy = K_FUNC::E_abs(ptB1[1]-ptA1[1]);
  dy = K_FUNC::E_min(dy, K_FUNC::E_abs(ptC1[1]-ptA1[1]));
  dy = K_FUNC::E_min(dy, K_FUNC::E_abs(ptC1[1]-ptB1[1]));
  dy = K_FUNC::E_min(dy, K_FUNC::E_abs(ptB2[1]-ptA2[1]));
  dy = K_FUNC::E_min(dy, K_FUNC::E_abs(ptC2[1]-ptA2[1]));
  dy = K_FUNC::E_min(dy, K_FUNC::E_abs(ptB2[1]-ptC2[1]));

  dz = K_FUNC::E_abs(ptB1[2]-ptA1[2]);
  dz = K_FUNC::E_min(dz, K_FUNC::E_abs(ptC1[2]-ptA1[2]));
  dz = K_FUNC::E_min(dz, K_FUNC::E_abs(ptC1[2]-ptB1[2]));
  dz = K_FUNC::E_min(dz, K_FUNC::E_abs(ptB2[2]-ptA2[2]));
  dz = K_FUNC::E_min(dz, K_FUNC::E_abs(ptC2[2]-ptA2[2]));
  dz = K_FUNC::E_min(dz, K_FUNC::E_abs(ptB2[2]-ptC2[2]));
  lref = K_FUNC::E_min(dx,dy); lref = K_FUNC::E_min(lref,dz);

  // dbx
  lref = 1.;
  // fin dbx
  E_Float tol = lref*eps;
  E_Float tol2 = tol*tol;
  //1- calcul de l equation de plan du triangle T1
  E_Float N1[3];
  N1[0] = (ptB1[1]-ptA1[1])*(ptC1[2]-ptA1[2])-(ptC1[1]-ptA1[1])*(ptB1[2]-ptA1[2]);
  N1[1] = (ptB1[2]-ptA1[2])*(ptC1[0]-ptA1[0])-(ptC1[2]-ptA1[2])*(ptB1[0]-ptA1[0]);
  N1[2] = (ptB1[0]-ptA1[0])*(ptC1[1]-ptA1[1])-(ptC1[0]-ptA1[0])*(ptB1[1]-ptA1[1]);
  E_Float d1 = - (N1[0]*ptA1[0] + N1[1]*ptA1[1] + N1[2]*ptA1[2]); 
  E_Float resA2 = N1[0]*ptA2[0] + N1[1]*ptA2[1] + N1[2]*ptA2[2]+d1;
  E_Float resB2 = N1[0]*ptB2[0] + N1[1]*ptB2[1] + N1[2]*ptB2[2]+d1;
  E_Float resC2 = N1[0]*ptC2[0] + N1[1]*ptC2[1] + N1[2]*ptC2[2]+d1;

  if (resA2 > tol2 && resB2 > tol2 && resC2 > tol2) return 0;
  if (resA2 <-tol2 && resB2 <-tol2 && resC2 <-tol2) return 0;
  //3- calcul de l equation de plan de T2
  E_Float N2[3];
  N2[0] = (ptB2[1]-ptA2[1])*(ptC2[2]-ptA2[2])-(ptC2[1]-ptA2[1])*(ptB2[2]-ptA2[2]);
  N2[1] = (ptB2[2]-ptA2[2])*(ptC2[0]-ptA2[0])-(ptC2[2]-ptA2[2])*(ptB2[0]-ptA2[0]);
  N2[2] = (ptB2[0]-ptA2[0])*(ptC2[1]-ptA2[1])-(ptC2[0]-ptA2[0])*(ptB2[1]-ptA2[1]);
  E_Float d2 = -(N2[0] * ptA2[0] + N2[1]*ptA2[1] + N2[2]*ptA2[2]); 

  //4- rejette si tous les pts de T1 sont du meme cote
  E_Float resA1 = N2[0]*ptA1[0] + N2[1]*ptA1[1] + N2[2]*ptA1[2]+d2;
  E_Float resB1 = N2[0]*ptB1[0] + N2[1]*ptB1[1] + N2[2]*ptB1[2]+d2;
  E_Float resC1 = N2[0]*ptC1[0] + N2[1]*ptC1[1] + N2[2]*ptC1[2]+d2;
  if (resA1 > tol2 && resB1 > tol2 && resC1 > tol2) return 0;
  if (resA1 <-tol2 && resB1 <-tol2 && resC1 <-tol2) return 0;
 
  // triangles coplanaires ? 
  E_Int nc1 = 0;
  if (K_FUNC::E_abs(resA1) < tol2) nc1++;
  if (K_FUNC::E_abs(resB1) < tol2) nc1++;
  if (K_FUNC::E_abs(resC1) < tol2) nc1++;
  E_Int nc2 = 0;
  if (K_FUNC::E_abs(resA2)<tol2) nc2++;
  if (K_FUNC::E_abs(resB2)<tol2) nc2++;
  if (K_FUNC::E_abs(resC2)<tol2) nc2++;
  // coplanaires : est-ce que les 2 triangles s intersectent ou non ?
  if (nc1 == 3) return testCoplanarTrianglesIntersection(ptA2, ptB2, ptC2, ptA1, ptB1, ptC1, N2, N1,tol);
  if (nc2 == 3) return testCoplanarTrianglesIntersection(ptA1, ptB1, ptC1, ptA2, ptB2, ptC2, N1, N2,tol);
  // un sommet commun 
  if (nc1 == 1) return testCommonVertexIntersection(ptA1, ptB1, ptC1, ptA2, ptB2, ptC2, 
                                                    resA1, resB1, resC1, resA2, resB2, resC2, tol);
  if (nc2 == 1) return testCommonVertexIntersection(ptA2, ptB2, ptC1, ptA1, ptB1, ptC1, 
                                                      resA2, resB2, resC2, resA1, resB1, resC1, tol); 
  // une arete d un triangle dans le plan de l autre triangle 
  if (nc1 == 2) return testCommonEdgeIntersection(ptA2, ptB2, ptC2, ptA1, ptB1, ptC1, 
                                                    resA2, resB2, resC2, resA1, resB1, resC1, N2, tol);
  if (nc2 == 2) return testCommonEdgeIntersection(ptA1, ptB1, ptC1, ptA2, ptB2, ptC2, 
                                                    resA1, resB1, resC1, resA2, resB2, resC2, N1, tol);

  //5- Non coplanaires : calcule la ligne d intersection + projection
               
  //   calcule les intervalles pour chq triangle
  E_Float D[3]; // vecteur directeur de la ligne 
  D[0] = N1[1]*N2[2]-N1[2]*N2[1];
  D[1] = N1[2]*N2[0]-N1[0]*N2[2];
  D[2] = N1[0]*N2[1]-N1[1]*N2[0];
  E_Int dir = 0; 
  if (K_FUNC::E_abs(D[0])>K_FUNC::E_abs(D[1]) && K_FUNC::E_abs(D[0])>K_FUNC::E_abs(D[2])) dir = 1;
  else if (K_FUNC::E_abs(D[1])>K_FUNC::E_abs(D[0]) && K_FUNC::E_abs(D[1])>K_FUNC::E_abs(D[2])) dir = 2;
  else dir = 3;
  E_Float pA1, pB1, pC1, pA2, pB2, pC2;
  switch (dir) 
  {
    case 1:
      pA1 = D[0]*ptA1[0]; pB1 = D[0]*ptB1[0]; pC1 = D[0]*ptC1[0]; 
      pA2 = D[0]*ptA2[0]; pB2 = D[0]*ptB2[0]; pC2 = D[0]*ptC2[0]; 
      break;
    case 2:
      pA1 = D[1]*ptA1[1]; pB1 = D[1]*ptB1[1]; pC1 = D[1]*ptC1[1]; 
      pA2 = D[1]*ptA2[1]; pB2 = D[1]*ptB2[1]; pC2 = D[1]*ptC2[1]; 
      break;
    case 3:
      pA1 = D[2]*ptA1[2]; pB1 = D[2]*ptB1[2]; pC1 = D[2]*ptC1[2]; 
      pA2 = D[2]*ptA2[2]; pB2 = D[2]*ptB2[2]; pC2 = D[2]*ptC2[2]; 
      break;
  }
  //detection des pts de T1 du meme cote
  E_Float t1T1, t2T1, t1T2, t2T2; 
  if ( K_FUNC::E_sign(resA1) == K_FUNC::E_sign(resB1) ) //C1 de l autre cote
  {t1T1 = pA1 + (pC1-pA1)*resA1/(resA1-resC1); t2T1 = pB1 + (pC1-pB1)*resB1/(resB1-resC1);}
  else if (K_FUNC::E_sign(resA1) == K_FUNC::E_sign(resC1) )// B1 de l autre cote
  {t1T1 = pA1 + (pB1-pA1)*resA1/(resA1-resB1); t2T1 = pC1 + (pB1-pC1)*resC1/(resC1-resB1);}
  else //A1 de l autre cote
  {t1T1 = pB1 + (pA1-pB1)*resB1/(resB1-resA1); t2T1 = pC1 + (pA1-pC1)*resC1/(resC1-resA1);}

  //idem pour T2
  if ( K_FUNC::E_sign(resA2) == K_FUNC::E_sign(resB2) ) //C2 de l autre cote
  {t1T2 = pA2 + (pC2-pA2)*resA2/(resA2-resC2); t2T2 = pB2 + (pC2-pB2)*resB2/(resB2-resC2);}
  else if (K_FUNC::E_sign(resA2) == K_FUNC::E_sign(resC2) )// B2 de l autre cote
  {t1T2 = pA2 + (pB2-pA2)*resA2/(resA2-resB2); t2T2 = pC2 + (pB2-pC2)*resC2/(resC2-resB2);}
  else //A2 de l autre cote
  {t1T2 = pB2 + (pA2-pB2)*resB2/(resB2-resA2); t2T2 = pC2 + (pA2-pC2)*resC2/(resC2-resA2);}

  //6- intersection des triangles sur la ligne
  // intervalle interne au triangle
  if ( K_FUNC::E_min(t1T1,t2T1) > K_FUNC::E_max(t1T2,t2T2) + tol ||
       K_FUNC::E_max(t1T1,t2T1) < K_FUNC::E_min(t1T2,t2T2) - tol ) return 0;
  else return -1;
}

//=============================================================================
E_Int K_COMPGEOM::getTypeOfSegmentIntersection(E_Float* P1T1, E_Float* P2T1,
                                               E_Float* P1T2, E_Float* P2T2,
                                               E_Float eps)
{
  E_Float pi0[3]; E_Float pi1[3];
  E_Int res = intersect2Segments(P1T1, P2T1, P1T2, P2T2, pi0, pi1);
  E_Float eps2 = eps*eps;
  E_Float dist, dist1, dist2;

  if (res == 1) // intersection en 1 pt
  {
    if (inSegment(pi0, P1T1, P2T1) == 1 && inSegment(pi0, P1T2, P2T2) == 1) 
    {
      //sommet ? 
      dist = (P1T1[0]-pi0[0])*(P1T1[0]-pi0[0]) + (P1T1[1]-pi0[1])*(P1T1[1]-pi0[1]) + (P1T1[2]-pi0[2])*(P1T1[2]-pi0[2]);
      if (K_FUNC::fEqualZero(dist, eps2) == true ) return 1;// intersecte sur sommet1
      
      dist = (P2T1[0]-pi0[0])*(P2T1[0]-pi0[0]) + (P2T1[1]-pi0[1])*(P2T1[1]-pi0[1]) + (P2T1[2]-pi0[2])*(P2T1[2]-pi0[2]);
      if (K_FUNC::fEqualZero(dist, eps2) == true ) return 1;// intersecte sur sommet2
      
      dist = (P1T2[0]-pi0[0])*(P1T2[0]-pi0[0]) + (P1T2[1]-pi0[1])*(P1T2[1]-pi0[1]) + (P1T2[2]-pi0[2])*(P1T2[2]-pi0[2]);
      if (K_FUNC::fEqualZero(dist, eps2) == true ) return 1;// intersecte sur sommet1
    
      dist = (P2T2[0]-pi0[0])*(P2T2[0]-pi0[0]) + (P2T2[1]-pi0[1])*(P2T2[1]-pi0[1]) + (P2T2[2]-pi0[2])*(P2T2[2]-pi0[2]);
      if (K_FUNC::fEqualZero(dist, eps2) == true ) return 1;// intersecte sur sommet2
      return -1;// intersection sur un pt interieur
    }  
    return 0;
  }
  else if ( res == 2 ) // intersection sur un segment
  {
    //arete ou segment interne ?
    // test [P1T1,P2T1] - [pi0,pi1]
    dist1 = (P1T1[0]-pi0[0])*(P1T1[0]-pi0[0]) + (P1T1[1]-pi0[1])*(P1T1[1]-pi0[1]) + (P1T1[2]-pi0[2])*(P1T1[2]-pi0[2]);
    dist2 = (P2T1[0]-pi1[0])*(P2T1[0]-pi1[0]) + (P2T1[1]-pi1[1])*(P2T1[1]-pi1[1]) + (P2T1[2]-pi1[2])*(P2T1[2]-pi1[2]); 
    if (K_FUNC::fEqualZero(dist1, eps2) == true && K_FUNC::fEqualZero(dist2, eps2) == true) return 2;
    // test [P1T1,P2T1] - [pi1,pi0] 
    dist1 = (P2T1[0]-pi0[0])*(P2T1[0]-pi0[0]) + (P2T1[1]-pi0[1])*(P2T1[1]-pi0[1]) + (P2T1[2]-pi0[2])*(P2T1[2]-pi0[2]);
    dist2 = (P1T1[0]-pi1[0])*(P1T1[0]-pi1[0]) + (P1T1[1]-pi1[1])*(P1T1[1]-pi1[1]) + (P1T1[2]-pi1[2])*(P1T1[2]-pi1[2]); 
    if (K_FUNC::fEqualZero(dist1, eps2) == true && K_FUNC::fEqualZero(dist2, eps2) == true) return 2;
    // test [P1T2,P2T2] - [pi0,pi1]
    dist1 = (P1T2[0]-pi0[0])*(P1T2[0]-pi0[0]) + (P1T2[1]-pi0[1])*(P1T2[1]-pi0[1]) + (P1T2[2]-pi0[2])*(P1T2[2]-pi0[2]);
    dist2 = (P2T2[0]-pi1[0])*(P2T2[0]-pi1[0]) + (P2T2[1]-pi1[1])*(P2T2[1]-pi1[1]) + (P2T2[2]-pi1[2])*(P2T2[2]-pi1[2]); 
    if (K_FUNC::fEqualZero(dist1, eps2) == true && K_FUNC::fEqualZero(dist2, eps2) == true) return 2;
    // test [P1T2,P2T2] - [pi1,pi0] 
    dist1 = (P2T2[0]-pi0[0])*(P2T2[0]-pi0[0]) + (P2T2[1]-pi0[1])*(P2T2[1]-pi0[1]) + (P2T2[2]-pi0[2])*(P2T2[2]-pi0[2]);
    dist2 = (P1T2[0]-pi1[0])*(P1T2[0]-pi1[0]) + (P1T2[1]-pi1[1])*(P1T2[1]-pi1[1]) + (P1T2[2]-pi1[2])*(P1T2[2]-pi1[2]); 
    if (K_FUNC::fEqualZero(dist1, eps2) == true && K_FUNC::fEqualZero(dist2, eps2) == true) return 2;
    
    //segment interne 
    if ( inSegment(pi0, P1T1, P2T1) == 1 && inSegment(pi1, P1T1, P2T1) == 1 ) 
    {
      if ( inSegment(pi0, P1T2, P2T2) == 1 && inSegment(pi1, P1T2, P2T2) == 1 ) 
        return -2;
    }
    return 0;
  }
  return 0;
}
//=============================================================================
/* Teste si 2 triangles coplanaires s'intersectent 
   Retourne  0 si pas d intersection 
   Retourne -1 si intersection en au moins un point different d un sommet 
   Retourne  1 si intersection en un sommet 
   Retourne  2 si intersection en une arete du triangle
   Retourne 10 si un triangle est contenu dans l autre
*/
//=============================================================================
E_Int K_COMPGEOM::testCoplanarTrianglesIntersection(
  E_Float* ptA1, E_Float* ptB1, E_Float* ptC1,
  E_Float* ptA2, E_Float* ptB2, E_Float* ptC2,
  E_Float* N1, E_Float* N2, E_Float eps)
{
  E_Float eps2 = eps*eps;
  E_Float uA1, uB1, uC1, vA1, vB1, vC1, wA1, wB1, wC1;
  E_Float uA2, uB2, uC2, vA2, vB2, vC2, wA2, wB2, wC2;
  //projection des points sur l axe ou l aire est maximisee
  E_Float nx = K_FUNC::E_abs(N1[0]);
  E_Float ny = K_FUNC::E_abs(N1[1]);
  E_Float nz = K_FUNC::E_abs(N1[2]);

  if ( nx > ny-eps2 && nx >nz-eps2 ) // projete sur plan x
  {
    uA1 = ptA1[1]; vA1 = ptA1[2]; wA1 = ptA1[0]; 
    uB1 = ptB1[1]; vB1 = ptB1[2]; wB1 = ptB1[0]; 
    uC1 = ptC1[1]; vC1 = ptC1[2]; wC1 = ptC1[0];
    uA2 = ptA2[1]; vA2 = ptA2[2]; wA2 = ptA2[0]; 
    uB2 = ptB2[1]; vB2 = ptB2[2]; wB2 = ptB2[0];
    uC2 = ptC2[1]; vC2 = ptC2[2]; wC2 = ptC2[0];
  }
  else if ( ny > nx-eps2 && ny >nz-eps2 )// projete sur plan y
  {
    uA1 = ptA1[2]; vA1 = ptA1[0]; wA1 = ptA1[1]; 
    uB1 = ptB1[2]; vB1 = ptB1[0]; wB1 = ptB1[1]; 
    uC1 = ptC1[2]; vC1 = ptC1[0]; wC1 = ptC1[1];
    uA2 = ptA2[2]; vA2 = ptA2[0]; wA2 = ptA2[1]; 
    uB2 = ptB2[2]; vB2 = ptB2[0]; wB2 = ptB2[1]; 
    uC2 = ptC2[2]; vC2 = ptC2[0]; wC2 = ptC2[1];
  }
  else if ( nz > nx-eps2 && nz >ny-eps2 )// projete sur plan z
  {
    uA1 = ptA1[0]; vA1 = ptA1[1]; wA1 = ptA1[2]; 
    uB1 = ptB1[0]; vB1 = ptB1[1]; wB1 = ptB1[2]; 
    uC1 = ptC1[0]; vC1 = ptC1[1]; wC1 = ptC1[2];
    uA2 = ptA2[0]; vA2 = ptA2[1]; wA2 = ptA2[2];  
    uB2 = ptB2[0]; vB2 = ptB2[1]; wB2 = ptB2[2]; 
    uC2 = ptC2[0]; vC2 = ptC2[1]; wC2 = ptC2[2];
  }
  else return 0;

  /* 1er passage : test des 3 edges de T1 avec les 3 edges de T2
     des que l un intersecte -> sortie */
  E_Float P1T1[3]; E_Float P2T1[3];E_Float P1T2[3]; E_Float P2T2[3]; E_Float P3T2[3];
  E_Int okt[9]; 
  //arete A1B1 avec A2B2
  P1T1[0] = uA1; P1T1[1] = vA1; P1T1[2] = wA1; 
  P2T1[0] = uB1; P2T1[1] = vB1; P2T1[2] = wB1; 
  P1T2[0] = uA2; P1T2[1] = vA2; P1T2[2] = wA2;
  P2T2[0] = uB2; P2T2[1] = vB2; P2T2[2] = wB2; 
  okt[0] = getTypeOfSegmentIntersection(P1T1, P2T1, P1T2, P2T2, eps);

  //arete A1B1 avec A2C2
  P2T2[0] = uC2; P2T2[1] = vC2; P2T2[2] = wC2; 
  okt[1] = getTypeOfSegmentIntersection(P1T1, P2T1, P1T2, P2T2, eps);

  //arete A1B1 avec B2C2
  P1T2[0] = uB2; P1T2[1] = vB2; P1T2[2] = wB2; 
  okt[2] = getTypeOfSegmentIntersection(P1T1, P2T1, P1T2, P2T2, eps);

  //arete A1C1 avec A2B2
  P2T1[0] = uC1; P2T1[1] = vC1; P2T1[2] = wC1; 
  P1T2[0] = uA2; P1T2[1] = vA2; P1T2[2] = wA2; 
  P2T2[0] = uB2; P2T2[1] = vB2; P2T2[2] = wB2;
  okt[3] = getTypeOfSegmentIntersection(P1T1, P2T1, P1T2, P2T2, eps);

  //arete A1C1 avec A2C2
  P2T2[0] = uC2; P2T2[1] = vC2; P2T2[2] = wC2; 
  okt[4] = getTypeOfSegmentIntersection(P1T1, P2T1, P1T2, P2T2, eps);

  //arete A1C1 avec B2C2
  P1T2[0] = uB2; P1T2[1] = vB2; P1T2[2] = wB2; 
  okt[5] = getTypeOfSegmentIntersection(P1T1, P2T1, P1T2, P2T2, eps);

  //arete B1C1 avec A2B2
  P1T1[0] = uB1; P1T1[1] = vB1; P1T1[2] = wB1; 
  P1T2[0] = uA2; P1T2[1] = vA2; P1T2[2] = wA2; 
  P2T2[0] = uB2; P2T2[1] = vB2; P2T2[2] = wB2;
  okt[6] = getTypeOfSegmentIntersection(P1T1, P2T1, P1T2, P2T2, eps);

  //arete A1B1 avec A2C2
  P2T2[0] = uC2; P2T2[1] = vC2; P2T2[2] = wC2; 
  okt[7] = getTypeOfSegmentIntersection(P1T1, P2T1, P1T2, P2T2, eps);

  //arete A1B1 avec B2C2
  P1T2[0] = uB2; P1T2[1] = vB2; P1T2[2] = wB2; 
  okt[8] = getTypeOfSegmentIntersection(P1T1, P2T1, P1T2, P2T2, eps);

  //cas ou un pt interieur est intersection
  for (E_Int i = 0; i < 9; i++) { if (okt[i] == -1) return -1; }
  for (E_Int i = 0; i < 9; i++)
  { if (okt[i] != 0) return okt[i]; }// cas possibles 1, 2, -2

  /* 2e passage : test si T1 interieur a T2 */
  E_Float* Ball = NULL; E_Float* BB = NULL;
  P1T2[0] = uA2; P1T2[1] = vA2; P1T2[2] = wA2; 
  P2T2[0] = uB2; P2T2[1] = vB2; P2T2[2] = wB2; 
  P3T2[0] = uC2; P3T2[1] = vC2; P3T2[2] = wC2; 
  P1T1[0] = uA1; P1T1[1] = vA1; P1T1[2] = wA1;
  if (pointInTriangle2D(P1T2, P2T2, P3T2, P1T1, Ball, BB) == 1) return 1;
                       
  /* 3e passage : test si T2 interieur a T1 */
  P1T2[0] = uA1; P1T2[1] = vA1; P1T2[2] = wA1; 
  P2T2[0] = uB1; P2T2[1] = vB1; P2T2[2] = wB1; 
  P3T2[0] = uC1; P3T2[1] = vC1; P3T2[2] = wC1; 
  P1T1[0] = uA2; P1T1[1] = vA2; P1T1[2] = wA2;
  if (pointInTriangle2D(P1T2, P2T2, P3T2, P1T1, Ball,BB) == 1) return 1; 
  return 0;
}

//=============================================================================
E_Int K_COMPGEOM::testCommonVertexIntersection(
  E_Float* ptA1, E_Float* ptB1, E_Float* ptC1,
  E_Float* ptA2, E_Float* ptB2, E_Float* ptC2,
  E_Float resA1, E_Float resB1, E_Float resC1, 
  E_Float resA2, E_Float resB2, E_Float resC2, 
  E_Float eps)
{
  E_Float eps2 = eps*eps;
  E_Float dist1;
  if (K_FUNC::fEqualZero(resA1,eps) == true) 
  {
    dist1 = (ptA2[0]-ptA1[0])*(ptA2[0]-ptA1[0]) + (ptA2[1]-ptA1[1])*(ptA2[1]-ptA1[1])+ (ptA2[2]-ptA1[2])*(ptA2[2]-ptA1[2]);
    if (K_FUNC::fEqualZero(dist1,eps2) == true) return 1;
    dist1 = (ptB2[0]-ptA1[0])*(ptB2[0]-ptA1[0]) + (ptB2[1]-ptA1[1])*(ptB2[1]-ptA1[1])+ (ptB2[2]-ptA1[2])*(ptB2[2]-ptA1[2]);
    if (K_FUNC::fEqualZero(dist1,eps2) == true) return 1;
    dist1 = (ptC2[0]-ptA1[0])*(ptC2[0]-ptA1[0]) + (ptC2[1]-ptA1[1])*(ptC2[1]-ptA1[1])+ (ptC2[2]-ptA1[2])*(ptC2[2]-ptA1[2]);
    if (K_FUNC::fEqualZero(dist1,eps2) == true) return 1;
  }
  else if (K_FUNC::fEqualZero(resB1,eps) == true) 
  {
    dist1 = (ptA2[0]-ptB1[0])*(ptA2[0]-ptB1[0]) + (ptA2[1]-ptB1[1])*(ptA2[1]-ptB1[1])+ (ptA2[2]-ptB1[2])*(ptA2[2]-ptB1[2]);
    if ( K_FUNC::fEqualZero(dist1,eps2) == true ) return 1;
    dist1 = (ptB2[0]-ptB1[0])*(ptB2[0]-ptB1[0]) + (ptB2[1]-ptB1[1])*(ptB2[1]-ptB1[1])+ (ptB2[2]-ptB1[2])*(ptB2[2]-ptB1[2]);
    if ( K_FUNC::fEqualZero(dist1,eps2) == true ) return 1;
    dist1 = (ptC2[0]-ptB1[0])*(ptC2[0]-ptB1[0]) + (ptC2[1]-ptB1[1])*(ptC2[1]-ptB1[1])+ (ptC2[2]-ptB1[2])*(ptC2[2]-ptB1[2]);
    if ( K_FUNC::fEqualZero(dist1,eps2) == true ) return 1;
  }
  else if (K_FUNC::fEqualZero(resC1,eps) == true) 
  {
    dist1 = (ptA2[0]-ptC1[0])*(ptA2[0]-ptC1[0]) + (ptA2[1]-ptC1[1])*(ptA2[1]-ptC1[1])+ (ptA2[2]-ptC1[2])*(ptA2[2]-ptC1[2]);
    if (K_FUNC::fEqualZero(dist1,eps2) == true) return 1;
    dist1 = (ptB2[0]-ptC1[0])*(ptB2[0]-ptC1[0]) + (ptB2[1]-ptC1[1])*(ptB2[1]-ptC1[1])+ (ptB2[2]-ptC1[2])*(ptB2[2]-ptC1[2]);
    if (K_FUNC::fEqualZero(dist1,eps2) == true) return 1;
    dist1 = (ptC2[0]-ptC1[0])*(ptC2[0]-ptC1[0]) + (ptC2[1]-ptC1[1])*(ptC2[1]-ptC1[1])+ (ptC2[2]-ptC1[2])*(ptC2[2]-ptC1[2]);
    if (K_FUNC::fEqualZero(dist1,eps2) == true) return 1;
  }
  return 0;
}
//=============================================================================
/* 4 possibilites : coincidence des aretes, intersection sur un segment 
   interieur, en un point ou pas d'intersection */
//=============================================================================
E_Int K_COMPGEOM::testCommonEdgeIntersection(
  E_Float* ptA1, E_Float* ptB1, E_Float* ptC1,
  E_Float* ptA2, E_Float* ptB2, E_Float* ptC2,
  E_Float resA1, E_Float resB1, E_Float resC1, 
  E_Float resA2, E_Float resB2, E_Float resC2, 
  E_Float* N1, E_Float eps)
{
  E_Float eps2 = eps*eps;
  E_Int ok[3]; E_Int ok0 = 0;
  //projection des points sur l axe ou l aire est maximisee
  E_Float uA1[3], uB1[3], uC1[3];
  E_Float uA2[3], uB2[3], uC2[3];
  E_Float nx = K_FUNC::E_abs(N1[0]);
  E_Float ny = K_FUNC::E_abs(N1[1]);
  E_Float nz = K_FUNC::E_abs(N1[2]);
  if (nx > ny-eps2 && nx > nz-eps2) // projete sur plan x
  {
    uA1[0] = ptA1[1]; uA1[1] = ptA1[2]; uA1[2] = ptA1[0]; 
    uB1[0] = ptB1[1]; uB1[1] = ptB1[2]; uB1[2] = ptB1[0]; 
    uC1[0] = ptC1[1]; uC1[1] = ptC1[2]; uC1[2] = ptC1[0];
    uA2[0] = ptA2[1]; uA2[1] = ptA2[2]; uA2[2] = ptA2[0]; 
    uB2[0] = ptB2[1]; uB2[1] = ptB2[2]; uB2[2] = ptB2[0];
    uC2[0] = ptC2[1]; uC2[1] = ptC2[2]; uC2[2] = ptC2[0];
  }
  else if (ny > nx-eps2 && ny > nz-eps2) // projete sur plan y
  {
    uA1[0] = ptA1[2]; uA1[1] = ptA1[0]; uA1[2] = ptA1[1]; 
    uB1[0] = ptB1[2]; uB1[1] = ptB1[0]; uB1[2] = ptB1[1]; 
    uC1[0] = ptC1[2]; uC1[1] = ptC1[0]; uC1[2] = ptC1[1];
    uA2[0] = ptA2[2]; uA2[1] = ptA2[0]; uA2[2] = ptA2[1]; 
    uB2[0] = ptB2[2]; uB2[1] = ptB2[0]; uB2[2] = ptB2[1]; 
    uC2[0] = ptC2[2]; uC2[1] = ptC2[0]; uC2[2] = ptC2[1];
  }
  else if (nz > nx-eps2 && nz >ny-eps2) // projete sur plan z
  {
    uA1[0] = ptA1[0]; uA1[1] = ptA1[1]; uA1[2] = ptA1[2]; 
    uB1[0] = ptB1[0]; uB1[1] = ptB1[1]; uB1[2] = ptB1[2]; 
    uC1[0] = ptC1[0]; uC1[1] = ptC1[1]; uC1[2] = ptC1[2];
    uA2[0] = ptA2[0]; uA2[1] = ptA2[1]; uA2[2] = ptA2[2];  
    uB2[0] = ptB2[0]; uB2[1] = ptB2[1]; uB2[2] = ptB2[2]; 
    uC2[0] = ptC2[0]; uC2[1] = ptC2[1]; uC2[2] = ptC2[2];
  }
  //1- arete B2C2 dans le plan de T1
  if (K_FUNC::E_abs(resB2)<eps && K_FUNC::E_abs(resC2) < eps) 
  {
    ok[0] = getTypeOfSegmentIntersection(uB1, uC1, uB2, uC2, eps);
    ok[1] = getTypeOfSegmentIntersection(uA1, uC1, uB2, uC2, eps);
    ok[2] = getTypeOfSegmentIntersection(uA1, uB1, uB2, uC2, eps);
    if (K_FUNC::E_abs(ok[0]) > K_FUNC::E_abs(ok[1])) ok0 = ok[0];
    else ok0 = ok[1];
    if (K_FUNC::E_abs(ok[2]) > K_FUNC::E_abs(ok0)) ok0 = ok[2];
  }
  //2- arete A2C2 dans le plan de T1
  if (K_FUNC::E_abs(resA2)<eps && K_FUNC::E_abs(resC2) < eps) 
  {
    ok[0] = getTypeOfSegmentIntersection(uB1, uC1, uA2, uC2, eps);
    ok[1] = getTypeOfSegmentIntersection(uA1, uC1, uA2, uC2, eps);
    ok[2] = getTypeOfSegmentIntersection(uA1, uB1, uA2, uC2, eps);
    if (K_FUNC::E_abs(ok[0]) > K_FUNC::E_abs(ok[1])) ok0 = ok[0];
    else ok0 = ok[1];
    if (K_FUNC::E_abs(ok[2]) > K_FUNC::E_abs(ok0 )) ok0 = ok[2];
  }
  //3- arete A2B2 dans le plan de T1
  if (K_FUNC::E_abs(resB2)<eps && K_FUNC::E_abs(resA2) < eps) 
  {
    ok[0] = getTypeOfSegmentIntersection(uB1, uC1, uB2, uA2, eps);
    ok[1] = getTypeOfSegmentIntersection(uA1, uC1, uB2, uA2, eps);
    ok[2] = getTypeOfSegmentIntersection(uA1, uB1, uB2, uA2, eps);
    if (K_FUNC::E_abs(ok[0]) > K_FUNC::E_abs(ok[1])) ok0 = ok[0];
    else ok0 = ok[1];
    if (K_FUNC::E_abs(ok[2]) > K_FUNC::E_abs(ok0 )) ok0 = ok[2];
  }
  return ok0;
}

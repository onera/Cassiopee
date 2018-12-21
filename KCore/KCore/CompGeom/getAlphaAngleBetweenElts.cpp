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

#include "compGeom.h"
#include <stdio.h>

using namespace std;

//=============================================================================
/* Return the alpha angle between the normals of two triangles.
   Return an angle in [0,360].
   Pour que l'angle soit correct a 180 deg, il faut que les triangles
   soient numerotes dans le meme sens.
   Return -1000 if one of the triangles is degenerated or if triangles 
   do not share a common edge. */
//=============================================================================
E_Float K_COMPGEOM::getAlphaAngleBetweenTriangles(
  E_Float* ptA1, E_Float* ptB1, E_Float* ptC1,
  E_Float* ptA2, E_Float* ptB2, E_Float* ptC2)
{
  E_Float eps = 1.e-12;
  E_Float tolps = 1.e-5; // tolerance sur le produit scalaire
  // detection des sommets opposes 
  E_Int foundA1=0; E_Int foundB1=0; E_Int foundC1=0;
  E_Int c = 0; 
  E_Float dx = ptA1[0]-ptA2[0]; 
  E_Float dy = ptA1[1]-ptA2[1]; 
  E_Float dz = ptA1[2]-ptA2[2];
  if (K_FUNC::fEqualZero(dx,eps) == true && 
      K_FUNC::fEqualZero(dy,eps) == true && 
      K_FUNC::fEqualZero(dz,eps) == true) 
  {foundA1 = 1; c++; goto B1;}
  dx = ptA1[0]-ptB2[0]; dy = ptA1[1]-ptB2[1]; dz = ptA1[2]-ptB2[2]; 
  if (K_FUNC::fEqualZero(dx,eps) == true && 
      K_FUNC::fEqualZero(dy,eps) == true && 
      K_FUNC::fEqualZero(dz,eps) == true) 
  {foundA1 = 2; c++;  goto B1;}
  dx = ptA1[0]-ptC2[0]; dy = ptA1[1]-ptC2[1]; dz = ptA1[2]-ptC2[2]; 
  if (K_FUNC::fEqualZero(dx,eps) == true && 
      K_FUNC::fEqualZero(dy,eps) == true && 
      K_FUNC::fEqualZero(dz,eps) == true) 
  {foundA1 = 3; c++; goto B1;}
  B1:;
  dx = ptB1[0]-ptA2[0]; dy = ptB1[1]-ptA2[1]; dz = ptB1[2]-ptA2[2]; 
  if (K_FUNC::fEqualZero(dx,eps) == true && 
      K_FUNC::fEqualZero(dy,eps) == true && 
      K_FUNC::fEqualZero(dz,eps) == true) 
  {foundB1 = 1; c++;  goto C1;}
  dx = ptB1[0]-ptB2[0]; dy = ptB1[1]-ptB2[1]; dz = ptB1[2]-ptB2[2]; 
  if (K_FUNC::fEqualZero(dx,eps) == true && 
      K_FUNC::fEqualZero(dy,eps) == true && 
      K_FUNC::fEqualZero(dz,eps) == true) 
  {foundB1 = 2; c++; goto C1;}
  dx = ptB1[0]-ptC2[0]; dy = ptB1[1]-ptC2[1]; dz = ptB1[2]-ptC2[2]; 
  if (K_FUNC::fEqualZero(dx,eps) == true && 
      K_FUNC::fEqualZero(dy,eps) == true && 
      K_FUNC::fEqualZero(dz,eps) == true) 
  {foundB1 = 3; c++; goto C1;}
  C1:;
  dx = ptC1[0]-ptA2[0]; dy = ptC1[1]-ptA2[1]; dz = ptC1[2]-ptA2[2]; 
  if (K_FUNC::fEqualZero(dx,eps) == true && 
      K_FUNC::fEqualZero(dy,eps) == true && 
      K_FUNC::fEqualZero(dz,eps) == true) 
  {foundC1 = 1; c++; goto fin;}
  dx = ptC1[0]-ptB2[0]; dy = ptC1[1]-ptB2[1]; dz = ptC1[2]-ptB2[2]; 
  if (K_FUNC::fEqualZero(dx,eps) == true && 
      K_FUNC::fEqualZero(dy,eps) == true && 
      K_FUNC::fEqualZero(dz,eps) == true) 
  {foundC1 = 2; c++; goto fin;}
  dx = ptC1[0]-ptC2[0]; dy = ptC1[1]-ptC2[1]; dz = ptC1[2]-ptC2[2]; 
  if (K_FUNC::fEqualZero(dx,eps) == true && 
      K_FUNC::fEqualZero(dy,eps) == true && 
      K_FUNC::fEqualZero(dz,eps) == true) 
  {foundC1 = 3; c++; goto fin;}
  fin:;
  if (c < 2) {return -1000.;}
  // Calcul les vecteurs t1 et t2
  E_Float t1[3]; E_Float t2[3];
  t1[0]=0.;t1[1]=0.;t1[2]=0.;
  t2[0]=0.;t2[1]=0.;t2[2]=0.;
  if (foundA1 == 0)
  {
    t1[0] = K_CONST::ONE_HALF*(ptB1[0]+ptC1[0])-ptA1[0];
    t1[1] = K_CONST::ONE_HALF*(ptB1[1]+ptC1[1])-ptA1[1];
    t1[2] = K_CONST::ONE_HALF*(ptB1[2]+ptC1[2])-ptA1[2];
  }
  if (foundB1 == 0) 
  {
    t1[0] = K_CONST::ONE_HALF*(ptC1[0]+ptA1[0])-ptB1[0];
    t1[1] = K_CONST::ONE_HALF*(ptC1[1]+ptA1[1])-ptB1[1];
    t1[2] = K_CONST::ONE_HALF*(ptC1[2]+ptA1[2])-ptB1[2];
  }
  if (foundC1 == 0) 
  {
    t1[0] = K_CONST::ONE_HALF*(ptB1[0]+ptA1[0])-ptC1[0];
    t1[1] = K_CONST::ONE_HALF*(ptB1[1]+ptA1[1])-ptC1[1];
    t1[2] = K_CONST::ONE_HALF*(ptB1[2]+ptA1[2])-ptC1[2];
  }
  if (foundA1 != 1 && foundB1 != 1 && foundC1 != 1)
  { 
    t2[0] = K_CONST::ONE_HALF*(ptB2[0]+ptC2[0])-ptA2[0];
    t2[1] = K_CONST::ONE_HALF*(ptB2[1]+ptC2[1])-ptA2[1];
    t2[2] = K_CONST::ONE_HALF*(ptB2[2]+ptC2[2])-ptA2[2];
  }
  else if (foundA1 != 2 && foundB1 != 2 && foundC1 != 2) 
  { 
    t2[0] = K_CONST::ONE_HALF*(ptA2[0]+ptC2[0])-ptB2[0];
    t2[1] = K_CONST::ONE_HALF*(ptA2[1]+ptC2[1])-ptB2[1];
    t2[2] = K_CONST::ONE_HALF*(ptA2[2]+ptC2[2])-ptB2[2];
  }
  else if (foundA1 != 3 && foundB1 != 3 && foundC1 != 3) 
  { 
    t2[0] = K_CONST::ONE_HALF*(ptB2[0]+ptA2[0])-ptC2[0];
    t2[1] = K_CONST::ONE_HALF*(ptB2[1]+ptA2[1])-ptC2[1];
    t2[2] = K_CONST::ONE_HALF*(ptB2[2]+ptA2[2])-ptC2[2];
  }
  
  // 1er triangle
  E_Float l1x = ptA1[0]-ptB1[0];
  E_Float l1y = ptA1[1]-ptB1[1];
  E_Float l1z = ptA1[2]-ptB1[2];
  E_Float l2x = ptA1[0]-ptC1[0];
  E_Float l2y = ptA1[1]-ptC1[1];
  E_Float l2z = ptA1[2]-ptC1[2];
    
  E_Float surfx = (l1y*l2z-l1z*l2y);
  E_Float surfy = (l1z*l2x-l1x*l2z);
  E_Float surfz = (l1x*l2y-l1y*l2x);
  E_Float surf1 = sqrt(surfx*surfx+surfy*surfy+surfz*surfz);
  E_Float sx1 = surfx; E_Float sy1 = surfy; E_Float sz1 = surfz;

  // 2nd triangle 
  l1x = ptA2[0]-ptB2[0];
  l1y = ptA2[1]-ptB2[1];
  l1z = ptA2[2]-ptB2[2];
  l2x = ptA2[0]-ptC2[0];
  l2y = ptA2[1]-ptC2[1];
  l2z = ptA2[2]-ptC2[2];
  
  surfx = (l1y*l2z-l1z*l2y);
  surfy = (l1z*l2x-l1x*l2z);
  surfz = (l1x*l2y-l1y*l2x);
  E_Float surf2 = sqrt(surfx*surfx+surfy*surfy+surfz*surfz);
  
  E_Float sx2 = surfx;
  E_Float sy2 = surfy;
  E_Float sz2 = surfz;

  if (K_FUNC::fEqualZero(surf1, 1.e-16) == true || 
      K_FUNC::fEqualZero(surf2, 1.e-16) == true)
  {
    return -1000.;
  }
  E_Float alp0 = 180./K_CONST::E_PI;
  E_Float inv1 = 1./surf1; E_Float inv2 = 1./surf2;
  E_Float n1[3]; E_Float n2[3];
  n1[0] = sx1*inv1; n1[1] = sy1*inv1; n1[2] = sz1*inv1;
  n2[0] = sx2*inv2; n2[1] = sy2*inv2; n2[2] = sz2*inv2;
  // angle entre les deux normales
  E_Float ps = K_FUNC::dot<3>(n1, n2);
  if (ps >= 1.-tolps) return 180.;
  if (ps <= -1.+tolps) return 180.;
  E_Float alp1 = acos(ps);
  E_Float alpha1 = alp1*alp0;

  // Angle entre la normale n1 avec un vecteur tangent de T2 
  E_Float normt1 = K_FUNC::normalize<3>(t1);
  if (K_FUNC::fEqualZero(normt1,eps) == true) return -1000.;
  E_Float normt2 = K_FUNC::normalize<3>(t2);
  if (K_FUNC::fEqualZero(normt2,eps) == true) return -1000.;

  E_Float ps1 = K_FUNC::dot<3>(n2, t1);
  E_Float ps2 = K_FUNC::dot<3>(n1, t2);
 
  if (ps > tolps) 
  { 
    if (ps1 >= tolps && ps2 >= tolps) alpha1 = 180.+alpha1;
    else if (ps1 < tolps && ps2 < tolps) alpha1 = 180.-alpha1;
    //else { return -1000.; }
    else alpha1 = 180.-alpha1; // a 180 deg pres
  }
  else
  {
    if (ps1 > tolps && ps2 > tolps) alpha1 = 180.+alpha1;
    else if (ps1 <= tolps && ps2 <= tolps) alpha1 = 180.-alpha1;
    //else { return -1000.; }
    else alpha1 = 180.-alpha1; // a 180 degre pres
  }
  return alpha1;
}

//=============================================================================
/* Return the alpha angle between the normals of two quads.
   Return an angle in [0,360]
   Return -1000 if one of the quads is degenerated or if quads 
   do not share a common edge. */
//=============================================================================
E_Float K_COMPGEOM::getAlphaAngleBetweenQuads(
  E_Float* ptA1, E_Float* ptB1, E_Float* ptC1, E_Float* ptD1, 
  E_Float* ptA2, E_Float* ptB2, E_Float* ptC2, E_Float* ptD2)
{
  E_Float eps = 1.e-12;
  E_Float tolps = 1.e-5; // tolerance sur le produit scalaire

  // detection des sommets opposes 
  E_Int foundA1 = 0; E_Int foundB1 = 0; E_Int foundC1 = 0; E_Int foundD1 = 0;
  E_Int c = 0; 
  E_Float dx = ptA1[0]-ptA2[0]; 
  E_Float dy = ptA1[1]-ptA2[1]; 
  E_Float dz = ptA1[2]-ptA2[2]; 
  if ( K_FUNC::fEqualZero(dx,eps) == true && 
       K_FUNC::fEqualZero(dy,eps) == true && 
       K_FUNC::fEqualZero(dz,eps) == true) 
  {foundA1 = 1; c++; goto B1;}
  dx = ptA1[0]-ptB2[0]; dy = ptA1[1]-ptB2[1]; dz = ptA1[2]-ptB2[2]; 
  if ( K_FUNC::fEqualZero(dx,eps) == true && 
       K_FUNC::fEqualZero(dy,eps) == true && 
       K_FUNC::fEqualZero(dz,eps) == true) 
  {foundA1 = 2; c++;  goto B1;}
  dx = ptA1[0]-ptC2[0]; dy = ptA1[1]-ptC2[1]; dz = ptA1[2]-ptC2[2]; 
  if ( K_FUNC::fEqualZero(dx,eps) == true && 
       K_FUNC::fEqualZero(dy,eps) == true && 
       K_FUNC::fEqualZero(dz,eps) == true) 
  {foundA1 = 3; c++; goto B1;}
  dx = ptA1[0]-ptD2[0]; dy = ptA1[1]-ptD2[1]; dz = ptA1[2]-ptD2[2]; 
  if ( K_FUNC::fEqualZero(dx,eps) == true && 
       K_FUNC::fEqualZero(dy,eps) == true && 
       K_FUNC::fEqualZero(dz,eps) == true) 
  {foundA1 = 4; c++; goto B1;}
  B1:;

  dx = ptB1[0]-ptA2[0]; dy = ptB1[1]-ptA2[1]; dz = ptB1[2]-ptA2[2]; 
  if ( K_FUNC::fEqualZero(dx,eps) == true && 
       K_FUNC::fEqualZero(dy,eps) == true && 
       K_FUNC::fEqualZero(dz,eps) == true) 
  {foundB1 = 1; c++;  goto C1;}
  dx = ptB1[0]-ptB2[0]; dy = ptB1[1]-ptB2[1]; dz = ptB1[2]-ptB2[2]; 
  if ( K_FUNC::fEqualZero(dx,eps) == true && 
       K_FUNC::fEqualZero(dy,eps) == true && 
       K_FUNC::fEqualZero(dz,eps) == true)  
  {foundB1 = 2; c++; goto C1;}
  dx = ptB1[0]-ptC2[0]; dy = ptB1[1]-ptC2[1]; dz = ptB1[2]-ptC2[2]; 
  if ( K_FUNC::fEqualZero(dx,eps) == true && 
       K_FUNC::fEqualZero(dy,eps) == true && 
       K_FUNC::fEqualZero(dz,eps) == true) 
  {foundB1 = 3; c++; goto C1;}
  dx = ptB1[0]-ptD2[0]; dy = ptB1[1]-ptD2[1]; dz = ptB1[2]-ptD2[2]; 
  if ( K_FUNC::fEqualZero(dx,eps) == true && 
       K_FUNC::fEqualZero(dy,eps) == true && 
       K_FUNC::fEqualZero(dz,eps) == true) 
  {foundB1 = 4; c++; goto C1;}
  C1:;

  dx = ptC1[0]-ptA2[0]; dy = ptC1[1]-ptA2[1]; dz = ptC1[2]-ptA2[2]; 
  if ( K_FUNC::fEqualZero(dx,eps) == true && 
       K_FUNC::fEqualZero(dy,eps) == true && 
       K_FUNC::fEqualZero(dz,eps) == true) 
  {foundC1 = 1; c++; goto D1;}
  dx = ptC1[0]-ptB2[0]; dy = ptC1[1]-ptB2[1]; dz = ptC1[2]-ptB2[2]; 
  if ( K_FUNC::fEqualZero(dx,eps) == true && 
       K_FUNC::fEqualZero(dy,eps) == true && 
       K_FUNC::fEqualZero(dz,eps) == true) 
  {foundC1 = 2; c++; goto D1;}
  dx = ptC1[0]-ptC2[0]; dy = ptC1[1]-ptC2[1]; dz = ptC1[2]-ptC2[2]; 
  if ( K_FUNC::fEqualZero(dx,eps) == true && 
       K_FUNC::fEqualZero(dy,eps) == true && 
       K_FUNC::fEqualZero(dz,eps) == true) 
  {foundC1 = 3; c++; goto D1;}
  dx = ptC1[0]-ptD2[0]; dy = ptC1[1]-ptD2[1]; dz = ptC1[2]-ptD2[2]; 
  if ( K_FUNC::fEqualZero(dx,eps) == true && 
       K_FUNC::fEqualZero(dy,eps) == true && 
       K_FUNC::fEqualZero(dz,eps) == true) 
  {foundC1 = 4; c++; goto D1;}
  D1:;

  dx = ptD1[0]-ptA2[0]; dy = ptD1[1]-ptA2[1]; dz = ptD1[2]-ptA2[2]; 
  if ( K_FUNC::fEqualZero(dx,eps) == true && 
       K_FUNC::fEqualZero(dy,eps) == true && 
       K_FUNC::fEqualZero(dz,eps) == true) 
  {foundD1 = 1; c++; goto fin;}
  dx = ptD1[0]-ptB2[0]; dy = ptD1[1]-ptB2[1]; dz = ptD1[2]-ptB2[2]; 
  if ( K_FUNC::fEqualZero(dx,eps) == true && 
       K_FUNC::fEqualZero(dy,eps) == true && 
       K_FUNC::fEqualZero(dz,eps) == true) 
  {foundD1 = 2; c++; goto fin;}
  dx = ptD1[0]-ptC2[0]; dy = ptD1[1]-ptC2[1]; dz = ptD1[2]-ptC2[2]; 
  if ( K_FUNC::fEqualZero(dx,eps) == true && 
       K_FUNC::fEqualZero(dy,eps) == true && 
       K_FUNC::fEqualZero(dz,eps) == true) 
  {foundD1 = 3; c++; goto fin;}
  dx = ptD1[0]-ptD2[0]; dy = ptD1[1]-ptD2[1]; dz = ptD1[2]-ptD2[2]; 
  if ( K_FUNC::fEqualZero(dx,eps) == true && 
       K_FUNC::fEqualZero(dy,eps) == true && 
       K_FUNC::fEqualZero(dz,eps) == true) 
  {foundD1 = 4; c++; goto fin;}
  fin:;
  if (c < 2) {return -1000.;}

  E_Float t1[3]; E_Float t2[3];
  t1[0]=0.;t1[1]=0.;t1[2]=0.;
  t2[0]=0.;t2[1]=0.;t2[2]=0.;
  if (foundA1 == 0 && foundB1 == 0) 
  {
    t1[0] = K_CONST::ONE_HALF*((ptD1[0]+ptC1[0])-(ptA1[0]+ptB1[0]));
    t1[1] = K_CONST::ONE_HALF*((ptD1[1]+ptC1[1])-(ptA1[1]+ptB1[1]));
    t1[2] = K_CONST::ONE_HALF*((ptD1[2]+ptC1[2])-(ptA1[2]+ptB1[2]));
  }
  else if (foundA1 == 0 && foundC1 == 0) 
  {
    t1[0] = K_CONST::ONE_HALF*((ptD1[0]+ptB1[0])-(ptA1[0]+ptC1[0]));
    t1[1] = K_CONST::ONE_HALF*((ptD1[1]+ptB1[1])-(ptA1[1]+ptC1[1]));
    t1[2] = K_CONST::ONE_HALF*((ptD1[2]+ptB1[2])-(ptA1[2]+ptC1[2]));
  }
  else if (foundA1 == 0 && foundD1 == 0) 
  {
    t1[0] = K_CONST::ONE_HALF*((ptC1[0]+ptB1[0])-(ptA1[0]+ptD1[0]));
    t1[1] = K_CONST::ONE_HALF*((ptC1[1]+ptB1[1])-(ptA1[1]+ptD1[1]));
    t1[2] = K_CONST::ONE_HALF*((ptC1[2]+ptB1[2])-(ptA1[2]+ptD1[2]));
  }
  else if (foundB1 == 0 && foundC1 == 0) 
  {
    t1[0] = K_CONST::ONE_HALF*((ptA1[0]+ptD1[0])-(ptB1[0]+ptC1[0]));
    t1[1] = K_CONST::ONE_HALF*((ptA1[1]+ptD1[1])-(ptB1[1]+ptC1[1]));
    t1[2] = K_CONST::ONE_HALF*((ptA1[2]+ptD1[2])-(ptB1[2]+ptC1[2]));
  }
  else if (foundB1 == 0 && foundD1 == 0) 
  {
    t1[0] = K_CONST::ONE_HALF*((ptA1[0]+ptC1[0])-(ptD1[0]+ptB1[0]));
    t1[1] = K_CONST::ONE_HALF*((ptA1[1]+ptC1[1])-(ptD1[1]+ptB1[1]));
    t1[2] = K_CONST::ONE_HALF*((ptA1[2]+ptC1[2])-(ptD1[2]+ptB1[2]));
  }
  else if (foundC1 == 0 && foundD1 == 0) 
  {
    t1[0] = K_CONST::ONE_HALF*((ptA1[0]+ptB1[0])-(ptD1[0]+ptC1[0]));
    t1[1] = K_CONST::ONE_HALF*((ptA1[1]+ptB1[1])-(ptD1[1]+ptC1[1]));
    t1[2] = K_CONST::ONE_HALF*((ptA1[2]+ptB1[2])-(ptD1[2]+ptC1[2]));
  }
  if (foundA1 != 1 && foundB1 != 1 && foundC1 != 1 && foundD1 != 1 &&
      foundA1 != 2 && foundB1 != 2 && foundC1 != 2 && foundD1 != 2) 
  { 
    t2[0] = K_CONST::ONE_HALF*((ptD2[0]+ptC2[0])-(ptA2[0]+ptB2[0]));
    t2[1] = K_CONST::ONE_HALF*((ptD2[1]+ptC2[1])-(ptA2[1]+ptB2[1]));
    t2[2] = K_CONST::ONE_HALF*((ptD2[2]+ptC2[2])-(ptA2[2]+ptB2[2]));
  }
  else if (foundA1 != 3 && foundB1 != 3 && foundC1 != 3 && foundD1 != 3 &&
           foundA1 != 2 && foundB1 != 2 && foundC1 != 2 && foundD1 != 2) 
  { 
    t2[0] = K_CONST::ONE_HALF*((ptA2[0]+ptD2[0])-(ptB2[0]+ptC2[0]));
    t2[1] = K_CONST::ONE_HALF*((ptA2[1]+ptD2[1])-(ptB2[1]+ptC2[1]));
    t2[2] = K_CONST::ONE_HALF*((ptA2[2]+ptD2[2])-(ptB2[2]+ptC2[2]));
  }
  else if (foundA1 != 3 && foundB1 != 3 && foundC1 != 3 && foundD1 != 3 &&
           foundA1 != 4 && foundB1 != 4 && foundC1 != 4 && foundD1 != 4) 
  { 
    t2[0] = K_CONST::ONE_HALF*((ptA2[0]+ptB2[0])-(ptD2[0]+ptC2[0]));
    t2[1] = K_CONST::ONE_HALF*((ptA2[1]+ptB2[1])-(ptD2[1]+ptC2[1]));
    t2[2] = K_CONST::ONE_HALF*((ptA2[2]+ptB2[2])-(ptD2[2]+ptC2[2]));
  }
  else if (foundA1 != 1 && foundB1 != 1 && foundC1 != 1 && foundD1 != 1 &&
           foundA1 != 4 && foundB1 != 4 && foundC1 != 4 && foundD1 != 4) 
  { 
    t2[0] = K_CONST::ONE_HALF*((ptC2[0]+ptB2[0])-(ptD2[0]+ptA2[0]));
    t2[1] = K_CONST::ONE_HALF*((ptC2[1]+ptB2[1])-(ptD2[1]+ptA2[1]));
    t2[2] = K_CONST::ONE_HALF*((ptC2[2]+ptB2[2])-(ptD2[2]+ptA2[2]));
  }
  else if (foundA1 != 3 && foundB1 != 3 && foundC1 != 3 && foundD1 != 3 &&
           foundA1 != 1 && foundB1 != 1 && foundC1 != 1 && foundD1 != 1) 
  { 
    t2[0] = K_CONST::ONE_HALF*((ptD2[0]+ptB2[0])-(ptA2[0]+ptC2[0]));
    t2[1] = K_CONST::ONE_HALF*((ptD2[1]+ptB2[1])-(ptA2[1]+ptC2[1]));
    t2[2] = K_CONST::ONE_HALF*((ptD2[2]+ptB2[2])-(ptA2[2]+ptC2[2]));
  }
  else if (foundA1 != 2 && foundB1 != 2 && foundC1 != 2 && foundD1 != 2 &&
           foundA1 != 4 && foundB1 != 4 && foundC1 != 4 && foundD1 != 4) 
  { 
    t2[0] = K_CONST::ONE_HALF*((ptC2[0]+ptA2[0])-(ptD2[0]+ptB2[0]));
    t2[1] = K_CONST::ONE_HALF*((ptC2[1]+ptA2[1])-(ptD2[1]+ptB2[1]));
    t2[2] = K_CONST::ONE_HALF*((ptC2[2]+ptA2[2])-(ptD2[2]+ptB2[2]));
  }
  //1er quad
  // AB x AC
  E_Float l1x = ptA1[0]-ptB1[0];
  E_Float l1y = ptA1[1]-ptB1[1];
  E_Float l1z = ptA1[2]-ptB1[2];
  E_Float l2x = ptA1[0]-ptC1[0];
  E_Float l2y = ptA1[1]-ptC1[1];
  E_Float l2z = ptA1[2]-ptC1[2];  
  E_Float surf1x = (l1y*l2z-l1z*l2y);
  E_Float surf1y = (l1z*l2x-l1x*l2z);
  E_Float surf1z = (l1x*l2y-l1y*l2x);
  // AC x AD
  l1x = ptC1[0]-ptA1[0];
  l1y = ptC1[0]-ptA1[0];
  l1z = ptC1[0]-ptA1[0];
  l2x = ptD1[0]-ptA1[0];
  l2y = ptD1[0]-ptA1[0];
  l2z = ptD1[0]-ptA1[0];         
  E_Float surf2x = (l1y*l2z-l1z*l2y);
  E_Float surf2y = (l1z*l2x-l1x*l2z);
  E_Float surf2z = (l1x*l2y-l1y*l2x);
   
  E_Float sx1 = surf1x + surf2x;
  E_Float sy1 = surf1y + surf2y;
  E_Float sz1 = surf1z + surf2z;
  E_Float surf1 = sqrt(sx1*sx1+sy1*sy1+sz1*sz1);

  //2nd quad
  // AB x AC
  l1x = ptA2[0]-ptB2[0];
  l1y = ptA2[1]-ptB2[1];
  l1z = ptA2[2]-ptB2[2];
  l2x = ptA2[0]-ptC2[0];
  l2y = ptA2[1]-ptC2[1];
  l2z = ptA2[2]-ptC2[2];  
  surf1x = (l1y*l2z-l1z*l2y);
  surf1y = (l1z*l2x-l1x*l2z);
  surf1z = (l1x*l2y-l1y*l2x);
  // AC x AD
  l1x = ptC2[0]-ptA2[0];
  l1y = ptC2[0]-ptA2[0];
  l1z = ptC2[0]-ptA2[0];
  l2x = ptD2[0]-ptA2[0];
  l2y = ptD2[0]-ptA2[0];
  l2z = ptD2[0]-ptA2[0];         
  surf2x = (l1y*l2z-l1z*l2y);
  surf2y = (l1z*l2x-l1x*l2z);
  surf2z = (l1x*l2y-l1y*l2x);
   
  E_Float sx2 = surf1x + surf2x;
  E_Float sy2 = surf1y + surf2y;
  E_Float sz2 = surf1z + surf2z;
  E_Float surf2 = sqrt(sx2*sx2+sy2*sy2+sz2*sz2);

  if ( K_FUNC::fEqualZero(surf1, 1.e-16) == true || 
       K_FUNC::fEqualZero(surf2, 1.e-16) == true )
  {
    return -1000.;
  }
  E_Float alp0 = 180./K_CONST::E_PI;
  E_Float inv1 = 1./surf1; E_Float inv2 = 1./surf2;
  E_Float n1[3]; E_Float n2[3];
  n1[0] = sx1*inv1; n1[1] = sy1*inv1; n1[2] = sz1*inv1;
  n2[0] = sx2*inv2; n2[1] = sy2*inv2; n2[2] = sz2*inv2;
  //angle entre les deux normales
  E_Float ps = K_FUNC::dot<3>(n1, n2);
  if ( ps >= 1.-tolps) return 180.; 
  if ( ps <=-1.+tolps) return 360.;
  E_Float alp1 = acos(ps);
  E_Float alpha1 = alp1*alp0;

  //Angle entre la normale n1 avec un vecteur tangent de T2 avec 
  E_Float normt1 = K_FUNC::normalize<3>(t1);
  if ( K_FUNC::fEqualZero(normt1,eps) == true ) return -1000.;
  E_Float normt2 = K_FUNC::normalize<3>(t2);
  if ( K_FUNC::fEqualZero(normt2,eps) == true ) return -1000.;
  E_Float ps1 = K_FUNC::dot<3>(n2, t1);
  E_Float ps2 = K_FUNC::dot<3>(n1, t2);
  if (ps > tolps)
  { 
    if (ps1 >= tolps && ps2 >= tolps) alpha1 = 180.+alpha1;
    else if (ps1 < tolps && ps2 < tolps) alpha1 = 180.-alpha1;
    else {return -1000.;}
  }
  else if (ps < -tolps) //ps des normales < 0 
  {
    if (ps1 > tolps && ps2  > tolps) alpha1 = 180.+alpha1;
    else if ( ps1 <= tolps && ps2 <= tolps ) alpha1 = 180.-alpha1; 
    else {return -1000.;}
  }
  else // ps = 0 
  {
    if (ps1 > tolps && ps2  > tolps) alpha1 = 270.;
    else if (ps1 <= tolps && ps2 <= tolps) alpha1 = 90.; 
    else {printf("Warning: getAlphaAngleBetweenQuads: quads are not oriented in the same direction.\n"); return -1000.;}
  }
  return alpha1;
}
//=============================================================================
/* Return the alpha angle between two vectors ptA1ptB1 and ptA2ptB2
   dirVect must be (approximatively) the direction vector orthogonal to the 
   plane defined by the vectors
   Return -1000. if vectors do not share a common point. */
//=============================================================================
E_Float K_COMPGEOM::getAlphaAngleBetweenBars(E_Float* ptA1, E_Float* ptB1, 
                                             E_Float* ptA2, E_Float* ptB2,
                                             E_Float* dirVect)
{
  E_Float eps = 1.e-12;

  E_Float alp0 = 180./K_CONST::E_PI;
  E_Float t1[3]; E_Float t2[3];
  // detection du pt commun
  E_Float dx, dy, dz;
  
  dx = ptA1[0]-ptA2[0]; dy = ptA1[1]-ptA2[1]; dz = ptA1[2]-ptA2[2];
  if (K_FUNC::fEqualZero(dx,eps) == true && 
      K_FUNC::fEqualZero(dy,eps) == true && 
      K_FUNC::fEqualZero(dz,eps) == true)
  {
    t1[0] = ptB1[0]-ptA1[0]; t1[1] = ptB1[1]-ptA1[1]; t1[2] = ptB1[2]-ptA1[2];
    t2[0] = ptB2[0]-ptA2[0]; t2[1] = ptB2[1]-ptA2[1]; t2[2] = ptB2[2]-ptA2[2];
    goto fin;
  }
  dx = ptA1[0]-ptB2[0]; dy = ptA1[1]-ptB2[1]; dz = ptA1[2]-ptB2[2];
  if (K_FUNC::fEqualZero(dx,eps) == true && 
      K_FUNC::fEqualZero(dy,eps) == true && 
      K_FUNC::fEqualZero(dz,eps) == true)
  {
    t1[0] = ptB1[0]-ptA1[0]; t1[1] = ptB1[1]-ptA1[1]; t1[2] = ptB1[2]-ptA1[2];
    t2[0] = ptA2[0]-ptB2[0]; t2[1] = ptA2[1]-ptB2[1]; t2[2] = ptA2[2]-ptB2[2];
    goto fin;
  }  
  dx = ptB1[0]-ptA2[0]; dy = ptB1[1]-ptA2[1]; dz = ptB1[2]-ptA2[2];
  if (K_FUNC::fEqualZero(dx,eps) == true && 
      K_FUNC::fEqualZero(dy,eps) == true && 
      K_FUNC::fEqualZero(dz,eps) == true)
  {
    t1[0] = ptA1[0]-ptB1[0]; t1[1] = ptA1[1]-ptB1[1]; t1[2] = ptA1[2]-ptB1[2];
    t2[0] = ptB2[0]-ptA2[0]; t2[1] = ptB2[1]-ptA2[1]; t2[2] = ptB2[2]-ptA2[2];
    goto fin;
  }
  dx = ptB1[0]-ptB2[0]; dy = ptB1[1]-ptB2[1]; dz = ptB1[2]-ptB2[2];
  if (K_FUNC::fEqualZero(dx,eps) == true && 
      K_FUNC::fEqualZero(dy,eps) == true && 
      K_FUNC::fEqualZero(dz,eps) == true)
  {
    t1[0] = ptA1[0]-ptB1[0]; t1[1] = ptA1[1]-ptB1[1]; t1[2] = ptA1[2]-ptB1[2];
    t2[0] = ptA2[0]-ptB2[0]; t2[1] = ptA2[1]-ptB2[1]; t2[2] = ptA2[2]-ptB2[2];
    goto fin;
  }
  return -1000.;

  fin:;

  // Compute angle
  E_Float dx1 = t1[0]; E_Float dx2 = t2[0];
  E_Float dy1 = t1[1]; E_Float dy2 = t2[1];
  E_Float dz1 = t1[2]; E_Float dz2 = t2[2];
  E_Float n1 = dx1*dx1+dy1*dy1+dz1*dz1;
  E_Float n2 = dx2*dx2+dy2*dy2+dz2*dz2;
  n1 = K_FUNC::E_max(n1, 1.e-12);
  n2 = K_FUNC::E_max(n2, 1.e-12);

  E_Float inv = 1./sqrt(n1*n2);
  E_Float pv1 = dy1*dz2 - dz1*dy2;
  E_Float pv2 = dz1*dx2 - dx1*dz2;
  E_Float pv3 = dx1*dy2 - dy1*dx2;
  E_Float pv = pv1*dirVect[0] + pv2*dirVect[1] + pv3*dirVect[2]; 
  pv = pv * inv;
  E_Float ps = dx2*dx1 + dy2*dy1 + dz1*dz2;
  ps = ps * inv;
  
  // avec un meilleur conditionnement
  /*
  E_Float inv = 1./sqrt(n1);
  dx1 *= inv; dy1 *= inv; dz1 *= inv;
  inv = 1./sqrt(n2);
  dx2 *= inv; dy2 *= inv; dz2 *= inv;
  E_Float pv1 = dy1*dz2 - dz1*dy2;
  E_Float pv2 = dz1*dx2 - dx1*dz2;
  E_Float pv3 = dx1*dy2 - dy1*dx2;
  E_Float pv = pv1*dirVect[0] + pv2*dirVect[1] + pv3*dirVect[2]; 
  E_Float ps = dx2*dx1 + dy2*dy1 + dz1*dz2;
  */

  pv = K_FUNC::E_min(pv, 1.);
  pv = K_FUNC::E_max(pv,-1.);
  ps = K_FUNC::E_min(ps, 1.);
  ps = K_FUNC::E_max(ps,-1.);
  
  if (pv >= 0. && ps >= 0.) return asin(pv)*alp0; 
  else if (pv <= 0. && ps >= 0.) return 360. + asin(pv)*alp0;
  else if (pv >= 0. && ps <= 0.) return 180. - asin(pv)*alp0;
  else return 180. - asin(pv)*alp0;
}
//=============================================================================
/*  Return the alpha angle between 2 polygons p1 and p2.
    Return an angle in [0,360].
    Return -1000 if one of the polygons is degenerated or if polygons 
    do not share a common edge. 
    Polygons must be indexed in a rotating manner. */
//=============================================================================
E_Float K_COMPGEOM::getAlphaAngleBetweenPolygons(vector<E_Float*>& p1, 
                                                 vector<E_Float*>& p2)
{
  E_Int nv1 = p1.size();
  E_Int nv11 = nv1-1;
  E_Int nv2 = p2.size();
  E_Int nv21 = nv2-1;
  E_Float* pt1; E_Float* pt2;
  E_Float dx, dy, dz;
  E_Float eps = 1.e-12;
  E_Int i, j;

  // identifie un edge commun entre p1 et p2
  for (i = 0; i < nv1; i++)
  { 
    pt1 = p1[i];
    for (j = 0; j < nv2; j++)
    {
      pt2 = p2[j];
      dx = pt2[0]-pt1[0]; dy = pt2[1]-pt1[1]; dz = pt2[2]-pt1[2];
      if (dx*dx+dy*dy+dz*dz < eps)
      {
        // i-j matches
        goto next;
      }
    }
  }
  return -1000.;
  next: ;

#define CHECK pt1 = p1[i2]; pt2 = p2[j2];                       \
  dx = pt2[0]-pt1[0]; dy = pt2[1]-pt1[1]; dz = pt2[2]-pt1[2];   \
  if (dx*dx+dy*dy+dz*dz < eps) goto next2;
    
  E_Int i2, j2;
  // i+1, j+1
  if (i < nv11 && j < nv21)
  {
    i2 = i+1; j2 = j+1;
    CHECK;
  }
  else if (i == nv11 && j < nv21)
  {
    i2 = 0; j2 = j+1;
    CHECK;
  }
  else if (i < nv11 && j == nv21)
  {
    i2 = i+1; j2 = 0;
    CHECK;
  }
  else if (i == nv11 && j == nv21)
  {
    i2 = 0; j2 = 0;
    CHECK;
  }
  // i+1, j-1
  if (i < nv11 && j > 0)
  {
    i2 = i+1; j2 = j-1; 
    CHECK;
  }
  else if (i == nv11 && j > 0)
  {
    i2 = 0; j2 = j-1; 
    CHECK;
  }
  else if (i < nv11 && j == 0)
  {
    i2 = i+1; j2 = nv21; 
    CHECK;
  }
  else if (i == nv11 && j == 0)
  {
    i2 = 0; j2 = nv21; 
    CHECK;
  }
  
  // i-1, j-1
  if (i > 0 && j > 0)
  {
    i2 = i-1; j2 = j-1;
    CHECK;
  }
  else if (i == 0 && j > 0)
  {
    i2 = nv11; j2 = j-1;
    CHECK;
  }
  else if (i > 0 && j == 0)
  {
    i2 = i-1; j2 = nv21;
    CHECK;
  }
  else if (i == 0 && j == 0)
  {
    i2 = nv11; j2 = nv21;
    CHECK;
  }

  // i-1, j+1
  if (i > 0 && j < nv21)
  {
    i2 = i-1; j2 = j+1;
    CHECK;
  }
  else if (i == 0 && j < nv21)
  {
    i2 = nv11; j2 = j+1;
    CHECK;
  }
  else if (i > 0 && j == nv21)
  {
    i2 = i-1; j2 = 0;
    CHECK;
  }
  else if (i == 0 && j == nv21)
  {
    i2 = nv11; j2 = 0;
    CHECK;
  }
  
  return -1000.;
  next2: ;

  // les triangles doivent etre numerotes dans le meme sens
  E_Float* A1; E_Float* B1; E_Float* C1;
  E_Float* A2; E_Float* B2; E_Float* C2;
  A1 = p1[i]; B1 = p1[i2];
  A2 = p2[j2]; B2 = p2[j];

  E_Int i2p1 = i2+1;
  if (i2p1 == nv1) i2p1 = 0;
  E_Int i2m1 = i2-1;
  if (i2m1 == -1) i2m1 = nv11;

  if (i2p1 != i) C1 = p1[i2p1];
  else C1 = p1[i2m1];
 
  E_Int j2p1 = j+1;
  if (j2p1 == nv2) j2p1 = 0;
  E_Int j2m1 = j-1;
  if (j2m1 == -1) j2m1 = nv21;

  if (j2p1 != j2) C2 = p2[j2p1];
  else C2 = p2[j2m1];

  /*
  printf("A1 %f %f %f\n", A1[0], A1[1], A1[2]);
  printf("B1 %f %f %f\n", B1[0], B1[1], B1[2]);
  printf("C1 %f %f %f\n", C1[0], C1[1], C1[2]);
  printf("A2 %f %f %f\n", A2[0], A2[1], A2[2]);
  printf("B2 %f %f %f\n", B2[0], B2[1], B2[2]);
  printf("C2 %f %f %f\n", C2[0], C2[1], C2[2]);
  */

  return getAlphaAngleBetweenTriangles(A1, B1, C1, A2, B2, C2);
}

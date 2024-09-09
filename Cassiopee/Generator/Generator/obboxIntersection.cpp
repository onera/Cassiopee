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

# include "generator.h"
# include <iostream>     // for printing //cout
# include <math.h>       // for the maths

#define norm(x0,y0,z0,x1,y1,z1,result)\
 result = sqrt(pow(x1-x0,2)+pow(y1-y0,2)+pow(z1-z0,2));

using namespace K_FLD;
using namespace std;

//===========================================================================
/* Determine if 2 bounding boxes of 2 arrays intersects
   (array, copy) */
//===========================================================================
PyObject* K_GENERATOR::obboxIntersection(PyObject* self, PyObject* args)
{
  PyObject* array1; PyObject* array2;
  if (!PYPARSETUPLE_(args, OO_, &array1, &array2)) return NULL;

  // Check array1
  E_Int ni1, nj1, nk1;
  FldArrayF* f1; FldArrayI* cn1;
  char* varString1; char* eltType1;
  E_Int res1 = K_ARRAY::getFromArray(
    array1, varString1, f1, ni1, nj1, nk1, cn1, eltType1, true);
    
  if (res1 != 1 && res1 != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "obboxIntersection: array1 is invalid.");
    return NULL;
  }
  // Check array2
  E_Int ni2, nj2, nk2;
  FldArrayF* f2; FldArrayI* cn2;
  char* varString2; char* eltType2;
  E_Int res2 = K_ARRAY::getFromArray(
    array2, varString2, f2, ni2, nj2, nk2, cn2, eltType2, true);
  
  if (res2 != 1 && res2 != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "obboxIntersection: array2 is invalid.");
    return NULL;
    RELEASESHAREDB(res1, array1, f1, cn1); 
  }

  
  // determination de x,y et z ds les arrays
  E_Int posx1 = K_ARRAY::isCoordinateXPresent(varString1);
  E_Int posy1 = K_ARRAY::isCoordinateYPresent(varString1);
  E_Int posz1 = K_ARRAY::isCoordinateZPresent(varString1);
  E_Int posx2 = K_ARRAY::isCoordinateXPresent(varString2);
  E_Int posy2 = K_ARRAY::isCoordinateYPresent(varString2);
  E_Int posz2 = K_ARRAY::isCoordinateZPresent(varString2);
  if (posx1 == -1 || posy1 == -1 || posz1 == -1)
  {
    RELEASESHAREDB(res1, array1, f1, cn1); 
    RELEASESHAREDB(res2, array2, f2, cn2); 
    PyErr_SetString(PyExc_TypeError,
                    "obboxIntersection: can't find coordinates in array1.");
    return NULL;
  }
  if (posx2 == -1 || posy2 == -1 || posz2 == -1)
  {
    RELEASESHAREDB(res1, array1, f1, cn1); 
    RELEASESHAREDB(res2, array2, f2, cn2); 
    PyErr_SetString(PyExc_TypeError,
                    "obboxIntersection: can't find coordinates in array2.");
    return NULL;
  }
  posx1++; posy1++; posz1++; 
  posx2++; posy2++; posz2++;
  
  /* 
  Operations for OBB(A) 
  */
  E_Float* xt1 = f1->begin(1);
  E_Float* yt1 = f1->begin(2);
  E_Float* zt1 = f1->begin(3);
  FldArrayF a(3);

  // Calculates the OBB(A) radii
  norm(xt1[0],yt1[0],zt1[0],xt1[1],yt1[1],zt1[1],a[0]);
  norm(xt1[0],yt1[0],zt1[0],xt1[2],yt1[2],zt1[2],a[1]);
  norm(xt1[0],yt1[0],zt1[0],xt1[4],yt1[4],zt1[4],a[2]);
  for (E_Int i=0; i < 3; i++) {a[i] /= 2.;}

  FldArrayF A(3,3);
  A(0,1) = (xt1[1]-xt1[0])/(2.*a[0]);
  A(1,1) = (yt1[1]-yt1[0])/(2.*a[0]);
  A(2,1) = (zt1[1]-zt1[0])/(2.*a[0]);
  A(0,2) = (xt1[2]-xt1[0])/(2.*a[1]);
  A(1,2) = (yt1[2]-yt1[0])/(2.*a[1]);
  A(2,2) = (zt1[2]-zt1[0])/(2.*a[1]);
  A(0,3) = (xt1[4]-xt1[0])/(2.*a[2]);
  A(1,3) = (yt1[4]-yt1[0])/(2.*a[2]);
  A(2,3) = (zt1[4]-zt1[0])/(2.*a[2]);

  if (a[0] == 0.) // then problem is fully 2D
  {
    for (E_Int i=0; i < 3; i++)
    {
      if (A(i,2) == 0.) {A(i,1)=1.;} else {A(i,1)=0.;}
    }
  }
 
  // Center of A
  FldArrayF Ca(3);
  Ca[0] = (xt1[7]-xt1[0])/2. + xt1[0];
  Ca[1] = (yt1[7]-yt1[0])/2. + yt1[0];
  Ca[2] = (zt1[7]-zt1[0])/2. + zt1[0];

  /* 
     Operations for OBB(B) 
  */
  E_Float* xt2 = f2->begin(1);
  E_Float* yt2 = f2->begin(2);;
  E_Float* zt2 = f2->begin(3);
  FldArrayF b(3);

  // Calculates the OBB(B) radii
  norm(xt2[0],yt2[0],zt2[0],xt2[1],yt2[1],zt2[1],b[0]);
  norm(xt2[0],yt2[0],zt2[0],xt2[2],yt2[2],zt2[2],b[1]);
  norm(xt2[0],yt2[0],zt2[0],xt2[4],yt2[4],zt2[4],b[2]);
  for (E_Int i=0; i < 3; i++) {b[i] /= 2.;}

  FldArrayF B(3,3);
  B(0,1) = (xt2[1]-xt2[0])/(2.*b[0]);
  B(1,1) = (yt2[1]-yt2[0])/(2.*b[0]);
  B(2,1) = (zt2[1]-zt2[0])/(2.*b[0]);
  B(0,2) = (xt2[2]-xt2[0])/(2.*b[1]);
  B(1,2) = (yt2[2]-yt2[0])/(2.*b[1]);
  B(2,2) = (zt2[2]-zt2[0])/(2.*b[1]);
  B(0,3) = (xt2[4]-xt2[0])/(2.*b[2]);
  B(1,3) = (yt2[4]-yt2[0])/(2.*b[2]);
  B(2,3) = (zt2[4]-zt2[0])/(2.*b[2]);

  if (b[0] == 0.) // then problem is fully 2D
  {
    for (E_Int i=0; i < 3; i++)
    {
      if (B(i,2) == 0.) {B(i,1)=1.;} else {B(i,1)=0.;}
    }
  }

  // Center of B
  FldArrayF Cb(3);
  Cb[0] = (xt2[7]-xt2[0])/2. + xt2[0];
  Cb[1] = (yt2[7]-yt2[0])/2. + yt2[0];
  Cb[2] = (zt2[7]-zt2[0])/2. + zt2[0];

  /* 
     Rotates and translates the OBBs.
     This is done in order to change the reference frame from the absolute
     one to a OBB(A) relative one, so that A becomes the identity matrix
     and B the Rotation matrix relative to the reference frame A.
  */

  // rotation matrix
  FldArrayF R(3,3);
  FldArrayF aR(3,3);
  FldArrayF invA(3,3);
  FldArrayF Tr(3);
  FldArrayF CbT(3);
  FldArrayF CaT(3);
  K_LINEAR::inv3(A.begin(), invA.begin());
  K_LINEAR::prodm(3, 3, 3, invA.begin(), B.begin(), R.begin());
  K_LINEAR::prodv(3, 3, invA.begin(), Ca.begin(), CaT.begin());
  K_LINEAR::prodv(3, 3, invA.begin(), Cb.begin(), CbT.begin());
  
  // Absolute value of the rotation matrix and Epsilon correction
  E_Float aRij;
  for (E_Int j=1; j < 4; j++)
  {
    for (E_Int i=0; i < 3; i++)
    {
       aRij = K_FUNC::E_abs(R(i,j));
       if (aRij < 1.e-12) {aRij=0.; R(i,j)=0.;}
       aR(i,j) = aRij;
     }  
  }
  // Translation vector
  for (E_Int i=0; i < 3; i++)
  {
    Tr[i] = CbT[i] - CaT[i];
  }
  
  // Tests for existing separating axis:
  E_Int  isIntersect = 0;
  if( K_FUNC::E_abs(Tr[0]) > (a[0] + b[0]*aR(0,1) + b[1]*aR(0,2) + b[2]*aR(0,3)))
  {
   isIntersect = 0;
  //cout<<"L=A1 is a Separating Axis"<<endl;
  }
  else if( K_FUNC::E_abs(Tr[1]) > (a[1] + b[0]*aR(1,1) + b[1]*aR(1,2) + b[2]*aR(1,3)))
  {
   isIntersect = 0;
  //cout<<"L=A2 is a Separating Axis"<<endl;
  }
  else if(   K_FUNC::E_abs(Tr[2]) > (a[2] + b[0]*aR(2,1) + b[1]*aR(2,2) + b[2]*aR(2,3))  )
  {
   isIntersect = 0;
  //cout<<"L=A3 is a Separating Axis"<<endl;
  }
  else if(  K_FUNC::E_abs(Tr[0]*R(0,1) + Tr[1]*R(1,1) + Tr[2]*R(2,1)) > (a[0]*aR(0,1) + a[1]*aR(1,1) + a[2]*aR(2,1) + b[0])  )
  {
   isIntersect = 0;
  //cout<<"L=B1 is a Separating Axis"<<endl;
  }
  else if(  K_FUNC::E_abs(Tr[0]*R(0,2) + Tr[1]*R(1,2) + Tr[2]*R(2,2)) > (a[0]*aR(0,2) + a[1]*aR(1,2) + a[2]*aR(2,2) + b[1]) )
  {
   isIntersect = 0;
  //cout<<"L=B2 is a Separating Axis"<<endl;
  }
  else if(  K_FUNC::E_abs(Tr[0]*R(0,3) + Tr[1]*R(1,3) + Tr[2]*R(2,3)) > (a[0]*aR(0,3) + a[1]*aR(1,3) + a[2]*aR(2,3) + b[2]) )
  {
   isIntersect = 0;
  //cout<<"L=B3 is a Separating Axis"<<endl;
  }
  else if(  K_FUNC::E_abs(Tr[2]*R(1,1) - Tr[1]*R(2,1)) > (a[1]*aR(2,1) + a[2]*aR(1,1) + b[1]*aR(0,3) + b[2]*aR(0,2)) )
  {
   isIntersect = 0;
  //cout<<"A1xB1 is a Separating Axis"<<endl;
  }
  else if(   K_FUNC::E_abs(Tr[2]*R(1,2) - Tr[1]*R(2,2)) > (a[1]*aR(2,2) + a[2]*aR(1,2) + b[0]*aR(0,3) + b[2]*aR(0,1))  )
  {
   isIntersect = 0;
  //cout<<"A1xB2 is a Separating Axis"<<endl;
  }
  else if(   K_FUNC::E_abs(Tr[2]*R(1,3) - Tr[1]*R(2,3)) > (a[1]*aR(2,3) + a[2]*aR(1,3) + b[0]*aR(0,2) + b[1]*aR(0,1))  )
  {
   isIntersect = 0;
  //cout<<"A1xB3 is a Separating Axis"<<endl;
  }
  else if(  K_FUNC::E_abs(Tr[0]*R(2,1) - Tr[2]*R(0,1)) > (a[0]*aR(2,1) + a[2]*aR(0,1) + b[1]*aR(1,3) + b[2]*aR(1,2))  )
  {
   isIntersect = 0;
  //cout<<"A2xB1 is a Separating Axis"<<endl;
  }
  else if(  K_FUNC::E_abs(Tr[0]*R(2,2) - Tr[2]*R(0,2)) > (a[0]*aR(2,2) + a[2]*aR(0,2) + b[0]*aR(1,3) + b[2]*aR(1,1)) )
  {
   isIntersect = 0;
  //cout << "A2xB2 is a Separating Axis" << endl;
  }
  else if(  K_FUNC::E_abs(Tr[0]*R(2,3) - Tr[2]*R(0,3)) > (a[0]*aR(2,3) + a[2]*aR(0,3) + b[0]*aR(1,2) + b[1]*aR(1,1)) )
  {
   isIntersect = 0;
  //cout<<"A2xB3 is a Separating Axis"<<endl;
  }
  else if(  K_FUNC::E_abs(Tr[1]*R(0,1) - Tr[0]*R(1,1)) > (a[0]*aR(1,1) + a[1]*aR(0,1) + b[1]*aR(2,3) + b[2]*aR(2,2)) )
  {
   isIntersect = 0;
  //cout<<"A3xB1 is a Separating Axis"<<endl;
  }
  else if(  K_FUNC::E_abs(Tr[1]*R(0,2) - Tr[0]*R(1,2)) > (a[0]*aR(1,2) + a[1]*aR(0,2) + b[0]*aR(2,3) + b[2]*aR(2,1)) )
  {
   isIntersect = 0;
  //cout<<"A3xB2 is a Separating Axis"<<endl;
  }
  else if(  K_FUNC::E_abs(Tr[1]*R(0,3) - Tr[0]*R(1,3)) > (a[0]*aR(1,3) + a[1]*aR(0,3) + b[0]*aR(2,2) + b[1]*aR(2,1)) )
  {
   isIntersect = 0;
  //cout<<"A3xB3 is a Separating Axis"<<endl;
  }
  else 
  {
   isIntersect = 1;
  //cout << "OBBs are overlapping" << endl;
  }

  RELEASESHAREDB(res1, array1, f1, cn1); 
  RELEASESHAREDB(res2, array2, f2, cn2); 

  return Py_BuildValue(I_, isIntersect);
}

//===========================================================================
/* Determine if 2 bounding boxes of 2 arrays intersects
   (zone, in place) */
//===========================================================================
PyObject* K_GENERATOR::_obboxIntersectionZ(PyObject* self, PyObject* args)
{
  char* GridCoordinates; char* FlowSolutionNodes; char* FlowSolutionCenters;
  PyObject* zone1; PyObject* zone2;
  if (!PYPARSETUPLE_(args, OO_ SSS_,
                    &zone1, &zone2, &GridCoordinates,  &FlowSolutionNodes,
                    &FlowSolutionCenters)) return NULL;

  // Checks coordinates of zone 1
  vector<PyArrayObject*> hook1;
  E_Int im1, jm1, km1, cnSize1, cnNfld1;
  char* varString1; char* eltType1;
  vector<E_Float*> fields1; vector<E_Int> locs1;
  vector<E_Int*> cn1;
  K_PYTREE::getFromZone(zone1, 1, 0, varString1, fields1, locs1, im1, jm1, km1, 
                        cn1, cnSize1, cnNfld1, eltType1, hook1,
                        GridCoordinates, FlowSolutionNodes, FlowSolutionCenters);
  E_Int posx1 = K_ARRAY::isCoordinateXPresent(varString1);
  E_Int posy1 = K_ARRAY::isCoordinateYPresent(varString1);
  E_Int posz1 = K_ARRAY::isCoordinateZPresent(varString1);
  if (posx1 == -1 || posy1 == -1 || posz1 == -1)
  {
    delete [] varString1; delete [] eltType1;
    RELEASESHAREDZ(hook1, (char*)NULL, (char*)NULL);
    PyErr_SetString(PyExc_TypeError,
                    "obboxIntersection: cannot find coordinates in zone1.");
    return NULL;
  }

  // Checks coordinates of zone 2
  vector<PyArrayObject*> hook2;
  E_Int im2, jm2, km2, cnSize2, cnNfld2;
  char* varString2; char* eltType2;
  vector<E_Float*> fields2; vector<E_Int> locs2;
  vector<E_Int*> cn2;
  K_PYTREE::getFromZone(zone2, 1, 0, varString2, fields2, locs2, im2, jm2, km2, 
                        cn2, cnSize2, cnNfld2, eltType2, hook2,
                        GridCoordinates, FlowSolutionNodes, FlowSolutionCenters);
  E_Int posx2 = K_ARRAY::isCoordinateXPresent(varString2);
  E_Int posy2 = K_ARRAY::isCoordinateYPresent(varString2);
  E_Int posz2 = K_ARRAY::isCoordinateZPresent(varString2);
  if (posx2 == -1 || posy2 == -1 || posz2 == -1)
  {
    delete [] varString1; delete [] eltType1;
    delete [] varString2; delete [] eltType2;
    RELEASESHAREDZ(hook1, (char*)NULL, (char*)NULL);
    RELEASESHAREDZ(hook2, (char*)NULL, (char*)NULL);
    PyErr_SetString(PyExc_TypeError,
                    "obboxIntersection: cannot find coordinates in zone2.");
    return NULL;
  }

  /* 
  Operations for OBB(A) 
  */
  E_Float* xt1 = fields1[posx1]; 
  E_Float* yt1 = fields1[posy1]; 
  E_Float* zt1 = fields1[posz1]; 
  FldArrayF a(3);

  // Calculates the OBB(A) radii
  norm(xt1[0],yt1[0],zt1[0],xt1[1],yt1[1],zt1[1],a[0]);
  norm(xt1[0],yt1[0],zt1[0],xt1[2],yt1[2],zt1[2],a[1]);
  norm(xt1[0],yt1[0],zt1[0],xt1[4],yt1[4],zt1[4],a[2]);
  for (E_Int i=0; i < 3; i++) {a[i] /= 2.;}


  FldArrayF A(3,3);
  A(0,1) = (xt1[1]-xt1[0])/(2.*a[0]);
  A(1,1) = (yt1[1]-yt1[0])/(2.*a[0]);
  A(2,1) = (zt1[1]-zt1[0])/(2.*a[0]);
  A(0,2) = (xt1[2]-xt1[0])/(2.*a[1]);
  A(1,2) = (yt1[2]-yt1[0])/(2.*a[1]);
  A(2,2) = (zt1[2]-zt1[0])/(2.*a[1]);
  A(0,3) = (xt1[4]-xt1[0])/(2.*a[2]);
  A(1,3) = (yt1[4]-yt1[0])/(2.*a[2]);
  A(2,3) = (zt1[4]-zt1[0])/(2.*a[2]);

  if (a[0] == 0.) // then problem is fully 2D
  {
    for (E_Int i=0; i < 3; i++)
    {
      if (A(i,2) == 0.) {A(i,1)=1.;} else {A(i,1)=0.;}
    }
  }
 
  // Center of A
  FldArrayF Ca(3);
  Ca[0] = (xt1[7]-xt1[0])/2. + xt1[0];
  Ca[1] = (yt1[7]-yt1[0])/2. + yt1[0];
  Ca[2] = (zt1[7]-zt1[0])/2. + zt1[0];


  /* 
  Operations for OBB(B) 
  */
  E_Float* xt2 = fields2[posx2]; 
  E_Float* yt2 = fields2[posy2]; 
  E_Float* zt2 = fields2[posz2]; 
  FldArrayF b(3);

  // Calculates the OBB(B) radii
  norm(xt2[0],yt2[0],zt2[0],xt2[1],yt2[1],zt2[1],b[0]);
  norm(xt2[0],yt2[0],zt2[0],xt2[2],yt2[2],zt2[2],b[1]);
  norm(xt2[0],yt2[0],zt2[0],xt2[4],yt2[4],zt2[4],b[2]);
  for (E_Int i=0; i < 3; i++) {b[i] /= 2.;}

  FldArrayF B(3,3);
  B(0,1) = (xt2[1]-xt2[0])/(2.*b[0]);
  B(1,1) = (yt2[1]-yt2[0])/(2.*b[0]);
  B(2,1) = (zt2[1]-zt2[0])/(2.*b[0]);
  B(0,2) = (xt2[2]-xt2[0])/(2.*b[1]);
  B(1,2) = (yt2[2]-yt2[0])/(2.*b[1]);
  B(2,2) = (zt2[2]-zt2[0])/(2.*b[1]);
  B(0,3) = (xt2[4]-xt2[0])/(2.*b[2]);
  B(1,3) = (yt2[4]-yt2[0])/(2.*b[2]);
  B(2,3) = (zt2[4]-zt2[0])/(2.*b[2]);

  if (b[0] == 0.) // then problem is fully 2D
  {
    for (E_Int i=0; i < 3; i++)
    {
      if (B(i,2) == 0.) {B(i,1)=1.;} else {B(i,1)=0.;}
    }
  }

  // Center of B
  FldArrayF Cb(3);
  Cb[0] = (xt2[7]-xt2[0])/2. + xt2[0];
  Cb[1] = (yt2[7]-yt2[0])/2. + yt2[0];
  Cb[2] = (zt2[7]-zt2[0])/2. + zt2[0];
  

  /* 
     Rotates and translates the OBBs.
     This is done in order to change the reference frame from the absolute
     one to a OBB(A) relative one, so that A becomes the identity matrix
     and B the Rotation matrix relative to the reference frame A.
  */

  // rotation matrix
  FldArrayF R(3,3);
  FldArrayF aR(3,3);
  FldArrayF invA(3,3);
  FldArrayF Tr(3);
  FldArrayF CbT(3);
  FldArrayF CaT(3);
  K_LINEAR::inv3(A.begin(), invA.begin());
  K_LINEAR::prodm(3, 3, 3, invA.begin(), B.begin(), R.begin());
  K_LINEAR::prodv(3, 3, invA.begin(), Ca.begin(), CaT.begin());
  K_LINEAR::prodv(3, 3, invA.begin(), Cb.begin(), CbT.begin());
  
  // Absolute value of the rotation matrix and Epsilon correction
  E_Float aRij;
  for (E_Int j=1; j < 4; j++)
  {
    for (E_Int i=0; i < 3; i++)
    {
       aRij = K_FUNC::E_abs(R(i,j));
       if (aRij < 1.e-12) {aRij=0.; R(i,j)=0.;}
       aR(i,j) = aRij;
     }  
  }
  // Translation vector
  for (E_Int i=0; i < 3; i++)
  {
    Tr[i] = CbT[i] - CaT[i];
  }
  
  // Tests for existing separating axis:
  E_Int  isIntersect = 0;
  if( K_FUNC::E_abs(Tr[0]) > (a[0] + b[0]*aR(0,1) + b[1]*aR(0,2) + b[2]*aR(0,3)))
  {
   isIntersect = 0;
  //cout<<"L=A1 is a Separating Axis"<<endl;
  }
  else if( K_FUNC::E_abs(Tr[1]) > (a[1] + b[0]*aR(1,1) + b[1]*aR(1,2) + b[2]*aR(1,3)))
  {
   isIntersect = 0;
  //cout<<"L=A2 is a Separating Axis"<<endl;
  }
  else if(   K_FUNC::E_abs(Tr[2]) > (a[2] + b[0]*aR(2,1) + b[1]*aR(2,2) + b[2]*aR(2,3))  )
  {
   isIntersect = 0;
  //cout<<"L=A3 is a Separating Axis"<<endl;
  }
  else if(  K_FUNC::E_abs(Tr[0]*R(0,1) + Tr[1]*R(1,1) + Tr[2]*R(2,1)) > (a[0]*aR(0,1) + a[1]*aR(1,1) + a[2]*aR(2,1) + b[0])  )
  {
   isIntersect = 0;
  //cout<<"L=B1 is a Separating Axis"<<endl;
  }
  else if(  K_FUNC::E_abs(Tr[0]*R(0,2) + Tr[1]*R(1,2) + Tr[2]*R(2,2)) > (a[0]*aR(0,2) + a[1]*aR(1,2) + a[2]*aR(2,2) + b[1]) )
  {
   isIntersect = 0;
  //cout<<"L=B2 is a Separating Axis"<<endl;
  }
  else if(  K_FUNC::E_abs(Tr[0]*R(0,3) + Tr[1]*R(1,3) + Tr[2]*R(2,3)) > (a[0]*aR(0,3) + a[1]*aR(1,3) + a[2]*aR(2,3) + b[2]) )
  {
   isIntersect = 0;
  //cout<<"L=B3 is a Separating Axis"<<endl;
  }
  else if(  K_FUNC::E_abs(Tr[2]*R(1,1) - Tr[1]*R(2,1)) > (a[1]*aR(2,1) + a[2]*aR(1,1) + b[1]*aR(0,3) + b[2]*aR(0,2)) )
  {
   isIntersect = 0;
  //cout<<"A1xB1 is a Separating Axis"<<endl;
  }
  else if(   K_FUNC::E_abs(Tr[2]*R(1,2) - Tr[1]*R(2,2)) > (a[1]*aR(2,2) + a[2]*aR(1,2) + b[0]*aR(0,3) + b[2]*aR(0,1))  )
  {
   isIntersect = 0;
  //cout<<"A1xB2 is a Separating Axis"<<endl;
  }
  else if(   K_FUNC::E_abs(Tr[2]*R(1,3) - Tr[1]*R(2,3)) > (a[1]*aR(2,3) + a[2]*aR(1,3) + b[0]*aR(0,2) + b[1]*aR(0,1))  )
  {
   isIntersect = 0;
  //cout<<"A1xB3 is a Separating Axis"<<endl;
  }
  else if(  K_FUNC::E_abs(Tr[0]*R(2,1) - Tr[2]*R(0,1)) > (a[0]*aR(2,1) + a[2]*aR(0,1) + b[1]*aR(1,3) + b[2]*aR(1,2))  )
  {
   isIntersect = 0;
  //cout<<"A2xB1 is a Separating Axis"<<endl;
  }
  else if(  K_FUNC::E_abs(Tr[0]*R(2,2) - Tr[2]*R(0,2)) > (a[0]*aR(2,2) + a[2]*aR(0,2) + b[0]*aR(1,3) + b[2]*aR(1,1)) )
  {
   isIntersect = 0;
  //cout << "A2xB2 is a Separating Axis" << endl;
  }
  else if(  K_FUNC::E_abs(Tr[0]*R(2,3) - Tr[2]*R(0,3)) > (a[0]*aR(2,3) + a[2]*aR(0,3) + b[0]*aR(1,2) + b[1]*aR(1,1)) )
  {
   isIntersect = 0;
  //cout<<"A2xB3 is a Separating Axis"<<endl;
  }
  else if(  K_FUNC::E_abs(Tr[1]*R(0,1) - Tr[0]*R(1,1)) > (a[0]*aR(1,1) + a[1]*aR(0,1) + b[1]*aR(2,3) + b[2]*aR(2,2)) )
  {
   isIntersect = 0;
  //cout<<"A3xB1 is a Separating Axis"<<endl;
  }
  else if(  K_FUNC::E_abs(Tr[1]*R(0,2) - Tr[0]*R(1,2)) > (a[0]*aR(1,2) + a[1]*aR(0,2) + b[0]*aR(2,3) + b[2]*aR(2,1)) )
  {
   isIntersect = 0;
  //cout<<"A3xB2 is a Separating Axis"<<endl;
  }
  else if(  K_FUNC::E_abs(Tr[1]*R(0,3) - Tr[0]*R(1,3)) > (a[0]*aR(1,3) + a[1]*aR(0,3) + b[0]*aR(2,2) + b[1]*aR(2,1)) )
  {
   isIntersect = 0;
  //cout<<"A3xB3 is a Separating Axis"<<endl;
  }
  else 
  {
   isIntersect = 1;
  //cout << "OBBs are overlapping" << endl;
  }
  
  delete [] varString1; delete [] eltType1;
  delete [] varString2; delete [] eltType2;
  RELEASESHAREDZ(hook1, (char*)NULL, (char*)NULL);
  RELEASESHAREDZ(hook2, (char*)NULL, (char*)NULL);

  return Py_BuildValue(I_, isIntersect);
}

//===========================================================================
/* Determine if an Axis-Aligned Bounding Box (array1) intersects an
   Oriented Bounding Box (array2)
   (zone, in place) */
//===========================================================================
PyObject* K_GENERATOR::crossIntersection(PyObject* self, PyObject* args)
{
  PyObject* array1; PyObject* array2;
  if (!PYPARSETUPLE_(args, OO_, &array1, &array2)) return NULL;

  // Check array1
  E_Int ni1, nj1, nk1;
  FldArrayF* f1; FldArrayI* cn1;
  char* varString1; char* eltType1;
  E_Int res1 = K_ARRAY::getFromArray(
    array1, varString1, f1, ni1, nj1, nk1, cn1, eltType1, true);
    
  if (res1 != 1 && res1 != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "crossIntersection: array1 is invalid.");
    return NULL;
  }
  // Check array2
  E_Int ni2, nj2, nk2;
  FldArrayF* f2; FldArrayI* cn2;
  char* varString2; char* eltType2;
  E_Int res2 = K_ARRAY::getFromArray(
    array2, varString2, f2, ni2, nj2, nk2, cn2, eltType2, true);
  
  if (res2 != 1 && res2 != 2)
  {
    PyErr_SetString(PyExc_TypeError,
                    "crossIntersection: array2 is invalid.");
    return NULL;
    RELEASESHAREDB(res1, array1, f1, cn1); 
  }

  // determination de x,y et z ds les arrays
  E_Int posx1 = K_ARRAY::isCoordinateXPresent(varString1);
  E_Int posy1 = K_ARRAY::isCoordinateYPresent(varString1);
  E_Int posz1 = K_ARRAY::isCoordinateZPresent(varString1);
  E_Int posx2 = K_ARRAY::isCoordinateXPresent(varString2);
  E_Int posy2 = K_ARRAY::isCoordinateYPresent(varString2);
  E_Int posz2 = K_ARRAY::isCoordinateZPresent(varString2);
  if (posx1 == -1 || posy1 == -1 || posz1 == -1)
  {
    RELEASESHAREDB(res1, array1, f1, cn1); 
    RELEASESHAREDB(res2, array2, f2, cn2); 
    PyErr_SetString(PyExc_TypeError,
                    "crossIntersection: can't find coordinates in array1.");
    return NULL;
  }
  if (posx2 == -1 || posy2 == -1 || posz2 == -1)
  {
    RELEASESHAREDB(res1, array1, f1, cn1); 
    RELEASESHAREDB(res2, array2, f2, cn2); 
    PyErr_SetString(PyExc_TypeError,
                    "crossIntersection: can't find coordinates in array2.");
    return NULL;
  }
  posx1++; posy1++; posz1++; 
  posx2++; posy2++; posz2++;

  /* 
  Operations for OBB(A) 
  */
  E_Float* xt1 = f1->begin(1);
  E_Float* yt1 = f1->begin(2);
  E_Float* zt1 = f1->begin(3);
  FldArrayF a(3);

  // Calculates the OBB(A) radii
  norm(xt1[0],yt1[0],zt1[0],xt1[1],yt1[1],zt1[1],a[0]);
  norm(xt1[0],yt1[0],zt1[0],xt1[2],yt1[2],zt1[2],a[1]);
  norm(xt1[0],yt1[0],zt1[0],xt1[4],yt1[4],zt1[4],a[2]);
  for (E_Int i=0; i < 3; i++) {a[i] /= 2.;}

  FldArrayF A(3,3);
  A(0,1) = (xt1[1]-xt1[0])/(2.*a[0]);
  A(1,1) = (yt1[1]-yt1[0])/(2.*a[0]);
  A(2,1) = (zt1[1]-zt1[0])/(2.*a[0]);
  A(0,2) = (xt1[2]-xt1[0])/(2.*a[1]);
  A(1,2) = (yt1[2]-yt1[0])/(2.*a[1]);
  A(2,2) = (zt1[2]-zt1[0])/(2.*a[1]);
  A(0,3) = (xt1[4]-xt1[0])/(2.*a[2]);
  A(1,3) = (yt1[4]-yt1[0])/(2.*a[2]);
  A(2,3) = (zt1[4]-zt1[0])/(2.*a[2]);

  if (a[0] == 0.) // then problem is fully 2D
  {
    for (E_Int i=0; i < 3; i++)
    {
      if (A(i,2) == 0.) {A(i,1)=1.;} else {A(i,1)=0.;}
    }
  }

  /* 
     Operations for OBB(B) 
  */
  E_Float* xt2 = f2->begin(1);
  E_Float* yt2 = f2->begin(2);
  E_Float* zt2 = f2->begin(3);
  FldArrayF b(3);

  // Calculates the OBB(B) radii
  norm(xt2[0],yt2[0],zt2[0],xt2[1],yt2[1],zt2[1],b[0]);
  norm(xt2[0],yt2[0],zt2[0],xt2[2],yt2[2],zt2[2],b[1]);
  norm(xt2[0],yt2[0],zt2[0],xt2[4],yt2[4],zt2[4],b[2]);
  for (E_Int i=0; i < 3; i++) {b[i] /= 2.;}

  // Rotation matrix is exactly B
  FldArrayF R(3,3);
  R(0,1) = (xt2[1]-xt2[0])/(2.*b[0]);
  R(1,1) = (yt2[1]-yt2[0])/(2.*b[0]);
  R(2,1) = (zt2[1]-zt2[0])/(2.*b[0]);
  R(0,2) = (xt2[2]-xt2[0])/(2.*b[1]);
  R(1,2) = (yt2[2]-yt2[0])/(2.*b[1]);
  R(2,2) = (zt2[2]-zt2[0])/(2.*b[1]);
  R(0,3) = (xt2[4]-xt2[0])/(2.*b[2]);
  R(1,3) = (yt2[4]-yt2[0])/(2.*b[2]);
  R(2,3) = (zt2[4]-zt2[0])/(2.*b[2]);

  if (b[0] == 0.) // then problem is fully 2D
  {
    for (E_Int i=0; i < 3; i++)
    {
      if (R(i,2) == 0.) {R(i,1)=1.;} else {R(i,1)=0.;}
    }
  }

  // Translation vector is center of B minus center of A
  FldArrayF Tr(3);
  Tr[0] = (xt2[7]-xt2[0])/2. + xt2[0] - ((xt1[7]-xt1[0])/2. + xt1[0]);
  Tr[1] = (yt2[7]-yt2[0])/2. + yt2[0] - ((yt1[7]-yt1[0])/2. + yt1[0]);
  Tr[2] = (zt2[7]-zt2[0])/2. + zt2[0] - ((zt1[7]-zt1[0])/2. + zt1[0]);
  

  // Absolute value of the rotation matrix and Epsilon correction
  FldArrayF aR(3,3);
  E_Float aRij;
  for (E_Int j=1; j < 4; j++)
  {
    for (E_Int i=0; i < 3; i++)
    {
       aRij = K_FUNC::E_abs(R(i,j));
       if (aRij < 1.e-12) {aRij=0.; R(i,j)=0.;}
       aR(i,j) = aRij;
     }  
  }
  
  // Tests for existing separating axis:
  E_Int  isIntersect = 0;
  if( K_FUNC::E_abs(Tr[0]) > (a[0] + b[0]*aR(0,1) + b[1]*aR(0,2) + b[2]*aR(0,3)))
  {
   isIntersect = 0;
  //cout<<"L=A1 is a Separating Axis"<<endl;
  }
  else if( K_FUNC::E_abs(Tr[1]) > (a[1] + b[0]*aR(1,1) + b[1]*aR(1,2) + b[2]*aR(1,3)))
  {
   isIntersect = 0;
  //cout<<"L=A2 is a Separating Axis"<<endl;
  }
  else if(   K_FUNC::E_abs(Tr[2]) > (a[2] + b[0]*aR(2,1) + b[1]*aR(2,2) + b[2]*aR(2,3))  )
  {
   isIntersect = 0;
  //cout<<"L=A3 is a Separating Axis"<<endl;
  }
  else if(  K_FUNC::E_abs(Tr[0]*R(0,1) + Tr[1]*R(1,1) + Tr[2]*R(2,1)) > (a[0]*aR(0,1) + a[1]*aR(1,1) + a[2]*aR(2,1) + b[0])  )
  {
   isIntersect = 0;
  //cout<<"L=B1 is a Separating Axis"<<endl;
  }
  else if(  K_FUNC::E_abs(Tr[0]*R(0,2) + Tr[1]*R(1,2) + Tr[2]*R(2,2)) > (a[0]*aR(0,2) + a[1]*aR(1,2) + a[2]*aR(2,2) + b[1]) )
  {
   isIntersect = 0;
  //cout<<"L=B2 is a Separating Axis"<<endl;
  }
  else if(  K_FUNC::E_abs(Tr[0]*R(0,3) + Tr[1]*R(1,3) + Tr[2]*R(2,3)) > (a[0]*aR(0,3) + a[1]*aR(1,3) + a[2]*aR(2,3) + b[2]) )
  {
   isIntersect = 0;
  //cout<<"L=B3 is a Separating Axis"<<endl;
  }
  else if(  K_FUNC::E_abs(Tr[2]*R(1,1) - Tr[1]*R(2,1)) > (a[1]*aR(2,1) + a[2]*aR(1,1) + b[1]*aR(0,3) + b[2]*aR(0,2)) )
  {
   isIntersect = 0;
  //cout<<"A1xB1 is a Separating Axis"<<endl;
  }
  else if(   K_FUNC::E_abs(Tr[2]*R(1,2) - Tr[1]*R(2,2)) > (a[1]*aR(2,2) + a[2]*aR(1,2) + b[0]*aR(0,3) + b[2]*aR(0,1))  )
  {
   isIntersect = 0;
  //cout<<"A1xB2 is a Separating Axis"<<endl;
  }
  else if(   K_FUNC::E_abs(Tr[2]*R(1,3) - Tr[1]*R(2,3)) > (a[1]*aR(2,3) + a[2]*aR(1,3) + b[0]*aR(0,2) + b[1]*aR(0,1))  )
  {
   isIntersect = 0;
  //cout<<"A1xB3 is a Separating Axis"<<endl;
  }
  else if(  K_FUNC::E_abs(Tr[0]*R(2,1) - Tr[2]*R(0,1)) > (a[0]*aR(2,1) + a[2]*aR(0,1) + b[1]*aR(1,3) + b[2]*aR(1,2))  )
  {
   isIntersect = 0;
  //cout<<"A2xB1 is a Separating Axis"<<endl;
  }
  else if(  K_FUNC::E_abs(Tr[0]*R(2,2) - Tr[2]*R(0,2)) > (a[0]*aR(2,2) + a[2]*aR(0,2) + b[0]*aR(1,3) + b[2]*aR(1,1)) )
  {
   isIntersect = 0;
  //cout << "A2xB2 is a Separating Axis" << endl;
  }
  else if(  K_FUNC::E_abs(Tr[0]*R(2,3) - Tr[2]*R(0,3)) > (a[0]*aR(2,3) + a[2]*aR(0,3) + b[0]*aR(1,2) + b[1]*aR(1,1)) )
  {
   isIntersect = 0;
  //cout<<"A2xB3 is a Separating Axis"<<endl;
  }
  else if(  K_FUNC::E_abs(Tr[1]*R(0,1) - Tr[0]*R(1,1)) > (a[0]*aR(1,1) + a[1]*aR(0,1) + b[1]*aR(2,3) + b[2]*aR(2,2)) )
  {
   isIntersect = 0;
  //cout<<"A3xB1 is a Separating Axis"<<endl;
  }
  else if(  K_FUNC::E_abs(Tr[1]*R(0,2) - Tr[0]*R(1,2)) > (a[0]*aR(1,2) + a[1]*aR(0,2) + b[0]*aR(2,3) + b[2]*aR(2,1)) )
  {
   isIntersect = 0;
  //cout<<"A3xB2 is a Separating Axis"<<endl;
  }
  else if(  K_FUNC::E_abs(Tr[1]*R(0,3) - Tr[0]*R(1,3)) > (a[0]*aR(1,3) + a[1]*aR(0,3) + b[0]*aR(2,2) + b[1]*aR(2,1)) )
  {
   isIntersect = 0;
  //cout<<"A3xB3 is a Separating Axis"<<endl;
  }
  else 
  {
   isIntersect = 1;
  //cout << "OBBs are overlapping" << endl;
  }

  RELEASESHAREDB(res1, array1, f1, cn1); 
  RELEASESHAREDB(res2, array2, f2, cn2); 

  return Py_BuildValue(I_, isIntersect);
}

//===========================================================================
/* Determine if an Axis-Aligned Bounding Box (zone1) intersects an
   Oriented Bounding Box (zone2)
   (zone, in place) */
//===========================================================================
PyObject* K_GENERATOR::_crossIntersectionZ(PyObject* self, PyObject* args)
{
  char* GridCoordinates; char* FlowSolutionNodes; char* FlowSolutionCenters;
  PyObject* zone1; PyObject* zone2;
  if (!PYPARSETUPLE_(args, OO_ SSS_,
                    &zone1, &zone2, &GridCoordinates,  &FlowSolutionNodes,
                    &FlowSolutionCenters)) return NULL;

  // Checks coordinates of zone 1
  vector<PyArrayObject*> hook1;
  E_Int im1, jm1, km1, cnSize1, cnNfld1;
  char* varString1; char* eltType1;
  vector<E_Float*> fields1; vector<E_Int> locs1;
  vector<E_Int*> cn1;
  K_PYTREE::getFromZone(zone1, 1, 0, varString1, fields1, locs1, im1, jm1, km1, 
                        cn1, cnSize1, cnNfld1, eltType1, hook1,
                        GridCoordinates, FlowSolutionNodes, FlowSolutionCenters);
  E_Int posx1 = K_ARRAY::isCoordinateXPresent(varString1);
  E_Int posy1 = K_ARRAY::isCoordinateYPresent(varString1);
  E_Int posz1 = K_ARRAY::isCoordinateZPresent(varString1);
  if (posx1 == -1 || posy1 == -1 || posz1 == -1)
  {
    delete [] varString1; delete [] eltType1;
    RELEASESHAREDZ(hook1, (char*)NULL, (char*)NULL);
    PyErr_SetString(PyExc_TypeError,
                    "obboxIntersection: cannot find coordinates in zone1.");
    return NULL;
  }

  // Checks coordinates of zone 2
  vector<PyArrayObject*> hook2;
  E_Int im2, jm2, km2, cnSize2, cnNfld2;
  char* varString2; char* eltType2;
  vector<E_Float*> fields2; vector<E_Int> locs2;
  vector<E_Int*> cn2;
  K_PYTREE::getFromZone(zone2, 1, 0, varString2, fields2, locs2, im2, jm2, km2, 
                        cn2, cnSize2, cnNfld2, eltType2, hook2,
                        GridCoordinates, FlowSolutionNodes, FlowSolutionCenters);
  E_Int posx2 = K_ARRAY::isCoordinateXPresent(varString2);
  E_Int posy2 = K_ARRAY::isCoordinateYPresent(varString2);
  E_Int posz2 = K_ARRAY::isCoordinateZPresent(varString2);
  if (posx2 == -1 || posy2 == -1 || posz2 == -1)
  {
    delete [] varString1; delete [] eltType1;
    delete [] varString2; delete [] eltType2;
    RELEASESHAREDZ(hook1, (char*)NULL, (char*)NULL);
    RELEASESHAREDZ(hook2, (char*)NULL, (char*)NULL);
    PyErr_SetString(PyExc_TypeError,
                    "obboxIntersection: cannot find coordinates in zone2.");
    return NULL;
  }

  /* 
     Operations for OBB(A) 
  */
  E_Float* xt1 = fields1[posx1]; 
  E_Float* yt1 = fields1[posy1]; 
  E_Float* zt1 = fields1[posz1]; 
  FldArrayF a(3);

  // Calculates the OBB(A) radii
  norm(xt1[0],yt1[0],zt1[0],xt1[1],yt1[1],zt1[1],a[0]);
  norm(xt1[0],yt1[0],zt1[0],xt1[2],yt1[2],zt1[2],a[1]);
  norm(xt1[0],yt1[0],zt1[0],xt1[4],yt1[4],zt1[4],a[2]);
  for (E_Int i=0; i < 3; i++) {a[i] /= 2.;}

  FldArrayF A(3,3);
  A(0,1) = (xt1[1]-xt1[0])/(2.*a[0]);
  A(1,1) = (yt1[1]-yt1[0])/(2.*a[0]);
  A(2,1) = (zt1[1]-zt1[0])/(2.*a[0]);
  A(0,2) = (xt1[2]-xt1[0])/(2.*a[1]);
  A(1,2) = (yt1[2]-yt1[0])/(2.*a[1]);
  A(2,2) = (zt1[2]-zt1[0])/(2.*a[1]);
  A(0,3) = (xt1[4]-xt1[0])/(2.*a[2]);
  A(1,3) = (yt1[4]-yt1[0])/(2.*a[2]);
  A(2,3) = (zt1[4]-zt1[0])/(2.*a[2]);

  if (a[0] == 0.) // then problem is fully 2D
  {
    for (E_Int i=0; i < 3; i++)
    {
      if (A(i,2) == 0.) {A(i,1)=1.;} else {A(i,1)=0.;}
    }
  }
 

  /* 
     Operations for OBB(B) 
  */
  E_Float* xt2 = fields2[posx2]; 
  E_Float* yt2 = fields2[posy2]; 
  E_Float* zt2 = fields2[posz2]; 
  FldArrayF b(3);

  // Calculates the OBB(B) radii
  norm(xt2[0],yt2[0],zt2[0],xt2[1],yt2[1],zt2[1],b[0]);
  norm(xt2[0],yt2[0],zt2[0],xt2[2],yt2[2],zt2[2],b[1]);
  norm(xt2[0],yt2[0],zt2[0],xt2[4],yt2[4],zt2[4],b[2]);
  for (E_Int i=0; i < 3; i++) {b[i] /= 2.;}

  // Rotation matrix is exactly B
  FldArrayF R(3,3);
  R(0,1) = (xt2[1]-xt2[0])/(2.*b[0]);
  R(1,1) = (yt2[1]-yt2[0])/(2.*b[0]);
  R(2,1) = (zt2[1]-zt2[0])/(2.*b[0]);
  R(0,2) = (xt2[2]-xt2[0])/(2.*b[1]);
  R(1,2) = (yt2[2]-yt2[0])/(2.*b[1]);
  R(2,2) = (zt2[2]-zt2[0])/(2.*b[1]);
  R(0,3) = (xt2[4]-xt2[0])/(2.*b[2]);
  R(1,3) = (yt2[4]-yt2[0])/(2.*b[2]);
  R(2,3) = (zt2[4]-zt2[0])/(2.*b[2]);

  if (b[0] == 0.) // then problem is fully 2D
  {
    for (E_Int i=0; i < 3; i++)
    {
      if (R(i,2) == 0.) {R(i,1)=1.;} else {R(i,1)=0.;}
    }
  }

  // Translation vector is center of B minus center of A
  FldArrayF Tr(3);
  Tr[0] = (xt2[7]-xt2[0])/2. + xt2[0] - ((xt1[7]-xt1[0])/2. + xt1[0]);
  Tr[1] = (yt2[7]-yt2[0])/2. + yt2[0] - ((yt1[7]-yt1[0])/2. + yt1[0]);
  Tr[2] = (zt2[7]-zt2[0])/2. + zt2[0] - ((zt1[7]-zt1[0])/2. + zt1[0]);
  

  // Absolute value of the rotation matrix and Epsilon correction
  FldArrayF aR(3,3);
  E_Float aRij;
  for (E_Int j=1; j < 4; j++)
  {
    for (E_Int i=0; i < 3; i++)
    {
       aRij = K_FUNC::E_abs(R(i,j));
       if (aRij < 1.e-12) {aRij=0.; R(i,j)=0.;}
       aR(i,j) = aRij;
     }  
  }
 //cout <<  "Rotation matrix" << endl;
 //cout << R << endl;
 //cout <<  "Absolute value Rotation matrix" << endl;
 //cout << aR << endl;  
  
  // Tests for existing separating axis:
  E_Int  isIntersect = 0;
  if( K_FUNC::E_abs(Tr[0]) > (a[0] + b[0]*aR(0,1) + b[1]*aR(0,2) + b[2]*aR(0,3)))
  {
   isIntersect = 0;
  //cout<<"L=A1 is a Separating Axis"<<endl;
  }
  else if( K_FUNC::E_abs(Tr[1]) > (a[1] + b[0]*aR(1,1) + b[1]*aR(1,2) + b[2]*aR(1,3)))
  {
   isIntersect = 0;
  //cout<<"L=A2 is a Separating Axis"<<endl;
  }
  else if(   K_FUNC::E_abs(Tr[2]) > (a[2] + b[0]*aR(2,1) + b[1]*aR(2,2) + b[2]*aR(2,3))  )
  {
   isIntersect = 0;
  //cout<<"L=A3 is a Separating Axis"<<endl;
  }
  else if(  K_FUNC::E_abs(Tr[0]*R(0,1) + Tr[1]*R(1,1) + Tr[2]*R(2,1)) > (a[0]*aR(0,1) + a[1]*aR(1,1) + a[2]*aR(2,1) + b[0])  )
  {
   isIntersect = 0;
  //cout<<"L=B1 is a Separating Axis"<<endl;
  }
  else if(  K_FUNC::E_abs(Tr[0]*R(0,2) + Tr[1]*R(1,2) + Tr[2]*R(2,2)) > (a[0]*aR(0,2) + a[1]*aR(1,2) + a[2]*aR(2,2) + b[1]) )
  {
   isIntersect = 0;
  //cout<<"L=B2 is a Separating Axis"<<endl;
  }
  else if(  K_FUNC::E_abs(Tr[0]*R(0,3) + Tr[1]*R(1,3) + Tr[2]*R(2,3)) > (a[0]*aR(0,3) + a[1]*aR(1,3) + a[2]*aR(2,3) + b[2]) )
  {
   isIntersect = 0;
  //cout<<"L=B3 is a Separating Axis"<<endl;
  }
  else if(  K_FUNC::E_abs(Tr[2]*R(1,1) - Tr[1]*R(2,1)) > (a[1]*aR(2,1) + a[2]*aR(1,1) + b[1]*aR(0,3) + b[2]*aR(0,2)) )
  {
   isIntersect = 0;
  //cout<<"A1xB1 is a Separating Axis"<<endl;
  }
  else if(   K_FUNC::E_abs(Tr[2]*R(1,2) - Tr[1]*R(2,2)) > (a[1]*aR(2,2) + a[2]*aR(1,2) + b[0]*aR(0,3) + b[2]*aR(0,1))  )
  {
   isIntersect = 0;
  //cout<<"A1xB2 is a Separating Axis"<<endl;
  }
  else if(   K_FUNC::E_abs(Tr[2]*R(1,3) - Tr[1]*R(2,3)) > (a[1]*aR(2,3) + a[2]*aR(1,3) + b[0]*aR(0,2) + b[1]*aR(0,1))  )
  {
   isIntersect = 0;
  //cout<<"A1xB3 is a Separating Axis"<<endl;
  }
  else if(  K_FUNC::E_abs(Tr[0]*R(2,1) - Tr[2]*R(0,1)) > (a[0]*aR(2,1) + a[2]*aR(0,1) + b[1]*aR(1,3) + b[2]*aR(1,2))  )
  {
   isIntersect = 0;
  //cout<<"A2xB1 is a Separating Axis"<<endl;
  }
  else if(  K_FUNC::E_abs(Tr[0]*R(2,2) - Tr[2]*R(0,2)) > (a[0]*aR(2,2) + a[2]*aR(0,2) + b[0]*aR(1,3) + b[2]*aR(1,1)) )
  {
   isIntersect = 0;
  //cout << "A2xB2 is a Separating Axis" << endl;
  }
  else if(  K_FUNC::E_abs(Tr[0]*R(2,3) - Tr[2]*R(0,3)) > (a[0]*aR(2,3) + a[2]*aR(0,3) + b[0]*aR(1,2) + b[1]*aR(1,1)) )
  {
   isIntersect = 0;
  //cout<<"A2xB3 is a Separating Axis"<<endl;
  }
  else if(  K_FUNC::E_abs(Tr[1]*R(0,1) - Tr[0]*R(1,1)) > (a[0]*aR(1,1) + a[1]*aR(0,1) + b[1]*aR(2,3) + b[2]*aR(2,2)) )
  {
   isIntersect = 0;
  //cout<<"A3xB1 is a Separating Axis"<<endl;
  }
  else if(  K_FUNC::E_abs(Tr[1]*R(0,2) - Tr[0]*R(1,2)) > (a[0]*aR(1,2) + a[1]*aR(0,2) + b[0]*aR(2,3) + b[2]*aR(2,1)) )
  {
   isIntersect = 0;
  //cout<<"A3xB2 is a Separating Axis"<<endl;
  }
  else if(  K_FUNC::E_abs(Tr[1]*R(0,3) - Tr[0]*R(1,3)) > (a[0]*aR(1,3) + a[1]*aR(0,3) + b[0]*aR(2,2) + b[1]*aR(2,1)) )
  {
   isIntersect = 0;
  //cout<<"A3xB3 is a Separating Axis"<<endl;
  }
  else 
  {
   isIntersect = 1;
  //cout << "OBBs are overlapping" << endl;
  }

  delete [] varString1; delete [] eltType1;
  delete [] varString2; delete [] eltType2;
  RELEASESHAREDZ(hook1, (char*)NULL, (char*)NULL);
  RELEASESHAREDZ(hook2, (char*)NULL, (char*)NULL);

  return Py_BuildValue(I_, isIntersect);
}

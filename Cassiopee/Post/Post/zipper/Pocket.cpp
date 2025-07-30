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

#include "Pocket.h"
#include <vector>

using namespace std;
using namespace K_FUNC;
using namespace K_FLD;
using namespace K_CONST;

//=============================================================================
Pocket::Pocket(FldArrayF& field)
{
  _field = field; 
  _size = _field.getSize();  
  _nfld = field.getNfld();
}

//=============================================================================
Pocket::~Pocket()
{
}
//=============================================================================
/* Write the pocket*/
//=============================================================================
void Pocket::writeLine(char* fileName, E_Bool add)
{ 
  E_Int size1 = _size+1;
  FldArrayF field(size1, _nfld);
  
  for (E_Int eq = 1; eq <= _nfld; eq++)
  {
    for (E_Int i = 0; i < _size; i++)
      field(i, eq) = _field(i, eq);
    
    field(_size,eq) = _field(0,eq);
  }
  // Build connect by elements
  FldArrayI connect( _size, 3 );
  E_Int c = 0;
  for (E_Int i = 0; i < _size; i++)
  {     
    connect(c, 1) = i+1;
    connect(c, 2) = i+2;
    connect(c, 3) = i+2;
    c++;
  }
 
  //if (c > 0)
  //  K_IO::GenIO::getInstance()->tpwriteTriangles(fileName, field, 
  //                                               connect, add);
}
//============================================================================
/* Close pocket by building triangles. Return true if triangles are created */
//============================================================================
E_Bool Pocket::closePocket(FldArrayF& field, FldArrayI& connect)
{
  vector<TriangleZ*> triangles;
  E_Bool closed = false;
  FldArrayF fieldA(_nfld);
  FldArrayF fieldB(_nfld);
  FldArrayF fieldC(_nfld);
  FldArrayF fieldD(_nfld);

  // make triangles
  if ( _size < 3 )
    return false;
  
  else if ( _size == 3)
  {
    for (E_Int eq = 0; eq < _nfld; eq++)
    {
      fieldA[eq] = _field(0, eq+1);
      fieldB[eq] = _field(1, eq+1);
      fieldC[eq] = _field(2, eq+1);
    }
    TriangleZ* t = new TriangleZ(fieldA, fieldB, fieldC);
    triangles.push_back(t);
    closed = true;  
  }
  else 
    closed = computeTriangulation(triangles);  
  
  // Compute the triangle connectivity
  compConnectivity(triangles, field, connect);
  return closed;
}

//=============================================================================
/* Computes triangles starting from istart */
//=============================================================================
E_Bool Pocket::computeTriangulation(vector<TriangleZ*>& triangles)
{
  E_Int iA, iB, iC, iD;
  FldArrayF fieldA(_nfld);
  FldArrayF fieldB(_nfld);
  FldArrayF fieldC(_nfld);
  FldArrayF fieldD(_nfld);
  E_Bool closed = false;
  
  /*-----------------*/
  /* triangulation   */
  /*-----------------*/
  //precond : faire une string si des points sont alignes
  FldArrayF field1;
  FldArrayF field2;
  if (checkIfPtsAreAligned(field1, field2) == true )
  {
    closeSpecificPocket(field1, field2, triangles);
    return true;
  }
  else 
  {
    // INIT 
    iA = 0;
    iB = 1;
    iC = _size-2;
    iD = _size-1; 
    
    while ( iB != iC && iB < _size)
    {
      E_Bool isOK = compDelaunay(iA, iB, iC, iD, 
                                    _field, _field, 
                                    triangles);
      if ( isOK == false )
      {
        for (E_Int eq = 0; eq < _nfld; eq++)
        {
          fieldA[eq] = _field(iA, eq+1);
          fieldB[eq] = _field(iB, eq+1);
          fieldD[eq] = _field(iD, eq+1);
        }
        TriangleZ* t = new TriangleZ(fieldA, fieldB, fieldD);
        triangles.push_back(t);
        iA++;
        iB++;
      }
      if ( iB == iC)
      { 
        for (E_Int eq = 0; eq < _nfld; eq++)
        {
          fieldA[eq] = _field(iA, eq+1);
          fieldB[eq] = _field(iB, eq+1);
          fieldD[eq] = _field(iD, eq+1);
        } 
        TriangleZ* t = new TriangleZ(fieldA, fieldB, fieldD);
        triangles.push_back(t);
        closed = true;
      }
    }
  }
  return closed;  
}

//=============================================================================
/* Compute connectivity for triangles  */
//=============================================================================
void Pocket::compConnectivity( vector<TriangleZ*>& triangles,
                               FldArrayF& field,
                               FldArrayI& connect)
{
  E_Int nelts = triangles.size(); 
  connect.malloc(nelts,3);
  field.malloc(3*nelts, _nfld);
  
  for ( E_Int n = 0; n < nelts; n++ )
  {
    FldArrayF& field1 = triangles[n]->getField1();
    FldArrayF& field2 = triangles[n]->getField2();
    FldArrayF& field3 = triangles[n]->getField3();
    
    for (E_Int eq= 1 ; eq <= _nfld ; eq++)
    {  
      field(3*n  , eq) =  field1[eq-1];
      field(3*n+1, eq) =  field2[eq-1];
      field(3*n+2, eq) =  field3[eq-1];
    }
  }
  
  for ( E_Int n = 0 ; n < nelts; n++ )
  {
    connect(n,1) = 3*n+1;
    connect(n,2) = 3*n+2;
    connect(n,3) = 3*n+3;
  }
}

//=============================================================================
/* Compute triangulation for the quad PiPi+1PjPj+1 of the pocket */
//=============================================================================
E_Bool Pocket::compDelaunay( E_Int& iA, E_Int& iB,
                                E_Int& iC, E_Int& iD,
                                FldArrayF& field1, FldArrayF& field2,
                                vector<TriangleZ*>& triangles)
{  
  /* QUAD PiPi+1PjPj+1 = ABCD */
  FldArrayF fieldA(_nfld);
  FldArrayF fieldB(_nfld);
  FldArrayF fieldC(_nfld);
  FldArrayF fieldD(_nfld);
 
  for (E_Int eq = 0; eq < _nfld; eq++)
  {
    fieldA[eq] = field1(iA, eq+1);
    fieldB[eq] = field1(iB, eq+1);
    fieldC[eq] = field2(iC, eq+1);
    fieldD[eq] = field2(iD, eq+1);
  }
  
  //triangle ADC
  E_Int diagType = 1;
  TriangleZ* t1 = new TriangleZ(fieldA, fieldD, fieldC);
  //diagonal AC
  E_Float diagAC = compDistanceBetweenPoints(fieldA, fieldC);
  E_Bool isValid1 = checkTriangles( diagType, fieldA, fieldB,
                                       fieldC, fieldD);
  //triangle ABD
  diagType = 2;
  TriangleZ* t2 = new TriangleZ(fieldA, fieldD, fieldB);
  //diagonal AD
  E_Float diagBD = compDistanceBetweenPoints(fieldB, fieldD);
  E_Bool isValid2 = checkTriangles( diagType, fieldA, fieldB,
                                       fieldC, fieldD);
      
  if ( isValid1 == true && isValid2 == true )
  {
    //select the smaller diagonal
    if ( diagAC < diagBD ) // AC connection selected
    {
      triangles.push_back(t1);
      //eliminates D for new triangulation research
      iD = iC;
      iC = iC-1;
    }
    else// BD connection selected
    {
      triangles.push_back(t2);
      // eliminate A for new triangulation research 
      iA = iB;
      iB = iB + 1;
    }
  }
  else if ( isValid1 == true && isValid2 == false)
  {
    //(AC) connection valid
    triangles.push_back(t1);
    //eliminates D for new triangulation research
    iD = iC;
    iC = iC-1;
  }
  
  else if ( isValid1 == false && isValid2 == true ) 
  {
    // (BD connection valid)
    triangles.push_back(t2);
    // eliminate A for new triangulation research 
    iA = iB;
    iB = iB + 1;
  }
  
  else return false;
  return true;
}

//=============================================================================
/* Compute the  distance between ind1 and ind2  */
//=============================================================================
E_Float Pocket::compDistanceBetweenPoints(FldArrayF& field1, FldArrayF& field2)
{
  E_Float dx = field1[0] - field2[0];
  E_Float dy = field1[1] - field2[1];
  E_Float dz = field1[2] - field2[2];
  return  dx * dx + dy * dy + dz * dz;
}

//=============================================================================
/* Tell if the connection is valid  for triangulation
   Assume a quad ABCD. To test the validity of the triangle ACD, we must check
   if the point B is outside the circum radius of ACD.
*/
//=============================================================================
E_Bool Pocket::checkTriangles( E_Int diagType,
                                  FldArrayF& fieldA, FldArrayF& fieldB,
                                  FldArrayF& fieldC, FldArrayF& fieldD)
{
  // test if the quad ABCD is valid (convex)
  // if diagType = 1 (AC) common edge of triangles
  // triangles are ABC and ACD
  // n1 = (BA ^ BC)
  // n2 = (DC ^ DA)

  // if diagType = 2 (BD) common edge of triangles
  // triangles are ABD and BCD
  // n1 = (AB ^ AD)
  // n2 = (CD ^ CB)
  
  // quad is valid if n1.n2 > 0

  E_Float eps = 1.e-9;
  E_Float n1x, n1y, n1z; 
  E_Float n2x, n2y, n2z; // normals of triangles
  E_Float x1, x2, x3, x4;
  E_Float y1, y2, y3, y4;
  E_Float z1, z2, z3, z4;
  
  switch ( diagType )
  {
    case 1:
      //vector BA coords
      x1 = fieldA[0] - fieldB[0];
      y1 = fieldA[1] - fieldB[1];
      z1 = fieldA[2] - fieldB[2];
      
      //vector BC coords
      x2 = fieldC[0] - fieldB[0];
      y2 = fieldC[1] - fieldB[1];
      z2 = fieldC[2] - fieldB[2];

      //vector DC coords
      x3 = fieldC[0] - fieldD[0];
      y3 = fieldC[1] - fieldD[1];
      z3 = fieldC[2] - fieldD[2];

      //vector DA coords
      x4 = fieldA[0] - fieldD[0];
      y4 = fieldA[1] - fieldD[1];
      z4 = fieldA[2] - fieldD[2];
      break;
      
    case 2 :
      //vector AB coords
      x1 = fieldB[0] - fieldA[0];
      y1 = fieldB[1] - fieldA[1];
      z1 = fieldB[2] - fieldA[2];
      
      //vector AD coords
      x2 = fieldD[0] - fieldA[0];
      y2 = fieldD[1] - fieldA[1];
      z2 = fieldD[2] - fieldA[2];

      //vector CD coords
      x3 = fieldD[0] - fieldC[0];
      y3 = fieldD[1] - fieldC[1];
      z3 = fieldD[2] - fieldC[2];

      //vector CB coords
      x4 = fieldB[0] - fieldC[0];
      y4 = fieldB[1] - fieldC[1];
      z4 = fieldB[2] - fieldC[2];
      break;

    default:
      x1 = 0.;
      y1 = 0.;
      z1 = 0.;
      x2 = 0.;
      y2 = 0.;
      z2 = 0.;
      x3 = 0.;
      y3 = 0.;
      z3 = 0.;
      x4 = 0.;
      y4 = 0.;
      z4 = 0.;
  }

  n1x = y1 * z2 - y2 * z1;
  n1y = x2 * z1 - x1 * z2;
  n1z = x1 * y2 - x2 * y1;

  n2x = y3 * z4 - y4 * z3;
  n2y = x4 * z3 - x3 * z4;
  n2z = x3 * y4 - x4 * y3;

  E_Float norm1 = n1x*n1x + n1y*n1y + n1z*n1z;
  E_Float norm2 = n2x*n2x + n2y*n2y + n2z*n2z;
  E_Float norm = sqrt(norm1*norm2);
  E_Float ps = n1x * n2x + n1y * n2y + n1z * n2z ;
  ps = E_abs(ps)/norm;
  if ( ps > eps )
    return true ;

  else
    return false;
}
//=============================================================================
/* Check if a set of contiguous points are aligned */
//=============================================================================
E_Bool Pocket::checkIfPtsAreAligned(FldArrayF& field1, FldArrayF& field2)
{
  E_Int size = _field.getSize();
  E_Int size1, size2;
  E_Float dx1, dy1, dz1;
  E_Float dx2, dy2, dz2;
  E_Float n1, n2;
  E_Float ps, psmin, norm;
  E_Float eps = 1.e-5;
  psmin = 1.-eps;
  FldArrayI indices(size);
  FldArrayIS dejaVu(size);
  E_Int cnt = 0;
  E_Int im, ip;

  dejaVu.setAllValuesAt(0);
  for (E_Int ind = 0; ind < size; ind++)
  {
    im = ind-1;
    ip = ind+1; 
    if ( ind == 0 ) im = size-1;
    if ( ind == size-1) ip = 0;
   
    dx1 = _field(ind,1) - _field(im,1);
    dy1 = _field(ind,2) - _field(im,2);
    dz1 = _field(ind,3) - _field(im,3);

    dx2 = _field(ip,1) - _field(ind,1);
    dy2 = _field(ip,2) - _field(ind,2);
    dz2 = _field(ip,3) - _field(ind,3); 
    
    n1 = dx1 * dx1 + dy1 * dy1 + dz1 * dz1;
    n2 = dx2 * dx2 + dy2 * dy2 + dz2 * dz2;
    norm =sqrt(n1*n2);

    ps =  dx1*dx2 + dy1*dy2 + dz1*dz2;
    ps = ps/(norm);
    ps = E_abs(ps);
    if ( ps > psmin)
    {
      if (dejaVu[im] == 0)
      {
        indices[cnt] = im;
        dejaVu[im] = 1;
        cnt++;
      }
      if (dejaVu[ind] == 0)
      {
        indices[cnt] = ind;
        dejaVu[ind] = 1;
        cnt++;
      }
      if (dejaVu[ip] ==0)
      {
        indices[cnt] = ip;
        dejaVu[ip] = 1;
        cnt++;
      }
    }
    else 
    {         
      if ( cnt > 0 )
      {
        indices.reAlloc(cnt);        
        size1 = cnt;
        size2 = _size-size1;
        field1.malloc(size1,_nfld);
        field2.malloc(size2,_nfld);

        cnt = 0;
        for(E_Int i = 0; i < size1; i++)
          for(E_Int eq = 1; eq <= _nfld; eq++)
            field1(i,eq) = _field(indices[i],eq);
        
        E_Int istart = indices[0];
        E_Int iend = indices[size1-1];
        
        if ( istart < iend )
        {
          for (E_Int i = iend+1; i < size; i++)
          {
            for(E_Int eq = 1; eq <= _nfld; eq++)   
              field2(cnt,eq) = _field(i,eq);
            cnt++;
          }
          for (E_Int i = 0; i < istart; i++)
          {
            for(E_Int eq = 1; eq <= _nfld; eq++)   
              field2(cnt,eq) = _field(i,eq);
            cnt++;  
          }
        }
        else 
        {
          for (E_Int i = iend+1; i < istart; i++)
          {
            for(E_Int eq = 1; eq <= _nfld; eq++)   
              field2(cnt,eq) = _field(i,eq);
            cnt++;
          }
        }
        return true;
      }
    }
  }
  return false;
}

//=============================================================================
/* Close pocket when a set of points are aligned*/
//=============================================================================
void Pocket::closeSpecificPocket(FldArrayF& field1, FldArrayF& field2,
                                 vector<TriangleZ*>& triangles)
{
  const E_Int imax1 = field1.getSize();
  const E_Int imax2 = field2.getSize();
  
  E_Int iA, iB, iC, iD;
  FldArrayF fieldA(_nfld);
  FldArrayF fieldB(_nfld);
  FldArrayF fieldC(_nfld);
  FldArrayF fieldD(_nfld);
 
  /*-----------------*/
  /* triangulation   */
  /*-----------------*/
  iA = 0;
  iB = iA + 1;
  iD = imax2-1;
  iC = iD - 1;
  
  while ( iA < imax1-1 && iD > 0 )
  {    
    E_Bool isOK = compDelaunay(iA, iB, iC, iD, 
                                  field1, field2, triangles);
     if ( isOK == false )
     {
       for (E_Int eq = 0; eq < _nfld; eq++)
       {
         fieldA[eq] = field1(iA, eq+1);
         fieldB[eq] = field1(iB, eq+1);
         fieldD[eq] = field2(iD, eq+1);
        }
        TriangleZ* t = new TriangleZ(fieldA, fieldB, fieldD);
        triangles.push_back(t);
        iA++;
        iB++;
      }
  }

  // end of seg1. connect remaining pts of seg2 to end pt of seg1
  if ( (iA == imax1-1)  && ( iD > 0) )
  {
    for (E_Int eq = 1; eq<=_nfld; eq++)
      fieldA[eq-1] = field1(iA,eq);

    while (iD > 0)
    {
      iC = iD-1;
      for (E_Int eq = 1; eq<=_nfld; eq++)
      {
        fieldD[eq-1] = field2(iD,eq);
        fieldC[eq-1] = field2(iC,eq);
      }
      TriangleZ* t  = new TriangleZ( fieldA, fieldD, fieldC);
      triangles.push_back(t);
      iD--;
    }
  }
  //end of seg2. connect remaining pts of seg1 to end pt of seg1  
  else if ( (iA < imax1-1) && (iD == 0) )
  {
    for (E_Int eq = 1; eq<=_nfld; eq++)
      fieldD[eq-1] = field2(iD, eq);
    while (iA < imax1-1)
    {
      iB = iA+1;
      for (E_Int eq = 1; eq<=_nfld; eq++)
      {
        fieldA[eq-1] = field1(iA,eq);
        fieldB[eq-1] = field1(iB,eq);
      }
      TriangleZ* t  = new TriangleZ( fieldA, fieldB, fieldD);
      triangles.push_back(t);
      iA++;
    } 
  }
}


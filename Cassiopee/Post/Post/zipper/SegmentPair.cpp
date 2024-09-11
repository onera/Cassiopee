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

# include "SegmentPair.h"
# include <stdlib.h>

using namespace K_CONST;
using namespace K_FUNC;
using namespace std;
using namespace K_FLD;

//=============================================================================
SegmentPair::SegmentPair(CString* st1, CString* st2)
{
  _indseg1 = st1->getIndArray();
  _indseg2 = st2->getIndArray(); 
  _blk1 = st1->getBlock();
  _blk2 = st2->getBlock();
  _origSt1 = st1->getOriginalString();
  _origSt2 = st2->getOriginalString();

}
//=============================================================================
SegmentPair::~SegmentPair()
{
}

//=============================================================================
FldArrayI& SegmentPair::getIndArray1()
{
  return _indseg1;
}
//=============================================================================
FldArrayI& SegmentPair::getIndArray2()
{
  return _indseg2;
}
//=============================================================================
StructBlock* SegmentPair::getBlock1()
{
  return _blk1;
}
//=============================================================================
StructBlock* SegmentPair::getBlock2()
{
  return _blk2;
}
//=============================================================================
/* Compute the zipper for  the current segment pair */
//=============================================================================
E_Boolean SegmentPair::computeZipper(const E_Int nfieldTot,
                                     FldArrayF& field, FldArrayI& connect)
{
  vector<TriangleZ*> triangles;
  const E_Int imax1 = _indseg1.getSize();
  const E_Int imax2 = _indseg2.getSize();
  
  E_Int iprev1, iprev2, inext1, inext2;
  E_Int smallestseg;
  
  /*-----------------*/
  /* triangulation   */
  /*-----------------*/
  iprev1 = 0;
  iprev2 = 0;
  inext1 = iprev1 + 1;
  inext2 = iprev2 + 1;
  
  // NB: iprev1 & co are modified in compDelaunay
  while ( (iprev1 < imax1-1) && (iprev2 < imax2-1) )  
    compDelaunay( nfieldTot, iprev1, inext1, iprev2, inext2,
                  triangles);   
  
  // end of seg1. connect remaining pts of seg2 to end pt of seg1
  
  if ( (iprev1 == imax1-1)  && ( iprev2 < imax2-1) )
  {
    smallestseg = 1;
    connectRemainingPtsToEnd( nfieldTot, smallestseg,
                              imax1-1, iprev2, imax2-1,
                              triangles);
  }
  
  // end of seg2. connect remaining pts of seg1 to end pt of seg1
  else if ( (iprev2 == imax2-1) && (iprev1 < imax1-1) )
  {
    smallestseg = 2;
    connectRemainingPtsToEnd( nfieldTot, smallestseg,
                              imax2-1, iprev1, imax1-1,
                              triangles);
  }

  /*----------------------------------------*/
  /* Compute the connectivity for triangles */
  /*----------------------------------------*/
  compConnectivity(nfieldTot, triangles, field, connect);
  for (size_t i = 0; i < triangles.size(); i++) delete triangles[i];
  return true;  
}

//=============================================================================
/* Computes the  distance between ind1 of seg1 and ind2 of seg2  */
//=============================================================================
E_Float SegmentPair::compSegmentPtsDistance(E_Int ind1, E_Int ind2)
{
  FldArrayF& coord1 = _blk1->getCoord();
  FldArrayF& coord2 = _blk2->getCoord();
  
  E_Float dx = coord2(ind2, 1) - coord1(ind1, 1);
  E_Float dy = coord2(ind2, 2) - coord1(ind1, 2);
  E_Float dz = coord2(ind2, 3) - coord1(ind1, 3);
  return  dx * dx + dy * dy + dz * dz;
}

//=============================================================================
/* Computes triangulation for the quad PiPi+1Pj+1Pj   */
//=============================================================================
void SegmentPair::compDelaunay(const E_Int nfieldTot,
                               E_Int& iprev, E_Int& inext,
                               E_Int& jprev, E_Int& jnext,
                               vector<TriangleZ*>& triangles)
{  
  // Assume the quad [PiPi+1Pj+1Pj]=ABCD 

  /* ====================================================*/
  /*                IMPORTANT                            */
  /*   constructeur triangle:                            */
  /*   ordre A, D et C ou B pour la connectivité         */
  /* ====================================================*/

  FldArrayF& coord1 = _blk1->getCoord();
  FldArrayF& coord2 = _blk2->getCoord();
  FldArrayF& cfdField1 = _blk1->getCfdField();
  FldArrayF& cfdField2 = _blk2->getCfdField();

  FldArrayF fieldA(nfieldTot);
  FldArrayF fieldB(nfieldTot);
  FldArrayF fieldC(nfieldTot);
  FldArrayF fieldD(nfieldTot);

  E_Int indA = _indseg1[iprev];
  E_Int indB = _indseg1[inext];
  E_Int indC = _indseg2[jnext];
  E_Int indD = _indseg2[jprev];
  E_Float diagAC, diagBD;

  for (E_Int n = 0; n < 3; n++)
  {
    fieldA[n] = coord1(indA, n+1);
    fieldB[n] = coord1(indB, n+1);
    fieldC[n] = coord2(indC, n+1);
    fieldD[n] = coord2(indD, n+1);
  }
  
  for (E_Int n = 3; n < nfieldTot ; n++ )
  {
    fieldA[n] = cfdField1(indA,n-2);
    fieldB[n] = cfdField1(indB,n-2);
    fieldC[n] = cfdField2(indC,n-2);
    fieldD[n] = cfdField2(indD,n-2);
  }
  
  E_Int isValid = checkTriangles(fieldA, fieldB, fieldC, fieldD);
  TriangleZ* t1;
  TriangleZ* t2;
  switch (isValid)
  {
    case 0: // both diags ar valid . keep the smaller
      //triangle ADC
      t1 = new TriangleZ(fieldA, fieldD, fieldC);
      //diagonal AC
      diagAC = compSegmentPtsDistance(indA, indC);
           
      //triangle ABD
      t2 = new TriangleZ(fieldA, fieldD, fieldB);
      //diagonal AD
      diagBD = compSegmentPtsDistance(indB, indD);
           
      //select the smaller diagonal
      if ( diagAC < diagBD ) // AC connection selected
      {
        triangles.push_back(t1); delete t2;
        //eliminates D for new triangulation research
        jprev = jnext;
        jnext = jprev + 1;
      }
      else// BD connection selected
      {
        triangles.push_back(t2); delete t1;
        // eliminate A for new triangulation research 
        iprev = inext;
        inext = inext + 1;
      }
      break;
    case 1: // only diag AC is valid
      t1 = new TriangleZ(fieldA, fieldD, fieldC);                      
      triangles.push_back(t1);
      //eliminates D for new triangulation research
      jprev = jnext;
      jnext = jnext + 1;
      break;
    case 2: // BD connection valid
      t2 = new TriangleZ(fieldA, fieldD, fieldB);              
      triangles.push_back(t2);
      // eliminate A for new triangulation research 
      iprev = inext;
      inext = inext + 1;
      break;
    default:
      printf("Error: checkTriangles must return values from 0 to 2.\n");
      exit(0);
  }
}

//=============================================================================
/* Says if the connection is valid  for triangulation
   Computes surface of ABCD by diagonals:
   S1 = AC^AB + AC^AD 
   S2 = AB^BD + BD^BC 
   if both surfaces are equal, AC and BD are both valid : return 0
   if S1 < S2 only diagonal AC is valid : return 1
   if S1 > S2 only diagonal BD is valid : return 2 
*/
//=============================================================================
E_Int SegmentPair::checkTriangles(FldArrayF& fieldA, FldArrayF& fieldB,
                                  FldArrayF& fieldC, FldArrayF& fieldD)
{
  /* test if the quad ABCD is valid (not convex)
     Sabc = norm(AB^AC)
     Sacd = norm(AC^AD)*/
  E_Int projP = 0;
  E_Float xA = fieldA[0];
  E_Float yA = fieldA[1];
  E_Float zA = fieldA[2];

  E_Float xB = fieldB[0];
  E_Float yB = fieldB[1];
  E_Float zB = fieldB[2];
  
  E_Float xC = fieldC[0];
  E_Float yC = fieldC[1];
  E_Float zC = fieldC[2];

  E_Float xD = fieldD[0];
  E_Float yD = fieldD[1];
  E_Float zD = fieldD[2];

  E_Float n1x, n1y, n1z;
  E_Float n2x, n2y, n2z; // normals of triangles
  E_Float norm1, norm2;
  E_Float x1, y1, z1;    
  E_Float x2, y2, z2;  
  E_Float x3, y3, z3;     
  E_Float sum1, sum2;

  E_Float dist, distmin;
  E_Float xpA, ypA, zpA;
  E_Float xpB, ypB, zpB;
  E_Float xpC, ypC, zpC;
  E_Float xpD, ypD, zpD;
  distmin = -1.;

  /*---------------------------------------------------------*/
  /* project points A, B, C, D on planes BCD, ACD, ABD, ABC
     select smallest distance between point and its projection
     projP = 1 if plane is BCD
     projP = 2 if plane is ACD 
     projP = 3 if plane is ABD
     projP = 4 if plane is ABC */
  /*---------------------------------------------------------*/
  // project point A on triangle BCD->A'
  if ( projectPointOnPlane(xA, yA, zA, 
                           xB, yB, zB, xC, yC, zC, xD, yD, zD,
                           xpA, ypA, zpA, dist) == true)
  {
    distmin = dist;
    projP = 1; 
  }

  // project point B on triangle ACD->B'
  if (projectPointOnPlane(xB, yB, zB,
                          xA, yA, zA, xC, yC, zC, xD, yD, zD,
                          xpB, ypB, zpB, dist) == true)
  {
    if (dist < distmin )
    {
      distmin = dist;
      projP = 2;
    }
  }
 
  // project point C on triangle ABD->C'
  if( projectPointOnPlane(xC, yC, zC,
                          xA, yA, zA, xB, yB, zB, xD, yD, zD,
                          xpC, ypC, zpC, dist) == true)
  {
    if (dist < distmin)
    {
      distmin = dist;
      projP = 3;
    }
  }
  
  // project point D on triangle ABC->D'
  if (projectPointOnPlane(xD, yD, zD,
                          xA, yA, zA, xB, yB, zB, xC, yC, zC, 
                          xpD, ypD, zpD, dist) == true)
  {
    if (dist < distmin)
    {
      distmin = dist;
      projP = 4;
    } 
  }

  // build quad A'BCD or AB'CD or ABC'D or ABCD'
  switch (projP)
  {
    case 1://plane BCD
      xA = xpA;
      yA = ypA;
      zA = zpA;
      break;
      
    case 2://plane ACD
      xB = xpB;
      yB = ypB;
      zB = zpB;
      break;
    case 3://plane ABD
      xC = xpC;
      yC = ypC;
      zC = zpC;
      break;
    case 4 ://plane ABC
      xD = xpD;
      yD = ypD;
      zD = zpD;
      break;
    default :
      printf("Error: in checkTriangles: not a valid value for projP.\n");
      exit(0);
  }
  //----------------------------------------------------------------------
  //assume now that A'BCD (or AB'CD...depending  on projP) is called ABCD 
  //test diagonal AC
  x1 = xB-xA; //vector AB
  y1 = yB-yA;
  z1 = zB-zA;
  
  x2 = xC-xA;//vector AC
  y2 = yC-yA; 
  z2 = zC-zA; 
  
  x3 = xD-xA;//vector AD
  y3 = yD-yA;
  z3 = zD-zA;
  
  n1x = y1 * z2 - y2 * z1;
  n1y = x2 * z1 - x1 * z2;
  n1z = x1 * y2 - x2 * y1;

  n2x = y2 * z3 - y3 * z2;
  n2y = x3 * z2 - x2 * z3;
  n2z = x2 * y3 - x3 * y2;
 
  norm1 = n1x*n1x + n1y*n1y + n1z*n1z;//AB^AC
  norm2 = n2x*n2x + n2y*n2y + n2z*n2z;//AC^AD
  sum1 = norm1+norm2;
  
  //test diagonal BD :
  x2 = xD-xB;//vector BD
  y2 = yD-yB;
  z2 = zD-zB;

  x3 = xC-xB;//vector BC
  y3 = yC-yB;
  z3 = zC-zB;
 
  n1x = y1 * z2 - y2 * z1;
  n1y = x2 * z1 - x1 * z2;
  n1z = x1 * y2 - x2 * y1;

  n2x = y2 * z3 - y3 * z2;
  n2y = x3 * z2 - x2 * z3;
  n2z = x2 * y3 - x3 * y2;
 
  norm1 = n1x*n1x + n1y*n1y + n1z*n1z;//AB^BD
  norm2 = n2x*n2x + n2y*n2y + n2z*n2z;//BD^BC
  sum2 = norm1+norm2;
  
  E_Float delta = E_abs(sum1-sum2);
  
  if (delta < E_GEOM_CUTOFF)
    return 0;
  else 
  {
    if ( sum1 < sum2)
      return 1;
    else 
      return 2;
  }
}

//============================================================================
/* Connect the remaining pts of segment of indices indseg from the index iflag
   with the increment inc,  with the extremum pt iend of the other segment
   and store them in the triangulation vector
*/
//============================================================================
void SegmentPair::connectRemainingPtsToEnd(const E_Int nfieldTot,
                                           E_Int smallestseg, E_Int iend,
                                           E_Int iflag, E_Int imax, 
                                           vector<TriangleZ*>& triangles)
{
  // seg2 remaining points must be connected to seg1 last point
  // iend : the last index of the segment seg1
  // iflag is the starting index of the remaining pts of seg2 
  // imax is the last index of seg2
  
  //smallestseg says which segment is the segment to which
  //remaining pts of the  other segments are connected
  FldArrayI seg1, seg2;
  E_Int indA, indB, indC;
  StructBlock* blk1;
  StructBlock* blk2;
  
  switch ( smallestseg )
  {
    case 1:
      blk1 = _blk1;
      blk2 = _blk2;
      seg1.malloc( _indseg1.getSize());
      seg2.malloc( _indseg2.getSize());
      seg1 = _indseg1;
      seg2 = _indseg2;
      break;
      
    case 2 :
      blk1 = _blk2;
      blk2 = _blk1;
      seg1.malloc( _indseg2.getSize());
      seg2.malloc( _indseg1.getSize());
      seg1 = _indseg2;
      seg2 = _indseg1;
      break;

    default:
      blk1 = NULL;
      blk2 = NULL;
      printf("Error: cannot connect remaining points.\n");
  }
  
  FldArrayF& coord1 = blk1->getCoord();
  FldArrayF& coord2 = blk2->getCoord();
  FldArrayF& cfdField1 = blk1->getCfdField();
  FldArrayF& cfdField2 = blk2->getCfdField();
  
  FldArrayF fieldA(nfieldTot);
  FldArrayF fieldB(nfieldTot);
  FldArrayF fieldC(nfieldTot);
  E_Int nbCfdField = nfieldTot-3;
  
  indA = seg1[iend];
  
  for (E_Int i = iflag; i < imax; i++)
  {
    indB = seg2[i];
    indC = seg2[i+1];
    for (E_Int n = 0; n < 3; n++)
    {
      fieldA[n] = coord1(indA, n+1);
      fieldB[n] = coord2(indB, n+1);
      fieldC[n] = coord2(indC, n+1);
    }
  
    for (E_Int nb = 1; nb <= nbCfdField; nb++ )
    {
      fieldA[nb+2] = cfdField1(indA,nb);
      fieldB[nb+2] = cfdField2(indB,nb);
      fieldC[nb+2] = cfdField2(indC,nb);
    }
    
    TriangleZ* t  = new TriangleZ(fieldA, fieldB, fieldC);
    triangles.push_back(t);
  }
}

//=============================================================================
/* Computes connectivity for triangles  */
//=============================================================================
void SegmentPair::compConnectivity( const E_Int nfieldTot,
                                    vector<TriangleZ*>& triangles,
                                    FldArrayF& field,
                                    FldArrayI& connect)
{
  E_Int nelts = triangles.size(); 
  connect.malloc(nelts,3);
  field.malloc(3*nelts, nfieldTot);
  
  for ( E_Int n = 0; n < nelts; n++ )
  {
    FldArrayF& field1 = triangles[n]->getField1();
    FldArrayF& field2 = triangles[n]->getField2();
    FldArrayF& field3 = triangles[n]->getField3();
    
    for (E_Int nfld = 1 ; nfld <= nfieldTot ; nfld++)
    {  
      field(3*n  , nfld) =  field1[nfld-1];
      field(3*n+1, nfld) =  field2[nfld-1];
      field(3*n+2, nfld) =  field3[nfld-1];
    }
  }
  
  for ( E_Int n = 0 ; n < nelts; n++ )
  {
    connect(n,1) = 3*n+1;
    connect(n,2) = 3*n+2;
    connect(n,3) = 3*n+3;
  }
}
// //=============================================================================
// /* Reverse the indseg arrays . Used when the triangulation is not
//    good in a direction, then reverse the direction 
//    NOT USED FOR THE MOMENT*/
// //=============================================================================
// void SegmentPair::reverseSegmentIndArray()
// {
//   FldArrayI indorig1 = _indseg1;
//   FldArrayI indorig2 = _indseg2;
//   E_Int size1 = indorig1.getSize() - 1;
//   E_Int size2 = indorig2.getSize() - 1;
  
//   for (E_Int i = 0 ; i <= size1;i++)
//     _indseg1[i] = indorig1[size1-i];

//   for (E_Int i = 0; i <= size2 ; i++)
//     _indseg2[i] = indorig2[size2-i];  
// }


//======================================================================
/* Modify the flag array (used to make pockets with remaining points) */
//======================================================================
void SegmentPair::updateFlagForStrings(vector<SegmentPair*>& segPairs)
{
  E_Int i11, i21, i12, i22;
  E_Int ind1, ind2;
  E_Int max1 = _indseg1.getSize();
  E_Int max2 = _indseg2.getSize();

  FldArrayIS& flag1 = _origSt1->getFlag();
  FldArrayIS& flag2 = _origSt2->getFlag();

  FldArrayI& indArray1 = _origSt1->getIndArray();
  FldArrayI& indArray2 = _origSt2->getIndArray();

  E_Boolean isFound;
  E_Int indArray1Size = indArray1.getSize();
  E_Int indArray2Size = indArray2.getSize();

  /*-----------------*/
  /* Premiers points */
  /*-----------------*/
  ind1 = _indseg1[0];
  ind2 = _indseg2[0];
 
  // Recherche du no correspondant a ind1 dans la string initiale
  E_Int i = 0;
 
  while ( indArray1[i] != ind1)
    i++;
  i11 = i;
 
  // Recherche du no correspondant a ind2 dans la string initiale
  i = 0;
  while ( indArray2[i] != ind2 )
    i++;
  i21 = i;

  /* mise a jour des flag:  si les deux extremites matchent avec deux
     extremites d'un meme segmentPair,  flag = 2. 
     Si ind1 matche avec une extremite d'une paire v1 et ind2 avec une 
     paire v2, flag = 1 pour ind1 et ind2 */
  if ( flag1[i11] == 0 || flag2[i21] == 0)
  {
    isFound = searchForCorrespondingSegments(ind1, ind2, segPairs);
    if ( isFound == true )
    {
      flag1[i11] = 2;
      flag2[i21] = 2;
    }
    else 
    {
      flag1[i11] = 1;
      flag2[i21] = 1;
    }
  }
  
  /*---------------*/
  /*derniers points*/
  /*---------------*/
  ind1 = _indseg1[max1-1];
  ind2 = _indseg2[max2-1];
  
  // Recherche du no correspondant a ind1 dans la string initiale
  i = 0;
  while ( indArray1[i] != ind1 && i < indArray1Size )
    i++;
  if ( indArray1[i] == ind1 )
    i12 = i;
  else
  {
    printf("SegmentPair: ind1 not found in string.\n");
    exit(0);
  }
  
  // Recherche du no correspondant a ind2 dans la string initiale
  i = 0;
  while ( indArray2[i] != ind2 && i < indArray2Size )
    i++;
  if ( indArray2[i] == ind2 )
    i22 = i;
  else
  {
    printf("SegmentPair: ind2 not found in string.\n");
    exit(0);
  }
 
  /* mise a jour des flag:  si les deux extremites matchent avec deux
     extremites d'un meme segmentPair,  flag = 2. 
     Si ind1 matche avec une extremite d'une paire v1 et ind2 avec une 
     paire v2, flag = 1 pour ind1 et ind2 */
  if (flag1[i12] == 0 || flag2[i22] == 0 )
  {
    isFound = searchForCorrespondingSegments(ind1, ind2, segPairs);
  
    if ( isFound == true )
    {
      flag1[i12] = 2;
      flag2[i22] = 2;
    }
    else 
    {
      flag1[i12] = 1;
      flag2[i22] = 1;
    }
  }

  /*-------------------*/
  /* points interieurs */
  /*-------------------*/
  
  //2 points : a singlesegment == a segpair element may occur
  //flag it
  if ( max1 == 2 )
  {
    if (flag1[i11] == 1 && flag1[i12] == 1) // min et max du segment 1
    {
      flag1[i11] = -1;
      flag1[i12] = -1;
    }
  }
  else 
  {
    for (E_Int k = 1 ; k < max1-1; k++)
    {
      for (E_Int j = 0; j < indArray1Size; j++)
      {
        if ( indArray1[j] == _indseg1[k] )
          flag1[j] = 2;
      }
    }
  }

  if ( max2 == 2 )
  {
    if (flag2[i21] == 1 && flag2[i22] == 1)
    {
      flag2[i21] = -1;
      flag2[i22] = -1;
    }
  }
  else 
  {
    for (E_Int k = 1 ; k < max2-1; k++)
    {
      for (E_Int j = 0; j < indArray2Size; j++)
      {
        if ( indArray2[j] == _indseg2[k] )
          flag2[j] = 2;
      }
    }
  }
}

//=============================================================================
/* Given extremities of a segmentPair, search a corresponding segment pair
   ind1 comes from segment1, ind2 comes from segment2.
   Return true if matching segments for ind1 and ind2 derive from the same
   pair of segments 
*/
//=============================================================================
E_Boolean 
SegmentPair::searchForCorrespondingSegments(E_Int ind1, E_Int ind2, 
                                            vector<SegmentPair*>& segPairs)
{  
  E_Float matchTol = _blk1->getMatchTol();
  FldArrayF& coord1 = _blk1->getCoord();
  FldArrayF& coord2 = _blk2->getCoord();
  E_Float x1 = coord1(ind1, 1);
  E_Float y1 = coord1(ind1, 2);
  E_Float z1 = coord1(ind1, 3);
  
  E_Float x2 = coord2(ind2, 1);
  E_Float y2 = coord2(ind2, 2);
  E_Float z2 = coord2(ind2, 3);
  
  E_Float dx, dy, dz;
  E_Int ind31, ind32, ind41, ind42;
  E_Int segPairsSize =  segPairs.size();
  for (E_Int v = 0; v < segPairsSize; v++)
  {
    if ( segPairs[v] != this)
    {
      FldArrayI& indseg3 = segPairs[v]->getIndArray1();
      FldArrayI& indseg4 = segPairs[v]->getIndArray2();
      E_Int max3 = indseg3.getSize();
      E_Int max4 = indseg4.getSize();
      FldArrayF& coord3 = segPairs[v]->getBlock1()->getCoord();
      FldArrayF& coord4 = segPairs[v]->getBlock2()->getCoord();
      
      ind31 = indseg3[0];
      ind41 = indseg4[0];
      ind32 = indseg3[max3-1];
      ind42 = indseg4[max4-1];
 
      /*------------------------------*/
      /* test ind1 et premier segment */
      /*------------------------------*/

      // test ind1 et indseg3[0]
      dx = coord3(ind31, 1) - x1;
      dy = coord3(ind31, 2) - y1;
      dz = coord3(ind31, 3) - z1;
      dx = E_abs(dx);
      dy = E_abs(dy);
      dz = E_abs(dz);
      if ( dx < matchTol && dy < matchTol && dz < matchTol )
      {
        dx = coord4(ind41, 1) - x2;
        dy = coord4(ind41, 2) - y2;
        dz = coord4(ind41, 3) - z2;
        dx = E_abs(dx);
        dy = E_abs(dy);
        dz = E_abs(dz);
            
        if ( dx < matchTol && dy < matchTol && dz < matchTol )
          return true;
        else return false;
      }
    
      // test ind1 et indseg3[max3-1]
      dx = coord3(ind32, 1) - x1;
      dy = coord3(ind32, 2) - y1;
      dz = coord3(ind32, 3) - z1;
      dx = E_abs(dx);
      dy = E_abs(dy);
      dz = E_abs(dz);
      
      if ( dx < matchTol && dy < matchTol && dz < matchTol )
      {
        dx = coord4(ind42, 1) - x2;
        dy = coord4(ind42, 2) - y2;
        dz = coord4(ind42, 3) - z2;
        dx = E_abs(dx);
        dy = E_abs(dy);
        dz = E_abs(dz);
        
        if ( dx < matchTol && dy < matchTol && dz < matchTol )
          return true;
        else return false;
      }
      
      // test ind1 et indseg4[0]
      dx = coord4(ind41, 1) - x1;
      dy = coord4(ind41, 2) - y1;
      dz = coord4(ind41, 3) - z1;
      dx = E_abs(dx);
      dy = E_abs(dy);
      dz = E_abs(dz);
      if ( dx < matchTol && dy < matchTol && dz < matchTol )
      {
        dx = coord3(ind31, 1) - x2;
        dy = coord3(ind31, 2) - y2;
        dz = coord3(ind31, 3) - z2;
        dx = E_abs(dx);
        dy = E_abs(dy);
        dz = E_abs(dz);
            
        if ( dx < matchTol && dy < matchTol && dz < matchTol )
          return true;
        else return false;
      }
    
      //test ind1 et indseg4[max4-1]
      dx = coord4(ind42, 1) - x1;
      dy = coord4(ind42, 2) - y1;
      dz = coord4(ind42, 3) - z1;
      dx = E_abs(dx);
      dy = E_abs(dy);
      dz = E_abs(dz);
      if ( dx < matchTol && dy < matchTol && dz < matchTol )
      {
        dx = coord3(ind32, 1) - x2;
        dy = coord3(ind32, 2) - y2;
        dz = coord3(ind32, 3) - z2;
        dx = E_abs(dx);
        dy = E_abs(dy);
        dz = E_abs(dz);
        
        if ( dx < matchTol && dy < matchTol && dz < matchTol )
          return true;
        else return false;
      }

      /*------------------------------*/
      /* test ind2 et premier segment */
      /*------------------------------*/
      // test ind2 et ind31
      dx = coord3(ind31, 1) - x2;
      dy = coord3(ind31, 2) - y2;
      dz = coord3(ind31, 3) - z2;
      dx = E_abs(dx);
      dy = E_abs(dy);
      dz = E_abs(dz);
      if ( dx < matchTol && dy < matchTol && dz < matchTol )
      {
        dx = coord4(ind41, 1) - x1;
        dy = coord4(ind41, 2) - y1;
        dz = coord4(ind41, 3) - z1;
        dx = E_abs(dx);
        dy = E_abs(dy);
        dz = E_abs(dz);
        
        if ( dx < matchTol && dy < matchTol && dz < matchTol )
          return true;
        else return false;
      }
     
      // test ind2 et ind3
      dx = coord3(ind32, 1) - x2;
      dy = coord3(ind32, 2) - y2;
      dz = coord3(ind32, 3) - z2;
      dx = E_abs(dx);
      dy = E_abs(dy);
      dz = E_abs(dz);
      if ( dx < matchTol && dy < matchTol && dz < matchTol )
      {
        dx = coord4(ind42, 1) - x1;
        dy = coord4(ind42, 2) - y1;
        dz = coord4(ind42, 3) - z1;
        dx = E_abs(dx);
        dy = E_abs(dy);
        dz = E_abs(dz);
        
        if ( dx < matchTol && dy < matchTol && dz < matchTol )
          return true;
        else return false;
      }
      
      // test ind2 et indseg4[0]
      dx = coord4(ind41, 1) - x2;
      dy = coord4(ind41, 2) - y2;
      dz = coord4(ind41, 3) - z2;
      dx = E_abs(dx);
      dy = E_abs(dy);
      dz = E_abs(dz);
      if ( dx < matchTol && dy < matchTol && dz < matchTol )
      {
        dx = coord3(ind31, 1) - x1;
        dy = coord3(ind31, 2) - y1;
        dz = coord3(ind31, 3) - z1;
        dx = E_abs(dx);
        dy = E_abs(dy);
        dz = E_abs(dz);
            
        if ( dx < matchTol && dy < matchTol && dz < matchTol )
          return true;
        else return false;
      }
    
      //test ind1 et indseg4[max4-1]
      dx = coord4(ind42, 1) - x2;
      dy = coord4(ind42, 2) - y2;
      dz = coord4(ind42, 3) - z2;
      dx = E_abs(dx);
      dy = E_abs(dy);
      dz = E_abs(dz);
      if ( dx < matchTol && dy < matchTol && dz < matchTol )
      {
        dx = coord3(ind32, 1) - x1;
        dy = coord3(ind32, 2) - y1;
        dz = coord3(ind32, 3) - z1;
        dx = E_abs(dx);
        dy = E_abs(dy);
        dz = E_abs(dz);
        
        if ( dx < matchTol && dy < matchTol && dz < matchTol )
          return true;
        else return false;
      }
    }
  }
  return false;
}
//=========================================================================
/* Add extremities of segment pairs to the list of single segments to 
   be treated */
//========================================================================
void SegmentPair::identifyLastEdges(vector<FldArrayF*>& lastEdges)
{ 
  E_Int ind1 = _indseg1[0];
  E_Int ind2 = _indseg2[0];
 
  FldArrayF& field1 = _blk1->getGlobalField();
  FldArrayF& field2 = _blk2->getGlobalField();
  E_Int nfld = field1.getNfld(); 
  E_Int max1 = _indseg1.getSize()-1;
  E_Int max2 = _indseg2.getSize()-1;

  FldArrayF* field = new FldArrayF(2,nfld);
  // Premiers points des deux segments
  for (E_Int eq = 1; eq <= nfld; eq++)
  {
    (*field)(0,eq) = field1(ind1,eq);
    (*field)(1,eq) = field2(ind2,eq);
  } 

  lastEdges.push_back(field);
 
  //Derniers points des deux segments
  field = new FldArrayF(2,nfld);
  ind1 = _indseg1[max1];
  ind2 = _indseg2[max2];
  for (E_Int eq = 1; eq <= nfld; eq++)
  {
    (*field)(0,eq) = field1(ind1,eq);
    (*field)(1,eq) = field2(ind2,eq);
  }
  lastEdges.push_back(field); 
}


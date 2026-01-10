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

# include "ZipLib.h"
# include <stdio.h>
# include <stdlib.h>

using namespace K_FUNC;
using namespace K_FLD;
using namespace std;
using namespace K_CONST;

//=============================================================================
/* Test if edge1 (coord1) and edge2 (coord2) coming from different
   triangles are matching */
//=============================================================================
E_Bool testIfEdgesAreMatching(E_Float matchTol, 
                                 FldArrayF& field1, FldArrayF& field2)
{
  // premier cote 
  E_Float x11 = field1(0,1);
  E_Float y11 = field1(0,2);
  E_Float z11 = field1(0,3);

  E_Float x12 = field1(1,1);
  E_Float y12 = field1(1,2);
  E_Float z12 = field1(1,3);

  // second cote
  E_Float x21 = field2(0,1);
  E_Float y21 = field2(0,2);
  E_Float z21 = field2(0,3);

  E_Float x22 = field2(1,1);
  E_Float y22 = field2(1,2);
  E_Float z22 = field2(1,3);

  E_Float dx1 = E_abs(x21-x11);
  E_Float dy1 = E_abs(y21-y11);
  E_Float dz1 = E_abs(z21-z11);

  E_Float dx2 = E_abs(x22-x12);
  E_Float dy2 = E_abs(y22-y12);
  E_Float dz2 = E_abs(z22-z12);
 
  if (dx1 < matchTol && dy1 < matchTol && dz1 < matchTol)
  {
    if ( dx2 < matchTol && dy2 < matchTol && dz2 < matchTol)
      return true;
    else return false;
  }

  // 2 eme test 
  dx1 = E_abs(x22-x11);
  dy1 = E_abs(y22-y11);
  dz1 = E_abs(z22-z11);
  dx2 = E_abs(x21-x12);
  dy2 = E_abs(y21-y12);
  dz2 = E_abs(z21-z12); 
   
  if (dx1 < matchTol && dy1 < matchTol && dz1 < matchTol )
  {
    if (dx2 < matchTol && dy2 < matchTol && dz2 < matchTol)
      return true;
    else return false;
  }
  return false;
}

//============================================================================
/*Test if pt1 (field1(i1,:)) and one extremity of field2 coming from different
  singleSegments are matching. If true return i2. if not found return -1  */
//============================================================================
E_Int testIfPointsAreMatching(E_Float matchTol, E_Int i1, 
                              FldArrayF& field1, FldArrayF& field2)
{
  E_Int max2 = field2.getSize()-1;
  // premier point du segment 1
  E_Float x1 = field1(i1,1);
  E_Float y1 = field1(i1,2);
  E_Float z1 = field1(i1,3);

  // extremites du segment 2
  E_Float x21 = field2(0,1);
  E_Float y21 = field2(0,2);
  E_Float z21 = field2(0,3);

  E_Float x22 = field2(max2,1);
  E_Float y22 = field2(max2,2);
  E_Float z22 = field2(max2,3);

  E_Float dx1 = E_abs(x21-x1);
  E_Float dy1 = E_abs(y21-y1);
  E_Float dz1 = E_abs(z21-z1);

  E_Float dx2 = E_abs(x22-x1);
  E_Float dy2 = E_abs(y22-y1);
  E_Float dz2 = E_abs(z22-z1);

  //test if i1 and first element of seg2 are matching
  if (dx1 < matchTol && dy1 < matchTol && dz1 < matchTol)
    return 0;
  
  //test if i1 and last element of seg2 are matching
  if (dx2 < matchTol && dy2 < matchTol && dz2 < matchTol )
    return max2;

  // not found
  return -1;
}

//==========================================================================
/* Write the single segments to be pocketted */
//==========================================================================
void writeSingleSegments(vector<SingleSegment*>& singleSegments)
{
  //E_Bool add = false;
  E_Int ssSize = singleSegments.size();
  for (E_Int v = 0; v < ssSize; v++)
  {
    FldArrayF& field = singleSegments[v]->getGlobalField();
    E_Int np = field.getSize();
    FldArrayI connect( np-1,3);
  
    E_Int c = 0;
    for (E_Int i = 0; i < np-1; i++)
    {     
      connect(c, 1) = i+1;
      connect(c, 2) = i+2;
      connect(c, 3) = i+2;
      c++;
    }
    if (c > 0)
    {
      //K_IO::GenIO::getInstance()->tpwriteTriangles("strings2.tp", field, 
      //                                             connect, add);
      //add = true;
    }
  }
}

//==========================================================================
/* Build pockets */
//==========================================================================
void buildPocket(E_Float matchTol,
                 vector<SingleSegment*>& singleSegments,
                 vector<Pocket*>& pockets)
{
  E_Int size = singleSegments.size();
  if (size == 0) return;
  
  E_Int cnt;
  E_Int nfld = (singleSegments[0])->getGlobalField().getNfld();

  E_Int max = 0;
  for (E_Int v = 0; v < size; v++)
    max += singleSegments[v]->getGlobalField().getSize();

  FldArrayIS dejaVu(size);
  dejaVu.setAllValuesAtNull();
  FldArrayIS tag(size);

  E_Bool closed;
  E_Int iprev, inext, iend;
  E_Int vprev, vend;
  E_Int v2;

  for (E_Int v1 = 0; v1 < size; v1++)
  { 
    if (dejaVu[v1] == 0)
    {
      cnt = 0;
      FldArrayF* globalTab = new FldArrayF(max, nfld);
      closed = false;
      vprev = v1;
      vend = v1;
      iprev = 0;
      iend = 0;
      tag.setAllValuesAt(0);
      
      while (closed == false)
      {
        FldArrayF& field1 = singleSegments[vprev]->getGlobalField();
        
        for (v2 = 0; v2 < size; v2++)
        {
          if (dejaVu[v2] == 0 && v2 != vprev)
          {  
            FldArrayF& field2 = singleSegments[v2]->getGlobalField(); 
            E_Int max2 = field2.getSize()-1;
            inext = testIfPointsAreMatching(matchTol, iprev, field1, field2);

            if (inext != -1)
            {                
              vprev = v2;
              if ( inext == 0)
                iprev = max2;
              else iprev = 0;
              if (cnt + max2 >= globalTab->getSize())
              { // no room for storing fields in global segment
                delete globalTab;
                goto nextPocket;
              }
              else
                addSegmentInGlobalTab(inext, field2, *globalTab, cnt);
              tag[vprev] = 1;
       
              if ( iprev == iend && vprev == vend)
              {                
                eraseDoublePts(cnt, nfld, *globalTab);
                Pocket* pocket = new Pocket(*globalTab);
                pockets.push_back(pocket);
                globalTab->reAlloc(0);
                delete globalTab;
                closed = true;
                for (E_Int v3 = 0; v3 < size; v3++)
                {
                  if (tag[v3] == 1) 
                  {
                    dejaVu[v3] = 1;
                    break;
                  }
                }
              }
              goto next;
            }
          }
        }
        next:;

        // if closing is impossible: do not take it into account
        // for the moment
      
        if (v2 == size)
        {
          globalTab->reAlloc(0);
          delete globalTab;
          goto nextPocket; 
        }
      }
    }
    nextPocket:;
  }
 
  for (E_Int i = 0; i < dejaVu.getSize();i++)
  {
    if (dejaVu[i] == 0)
      printf("Warning: buildPocket: segment number " SF_D_ " is not in a pocket.\n", i+1);  
  }
  for (E_Int i = 0; i < size; i++)
    delete singleSegments[i];
  singleSegments.clear();
}

//=============================================================================
/* Add segment of field, starting from istart in listOfPockets */
//=============================================================================
void addSegmentInGlobalTab(E_Int istart, FldArrayF& field,
                           FldArrayF& globalTab, E_Int& cnt)
{
  E_Int nfld = field.getNfld();
  E_Int max = field.getSize()-1;

  if (istart == 0)
  {
    for (E_Int i = 0; i <= max; i++)
    { 
      for (E_Int eq = 1; eq <= nfld; eq++)
        globalTab(cnt, eq) = field(i, eq);
      cnt++;
    }
  }
  else if (istart == max)
  {
    for (E_Int i = 0; i <= max; i++)
    {
      for (E_Int eq = 1; eq <= nfld; eq++)
        globalTab(cnt, eq) = field(max-i, eq);
      cnt++;
    }
  }
  else 
  {
    printf("Error : addSegmentInGlobalTab : not a valid value for istart.");
    printf("Must be 0 or max.\n");
    exit(0);
  }
}

//====================================================================
/* Given an array of size sizeIni*nfld erase mutiply defined pts
   and resize the array*/
//====================================================================
void eraseDoublePts(E_Int sizeIni, E_Int nfld,
                    FldArrayF& globalTab)
{
  E_Float dx, dy, dz;
  FldArrayF work(sizeIni, nfld);

  E_Int c = 0;
  for (E_Int i = 0; i < sizeIni; i++)
  {
    E_Bool found = false;
    for (E_Int j = 0; j < i; j++)
    {
      dx = globalTab(i,1) - globalTab(j,1);
      dy = globalTab(i,2) - globalTab(j,2);
      dz = globalTab(i,3) - globalTab(j,3);
      
      if ( fEqualZero(dx) && fEqualZero(dy) && fEqualZero(dz))
      {
        found = true;
        break;
      }
    }
    if (found == false)
    {
      for (E_Int eq = 1; eq <= nfld; eq++)
        work(c,eq) = globalTab(i,eq);
      c++;
    }
  }
  work.reAllocMat(c,nfld);
  globalTab = work;
}


//=========================================================================
/* project point P on plane ABC. Return the coordinates of the projection H
   If projection is not possible return False
   IN : (x,y,z) : coordinates of point to be projected
   IN : (x0,y0,z0),(x1,y1,z1) (x2,y2,z2) coordinates of 3 points of the plane
   OUT : xp, yp, zp coordinates of the projection of (x,y,z) on plane (ABC)
   OUT : distProj : distance between (x,y,z) and (xp,yp,zp)
*/
//=========================================================================
E_Bool projectPointOnPlane(E_Float x, E_Float y, E_Float z,
                              E_Float x0, E_Float y0, E_Float z0,
                              E_Float x1, E_Float y1, E_Float z1,
                              E_Float x2, E_Float y2, E_Float z2,
                              E_Float& xp, E_Float& yp, E_Float& zp,
                              E_Float& distProj)
{
  E_Float e0x, e0y, e0z, e1x, e1y, e1z; 
  E_Float lTx, lTy, lTz;
  E_Float nx, ny, nz, t, norm;

  // Normale au triangle ABC
  e0x = (x1-x0);
  e0y = (y1-y0);
  e0z = (z1-z0);
  
  e1x = (x2-x0);
  e1y = (y2-y0);
  e1z = (z2-z0);
  
  nx = e0y*e1z-e0z*e1y;
  ny = e0z*e1x-e0x*e1z;
  nz = e0x*e1y-e0y*e1x;
  
  t = nx*(x0-x)+ny*(y0-y)+nz*(z0-z);
  norm = nx*nx+ny*ny+nz*nz;
  norm = E_max(norm,E_MIN_FLOAT);
  t = t/norm;
  
  lTx = x+t*nx;
  lTy = y+t*ny;
  lTz = z+t*nz;
  
  distProj = (lTx-x)*(lTx-x)+(lTy-y)*(lTy-y)+(lTz-z)*(lTz-z);
  
  xp = lTx;
  yp = lTy;
  zp = lTz;

  return true;
}

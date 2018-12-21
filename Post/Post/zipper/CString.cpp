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

# include "CString.h"
# include <vector>
# include <stdlib.h>

using namespace std;
using namespace K_FUNC;
using namespace K_FLD;
using namespace K_CONST;

//=============================================================================
CString::CString(FldArrayI& ind, StructBlock* blk)
{
  _blk = blk;
  _ind = ind;
  _dejaVu.malloc(_ind.getSize());
  _dejaVu.setAllValuesAt(0);
  _matchInfo.malloc(2,2);
  _matchInfo.setAllValuesAt(-1);
  _flag.malloc(_ind.getSize());
  _flag.setAllValuesAt(0);
}
//=============================================================================
/* Build a CString with elements from istart to iend of another string of the
   same block.
   Details : istart is NOT necessary smaller than iend here */
//=============================================================================
CString::CString(CString* st, E_Int istart, E_Int iend)
{
  _origString = st;
  _blk = st->getBlock();
  _matchInfo.malloc(2,2);
  _matchInfo.setAllValuesAt(-1);
  _flag.malloc(_ind.getSize());
  _flag.setAllValuesAt(0);
 
  FldArrayI& indstg = st->getIndArray();
  E_Int isize = E_abs(iend-istart)+1;
  _ind.malloc(isize);
  _dejaVu.malloc(isize);
  _dejaVu.setAllValuesAt(0);

  if (istart <= iend )
  {
    for (E_Int i =0; i < isize;i++)
      _ind[i] = indstg[istart+i];
  }
  else 
  {
    for (E_Int i = 0; i < isize; i++)
      _ind[i] = indstg[istart-i];
  }
}

//=============================================================================
CString::~CString()
{
}

//=============================================================================
/* Get the dejaVu array    */
//=============================================================================
FldArrayIS& CString::getDejaVu()
{
  return _dejaVu;
}
//=============================================================================
/* Get flag array    */
//=============================================================================
FldArrayIS& CString::getFlag()
{
  return _flag;
}

//=============================================================================
/* Return matchInfo array*/
//=============================================================================
FldArrayI& CString::getMatchInfo()
{
  return _matchInfo;
}
//=============================================================================
/* Reset deja vu to false */
//=============================================================================
void CString::resetDejaVu()
{
  _dejaVu.setAllValuesAt(0);
}

//=============================================================================
/* Get the indices array of points belonging to the cstring  */
//=============================================================================
FldArrayI& CString::getIndArray()
{
  return _ind;
}

//=============================================================================
/* Get the original string from which derives the current string */
//=============================================================================
CString*& CString::getOriginalString()
{
  return _origString;
}
//=============================================================================
/* Get the block to which belongs the cstring */
//=============================================================================
StructBlock* CString::getBlock()
{
  return _blk;
}

//=============================================================================
/* Return the list of strings that come from potentially matching strings,
   i.e. strings that are attached to a block that can be an overlap block
   of the current block */
//=============================================================================
vector<CString*>& CString::getListOfMatchingStrings()
{
  return _listOfMatchingStrings;
}

//=============================================================================
/* Find the first non deja vu point in string.
 Return -1 if no deja vu point exists in string. */
//=============================================================================
E_Int CString::findFirstNonDejaVu()
{
  E_Int ind = 0;
  while ( ind < _ind.getSize() )
  {
    if ( _dejaVu[ind] > 0  || _dejaVu[ind] == -1)
      ind++; 
    else break;
  }
  if (ind == _ind.getSize())
    return -1;
  else
    return ind;
}
//=============================================================================
/* Find the first non deja vu point in string.
 Return -1 if no deja vu point exists in string. */
//=============================================================================
E_Int CString::findFirstNonFlagged()
{
  E_Int ind = 0;
  while ( ind < _ind.getSize() )
  {
    if ( _flag[ind] == 2 )  
      ind++; 
    else break;
  }
  if (ind == _ind.getSize())
    return -1;
  else
    return ind;
}
//=============================================================================
/* For a string, search for matching segments of strings and order them in
   strOut1 and strOut2.
   For each point, find the nearest string and point.
   Go on till changing of string or going through the last point
   (for circular strings).
*/
//=============================================================================
void CString::searchForMatchingStrings( vector<CString*>& strings,
                                        vector<CString*>& strOut1,
                                        vector<CString*>& strOut2)
{
  CString* str2;
  CString* newStr = NULL;
  E_Int i2;
  E_Int fs1, fe1;
  E_Int fs2, fe2;
 
  E_Int ind = findFirstNonDejaVu();

  E_Boolean found;
  E_Int is1, ie1, js1, je1;
  E_Int is2, ie2, js2, je2;
  E_Boolean isSelected1, isSelected2;

  while ( ind != -1 )
  { 
    vector<CString*> stringsL = _listOfMatchingStrings;
    if (stringsL.size() == 0) return;

    found = this->nearestPoint( ind, stringsL, i2, str2 );

    if (found == false)
    {      
      printf("Error: searchForMatchingString: failed in %d.\n", ind);
      exit(0);
    }
   
    // Try to match this to str2
    is1 = ind;
    is2 = i2;
    matchStrings(this, str2, stringsL, is1, ie1, is2, ie2);

    // a matching string is found. Reverse
    //------------------------------------
    vector<CString*> stringsL2 = str2->getListOfMatchingStrings();
    js2 = -1;
    js1 = -1;
    matchStrings(str2, this, stringsL2, js2, je2, js1, je1);

    if ( is1 != -1 && is2 != -1 && js2 != -1 && je2 != -1)
    {      
      isSelected1 = areSegmentsOverlapping(is1, ie1, js1, je1);
      isSelected2 = areSegmentsOverlapping(js2, je2, is2, ie2);

      if (isSelected1 == true && isSelected2 == true )
        selectClosestSegment( this, str2, 
                              is1, js1, ie1, je1,
                              is2, js2, ie2, je2,
                              fs1, fe1, fs2, fe2);
      else // both segments dont overlap
      {
        if (fs1 == fe1 && fs2 == fe2 )
        {
          _dejaVu[ind] = -1;
          goto end;
        }
        else 
        {         
          // both segments do not overlap:
          // keep the selected segment at first match (no reverse)
          
          if ( is2 != ie2)
          {
            fs2 = is2;
            fe2 = ie2;
            // nearest point from is2 and ie2
            fs1 = findClosestPoint(str2, this, fs2, is1, ie1);
            fe1 = findClosestPoint(str2, this, fe2, is1, ie1); 
          }
          else
          {
            _dejaVu[ind] = -1;
            goto end;
          }
        }
      }
    }
    else 
    {
      _dejaVu[ind] = -1;
      goto end;
    }
      
    newStr = new CString(this, fs1, fe1);
    strOut1.push_back(newStr);
    newStr = new CString(str2, fs2, fe2);
    strOut2.push_back(newStr);
    setDejaVu(this, fs1, fe1, str2, fs2, fe2);
    end:;
    ind = findFirstNonDejaVu();
  }
}
// //=============================================================================
// /* For a string, search for matching segments of strings and order them in
//    strOut1 and strOut2.
//    For each point, find the nearest string and point.
//    Go on till changing of string or going through the last point
//    (for circular strings).
// */
// //=============================================================================
// void CString::searchForMatchingStrings( vector<CString*>& strings,
//                                         vector<CString*>& strOut1,
//                                         vector<CString*>& strOut2)
// {
//   CString* str2;
//   CString* newStr = NULL;
//   E_Int i2;
//   E_Int fs1, fe1;
//   E_Int fs2, fe2;
 
//   E_Int ind = findFirstNonDejaVu();

//   E_Boolean found;
//   E_Int is1, ie1, js1, je1;
//   E_Int is2, ie2, js2, je2;
//   E_Boolean isSelected1, isSelected2;
  
//   while ( ind != -1 )
//   { 
//     vector<CString*> stringsL = _listOfMatchingStrings;
//     if (stringsL.size() == 0) return;

//     found = this->nearestPoint( ind, stringsL, i2, str2 );

//     if (found == false)
//     {      
//       cout << "Error : searchForMatchingString : failed in "<< ind << endl;
//       exit(0);
//     }
    
//     //dbx
//     E_Int s1, s2;
//     for (s1 = 0; s1 < strings.size(); s1++)
//     {
//       if (strings[s1] == this ) 
//         break;
//     }
//     for (s2= 0; s2 < strings.size(); s2++)
//     {
//       if (strings[s2] == str2 ) 
//         break;
//     } 
//     cout << "string "<< s1+1 << " matches with " << s2+1<<endl;
//     //end dbx

//     // Try to match this to str2
//     is1 = ind;
//     is2 = i2;
//     matchStrings(this, str2, stringsL, is1, ie1, is2, ie2);

 
//     // a matching string is found. Reverse
//     //------------------------------------
//     vector<CString*> stringsL2 = str2->getListOfMatchingStrings();
//     js2 = is2;//-1;
//     js1 = -1;
//     matchStrings(str2, this, stringsL2, js2, je2, js1, je1);

//     if ( is1 != -1 && is2 != -1 && js2 != -1 && je2 != -1)
//     {      
//       isSelected1 = areSegmentsOverlapping(is1, ie1, js1, je1);
//       isSelected2 = areSegmentsOverlapping(js2, je2, is2, ie2);

//       cout <<"isSelect:"<<isSelected1 << " " << isSelected2<< endl;
//       if (isSelected1 == true && isSelected2 == true )
//       {
//         selectClosestSegment( this, str2, 
//                               is1, js1, ie1, je1,
//                               is2, js2, ie2, je2,
//                               fs1, fe1, fs2, fe2);
//         if (fs1 == fe1 && fs2 == fe2 )
//         {
//           _dejaVu[ind] = -1;
//           goto end;
//         }
//       }
//       else 
//       {         
//         // both segments do not overlap:
//         // keep the selected segment at first match (no reverse)
//         if ( is2 != ie2)
//         {
//           fs1 = is1;
//           fe1 = ie1;
//           fs2 = is2;
//           fe2 = ie2;
//         }
//         else
//         {
//           _dejaVu[ind] = -1;
//           goto end;
//         }
//       }
//     }
//     else 
//     {
//       _dejaVu[ind] = -1;
//       goto end;
//     }

//     newStr = new CString(this, fs1, fe1);
//     strOut1.push_back(newStr);
//     newStr = new CString(str2, fs2, fe2);
//     strOut2.push_back(newStr);
//     setDejaVu(this, fs1, fe1, str2, fs2, fe2);
//     end:;
//     ind = findFirstNonDejaVu();
//   }  
 
// }

//=============================================================================
/*
  Given two strings s1 and s2.
  Return the index that matches s1 on s2.
  Can return -1 : no match possible.
*/
//=============================================================================
void CString::matchStrings(CString* s1, CString* s2,
                           vector<CString*>& strings,
                           E_Int& istart1, E_Int& iend1,
                           E_Int& istart2, E_Int& iend2)
{
  E_Boolean found;
  E_Int i2;
  CString* str;

  FldArrayIS& dejaVu1 = s1->getDejaVu();
  FldArrayI& ind1 = s1->getIndArray();
  E_Int beg;
  E_Int ie1, ie2, is1, is2;
  iend1 = -1;
  iend2 = -1;

  if (istart1 == -1 )//&& istart2 == -1)
    beg = 0;
  else 
    beg = istart1;
   
  for (E_Int i1 = beg; i1 < ind1.getSize(); i1++)
  {
    if (dejaVu1[i1] < 2)
    {
      found = s1->nearestPoint(i1, strings, i2, str);
      if ( found == false ) 
      {
	printf("Warning: matchStrings: no nearest point found.\n"); 
	return;
      }
      
      if ( str == s2 )
      {
        is1 = i1;
        is2 = i2;
        E_Boolean isOK = compMatching( s1, s2, strings,
                                       is1, is2, ie1, ie2);
        
        if (isOK == true )
        {
          istart1 = is1;
          istart2 = is2;
          iend1 = ie1;
          iend2 = ie2;
          return;
        }
        else 
        {
          istart1 = -1;
          istart2 = -1;
          iend1 = -1;
          iend2 = -1;
        }
      }
    }
  }
}

//=============================================================================
/* Once the starting matching pts and corresponding strings found, search for 
   the starting and ending numbers of extremities of matching segments of 
   strings. 
*/
//=============================================================================
E_Boolean CString::compMatching(CString* s1, CString* s2,
                                vector<CString*>& strings,
                                E_Int& istart1, E_Int& istart2,
                                E_Int& iend1, E_Int& iend2) 
{
  E_Boolean match = false;
  FldArrayIS& dejaVu1 = s1->getDejaVu();
  FldArrayI& ind1 = s1->getIndArray();

  E_Int i2;
  CString* str;  
  iend1 = -1;
  iend2 = -1;
  E_Boolean found;

  for ( E_Int i1 = istart1+1; i1 < ind1.getSize(); i1++ )
  {
    // the point in current string has not been tested yet
    if ( dejaVu1[i1] < 2 )
    {
      found = s1->nearestPoint( i1, strings, i2, str);
      if ( found == false)
        return false;
      
      if ( str == s2 )
      { 
        
        iend1 = i1;
        iend2 = i2;
        match = true;
      }   
      else // sortie comp matching
        return match;
    }
    
    else // dejaVu   
      return match;
  }
  return match;
}

//=============================================================================
/* Return the nearest point on another string.
   This routine takes the dejaVu field into account.
   This routine can fail, if no point is available in all the strings :false
*/
//=============================================================================
E_Boolean CString::nearestPoint(E_Int ifirst1, vector<CString*>& strings,
                                E_Int& ifirst2, CString*& stringOut)
{
  E_Float dmin = E_MAX_FLOAT;
  E_Float dloc;
  E_Float dx, dy, dz;
  FldArrayF& coord = _blk->getCoord();
 
  // Coordinates of ind = _ind[ifirst1]
  E_Int ind = _ind[ifirst1];
  E_Float x1 = coord(ind, 1);
  E_Float y1 = coord(ind, 2);
  E_Float z1 = coord(ind, 3); 
  E_Boolean found = false;
  E_Int stringsSize =strings.size();
  for ( E_Int str = 0; str < stringsSize; str++)
  {
    CString* otherString = strings[str];
    StructBlock* blk2 = otherString->getBlock();
    // NB : strings must not be on the same block because they come 
    // from overlapping...
    if ( blk2 != _blk )
    {
      FldArrayIS& dejaVu = otherString->getDejaVu();
      FldArrayF& coord2 = blk2->getCoord();
      FldArrayI& otherInd = otherString->getIndArray();
      E_Float* coord2x = coord2.begin(1);
      E_Float* coord2y = coord2.begin(2);
      E_Float* coord2z = coord2.begin(3);
      
      for ( E_Int i2 = 0; i2 < otherInd.getSize(); i2++ )
      {
        if ( dejaVu[i2] < 2 )
        {
          E_Int ind2 = otherInd[i2];
          dx = coord2x[ind2] - x1;
          dy = coord2y[ind2] - y1;
          dz = coord2z[ind2] - z1;
          dloc = dx * dx + dy * dy + dz * dz;
          if ( dloc < dmin  )
          {
            dmin = dloc;
            ifirst2 = i2;
            stringOut = otherString;
            found = true;
          }
        }
      }
    }
  }
  return found;
}


//=============================================================================
/* Set dejaVu indices for points that are linked in segments
   from istart to iend indices (istart and iend are included ) */
//=============================================================================
void CString::setDejaVu(CString* str1, E_Int is1, E_Int ie1,
                        CString* str2, E_Int is2, E_Int ie2)
{
  FldArrayIS& dejaVu1 = str1->getDejaVu();
  FldArrayIS& dejaVu2 = str2->getDejaVu();
  E_Int start, end;

  // Treatment of points of str1 
  //----------------------------
  for (E_Int i = is1+1; i < ie1; i++)
    dejaVu1[i] = 2;
  
  start = E_min(is2, ie2);
  end = E_max(is2,ie2);
  for (E_Int i = start+1; i < end; i++)
    dejaVu2[i] = 2;
 
  //correct for extrema
  dejaVu1[is1]++;
  dejaVu1[ie1]++;
  dejaVu2[is2]++;
  dejaVu2[ie2]++;
}

//=======================================================================
/* Select the best segment 
   Choose between two pts is1 and js1 from str1 matching with is2 and js2
   from str2. Select then the segments that are the closest to each other
   Return ifound1 = best(is1, js1) and ifound2 = best(is2, js2)        
*/
//=======================================================================
void CString::selectClosestSegment( CString* s1, CString* s2,
                                    E_Int is1, E_Int js1,
                                    E_Int ie1, E_Int je1,
                                    E_Int is2, E_Int js2,
                                    E_Int ie2, E_Int je2,
                                    E_Int& fs1, E_Int& fe1,
                                    E_Int& fs2, E_Int& fe2)
{
  E_Int ind1, ind2;
  E_Float dx, dy, dz;
  E_Float dist, distmin;
  FldArrayF& coord1 = s1->getBlock()->getCoord();
  FldArrayF& coord2 = s2->getBlock()->getCoord();
  FldArrayI& indArray1 = s1->getIndArray();
  FldArrayI& indArray2 = s2->getIndArray();

  /* Correspondance entre les points */
  E_Int i11 = E_min(js1,je1);
  E_Int i21 = E_min(js2,je2);
  E_Int i12 = E_max(js1,je1);
  E_Int i22 = E_max(js2,je2);

  /* min distance between (is1,i11) and (is2, i12)*/
  fs1 = is1;
  fs2 = is2;
  ind1 = indArray1[is1];
  ind2 = indArray2[is2];
  dx = coord1(ind1, 1) - coord2(ind2, 1);
  dy = coord1(ind1, 2) - coord2(ind2, 2);
  dz = coord1(ind1, 3) - coord2(ind2, 3);
  distmin = dx * dx + dy * dy + dz * dz;

  ind2 = indArray2[i21];
  dx = coord1(ind1, 1) - coord2(ind2, 1);
  dy = coord1(ind1, 2) - coord2(ind2, 2);
  dz = coord1(ind1, 3) - coord2(ind2, 3);
  dist = dx * dx + dy * dy + dz * dz;

  if (dist < distmin)
  {
    fs2 = i21;
    distmin = dist;
  }
  
  ind1 = indArray1[i11];
  ind2 = indArray2[is2];
  dx = coord1(ind1, 1) - coord2(ind2, 1);
  dy = coord1(ind1, 2) - coord2(ind2, 2);
  dz = coord1(ind1, 3) - coord2(ind2, 3);
  dist = dx * dx + dy * dy + dz * dz;

  if (dist < distmin)
  {
    fs1 = i11;
    fs2 = is2;
    distmin = dist;
  }
  
  ind2 = indArray2[i21];
  dx = coord1(ind1, 1) - coord2(ind2, 1);
  dy = coord1(ind1, 2) - coord2(ind2, 2);
  dz = coord1(ind1, 3) - coord2(ind2, 3);
  dist = dx * dx + dy * dy + dz * dz;

  if (dist < distmin)
  {
    fs1 = i11;
    fs2 = i21;
    distmin = dist;
  }

  /* min distance between (ie1,i12) and (ie2, i22)*/
  ind1 = indArray1[ie1];
  ind2 = indArray2[ie2];
  fe1 = ie1;
  fe2 = ie2;
  dx = coord1(ind1, 1) - coord2(ind2, 1);
  dy = coord1(ind1, 2) - coord2(ind2, 2);
  dz = coord1(ind1, 3) - coord2(ind2, 3);
  distmin = dx * dx + dy * dy + dz * dz;

  ind2 = indArray2[i22];
  dx = coord1(ind1, 1) - coord2(ind2, 1);
  dy = coord1(ind1, 2) - coord2(ind2, 2);
  dz = coord1(ind1, 3) - coord2(ind2, 3);
  dist = dx * dx + dy * dy + dz * dz;

  if (dist < distmin)
  {
    fe2 = i22;
    distmin = dist;
  }
  
  ind1 = indArray1[i12];
  ind2 = indArray2[ie2];
  dx = coord1(ind1, 1) - coord2(ind2, 1);
  dy = coord1(ind1, 2) - coord2(ind2, 2);
  dz = coord1(ind1, 3) - coord2(ind2, 3);
  dist = dx * dx + dy * dy + dz * dz;

  if (dist < distmin)
  {
    fe1 = i12;
    fe2 = ie2;
    distmin = dist;
  }
  
  ind2 = indArray2[i22];
  dx = coord1(ind1, 1) - coord2(ind2, 1);
  dy = coord1(ind1, 2) - coord2(ind2, 2);
  dz = coord1(ind1, 3) - coord2(ind2, 3);
  dist = dx * dx + dy * dy + dz * dz;

  if (dist < distmin)
  {
    fe1 = i12;
    fe2 = i22;
    distmin = dist;
  }  
}
//=======================================================================
void CString::writeCoordOfStrings(CString* str2, E_Int is1, E_Int ie1,
                                  E_Int is2, E_Int ie2)
{
  StructBlock* blk1 = this->getBlock();
  StructBlock* blk2 = str2->getBlock();
  FldArrayF& coord1 = blk1->getCoord();
  FldArrayF& coord2 = blk2->getCoord();
  
  E_Int ind11 = this->getIndArray()[is1];
  E_Int ind12 = this->getIndArray()[ie1];
  E_Float x11 = coord1(ind11, 1);
  E_Float x12 = coord1(ind12, 1);
  
  E_Int ind21 = str2->getIndArray()[is2];
  E_Int ind22 = str2->getIndArray()[ie2];
  E_Float x21 = coord2(ind21, 1);
  E_Float x22 = coord2(ind22, 1); 
  printf(" xmin1= %f xmax1 = %f\n", x11, x12);
  printf(" xmin2= %f xmax2 = %f\n", x21, x22);
}

//=======================================================================
/* Prepare the list of strings that come from 
   1- overlapping blocks of the current string' s block
   2- matching blks but current string and matching block's string must 
   not match at extremities */
//=======================================================================
void CString::selectStrings(vector<CString*>& strings, E_Int noStr)
{
  E_Float matchTol = _blk->getMatchTol();
  
  vector<StructBlock*>&
    listOfOverlappingBlks =_blk->getListOfOverlappingBlks();
  
  vector<StructBlock*>& 
    listOfMatchingBlks = _blk->getListOfMatchingBlks();
  
  vector<StructBlock*> otherList;
  E_Int stringsSize = strings.size();
  for (E_Int s = 0; s < stringsSize; s++)
  {
    if (s != noStr)
      this->storeMatchingStringsInfo(strings, noStr, s);
  }

  E_Int cnt = 0;
  FldArrayIS& iblank1 = _blk->getIBlankArray();
  E_Float x1, y1, z1;
  E_Float dx, dy, dz;
  E_Int ind1, ind2;
  E_Int im1 = _blk->getIm();
  E_Int jm1 = _blk->getJm();
  FldArrayF& coord1 = _blk->getCoord();
  E_Int i1, j1;

  for (E_Int l = 0; l < _ind.getSize(); l++)
  {
    ind1 = _ind[l];
    j1 = ind1 / im1 + 1;
    i1 = ind1 - (j1-1) * im1 + 1;
    if ( i1 == 1 || i1 == im1 || j1 == 1 || j1 == jm1 )
    {
      x1 = coord1(ind1, 1);
      y1 = coord1(ind1, 2);
      z1 = coord1(ind1, 3);
      
      for (unsigned int v = 0; v < listOfMatchingBlks.size(); v++)
      {
        FldArrayIS& iblank2 = listOfMatchingBlks[v]->getIBlankArray();
        E_Int im2 = listOfMatchingBlks[v]->getIm();
        E_Int jm2 = listOfMatchingBlks[v]->getJm();
        FldArrayF& coord2 = listOfMatchingBlks[v]->getCoord();
        
        //frontiere j = 1
        for (E_Int i2 = 1; i2 <= im2 ; i2++)
        {
          ind2 = (i2-1);
          dx = coord2(ind2, 1) - x1;
          dy = coord2(ind2, 2) - y1;
          dz = coord2(ind2, 3) - z1;
          dx = E_abs(dx);
          dy = E_abs(dy);
          dz = E_abs(dz); 
     
          if ( dx < matchTol && dy < matchTol && dz < matchTol )
          {
            if ( iblank1[ind1] != 0 && iblank2[ind2] == 0 )
            {
              otherList.push_back(listOfMatchingBlks[v]);
              goto nextBlk;
            }
            else cnt++; 
          }
        }
        // frontiere j = jm2
        for (E_Int i2 = 1 ; i2 <= im2 ;i2++)
        {
          ind2 = (i2-1) + (jm2-1) * im2;
          dx = coord2(ind2, 1) - x1;
          dy = coord2(ind2, 2) - y1;
          dz = coord2(ind2, 3) - z1;
          dx = E_abs(dx);
          dy = E_abs(dy);
          dz = E_abs(dz);
        
          if ( dx < matchTol && dy < matchTol && dz < matchTol )
          {
            if ( iblank1[ind1] != 0 && iblank2[ind2] == 0  )
            {
              otherList.push_back(listOfMatchingBlks[v]);
              goto nextBlk;
            }
            else cnt++;
          }
        }
        //frontiere i= 1
        for (E_Int j2 = 1 ; j2 <= jm2 ;j2++)
        {
          ind2 = (j2-1) * im2;
          dx = coord2(ind2, 1) - x1;
          dy = coord2(ind2, 2) - y1;
          dz = coord2(ind2, 3) - z1;
          dx = E_abs(dx);
          dy = E_abs(dy);
          dz = E_abs(dz);

          if ( dx < matchTol && dy < matchTol && dz < matchTol )
          {
            if ( iblank1[ind1] != 0 && iblank2[ind2] == 0  )
            {
              otherList.push_back(listOfMatchingBlks[v]);
              goto nextBlk;
            }
            else cnt++;
          }
        }
        //frontiere i = im2
        for (E_Int j2 = 1 ; j2 <= jm2 ;j2++)
        {
          ind2 = (im2-1) + (j2-1) * im2;
          dx = coord2(ind2, 1) - x1;
          dy = coord2(ind2, 2) - y1;
          dz = coord2(ind2, 3) - z1;
          dx = E_abs(dx);
          dy = E_abs(dy);
          dz = E_abs(dz); 
          if ( dx < matchTol && dy < matchTol && dz < matchTol )
          {
            if ( iblank1[ind1] != 0 && iblank2[ind2] == 0 )
            {
              otherList.push_back(listOfMatchingBlks[v]);
              goto nextBlk;
            }
            else cnt++;
          }
        }
        // deux coins peuvent matcher 
        if ( cnt == 1 )
          otherList.push_back(listOfMatchingBlks[v]);
        nextBlk:;
      }
    }
  }
  
  // go through the string's list.
  // if the string comes from an overlapping block, add it to the possible
  // matching strings list
  E_Int otherListSize = otherList.size();
  for (E_Int i = 0 ; i < stringsSize; i++)
  {  
    if ( i != noStr )
    {
      StructBlock* blks = strings[i]->getBlock();
      
      for (unsigned int j = 0 ; j < listOfOverlappingBlks.size() ; j++ )
        if ( blks == listOfOverlappingBlks[j] )
        { 
          if ( _matchInfo(0,1) != i && _matchInfo(1,1) != i )
            _listOfMatchingStrings.push_back(strings[i]);
          goto fin;
        }
      
      
      for (E_Int j2 = 0 ; j2 < otherListSize; j2++ )
      {
        if ( blks == otherList[j2] )
        {
          if ( _matchInfo(0,1) != i && _matchInfo(1,1) != i )
            _listOfMatchingStrings.push_back(strings[i]);   
          break;
        }
      }
    }
    fin:;
  }
}
//=======================================================================
/* Test if strings' extremities match */ 
//=======================================================================
void CString::storeMatchingStringsInfo(vector<CString*>& strings,
                                       E_Int noStr1,
                                       E_Int noStr2)
{
  CString* str1 = strings[noStr1];
  CString* str2 = strings[noStr2];
  E_Float matchTol = _blk->getMatchTol();
  FldArrayI& ind1Array = str1->getIndArray();
  FldArrayI& ind2Array = str2->getIndArray();

  FldArrayF& coord1 = str1->getBlock()->getCoord();  
  FldArrayF& coord2 = str2->getBlock()->getCoord();
  E_Int max1 = ind1Array.getSize();
  E_Int max2 = ind2Array.getSize();
  
  E_Int ind10 = ind1Array[0];
  E_Int ind1M = ind1Array[max1-1];
  E_Int ind20 = ind2Array[0];
  E_Int ind2M = ind2Array[max2-1];

  E_Float x10 = coord1(ind10, 1);
  E_Float y10 = coord1(ind10, 2);
  E_Float z10 = coord1(ind10, 3);
  E_Float x20 = coord2(ind20, 1);
  E_Float y20 = coord2(ind20, 2);
  E_Float z20 = coord2(ind20, 3);

  E_Float x1M = coord1(ind1M, 1);
  E_Float y1M = coord1(ind1M, 2);
  E_Float z1M = coord1(ind1M, 3);
  E_Float x2M = coord2(ind2M, 1);
  E_Float y2M = coord2(ind2M, 2);
  E_Float z2M = coord2(ind2M, 3);
  
  E_Float dx = x10-x20;
  E_Float dy = y10-y20; 
  E_Float dz = z10-z20;
  dx = E_abs(dx);
  dy = E_abs(dy);
  dz = E_abs(dz);

  // min extremity of str1 ->store in 0 column
  //==========================================
  if ( dx < matchTol && dy < matchTol && dz < matchTol )
  {
    _matchInfo(0,1) = noStr2; 
    _matchInfo(0,2) = 0; // min extremity of str2 ->0 
    return;
  }
  dx = x10-x2M;
  dy = y10-y2M; 
  dz = z10-z2M;
  dx = E_abs(dx);
  dy = E_abs(dy);
  dz = E_abs(dz);

  if ( dx < matchTol && dy < matchTol && dz < matchTol )
  {
    _matchInfo(0,1) = noStr2;
    _matchInfo(0,2) = 1; //max extremity of str2 -> 1
    return;
  }
  
  // max extremity of str1 ->store in 1 column
  //==========================================
  dx = x1M-x20;
  dy = y1M-y20; 
  dz = z1M-z20;
  dx = E_abs(dx);
  dy = E_abs(dy);
  dz = E_abs(dz);
  
  if ( dx < matchTol && dy < matchTol && dz < matchTol )
  {
    _matchInfo(1,1) = noStr2;
    _matchInfo(1,2) = 0; // min extremity of str2 ->0 
    return;
  }
  dx = x1M-x2M;
  dy = y1M-y2M; 
  dz = z1M-z2M;
  dx = E_abs(dx);
  dy = E_abs(dy);
  dz = E_abs(dz);

  if ( dx < matchTol && dy < matchTol && dz < matchTol )
  {  
    _matchInfo(1,1) = noStr2;
    _matchInfo(1,2) = 1; // max extremity of str2 ->max   
    return;
  }
  return;
}
//=======================================================================
/* Test if selected segments are overlapping */
//=======================================================================
E_Boolean CString::areSegmentsOverlapping(E_Int is, E_Int ie,
                                          E_Int js, E_Int je)
{
  E_Int imin, imax;
  E_Int jmin, jmax;

  imin = E_min(is, ie); // normalement c'est deja ordonne
  imax = E_max(is, ie); //idem
  jmin = E_min(js, je);
  jmax = E_max(js, je);
      
  if ( jmax >= imin && jmin <= imax )
    return true;
  else 
    return false;
}

//=============================================================================
/* Identify points that are not assembled in segments and put them 
   in a list.
   Cette routine construit remainingSeg :
   C'est une liste de tableaux.
   Chaque tableau contient les indices des points qui n'ont pas ete deja
   inseres dans un segment et qui se suivent sur la grille.
*/
//=============================================================================
void CString::compRemainingSegments(list<FldArrayI*>& remainingSeg)
{
  E_Int size = _ind.getSize(); 
  E_Int cnt = 0;
  FldArrayI pointsForPockets(size);
  FldArrayI flag(size);
  E_Boolean isCut;
  E_Int isize;
  E_Int size1 = 0;
  E_Int size2 = 0;
  list<FldArrayI*>::iterator itr;
  list<FldArrayI*>::iterator itr2;
  list<FldArrayI*>::iterator itr3;
  list<FldArrayI*>::iterator itr4;
  list<FldArrayI*> listOfFlags;

  /* Identify remaining points */
  for (E_Int i = 0; i < size; i++)
  {
    if ( _flag[i] < 2 )
    {
      pointsForPockets[cnt] = _ind[i];
      flag[cnt] = _flag[i];
      cnt++;
    }
  }
  pointsForPockets.reAlloc(cnt);
  flag.reAlloc(cnt);
  FldArrayB dejaVu(cnt);
  dejaVu.setAllValuesAt(false);

  /*--------------*/
  /* Make strings */
  /*--------------*/
  FldArrayI links(cnt, 4); // overdimensionned, we never know
  _blk->compLinks(pointsForPockets, links);
 
  E_Boolean goOn = true;

  while (goOn == true)
  {
    FldArrayI* string = new FldArrayI(cnt);
    isize = _blk->chainPoints(pointsForPockets, links, 
                              *string, dejaVu);

    if (isize == 0)
    {
      delete string;
      goOn = false;
    }
    else if ( isize == 1) //elimination d'un point unique
    {
      string->reAlloc(0);
      delete string;
      goOn = true;
    }
    else if (isize > 1) 
    {
      string->reAlloc(isize);
      remainingSeg.push_back(string);
    }
  }

  // Change index into mesh index
  for (itr = remainingSeg.begin();
       itr != remainingSeg.end(); itr++)
  {
    isize = (*itr)->getSize();
    FldArrayI* oneFlag = new FldArrayI(isize);
    for (E_Int i = 0; i < isize; i++)
    {
      E_Int a = (**itr)[i];
      (**itr)[i] = pointsForPockets[a];
      (*oneFlag)[i] = flag[a];
    }
    listOfFlags.push_back(oneFlag);
  }
  
  
  /* If a string contains 2 consecutive pts with flag =-1
     cut the string into 2 strings */
  itr = remainingSeg.begin();
  itr2 = listOfFlags.begin();
  while (itr != remainingSeg.end() && 
         itr2 != listOfFlags.end())
  {
    FldArrayI& oneFlag = **itr2;
    FldArrayI& oneSeg = **itr;
    isCut = false;
    isize = (*itr)->getSize();
    // test if string need to be cut :
    for (E_Int ind = 0; ind < isize-1; ind++)
    {
      if ( oneFlag[ind] == -1 && oneFlag[ind+1] == -1)
      {
        size1 = ind+1;
        size2 = isize-size1;
        isCut = true;
        break;
      }
    }
    if (isCut == true)
    {
      FldArrayI* flag1 = new FldArrayI(size1);
      FldArrayI* flag2 = new FldArrayI(size2);
      FldArrayI* seg1 = new FldArrayI(size1);
      FldArrayI* seg2 = new FldArrayI(size2);
      // tile the current segment into 2 parts
      for (E_Int ind = 0; ind < size1; ind++)
      {
        (*flag1)[ind] = oneFlag[ind];
        (*seg1)[ind] = oneSeg[ind];
      }
    
      for (E_Int ind = 0; ind < size2; ind++)
      {
        (*flag2)[ind] = oneFlag[ind+size1];
        (*seg2)[ind] = oneSeg[ind+size1];
      }
      
      // Erase segment from list
      itr3 = itr;
      itr4 = itr2;
      itr++;
      itr2++;
      delete *itr3;
      delete *itr4;
      remainingSeg.erase(itr3);
      listOfFlags.erase(itr4);
      remainingSeg.push_back(seg1);
      remainingSeg.push_back(seg2);
      listOfFlags.push_back(flag1);
      listOfFlags.push_back(flag2);
    }
    else 
    {
      itr++;
      itr2++;
    }
  }
  for(list<FldArrayI*>::iterator itr1 = listOfFlags.begin();
      itr1 != listOfFlags.end(); itr1++)
    delete *itr1;
  listOfFlags.clear();
}
//======================================================================
/* Given a point ibeg on string s1, search for nearest point on string 
   s2 between bounds min and max of string s2
   Return the corresponding index if found, -1 elsewhere 
*/
//======================================================================
E_Int CString::findClosestPoint(CString* s1, CString* s2, E_Int ibeg, 
                                E_Int min, E_Int max) 
{
  E_Int ifound = -1;
  FldArrayF& coord1 = s1->getBlock()->getCoord();
  FldArrayF& coord2 = s2->getBlock()->getCoord();
  
  // find nearest pt ifound on s2 to ibeg coming from string s1:
  FldArrayI& indArray1 = s1->getIndArray();
  FldArrayI& indArray2 = s2->getIndArray();
  E_Int ind1 = indArray1[ibeg];

  E_Float x1 = coord1(ind1,1);
  E_Float y1 = coord1(ind1,2);
  E_Float z1 = coord1(ind1,3);

  E_Int ind2;
  E_Float dx, dy, dz;
  E_Float distmin = E_MAX_FLOAT;
  E_Float dist;

  for (E_Int i = min; i <= max; i++)
  {
    ind2 = indArray2[i];
    dx = coord2(ind2,1)-x1;
    dy = coord2(ind2,2)-y1;
    dz = coord2(ind2,3)-z1;

    dist = dx*dx + dy*dy + dz*dz;
    if ( dist < distmin)
    {
      distmin = dist;
      ifound = i;
    }
  } 
  
  if (ifound ==  -1)
    printf("WARNING: no close point was found between bounds.\n");
  return ifound;
}

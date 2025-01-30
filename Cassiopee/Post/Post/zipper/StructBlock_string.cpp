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

# include "StructBlock.h"
# include <stdlib.h>

using namespace K_FUNC;
using namespace std;
using namespace K_FLD;
using namespace K_CONST;

//=============================================================================
/* Stringing.
   Ordered the point in input fldarray, such that 1D strings can be formed.
   Output in a list of arrays of indices.
*/
//=============================================================================
void StructBlock::stringing(list<FldArrayI*>& strings)
{
  FldArrayI& inIndices = _pointsForStrings;
  list<FldArrayI*> strings1;

  // Keep track of already treated points
  E_Int n = inIndices.getSize();
  FldArrayB dejaVu(n);
  dejaVu.setAllValuesAt(false);

  // Compute the links inside inIndices
  FldArrayI links(n, 4); // overdimensionned, we never know
  compLinks(inIndices, links);
  
  // Temporary list storage
  list<FldArrayI*>::iterator itr;
  
  // Find strings in inIndices
  E_Boolean goOn = true;

  while (goOn)
  {
    FldArrayI* string = new FldArrayI(n);
    E_Int size = chainPoints(inIndices, links, *string, dejaVu);
    if (size == 0)
    {
      delete string;
      goOn = false;
    }
    else
    {
      string->reAlloc(size);
      strings1.push_back(string);
    }
  }

  // Change index into mesh index
  for (itr = strings1.begin(); itr != strings1.end(); itr++)
  {
    for (E_Int i = 0; i < (*itr)->getSize(); i++)
    {
      E_Int a = (**itr)[i];
      (**itr)[i] = inIndices[a];
    }
  }

  /* Test if string obtained is a loop */
  for (itr = strings1.begin(); itr != strings1.end(); itr++)
  {
    E_Int size = (*itr)->getSize();
    E_Boolean isALoop = isStringALoop(**itr);
    if ( isALoop )
    {
      E_Int t1 = E_Int(size/4);
      E_Int t3 = E_Int(3*size/4);

      if (t3 < size && t3 > t1)
      {
        E_Int cnt = 0;
        E_Int size1 = t3-t1+1;
        E_Int size2 = size-size1+2;
        
        FldArrayI* s1 = new FldArrayI(size1);
        FldArrayI* s2 = new FldArrayI(size2);
        
        for (E_Int i = t1; i <= t3; i++)
        {
          (*s1)[cnt] = (**itr)[i];
          cnt++;
        }
        cnt = 0;
        for (E_Int i = t3; i < size; i++)
        {
          (*s2)[cnt] = (**itr)[i];
          cnt++;
        }
        for (E_Int i = 0; i <= t1; i++)
        {
          (*s2)[cnt] = (**itr)[i];
          cnt++;
        }
        strings.push_back(s1);
        strings.push_back(s2);
      }
      else 
      {
        E_Int size1 = E_Int(size/2)+1;
        E_Int size2 = size-size1+1;
        
        FldArrayI* s1 = new FldArrayI(size1);
        FldArrayI* s2 = new FldArrayI(size2);
        for (E_Int i = 0; i < size1; i++)
          (*s1)[i] = (**itr)[i];
        
        for (E_Int i = 0; i < size2; i++)
          (*s2)[i] = (**itr)[i+size1-1];

        strings.push_back(s1);
        strings.push_back(s2);
      }
      delete *itr;
    }
    else // False
      strings.push_back(*itr);
  }
}

//=============================================================================
/* Tell if the string is a closed string */
//=============================================================================
E_Boolean StructBlock::isStringALoop(FldArrayI& string)
{
  E_Int size = string.getSize();
  E_Int ind0 = string[0];
  E_Int ind1 = string[size-1];

  E_Float dx = _coord(ind1, 1) - _coord(ind0, 1);
  E_Float dy = _coord(ind1, 2) - _coord(ind0, 2);
  E_Float dz = _coord(ind1, 3) - _coord(ind0, 3);
  dx = E_abs(dx);
  dy = E_abs(dy);
  dz = E_abs(dz); 
  if ( dx < _matchTol && dy < _matchTol && dz < _matchTol )
    return true;
  
  E_Int j0 = ind0/_im + 1;
  E_Int i0 = ind0 - (j0-1)*_im + 1; 
  
  E_Int j1 = ind1/_im + 1;
  E_Int i1 = ind1 - (j1-1)*_im + 1;
  
  E_Int di = E_abs(i1-i0);
  E_Int dj = E_abs(j1-j0);
  if( ( di == 1 && dj == 0) || (di == 0 && dj == 1) )
    return true;

  return false;
}

//=============================================================================
// Return size : string completed.
// Return 0 : no starting point can be found in inIndices and all points
// exhausted.
//=============================================================================
E_Int StructBlock::chainPoints(FldArrayI& inIndices, FldArrayI& links, 
                               FldArrayI& string, FldArrayB& dejaVu)
{
  E_Int n = inIndices.getSize();

  // Try to find a starting point
  E_Boolean startFound = false;
  for (E_Int ir = 0; ir < n; ir++)
  {
    if (dejaVu[ir] == false &&
        (links(ir, 1) == -999 || links(ir, 2) == -999))
    {
      // It is a starting point
      string[0] = ir;
      dejaVu[ir] = true;
      startFound = true;
      break;
    }
  }
  
  if (startFound == false)
  {
    E_Boolean alive = false;
    for (E_Int ir = 0; ir < n; ir++)
    {
      if (dejaVu[ir] == false)
      {
        string[0] = ir;
        dejaVu[ir] = true;
        alive = true;
        break;
      }
    }  
    
    if (alive == false)
      return 0;
  }
    
  // Try to find following elements
  E_Int c = 0;
  E_Int start = string[c]; c++;
  E_Int l1, l2;
  E_Boolean goOn = true;
  
  while (goOn)
  {
    l1 = links(start, 1);
    l2 = links(start, 2);
  
    if (l1 != -999 && l1 != start && dejaVu[l1] == false)
    {
      string[c] = l1; c++;
      start = l1;
      dejaVu[start] = true;
    }
    else if (l2 != -999 && l2 != start && dejaVu[l2] == false)
    {
      string[c] = l2; c++;
      start = l2;
      dejaVu[start] = true;
    }
    else
      goOn = false;
  }  
  return c;
}

//=============================================================================
/* compLinks
   Compute the links in a array of indices.
   Two point are linked if they are neighbour on the grid.
   WARNING:
   2 gap pts might be linked by a cell line that connect 2 unblanked cells.
   _iblankCC allows for their identification.
   
   Output in a array links gathering which index is linked to which one.
*/
//============================================================================
void StructBlock::compLinks(FldArrayI& inIndices, FldArrayI& links)
{ 
  // Init links array
  links.setAllValuesAt(-999);

  E_Int i, j, k;
  E_Int ind;
  E_Int n = inIndices.getSize();
  E_Int im = _im;
  E_Int jm = _jm;
  
  for ( E_Int ir = 0; ir < n; ir++)
  {
    ind = inIndices[ir];
    k = ind / (im*jm);
    j = (ind - k*im*jm) / im;
    i = ind - j*im -k*im*jm;
    
    if  ( ( i >= 1 && i < im-1 ) && ( j >= 1 && j < jm-1) )
      compLinksForInteriorPts(inIndices, ir, i, j, k, links);
    else
      compLinksForBndPts(inIndices, ir, i, j, k, links);
  }
}

//=============================================================================
/* Identify the points to be stringed, namely Gap Boundary Points
   They are of 2 types :
 TypeI: unblanked interior point which has from 3 to 7 over 8 unblanked
        neighbours, boundary or corner pt if any one of its 3 and 2
        neighbours respectively are blanked,
 TypeII: Overlapping borders (imax, jmax...),
 TypeIII: Point d'un raccord coincident d'un bloc 1 correspondant a un point
 iblank */
//=============================================================================
void StructBlock::identifyGapBndPts(E_Int noBlk,
                                    vector<StructBlock*> vectOfBlocks)
{
  E_Int n = 0; // index of the stringed point in _pointsForStrings

  E_Int im = _im;
  E_Int jm = _jm;

  _pointsForStrings.malloc(im*jm);
 
  // 1-type I points
  //----------------
  // Interior points
  compInteriorPtsForStringing(im, jm, n);

  // Boundary points except corners
  compBndPtsForStringing(im, jm, n);

  // Four corners
  compCornersForStringing(im, jm, n);
 
  // 2-type II points
  //-----------------
  compOverlapBordersForString( n);
 
  // 3-type III points
  //------------------
  compMatchingBordersForString(noBlk, n, vectOfBlocks);
 
  _pointsForStrings.reAlloc(n);

  // Eliminate double points
  eliminateDoublePoints(_pointsForStrings);
}

//=============================================================================
void StructBlock::eliminateDoublePoints(FldArrayI& indices)
{
  E_Int n = indices.getSize();
  FldArrayI work(n);

  E_Int c = 0;
  for (E_Int i = 0; i < n; i++)
  {
    E_Boolean found = false;
    for (E_Int j = 0; j < i; j++)
    {
      if (indices[i] == indices[j])
      {
        found = true;
        break;
      }
    }
    if (found == false)
    {
      work[c] = indices[i]; c++;
    }
  }
  work.reAlloc(c);
  indices = work; // copy
}

//=============================================================================
/* Compute the interior points to be stringed -
   To be valid, the interior pt must have from 3 to 7 unblanked neighbours
*/
//=============================================================================
void StructBlock::compInteriorPtsForStringing(E_Int im, E_Int jm, E_Int& n)
{
  FldArrayI ind(9);
  E_Int cnt;
 
  for (E_Int j = 2; j <= jm-1; j++)
    for (E_Int i = 2; i <= im-1; i++)
    {
      ind[0] = (i-1) + (j-1) * im ; // i,j
     
      // the current point must be unblanked
      if ( _iblank[ind[0]] != 0 )
      {
        ind[1] = (i-2) + (j-2) * im ; // i-1,j-1
        ind[2] = (i-1) + (j-2) * im ; // i,j-1
        ind[3] =  i    + (j-2) * im ; // i+1,j-1
        ind[4] =  i    + (j-1) * im ; // i+1,j
        ind[5] =  i    +  j    * im ; // i+1,j+1
        ind[6] = (i-1) +  j    * im ; // i,j+1
        ind[7] = (i-2) +  j    * im ; // i-1,j+1
        ind[8] = (i-2) + (j-1) * im ; // i-1,j
        
        // Count the number of neighbours that are unblanked
        cnt = 0;
        for (E_Int m = 1; m< ind.getSize() ; m++ )
          if ( _iblank[ind[m]] != 0 )
            cnt++;
        
        // If the interior point is valid for stringing 
        // to be valid: at least 3 and no more than 7 neighbours are unblanked
        if ( ( cnt >= 3 ) && ( cnt < 8 ) ) 
        {
          _pointsForStrings[n] = ind[0];
          n++;
        }
      }
    }
}

//=============================================================================
/* Compute the boundary points to be stringed-
   To be valid, the boundary pt must have at least one blanked neighbour */
//=============================================================================
void StructBlock::compBndPtsForStringing(E_Int im, E_Int jm, E_Int& n)
{
  FldArrayI ind(4);
  E_Int cnt;
  FldArrayF& cfd = getCfdField();
  E_Int cellNF;
  
  // I = cste boundaries
  for (E_Int j = 2; j <= jm-1; j++)
  {
    // i =1 border 
    ind[0] = (j-1) * im; // 1,j
    cellNF = E_Int(cfd(ind[0], _nfld));
    if ( _iblank[ind[0]] != 0 && cellNF != 2)
    {
      ind[1] =(j-2) * im ; // 1,j-1
      ind[2] =(j-1) * im + 1 ; // 2,j
      ind[3] = j    * im ; // 1,j+1
      // Count the number of neighbours that are unblanked
      cnt = 0;
      for (E_Int i = 1; i < ind.getSize() ; i++ )
        if ( _iblank[ind[i]] == 0 ) // criterium: blanked 
          cnt++;
      // if the bnd point is valid for stringing 
      // to be valid: at least 1 neighbour blanked
      if (  cnt > 0 ) 
      {
        _pointsForStrings[n] = ind[0];
        n++;
      }
    }
    
    // i = im border
    ind[0] = (im-1) + (j-1) * im; // im,j
    cellNF = E_Int(cfd(ind[0], _nfld));
    if ( _iblank[ind[0]] != 0 && cellNF != 2)
    {
      ind[1] = (im-1) + (j-2) * im ; // im,j-1
      ind[2] = (im-2) + (j-1) * im ; // im-1,j
      ind[3] = (im-1) +  j    * im ; // im,j+1
      //count the number of neighbours that are unblanked
      cnt = 0;
      for (E_Int i = 1; i < ind.getSize() ; i++ )
        if ( _iblank[ind[i]] == 0 ) // criterium: blanked 
          cnt++;
      // if the bnd point is valid for stringing 
      // to be valid: at least 1 neighbour blanked
      if (  cnt >0 ) 
      {
        _pointsForStrings[n] = ind[0];
        n++;
      }
    }
  }
  
  // J = cst boundaries
  for (E_Int i = 2; i <= im-1; i++)
  {
    // j = 1 border 
    ind[0] = (i-1); // i,1
    cellNF = E_Int(cfd(ind[0], _nfld));
    if ( _iblank[ind[0]] != 0 && cellNF != 2 )
    { 
      ind[1] = i-2;         // i-1,1
      ind[2] = i;           // i+1,1
      ind[3] = (i-1) + im ; // i,2
      //count the number of neighbours that are unblanked
      cnt = 0;
      for (E_Int l = 1; l < ind.getSize() ; l++ )
        if ( _iblank[ind[l]] == 0 ) // criterium: blanked 
          cnt++;
      // if the bnd point is valid for stringing 
      // to be valid: at least 1 neighbour blanked
      if (  cnt > 0 ) 
      {
        _pointsForStrings[n] = ind[0];
        n++;
      }
    }
    
    // j = jm border 
    ind[0] = (i-1) + (jm-1) * im; // i,jm
    cellNF = E_Int(cfd(ind[0], _nfld));
    if ( _iblank[ind[0]] != 0 && cellNF != 2 )
    { 
      ind[1] = i-2 + (jm-1) * im; // i-1,jm
      ind[2] = i   + (jm-1) * im; // i+1,jm
      ind[3] = i-1 + (jm-2) * im; // i,jm-1
      //count the number of neighbours that are unblanked
      cnt = 0;
      for (E_Int l = 1; l < ind.getSize() ; l++ )
        if ( _iblank[ind[l]] == 0 ) // criterium: blanked 
          cnt++;
      // if the bnd point is valid for stringing 
      // to be valid: at least 1 neighbour blanked
      if (  cnt > 0 ) 
      {
        _pointsForStrings[n] = ind[0];
        n++;
      }
    }
  }
}

//=============================================================================
/* Compute the boundary points to be stringed-
   To be valid, the boundary pt must have at least one blanked neighbour */
//=============================================================================
void StructBlock::compCornersForStringing(E_Int im, E_Int jm, E_Int& n)
{
  E_Int ind1, ind2, ind3, cellNF;
  FldArrayF& cfd = getCfdField();
  
  // Corner 1,1
  ind1 = 0;
  cellNF = E_Int(cfd(ind1, _nfld));
  if ( _iblank[ind1] != 0 && cellNF != 2)
  {
    ind2 = 1;//(2,1)
    ind3 = im;
    if ( (_iblank[ind2] == 0) || (_iblank[ind3] == 0 ) )
     {
       _pointsForStrings[n] = ind1;
       n++;
     }
  }
  
  //corner im,1
  ind1 = im-1;
  cellNF = E_Int(cfd(ind1, _nfld));
  if ( _iblank[ind1] != 0 && cellNF != 2)
  {
    ind2 = im-2;   //(im-1,1)
    ind3 = 2*im-1; //(im,2)
    if ( (_iblank[ind2] == 0) || (_iblank[ind3] == 0 ) )
     {
       _pointsForStrings[n] = ind1;
       n++;
     }
  }
  
  // Corner 1,jm
  ind1 = (jm-1) * im;
  cellNF = E_Int(cfd(ind1, _nfld));
  if ( _iblank[ind1] != 0 && cellNF != 2)
  {
    ind2 = (jm-2) * im;    //(1,jm-1)
    ind3 = (jm-1) * im + 1;//(2,jm) 
    if ( (_iblank[ind2] == 0) || (_iblank[ind3] == 0 ) )
     {
       _pointsForStrings[n] = ind1;
       n++;
     }
  }
  
  //corner im,jm
  ind1 = im-1 + (jm-1) * im;
  cellNF = E_Int(cfd(ind1, _nfld));
  if ( _iblank[ind1] != 0 && cellNF != 2)
  {
    ind2 = (im-2) + (jm-1) * im; //(im-1, jm)
    ind3 = (im-1) + (jm-2) * im; //(im,jm-1)
    if ( (_iblank[ind2] == 0) || (_iblank[ind3] == 0 ) )
     {
       _pointsForStrings[n] = ind1;
       n++;
     }
  }
}
//=============================================================================
/* Identify the overlapping borders */
//=============================================================================
void StructBlock::compOverlapBordersForString(E_Int& n)
{
  FldArrayF& cfd = _cfdField;
  E_Int ind, cellNF;
  E_Int indCell1, indCell2;
  E_Int im = _im;
  E_Int jm = _jm;

  E_Int im1 = im-1;
  E_Int jm1 = jm-1;
  
  // CORNERS
  //--------
  // 1st corner (1,1)
  ind = 0;
  cellNF = E_Int(cfd(ind, _nfld));
  if (cellNF == 2 && _iblank[ind] != 0)
  {
    indCell1 = 0;
    testForBlankedBorderNode(ind, indCell1, indCell1, n);
  }
  
  // 2nd corner (1,jm)
  ind = (jm-1) * im;
  cellNF = E_Int(cfd(ind, _nfld));
  if (cellNF == 2 && _iblank[ind] != 0)
  {
    indCell1 = (jm1-1) * im1;
    testForBlankedBorderNode(ind , indCell1, indCell1, n);
  }                        
  
  // 3rd corner (im,1)
  ind = im-1;
  cellNF = E_Int(cfd(ind, _nfld));
  if (cellNF == 2 && _iblank[ind] != 0)
  {
    indCell1 = im1-1;
    testForBlankedBorderNode(ind , indCell1, indCell1, n);
  }
  
  // last corner (im,jm)
  ind = (im-1) + (jm-1) * im;
  cellNF = E_Int(cfd(ind, _nfld));
  if (cellNF == 2 && _iblank[ind] != 0)
  {
    indCell1 = (im1-1) + (jm1-1) * im1;
    testForBlankedBorderNode(ind , indCell1, indCell1, n);
  }
  
  // INTERIOR NODES 
  //----------------
  // 1- Boundary i = 1
  //------------------
  for ( E_Int j = 2 ; j <= jm1 ; j++ )
  {
    ind = (j-1) * im;
    cellNF = E_Int(cfd(ind, _nfld));
    if ( cellNF == 2 && _iblank[ind] != 0)
    {
      indCell1 = (j-2) * im1;
      indCell2 = (j-1) * im1; 
      testForBlankedBorderNode(ind , indCell1, indCell2, n);
    }
  }
  // 2- Boundary i = imax
  //----------------------
  for ( E_Int j = 2; j <= jm1 ; j++)
  {
    ind = (im-1) + (j-1) * im;
    cellNF = E_Int(cfd(ind, _nfld));
    if ( cellNF == 2 && _iblank[ind] != 0)
    {
      indCell1 = (im1-1) +  (j-2)*im1;
      indCell2 = (im1-1) +  (j-1)*im1;
      testForBlankedBorderNode(ind , indCell1, indCell2, n);
    }
  }
  
  // 3- Boundary j = 1
  //---------------------
  for ( E_Int i = 2; i <= im1 ; i++)
  {
    ind = i-1;
    cellNF = E_Int(cfd(ind, _nfld));
    if ( cellNF == 2 && _iblank[ind] != 0)
    {
      indCell1 = i-2;
      indCell2 = i-1;
      testForBlankedBorderNode(ind , indCell1, indCell2, n);
    }
  }
  
  // 4- Boundary j = jm
  //----------------------
  for ( E_Int i = 2; i <= im1 ; i++)
  {
    ind = (i-1) + (jm-1) * im;
    cellNF = E_Int(cfd(ind, _nfld));
    if ( cellNF == 2 && _iblank[ind] != 0)
    {
      indCell1 = (i-2) + (jm1-1) * im1;
      indCell2 = (i-1) + (jm1-1) * im1;
      testForBlankedBorderNode(ind , indCell1, indCell2, n);
    }
  }
}

//===========================================================================
/* Test the Validity of the link between ind1 and ind2 in the direction dir
   If they are connected by a line that is inside the mesh resulting of
   the blanking then return false. */
//===========================================================================
E_Boolean StructBlock::testValidityOfLink(E_Int i, E_Int j, E_Int k,
                                          E_Int dir, E_Int sens)
{
  E_Int im = _im;
  E_Int jm = _jm;
  E_Int im1 = im-1;
  E_Int jm1 = jm-1;

  // indices of the cell centers
  E_Int indCell1, indCell2;
 
  switch ( dir )
  {
    case 1:
      
      i = E_min(i, i+sens); 
      if ( j >= 1 && j < jm1 )
      {
        indCell1 = i + (j-1)*im1;
        indCell2 = i + j*im1 ;
      }
      else if ( j == 0)
      {
        indCell1 = i;
        indCell2 = indCell1;
      }
      else// j= jm1 
      {
        indCell1 =  i + (jm1-1)*im1;
        indCell2 = indCell1;
      }
      
      if ( _iblankCC[indCell1] == 1 && _iblankCC[indCell2] == 1 )
        return false;
      else
        return true;
      
    case 2:

      j = E_min(j, j+sens);
      if ( i >= 1 && i < im-1 )
      {
        indCell1 = (i-1) + j*im1;
        indCell2 = i + j*im1;
      }
      else if ( i == 0)
      {
        indCell1 = j*im1 ;
        indCell2 = indCell1;
      }
      else// i= im1 
      {
        indCell1 = (im1-1) + j*im1;
        indCell2 = indCell1;
      }
       
      if ( _iblankCC[indCell1] == 1 && _iblankCC[indCell2] == 1 )
        return false;
      else
        return true;
   
    default:
      printf("Error: not a valid direction for linking.\n");
      exit(0);
      break;
  }
}

//=============================================================================
/* Compute Links for bnd Pts */
//=============================================================================
void StructBlock::compLinksForBndPts(FldArrayI& inIndices, E_Int ir,
                                     E_Int i, E_Int j, E_Int k,
                                     FldArrayI& links)
{
  E_Int n = inIndices.getSize();
  E_Int im = _im;
  E_Int jm = _jm;
  E_Int ind2, ind3;

  E_Boolean isValid;
  E_Int dir, sens;
  
  E_Int ind = inIndices[ir];

  E_Int l = 1;

  if ( i == 0 || i == im-1 )
  { 
    if (j + 1 < jm)
    {
      ind2 = ind+im;
      for (E_Int jr = 0; jr < n; jr++)
      {
        ind3 = inIndices[jr];
        if (ind3 == ind2)
        {
          links(ir, l) = jr;
          l++;
        }    
      }
    }

    if (j - 1 >= 0)
    {
      ind2 = ind-im;
      for (E_Int jr = 0; jr < n; jr++)
      {
        ind3 = inIndices[jr];
        if (ind3 == ind2)
        {
          links(ir, l) = jr;
          l++;
        }
      }
    }

    //last points
    if (i == 0 )
    {
      ind2 = ind+1;
      dir = 1 ;
      sens = 1;
    }
    else
    {
      ind2 = ind-1;
      dir = 1;
      sens = -1;
    }
    
    for (E_Int jr = 0; jr < n; jr++)
    {
      ind3 = inIndices[jr];
      if (ind3 == ind2)
      { 
        isValid = testValidityOfLink( i, j, k, dir, sens);
        if ( isValid == true )
        {
          links(ir, l) = jr;
          l++;
        }
      }
    }
  }
  
  if ( j == 0 || j == jm-1)
  {
    if (i+1 < im )
    {
      ind2 = ind+1;
      
      for (E_Int jr = 0; jr < n; jr++)
      {
        ind3 = inIndices[jr];
        if (ind3 == ind2)
        {
          links(ir, l) = jr;
          l++;
        }
      }
    }

    if (i - 1 >= 0)
    {
      ind2 = ind-1;
      for (E_Int jr = 0; jr < n; jr++)
      {
        ind3 = inIndices[jr];
        if (ind3 == ind2)
        {
          links(ir, l) = jr;
          l++;
        }
      }
    }

    //last point
    //-----------
    if ( j == 0 )
    {
      ind2 = ind + im;
      dir = 2;
      sens = 1;
    }
    else
    {
      ind2 = ind - im;
      dir = 2;
      sens = -1;
    }
    
    for (E_Int jr = 0; jr < n; jr++)
    {
      ind3 = inIndices[jr];
      if (ind3 == ind2)
      {
        isValid = testValidityOfLink( i, j, k, dir, sens);
        if ( isValid == true)
        {
          links(ir, l) = jr;
          l++;
        }
      }
    }
  }
}

//=============================================================================
/* Compute Links for Interior Pts */
//=============================================================================
void StructBlock::compLinksForInteriorPts(FldArrayI& inIndices, E_Int ir,
                                          E_Int i, E_Int j, E_Int k,
                                          FldArrayI& links)
{
  E_Int n = inIndices.getSize();
  E_Int im = _im;
  E_Int jm = _jm;
  E_Int ind2, ind3;

  E_Boolean isValid;
  E_Int dir; // direction of the line: 1=i, 2=j
  E_Int sens;
  E_Int ind = inIndices[ir];
  E_Int l = 1;
  
  if (i + 1 < im)
  {
    ind2 = ind+1;
    dir = 1;
    sens = 1;
      
    for (E_Int jr = 0; jr < n; jr++)
    {
      ind3 = inIndices[jr];
      if (ind3 == ind2)
      {
        isValid = testValidityOfLink( i, j, k, dir, sens);
        if ( isValid == true)
        {
          links(ir, l) = jr;
          l++;
        }
      }
    }
  }

  if (i - 1 >= 0)
  {
    ind2 = ind-1;
    dir = 1;
    sens = -1;
    for (E_Int jr = 0; jr < n; jr++)
    {
      ind3 = inIndices[jr];
      if (ind3 == ind2)
      {
        isValid = testValidityOfLink( i, j, k, dir, sens);
        if ( isValid == true )
        {
          links(ir, l) = jr;
          l++;
        }
      }
    }
  }
  
  if (j + 1 < jm)
  {
    ind2 = ind+im;
    dir = 2;
    sens = 1;
    for (E_Int jr = 0; jr < n; jr++)
    {
      ind3 = inIndices[jr];
      if (ind3 == ind2)
      {
        isValid = testValidityOfLink(i, j, k, dir, sens);
        if ( isValid == true )
        {
          links(ir, l) = jr;
          l++;
        }
      }    
    }
  }

  if (j - 1 >= 0)
  {
    ind2 = ind-im;
    dir = 2;
    sens = -1;
    for (E_Int jr = 0; jr < n; jr++)
    {
      ind3 = inIndices[jr];
      if (ind3 == ind2)
      {
        isValid = testValidityOfLink(i, j, k, dir, sens);
        if ( isValid == true )
        {
          links(ir, l) = jr;
          l++;
        }
      }
    }
  }  
}

//============================================================================
/* Test if the border node of index ind has to be set in a string or not    */
//============================================================================
void StructBlock::testForBlankedBorderNode(E_Int ind,
                                           E_Int indCell1,
                                           E_Int indCell2,
                                           E_Int& n )
{
  if ( _iblankCC[indCell1] != 0 || _iblankCC[indCell2] != 0 )
  {
    _pointsForStrings[n] = ind;
    n++;
  }
}

//===========================================================================
/* Identify pts for strings that are on a matching  boundary and whose
   matching pt is iblanked
*/
//===========================================================================
void
StructBlock::compMatchingBordersForString( E_Int noBlk, E_Int& n,
                                           vector<StructBlock*>& vectOfBlks)
{ 
  E_Int ind, indv;
  E_Float x, y, z;
  FldArrayI indBor;
  E_Int vectOfBlksSize = vectOfBlks.size();

  for (E_Int dir = 1; dir <= 4 ; dir++)
  {
    this->createArrayOfBorders(dir, indBor);
    
    for (E_Int i = 0; i < indBor.getSize(); i++)
    {
      ind = indBor[i];
      
      if (_iblank[ind] == 1)
      {
        x = _coord(ind, 1);
        y = _coord(ind, 2);
        z = _coord(ind, 3);
      
        E_Boolean match = false;
        E_Boolean validMatch = false;

        for (E_Int v = 0; v < vectOfBlksSize; v++)
        {
          if (v != _id)
          {
            FldArrayIS& iblank2 = vectOfBlks[v]->getIBlankArray();
         
            E_Boolean test =  
              searchForMatchingBnd2(x, y, z, vectOfBlks[v], indv);

            if ( test == true )
            {
              match = true;

              if ( iblank2[indv] == 1 )
              {
                validMatch = true;
                break;
              }
            }
          } 
        }

        if (match == true && validMatch == false)
        {
          _pointsForStrings[n] = ind; 
          n++;
        }
      }
    }
  }
}


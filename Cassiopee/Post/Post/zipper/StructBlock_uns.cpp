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

#include "StructBlock.h"

using namespace K_FUNC;
using namespace K_FLD;
using namespace K_CONST;

//=============================================================================
FldArrayI& StructBlock::getIBnd()
{
  return _ibnd;
}
//=============================================================================
E_Int StructBlock::getIBndSize()
{
  return _ibnd.getSize();
}
//=============================================================================
FldArrayI& StructBlock::getIDg()
{
  return _idg;
}
//=============================================================================
E_Int StructBlock::getIDgSize()
{
  return _idg.getSize();
}
//=============================================================================
FldArrayI& StructBlock::getUnsConnectEN()
{
  return _unsConnectEN;
}
//=============================================================================
E_Int StructBlock::getUnsConnectENSize()
{
  return _unsConnectEN.getSize();
}

//=============================================================================
/* Identify the points that are at the boundary of a StructBlock after blanking
   cells, namely Boundary Points
   They are of 2 types :
 Interior : unblanked interior point which has from 3 to 7 over 8 unblanked
   neighbours = unblanked interior points for stringing
 Exterior : unblanked exterior points
*/
//=============================================================================
void StructBlock::compIBnd()
{  
  E_Int im = _im;
  E_Int jm = _jm;
  E_Int n = 0 ;
  
  // Interior points
  //-----------------
  compIBndInt(im, jm, n);
  
  // Exterior points
  //----------------
  compIBndExt(im, jm, n);

  _ibnd.reAlloc(n);
}

//=============================================================================
/* Compute the interior points to be identified as Boundary Points
   To be a Boundary Point, the interior pt must be a Gap Point
*/
//=============================================================================
void StructBlock::compIBndInt(E_Int im, E_Int jm, E_Int& n)
{
  FldArrayI ind(9);
  E_Int cnt;
  for (E_Int j = 1; j < jm-1; j++)
    for (E_Int i = 1; i < im-1; i++)
    {
      ind[0] = i + j * im; // i,j
     
      // the current point must be unblanked
      if ( _iblank[ind[0]] != 0 )
      {
        ind[1] = (i-1) + (j-1) * im ; // i-1, j-1
        ind[2] =  i    + (j-1) * im ; // i,   j-1
        ind[3] = (i+1) + (j-1) * im ; // i+1, j-1
        ind[4] = (i+1) +  j    * im ; // i+1, j
        ind[5] = (i+1) + (j+1) * im ; // i+1, j+1
        ind[6] =  i    + (j+1) * im ; // i,   j+1
        ind[7] = (i-1) + (j+1) * im ; // i-1, j+1
        ind[8] = (i-1) +  j    * im ; // i-1, j
        
        // Count the number of neighbours that are unblanked
        cnt = 0;
        for (E_Int m = 1; m < ind.getSize(); m++ )
        {
          if ( _iblank[ind[m]] != 0 )
          {
            cnt++;
          }
        }
          
        // If the interior point is valid for stringing 
        // to be valid: at least 3 and no more than 7 neighbours are unblanked
        if ( ( cnt >= 3 ) && ( cnt < 8 ) ) 
        {
          _ibnd[n] = ind[0];
          n++;
        }
      }
    }
}

//=============================================================================
/* Compute the boundary points to be identified as Boundary Points
   To be a Boundary Point, the boundary point must have at least one blanked
   neighbor cell
*/
//=============================================================================
void StructBlock::compIBndExt(E_Int im, E_Int jm, E_Int& n)
{
  E_Int ind;

  // I = cste boundaries
  for (E_Int j = 0; j < jm; j++)
  {
    // i = 0 border 
    ind = j * im ; // 1,j
    if ( _iblank[ind] != 0)
    {
      _ibnd[n] = ind;
      n++;
    }
    
    // i = im-1 border
    ind = (im-1) + j * im ;// im,j
    if ( _iblank[ind] != 0)
    {
      _ibnd[n] = ind;
      n++;
    }
  }
  
  // J = cst boundaries
  for (E_Int i = 1; i < im-1; i++)
  {
    // j = 0 border 
    ind = i ; // i,0
    if ( _iblank[ind] != 0)
    { 
      _ibnd[n] = ind;
      n++;
    }
    
    // j = jm-1 border 
    ind = (jm-1) * im + i ; // i,jm
    if ( _iblank[ind] != 0)
    { 
      _ibnd[n] = ind;
      n++;
    }
  }
}
//=============================================================================
/* Update array _ibnd of StructBlock :
   if a node is in _ibnd AND is degenerated AND is not the reference node
   then its number is suppressed from _ibnd
   Then resize _ibnd.
*/
//=============================================================================
void StructBlock::updateIBnd()
{
  E_Int n = _ibnd.getSize();
  E_Int n2 = 0 ;
  for (E_Int i = 0 ; i < n ; i++)
  {
    if (_idg[_ibnd[i]] == _ibnd[i])
    {
      _ibnd[n2] = _ibnd[i] ;
      n2++;
    }    
  }
  _ibnd.reAlloc(n2);
}

//=============================================================================
/* Increment _ibnd in order to include it in global ibndG */
//=============================================================================
void StructBlock::incrIBnd(E_Int incr)
{
  E_Int n = _ibnd.getSize();
  for (E_Int i = 0; i < n; i++)
    _ibnd[i] = _ibnd[i] + incr;
}

//=============================================================================
/* Compute degenerations on the boundary of a StructBlock :
   if the distance between 2 nodes is lower than a fixed tolerance
   the nodes are considered as a unique node. */
//=============================================================================
void StructBlock::compIDg()
{
  E_Int im = _im;
  E_Int jm = _jm;

  E_Int ind1, ind2;
  E_Float x21, y21, z21;
  E_Float l;

  // i = 0 border 
  for (E_Int j = 0; j < jm-1; j++)
  {
    ind1 = j ;
    ind2 = j+1 ;
    x21 = _coord(ind2, 1) - _coord(ind1, 1);
    y21 = _coord(ind2, 2) - _coord(ind1, 2);
    z21 = _coord(ind2, 3) - _coord(ind1, 3);

    l = x21 * x21 + y21 * y21 + z21 * z21;
    l = sqrt(l);
    if (fEqualZero(l) == true)
      _idg[ind2] = _idg[ind1];
  }
  
  // j = 0 border
  for (E_Int i = 0; i < im-1; i++)
  {
    ind1 = i * jm;
    ind2 = ind1+jm; // (i+1)*jm
    x21 = _coord(ind2, 1) - _coord(ind1, 1);
    y21 = _coord(ind2, 2) - _coord(ind1, 2);
    z21 = _coord(ind2, 3) - _coord(ind1, 3);

    l = x21 * x21 + y21 * y21 + z21 * z21;
    l = sqrt(l);
    if (fEqualZero(l) == true)
      _idg[ind2] = _idg[ind1];
  }

  // i = im-1 border
  for (E_Int j = 0; j < jm-1; j++)
  {
    ind1 = (im-1) + j * im ;
    ind2 = ind1 + im; //(im-1) + (j+1) * im
    x21 = _coord(ind2, 1) - _coord(ind1, 1);
    y21 = _coord(ind2, 2) - _coord(ind1, 2);
    z21 = _coord(ind2, 3) - _coord(ind1, 3);

    l = x21 * x21 + y21 * y21 + z21 * z21;
    l = sqrt(l);
    if (fEqualZero(l) == true)
      _idg[ind2] = _idg[ind1];
  }
 
  // j = jm-1 border
  for (E_Int i = 0; i < im-1; i++)
  {
    ind1 = (jm-1) * im + i;
    ind2 =  ind1 + 1; //(jm-1) * im + i+1 ;
    x21 = _coord(ind2, 1) - _coord(ind1, 1);
    y21 = _coord(ind2, 2) - _coord(ind1, 2);
    z21 = _coord(ind2, 3) - _coord(ind1, 3);

    l = x21 * x21 + y21 * y21 + z21 * z21;
    l = sqrt(l);
    if (fEqualZero(l) == true)
      _idg[ind2] = _idg[ind1];
  }
}

//=============================================================================
/* Increment _idg in order to include it in global idgG
*/
//=============================================================================
void StructBlock::incrIDg(E_Int incr)
{
  E_Int n = _idg.getSize();
  for (E_Int i = 0; i < n; i++)
  {
    _idg[i] = _idg[i] + incr;
  }
}

//=============================================================================
/* Splits quadrilaterals whose iblankCC=1 of a block in two triangles.
   Quadrilaterals are cut on the smallest diagonal.
   Construct the connectivity Elements -> Nodes of a StructBlock.
*/
//=============================================================================
void StructBlock::compUnsConnectEN()
{
  E_Int ind1, ind2, ind3, ind4;
  E_Float x31, y31, z31, x42, y42, z42; 

  E_Int im = _im;
  E_Int im1 = im-1;
  E_Int jm = _jm;
  E_Int jm1 = jm-1;

  E_Int im1jm1 = im1*jm1;
  _unsConnectEN.malloc(2*im1jm1, 3);

  // Build connect by elements
  E_Int c = 0;
  for (E_Int j = 0; j < jm1; j++)
    for (E_Int i = 0; i < im1; i++)
    {
      if (_iblankCC[i+j*im1] == 1)
      {
        ind1 = i   +  j    * im;// i,j
        ind2 = i+1 +  j    * im;//i+1,j
        ind3 = i+1 + (j+1) * im;//i+1,j+1
        ind4 = i   + (j+1) * im;//i,j+1
        
        x31 = _coord(ind3, 1) - _coord(ind1, 1);
        y31 = _coord(ind3, 2) - _coord(ind1, 2);
        z31 = _coord(ind3, 3) - _coord(ind1, 3);
     
        x42 = _coord(ind4, 1) - _coord(ind2, 1);
        y42 = _coord(ind4, 2) - _coord(ind2, 2);
        z42 = _coord(ind4, 3) - _coord(ind2, 3);
        
        // First case : the cell is a triangle (degenerated quadrilateral)
        if ( _idg[ind1] == _idg[ind2] ) 
        {
           _unsConnectEN(c,1) = _idg[ind1];
           _unsConnectEN(c,2) = ind3;
           _unsConnectEN(c,3) = ind4;
           c++;    
        }
        else if ( _idg[ind1] == _idg[ind3] ) 
        {
          _unsConnectEN(c,1) = _idg[ind1];
          _unsConnectEN(c,2) = ind2;
          _unsConnectEN(c,3) = ind4;
          c++;
        }
        else if ( _idg[ind2] == _idg[ind4] ) 
        {
          _unsConnectEN(c,1) = ind1;
          _unsConnectEN(c,2) = _idg[ind2];
          _unsConnectEN(c,3) = ind3;
          c++;
        }
        else if ( _idg[ind3] == _idg[ind4] ) 
        {
          _unsConnectEN(c,1) = ind1;
          _unsConnectEN(c,2) = ind2;
          _unsConnectEN(c,3) = _idg[ind3];
          c++;
        }
        
        // Second case : the cell is a quadrilateral
        else
        {
          E_Float d1 = sqrt( x31*x31 + y31*y31 + z31*z31);
          E_Float d2 = sqrt( x42*x42 + y42*y42 + z42*z42);
          
          if (d1 < d2)
          {
            _unsConnectEN(c, 1) = ind1;
            _unsConnectEN(c, 2) = ind2;
            _unsConnectEN(c, 3) = ind3;
            c++;
            _unsConnectEN(c, 1) = ind1;
            _unsConnectEN(c, 2) = ind3;
            _unsConnectEN(c, 3) = ind4;
            c++;
          }
          else
          {
            _unsConnectEN(c, 1) = ind1;
            _unsConnectEN(c, 2) = ind2;
            _unsConnectEN(c, 3) = ind4;
            c++;
            _unsConnectEN(c, 1) = ind3;
            _unsConnectEN(c, 2) = ind2;
            _unsConnectEN(c, 3) = ind4;
            c++;
          }
        }
      }
    }
  _unsConnectEN.reAllocMat(c, 3);
}


//=============================================================================
/* Increment unsConnectEN in order to include it in global idgG
*/
//=============================================================================
void StructBlock::incrUnsConnectEN(E_Int incr)
{
  E_Int n = _unsConnectEN.getSize();
  for (E_Int j = 1 ; j < 4 ; j++)
    for (E_Int i = 0 ; i < n ; i++)
      _unsConnectEN(i,j) = _unsConnectEN(i,j) + incr;
}

//=============================================================================
/* Compute the boundary nodes of the StructBlock, the degenerations,
   and the triangular connectivity
*/
//=============================================================================
void StructBlock::compIBndIDgUnsConnectEN()
{
  compIBnd(); // _ibnd
  compIDg(); // _idg
  updateIBnd(); // _ibnd
  compUnsConnectEN(); // _unsConnectEN
}


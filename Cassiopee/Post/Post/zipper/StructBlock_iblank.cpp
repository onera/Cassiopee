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

#include "StructBlock.h"
#include "Intersection.h"

using namespace std;
using namespace K_FLD;
using namespace K_CONST;

//=============================================================================
/* Compute the iblank array for current block in the overlapping zone
   in: noBlk2 : number of the overlapping blk
*/
//=============================================================================
void StructBlock::compIBlank(vector<StructBlock*>& vectOfBlks, E_Int noBlk2)
{
  // Blanking Interpolable Pts
  E_Int nbI1 = _indirection[noBlk2];
  E_Int nbTotPrev1 = 0;
  E_Int ind, ind2;
  
  FldArrayIS& iblank2 = vectOfBlks[noBlk2]->getIBlankArray();
  
  for ( E_Int i = 0; i < noBlk2 ; i++)
    nbTotPrev1 = nbTotPrev1 + _indirection[i];
 
  for (E_Int n = nbTotPrev1; n < nbTotPrev1 + nbI1; n++)       
  {
    ind = _interpolableNode[n];
    ind2 = _interpolationNode[n];
    if ( ind != -1 && iblank2[ind2] != 0)
      _iblank[ind] = 0;

 //    if( _id==2 && noBlk2 ==4)
//       cout << "2,4 : ind "<< ind << " "<< ind2<<endl;
//     else if ( _id == 4 && noBlk2 == 2)
//       cout << "4,2 : ind "<< ind << " "<< ind2<<endl;
  }
}

//=============================================================================
/* Compute the iblankCC array: iblankCC = 0 for blanked cell center,
                                           1 otherwise
*/
//=============================================================================
void StructBlock::compIBlankCC()
{
  E_Int im1 = _im-1;
  E_Int jm1 = _jm-1;
  E_Int indCell;
  E_Int ind1, ind2, ind3, ind4;

  
  for (E_Int j = 0; j < jm1; j++)
    for (E_Int i = 0; i < im1; i++)
    {
      indCell = i + j*im1;
      ind1 = i + j*_im;
      ind2 = i+1 + j*_im;
      ind3 = i + (j+1)*_im;
      ind4 = i+1 + (j+1)*_im;
      
      if (_iblank[ind1]*_iblank[ind2]*_iblank[ind3]*_iblank[ind4] != 1)
        _iblankCC[indCell] = 0;
    }
}

//============================================================================
/* Test if the border node of index ind has to be set in a string or not    */
//============================================================================
void StructBlock::updateIBlankArray()
{
  E_Int im1 = _im-1;
  E_Int jm1 = _jm-1;
  
  E_Int indCell1, indCell2, indCell3, indCell4;
  E_Int ind;

  //interior points
  //---------------
  for (E_Int j = 2; j <= jm1; j++)
    for (E_Int i = 2; i <= im1; i++)
    {
      ind = (i-1) + (j-1) * _im; // index of the node (i,j)
      indCell1 = i-2 + (j-2) * im1; //index of the cell center (i-1,j-1)
      indCell2 = i-1 + (j-2) * im1; //index of the cell center (i,j-1)
      indCell3 = i-1 + (j-1) * im1; //index of the cell center (i,j)
      indCell4 = i-2 + (j-1) * im1;//index of the cell center (i-1,j)
      
      if ( _iblankCC[indCell1] == 0 && _iblankCC[indCell2] == 0 &&
           _iblankCC[indCell3] == 0 && _iblankCC[indCell4] == 0 )
        _iblank[ind] = 0;
    }

  // Corners
  //---------
  // 1st corner (1,1)
  ind = 0;
  indCell1 = 0;
  if ( _iblankCC[indCell1] == 0)
    _iblank[ind] = 0;

  // 2nd corner (1,jm)
  ind = (_jm-1) * _im;
  indCell1 = (jm1-1) * im1;
  if ( _iblankCC[indCell1] == 0)
    _iblank[ind] = 0;

  // 3rd corner (im,1)
  ind = _im-1;
  indCell1 = im1-1;
  if ( _iblankCC[indCell1] == 0)
    _iblank[ind] = 0;

  // last corner (im,jm)
  ind = im1 + jm1 * _im;
  indCell1 =  im1-1 + (jm1-1) * im1;
  if ( _iblankCC[indCell1] == 0)
    _iblank[ind] = 0;

  // bnd except corners

  // 1- Boundary i = 1
  //------------------
  for ( E_Int j = 2; j <= jm1; j++ )
  {
    ind = (j-1) * _im; // index of the node (i,j)
    indCell1 = (j-2) * im1; 
    indCell2 = (j-1) * im1; 
    if ( _iblankCC[indCell1] == 0 && _iblankCC[indCell2] == 0 )
      _iblank[ind] = 0;
  }

  // 2- Boundary i = imax
  //----------------------
  for ( E_Int j = 2; j <= jm1; j++)
  {
    ind = im1 + (j-1) * _im;
    indCell1 = im1-1 + (j-2) * im1;
    indCell2 = im1-1 + (j-1) * im1;
    if ( _iblankCC[indCell1] == 0 && _iblankCC[indCell2] == 0 )
      _iblank[ind] = 0;
  }

  // 3- Boundary j = 1
  //---------------------
  for ( E_Int i = 2; i <= im1; i++)
  {
    ind = i-1;
    indCell1 = i-2; 
    indCell2 = i-1;
    if ( _iblankCC[indCell1] == 0 && _iblankCC[indCell2] == 0 )
      _iblank[ind] = 0;
  }
  
  // 4- Boundary j = jm
  //----------------------
  for ( E_Int i = 2; i <= im1; i++)
  {
    ind = (i-1) + jm1 * _im;
    indCell1 = (i-2) + (jm1-1) * im1;
    indCell2 = (i-1) + (jm1-1) * im1;
    if ( _iblankCC[indCell1] == 0 && _iblankCC[indCell2] == 0 )
      _iblank[ind] = 0;
  }
}

//=============================================================================
/* Test the bounding box intersection of blk1 and blk2 */
//=============================================================================
E_Boolean StructBlock::testBBIntersection(StructBlock& block1,
                                          StructBlock& block2)
{
  E_Float xmax1, ymax1, zmax1, xmin1, ymin1, zmin1;
  E_Float xmin2, xmax2, ymin2, ymax2, zmin2, zmax2;
  
  block1.getBoundingBox(xmax1, ymax1, zmax1, xmin1, ymin1, zmin1);
  block2.getBoundingBox(xmax2, ymax2, zmax2, xmin2, ymin2, zmin2);

  if ( xmin1  <=  xmax2 && xmax1  >=  xmin2 &&
       ymin1  <=  ymax2 && ymax1  >=  ymin2 &&
       zmin1  <=  zmax2 && zmax1  >=  zmin2 )
    return true;
  
  else
    return false;
}

//=============================================================================
/* Compute interpolation information for current block :
   _interpolatedNode lists the interpolable nodes 
   (may be counted several times)
   _interpolationNode lists the interpolation nodes for interpolable node
   _indirection: says how many interpolable/interpolation pairs there are for
   each domain.

   Un point (noeud) est interpolable ssi :
   - il n'est pas en raccord coincident avec un autre bloc,
   - il est contenu dans la bounding box d'une cellule d'un autre bloc,
   - il se projete a une distance inferieur a ... dans cette cellule.
*/
//=============================================================================
void 
StructBlock::compInterpolationInformation(vector<StructBlock*>& vectOfBlks)
{
  E_Float tolx = _overlapTol;
  E_Float toly = _overlapTol;
  E_Float tolz = _overlapTol;
  E_Int cnt = 0;
 
  E_Int ind;
  E_Boolean test;
  E_Int i,j;
  E_Float x, y, z;
  //E_Float xmax1, ymax1, zmax1, xmin1, ymin1, zmin1;
  E_Float xmax2, ymax2, zmax2, xmin2, ymin2, zmin2;
  E_Int npts2;
  E_Int cellNF;
  E_Int nbOfBlks = vectOfBlks.size();
  
  _indirection.malloc( nbOfBlks);
  _indirection.setAllValuesAt(0);
  _interpolableNode.malloc(nbOfBlks*_nMeshPts);
  _interpolationNode.malloc(nbOfBlks*_nMeshPts);
  _interpolableNode.setAllValuesAtNull();
  _interpolationNode.setAllValuesAtNull();
  E_Int ovSize =  _overlappingBlks.size();
  for (E_Int v = 0; v < ovSize; v++)
  {
    StructBlock* blk2 = _overlappingBlks[v];
    npts2 =  blk2->getNbOfMeshPts();
    FldArrayF& cfd2 = blk2->getCfdField();
//     FldArrayF& coord2 = blk2->getCoord();

    for (ind = 0; ind < _nMeshPts; ind++)
    {
      x = _coord(ind,1);
      y = _coord(ind,2);
      z = _coord(ind,3);
      //xmax1 = _bbCell(ind,1);
      //ymax1 = _bbCell(ind,2);
      //zmax1 = _bbCell(ind,3);
      //xmin1 = _bbCell(ind,4);
      //ymin1 = _bbCell(ind,5);
      //zmin1 = _bbCell(ind,6);

      j = ind/_im + 1;
      i = ind - (j-1)*_im + 1;
      for (E_Int indi = 0; indi < npts2; indi++)
      {
        blk2->getBoundingBoxOfCell( indi, xmax2, ymax2, zmax2,
                                    xmin2, ymin2, zmin2);
        //E_Float xi = coord2(indi,1);
        //E_Float yi = coord2(indi,2);
        //E_Float zi = coord2(indi,3);
  
        // test if the center of the cell intersect s the bbox of the other cell
        if ( x  <=  xmax2 + tolx && x  >=  xmin2 - tolx &&
             y  <=  ymax2 + toly && y  >=  ymin2 - toly &&
             z  <=  zmax2 + tolz && z  >=  zmin2 - tolz)
        {
          cellNF = E_Int(cfd2(indi, _nfld));
         
          if (i == 1 || i == _im || j == 1 || j == _jm )
            test = searchForMatchingBnd(x, y, z, blk2);
          else test = false;
         
          if (cellNF != 0 && !test)
          { 
            if (testIntersectionOfCells(this, blk2, ind, indi, _overlapTol))
            {              
              _interpolableNode[cnt] = ind;
              _interpolationNode[cnt] = indi;
              _indirection[blk2->getBlockId()]++;
              cnt++;
              break;
            }
          }
        }
      }
    }
  }
  _interpolableNode.resize(cnt);
  _interpolationNode.resize(cnt);

}

//=============================================================================
// Construit la liste des blocs recouvrants le bloc courant : _overlappingBlks.
// Construit la liste des blocs en raccord coincident avec le bloc noBlk :
// _matchingBlks.
// 
// Un bloc j est en recouvrement avec le bloc noBlk ssi :
// - un des deux blocs au moins est chimere
// - leur bounding box s'intersecte
// - ils n'ont pas en commun uniquement un raccord coincident.
//
// Un bloc est dir en raccord coincident
// ssi ils ont en commun uniquement un raccord coincident.

// cellN must be last field of solution
//=============================================================================
void StructBlock::selectTypeOfBlks(vector<StructBlock*>& vectOfBlks)
{
  E_Int cellNF;
  E_Int noBlk = _id;
  assert(vectOfBlks[noBlk] == this);

//   cout << "For block "<<_id<<" : "<<endl;

  // Find if this block has interpolated or solid points
  E_Boolean isChimera = false;
  for (E_Int i = 0; i < _nMeshPts; i++ )
  {
    cellNF = E_Int(_cfdField(i, _nfld));
    if ( cellNF == 2  || cellNF == 0 )
    {
      isChimera = true;
      break;
    }
  }
  
  // Order vectOfBlks in :
  // _overlappingBlks : blocks that overlap with this
  // _matchingBlks : blocks with matching joins with this
  E_Int vectOfBlksSize = vectOfBlks.size();
  for ( E_Int n = 0; n < vectOfBlksSize; n++)
  {
    if (n != noBlk && testBBIntersection( *this, *vectOfBlks[n]) == true )
    {  
      if ( testMatchingBlks(this, vectOfBlks[n]) == false )
      {
        if ( isChimera == true )
          _overlappingBlks.push_back(vectOfBlks[n]);
        
        else 
        {
          // Is block n a chimera block
            FldArrayF& cfd2 = vectOfBlks[n]->getCfdField();
            for (E_Int i = 0; i < cfd2.getSize(); i++ )
            {
              cellNF = E_Int(cfd2(i, _nfld));
              if ( cellNF == 2 || cellNF == 0 )
              {
                _overlappingBlks.push_back(vectOfBlks[n]);
                break;
              }
            }
        }
      }
      
      else
      { 
        _overlappingBlks.push_back(vectOfBlks[n]);
        //_matchingBlks.push_back(vectOfBlks[n]);
      }
    }
  }

 //  for (E_Int l = 0; l < _overlappingBlks.size(); l++)
//   {
//     cout << "overlapped by "<<_overlappingBlks[l]->getBlockId()<<endl;
//   }
}

//=============================================================================
/* Compute the list of blks that have a bnd matching with one bnd of current
   blk. */
//=============================================================================
vector<StructBlock*>
StructBlock::compListOfMatchingBlks(vector<StructBlock*>& vectOfBlks)
{
  vector<StructBlock*> testBlks;
  E_Int im1 = _im;
  E_Int jm1 = _jm;
  FldArrayF& cfd1 =  _cfdField;
  E_Int ind1, cellNF;
  E_Float x1, y1, z1;
  E_Int vectOfBlksSize = vectOfBlks.size();
  for (E_Int v = 0 ; v < vectOfBlksSize; v++ )
  {
    if ( vectOfBlks[v] != this &&
         testBBIntersection( *this, *vectOfBlks[v]) == true )
    {
      // test frontiere j = 1
      for (E_Int i1 = 1; i1<= im1; i1++ )
      {
        ind1 = i1-1;
        cellNF = E_Int(cfd1( ind1, _nfld));
        if ( cellNF == 2 && i1 != 1 && i1 != im1)
          break;
        
        x1 = _coord(ind1, 1);
        y1 = _coord(ind1, 2);
        z1 = _coord(ind1, 3);
        
        if ( isAMatchingBnd(x1, y1, z1, vectOfBlks[v]) == true )
        {
          testBlks.push_back(vectOfBlks[v]);
          break;
        }
      }  
      
      // test frontiere j = jmax
      for (E_Int i1 = 1; i1 <= im1; i1++ )
      {
        ind1 = (i1-1) + (jm1-1) * im1;
        cellNF = E_Int(cfd1( ind1, _nfld));
        if ( cellNF == 2 && i1 != 1 && i1 != im1)
          break;
        
        x1 = _coord(ind1, 1);
        y1 = _coord(ind1, 2);
        z1 = _coord(ind1, 3);
      
        if ( isAMatchingBnd(x1, y1, z1, vectOfBlks[v]) == true )
        {
          testBlks.push_back(vectOfBlks[v]);
          break;
        }
      }
      
      // test frontiere i = 1
      for (E_Int j1 = 1; j1<= jm1; j1++ )
      {
        ind1 = (j1-1) * im1;
        cellNF = E_Int(cfd1( ind1, _nfld));
        if ( cellNF == 2 && j1 != 1 && j1 != jm1)
          break;
        
        x1 = _coord(ind1, 1);
        y1 = _coord(ind1, 2);
        z1 = _coord(ind1, 3);
        if ( isAMatchingBnd(x1, y1, z1, vectOfBlks[v]) == true )
        {
          testBlks.push_back(vectOfBlks[v]);
          break;
        }
      }
  
      // test frontiere i = imax
      for (E_Int j1 = 1; j1<= jm1; j1++ )
      {
        ind1 = (im1-1) + (j1-1) * im1;
        cellNF = E_Int(cfd1( ind1, _nfld));
        if ( cellNF == 2 && j1 != 1 && j1 != jm1)
          break;
        
        x1 = _coord(ind1, 1);
        y1 = _coord(ind1, 2);
        z1 = _coord(ind1, 3);
        if ( isAMatchingBnd(x1, y1, z1, vectOfBlks[v]) == true )
        {
          testBlks.push_back(vectOfBlks[v]);
          break;
        }
      }
    }
  }
  eliminateDoubleElts(testBlks);
  return testBlks;
}

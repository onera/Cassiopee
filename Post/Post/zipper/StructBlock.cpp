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

# include "StructBlock.h"
# include <stdlib.h>

using namespace std;
using namespace K_FLD;
using namespace K_FUNC;
using namespace K_CONST;

extern "C"
{
  void k6boundboxofstructcell_(E_Int* ind, E_Int& ni, 
                               E_Float* x, E_Float* y, E_Float* z,
                               E_Float& xmin, E_Float& xmax, 
                               E_Float& ymin, E_Float& ymax, 
                               E_Float& zmin, E_Float& zmax);
  void k6boundbox_(const E_Int& im, const E_Int& jm, const E_Int& km, 
                   const E_Float* x, const E_Float* y, const E_Float* z,
                   E_Float& xmax, E_Float& ymax, E_Float& zmax, 
                   E_Float& xmin, E_Float& ymin, E_Float& zmin);
}

//=============================================================================
/* Constructor from a field.
   field(ind, 1-3) : mesh coordinates
   field(ind, 4-nfld) : solution on mesh
   field(ind, nfld) : cell nature field.
*/
//=============================================================================
StructBlock::StructBlock(E_Int id, E_Int ni, E_Int nj, E_Int nk,
                         E_Int posx, E_Int posy, E_Int posz,
                         E_Float overlapTol, E_Float matchTol,
                         FldArrayF& field) :
  _cfdField(field, 4, field.getNfld())
{
  E_Int i, cellNF;

  _id = id;
  _nCfdFields = field.getNfld()-3;
  _nFields = field.getNfld();  
  _nfld = _cfdField.getNfld();
  
  //build coordinates
  FldArrayF coord(field.getSize(), 3);
  for (E_Int ind = 0; ind < field.getSize(); ind++)
  {
    coord(ind,1) = field(ind,posx);
    coord(ind,2) = field(ind,posy);
    coord(ind,3) = field(ind,posz);
  }

  _im = ni;
  _jm = nj;
  _km = nk;
  _coord = coord;

  _nMeshPts = ni*nj*nk;
  k6boundbox_(_im, _jm, _km, 
              _coord.begin(1), _coord.begin(2), _coord.begin(3),
              _xmax, _ymax, _zmax, _xmin, _ymin, _zmin);
  
  _globalField.malloc(_nMeshPts, 3+_cfdField.getNfld());

  for (i = 0; i < _nMeshPts; i++)
  {
    _globalField(i, 1) = _coord(i, 1);
    _globalField(i, 2) = _coord(i, 2);
    _globalField(i, 3) = _coord(i, 3);

    for (E_Int j = 0; j < _cfdField.getNfld(); j++)
      _globalField(i, j+4) = _cfdField(i, j+1);
  }
  
  _bbCell.malloc( _nMeshPts, 6 );
  E_Float xmax, xmin, ymax, ymin, zmax, zmin;
  for (i = 0; i < _nMeshPts; i++)
  {
    boundingBoxOfCell(i,
                      xmax, ymax, zmax,
                      xmin, ymin, zmin);
    _bbCell(i,1) = xmax;
    _bbCell(i,2) = ymax;
    _bbCell(i,3) = zmax;
    _bbCell(i,4) = xmin;
    _bbCell(i,5) = ymin;
    _bbCell(i,6) = zmin;
  }  
  _iblank.malloc(_nMeshPts);
  _iblank.setAllValuesAt(1);
  E_Int nCells = (ni-1)*(nj-1);
  _iblankCC.malloc(nCells);
  _iblankCC.setAllValuesAt(1);

  // Update iblank to 0 for masked pts
  for (E_Int ind = 0; ind < _nMeshPts ; ind++)
  {
    cellNF = E_Int(_cfdField(ind, _nfld));
    if ( cellNF == 0 )
      _iblank[ind] = 0;
  }
  
  _ibnd.malloc(field.getSize());
  _idg.malloc(field.getSize());

  // Init idg to the index of the point 
  for (i = 0 ; i < field.getSize() ; i++)
    _idg[i] = i;

  // Tolerances
  _matchTol = matchTol;//1.e-6;
  _overlapTol = overlapTol;//5.e-5;
}

//=============================================================================
// Destructor
//=============================================================================
StructBlock::~StructBlock()
{
}

//=============================================================================
E_Int StructBlock::getIm()
{
  return _im;
}
//=============================================================================
E_Int StructBlock::getJm()
{
  return _jm;
}
//=============================================================================
E_Int StructBlock::getKm()
{
  return _km;
}
//=============================================================================
FldArrayF& StructBlock::getCoord()
{
  return _coord;
}

//=============================================================================
E_Float StructBlock::getMatchTol()
{
  return _matchTol;
}
//=============================================================================
FldArrayF& StructBlock::getCfdField()
{
  return _cfdField;
}
//=============================================================================
FldArrayF& StructBlock::getGlobalField()
{
  return _globalField;
}

//=============================================================================
E_Int StructBlock::getNbOfCfdFields()
{
  return _nCfdFields;
}
//=============================================================================
E_Int StructBlock::getNbOfFields()
{
  return _nFields;
}
//=============================================================================
E_Int StructBlock::getNbOfMeshPts()
{
  return _nMeshPts;
}

//=============================================================================
/* Return the bounding box of mesh. */
//=============================================================================
void StructBlock::getBoundingBox(E_Float& xmax, E_Float& ymax, E_Float& zmax,
                                 E_Float& xmin, E_Float& ymin, E_Float& zmin)
{
  xmax = _xmax;
  ymax = _ymax;
  zmax = _zmax;
  xmin = _xmin;
  ymin = _ymin;
  zmin = _zmin;
}

//=============================================================================
/* Return the bounding box of cell ind of mesh. */
//=============================================================================
void StructBlock::getBoundingBoxOfCell(
  E_Int ind,
  E_Float& xmax, E_Float& ymax, E_Float& zmax,
  E_Float& xmin, E_Float& ymin, E_Float& zmin)
{
  xmax = _bbCell(ind, 1);
  ymax = _bbCell(ind, 2);
  zmax = _bbCell(ind, 3);
  xmin = _bbCell(ind, 4);
  ymin = _bbCell(ind, 5);
  zmin = _bbCell(ind, 6);
}
//=============================================================================
list<FldArrayI*>& StructBlock::getStrings()
{
  return _strings;
}
//=============================================================================
FldArrayI& StructBlock::getPointsForStrings()
{
  return _pointsForStrings;
}

//=============================================================================
FldArrayIS& StructBlock::getIBlankArray()
{
  return _iblank;
}

//=============================================================================
FldArrayI& StructBlock::getIndirectionArray()
{
  return _indirection;
}
//=============================================================================
FldArrayI& StructBlock::getInterpolableNodeArray()
{
  return _interpolableNode;
}

//=============================================================================
/* Return the number of interpolable pts with the block noBlk */
//=============================================================================
E_Int StructBlock::getNbOfInterpolablePtsForBlk(E_Int noBlk)
{
  return _indirection[noBlk];
}

//=============================================================================
/* Return the list of matching blocks of current block */
//=============================================================================
vector<StructBlock*>& StructBlock::getListOfMatchingBlks()
{
  return _matchingBlks;
}

//=============================================================================
/* Return the list of overlapping blocks of current block */
//=============================================================================
vector<StructBlock*>& StructBlock::getListOfOverlappingBlks()
{
  return _overlappingBlks;
}

//=============================================================================
/* Clean lists of blocks. */
//=============================================================================
void StructBlock::cleanListsOfBlks(vector<StructBlock*>& vectOfBlks)
{
  _overlappingBlks.clear();
  _matchingBlks.clear();
 
  _matchingBlks = compListOfMatchingBlks( vectOfBlks );
}

//=============================================================================
/* Update the list of overlapping blocks. Now overlapping blks are couple of
   blocks where computation of blanking is done.
 */
//=============================================================================
void StructBlock::addOverlappingBlk(StructBlock* overlappingBlk)
{
  _overlappingBlks.push_back( overlappingBlk );
}

//=============================================================================
/* Given a node of coordinates (x,y,z) of the current block,
   compute its distance to each cell vertex of the cell of ind indNode2 of 
   block noBlk2.
   Return the minimum distance between the node indNode1 and the cell vertices
   NOT USED*/
//=============================================================================
void StructBlock::computeDistanceOfNodeToCell(
  vector<StructBlock*>& vectOfBlks,
  E_Int ind,
  E_Float x1, E_Float y1, E_Float z1,
  E_Int noBlk2, E_Int indNode2,
  E_Float& dist)
{ 
  E_Int im2 = vectOfBlks[noBlk2]->getIm();
  E_Int jm2 = vectOfBlks[noBlk2]->getJm();
  FldArrayF& coord2 = vectOfBlks[noBlk2]->getCoord();

  E_Int j2 = indNode2 / im2 + 1;
  E_Int i2 = indNode2 - (j2-1) * im2 + 1;

  E_Int indi;
  FldArrayI i3(4);
  FldArrayI j3(4);
  
  E_Float x2, y2, z2;
  E_Float dx, dy, dz ;
  E_Float distmin;
  E_Int inci = 1;
  E_Int incj = 1;
  
  if ( i2 == im2)
    inci = -1;
  if ( j2 == jm2)
    incj = -1;
  
  i3[0] = i2;
  j3[0] = j2;
  i3[1] = i2 + inci;
  j3[1] = j2;
  i3[2] = i2 + inci;
  j3[2] = j2 + incj;
  i3[3] = i2;
  j3[3] = j2 + incj;

  distmin = E_MAX_FLOAT;
  
  for (E_Int m = 0; m < 4; m++)
  {
    indi = (i3[m]-1) + (j3[m]-1) * im2;
    
    x2 = coord2(indi, 1);
    y2 = coord2(indi, 2);
    z2 = coord2(indi, 3);
    dx = x2-x1;
    dy = y2-y1;
    dz = z2-z1;
          
    dist = dx*dx + dy*dy + dz*dz;
    dist = sqrt(dist);
    
    if ( dist < distmin )
      distmin = dist;
  }
}

//=============================================================================
/*
  Write a block and its solution in a unstructured tecplot file,
  taking iblank into account.
*/
//=============================================================================
void StructBlock::write(char* fileName, E_Boolean add)
{
  E_Int im = _im;
  E_Int im1 = im-1;
  E_Int jm = _jm;
  E_Int jm1 = jm-1;
  E_Int km = _km;
  E_Int km1;
  if (km == 1)
    km1 = 1;
  else
    km1 = km-1;
  E_Int imjmkm = im*jm*km;
  E_Int ind1, ind2, ind3, ind4;
  
  // Build field
  FldArrayF field(imjmkm, 3+_nfld);

  for (E_Int i = 0; i < imjmkm; i++)
  {
    field(i, 1) = _coord(i, 1);
    field(i, 2) = _coord(i, 2);
    field(i, 3) = _coord(i, 3);

    for (E_Int j = 0; j < _nfld; j++)
      field(i, j+4) = _cfdField(i, j+1);
  }
  
  // Build connect by elements
  FldArrayI connect( im1*jm1*km1, 4 );
  
  E_Int c = 0;
  for (E_Int k = 0; k < km1; k++)
    for (E_Int j = 0; j < jm1; j++)
      for (E_Int i = 0; i < im1; i++)
      {
        ind1 = i + j*im + k*im*jm;
        ind2 = i+1 + j*im + k*im*jm;
        ind3 = i + (j+1)*im + k*im*jm;
        ind4 = i+1 + (j+1)*im + k*im*jm;
        
        if (_iblank[ind1]*_iblank[ind2]*_iblank[ind3]*_iblank[ind4] == 1)
        {
          connect(c, 1) = ind1+1;
          connect(c, 2) = ind2+1;
          connect(c, 3) = ind4+1;
          connect(c, 4) = ind3+1;
          c++;
        }
      }
  connect.reAllocMat(c, 4);

  //if (c > 0)
  //  K_IO::GenIO::getInstance()->tpwriteQuads(
  //    fileName, field, connect, add);
  
}

//=============================================================================
void StructBlock::writeLine(FldArrayI& indices,
                            char* fileName, E_Boolean add)
{
  // Note : lines are defined by degenerated triangles
  E_Int ind;
  
  // Build field
  E_Int np = indices.getSize();
  
  FldArrayF field(np, 3+_nfld);

  for (E_Int i = 0; i < np; i++)
  {
    ind = indices[i];
    field(i, 1) = _coord(ind, 1);
    field(i, 2) = _coord(ind, 2);
    field(i, 3) = _coord(ind, 3);

    for (E_Int j = 0; j < _nfld; j++)
      field(i, j+4) = _cfdField(ind, j+1);
  }
  
  // Build connect by elements
  FldArrayI connect( np-1, 3 );
  
  E_Int c = 0;
  for (E_Int i = 0; i < np-1; i++)
  {     
    connect(c, 1) = i+1;
    connect(c, 2) = i+2;
    connect(c, 3) = i+2;
    c++;
  }

  //if (c > 0)
  //  K_IO::GenIO::getInstance()->tpwriteTriangles(
  //    fileName, field, connect, add);
}

//=============================================================================
/* Test if blk1 and blk2 have a matching boundary condition */
//=============================================================================
E_Boolean StructBlock::testMatchingBlks( StructBlock* blk1,
                                         StructBlock* blk2)
{
  E_Int im1 = blk1->getIm();
  E_Int jm1 = blk1->getJm();
  FldArrayF& coord1 = blk1->getCoord();
  FldArrayF& cfd1 =  blk1->getCfdField();
  E_Int ind1, cellNF;
  E_Float x1, y1, z1;

  // Compteur pour les raccords coincidents
  E_Int cnt = 0;

  // Test frontiere j = 1
  for (E_Int i1 = 1; i1 <= im1; i1++ )
  {
    ind1 = i1-1;
    cellNF = E_Int(cfd1( ind1, _nfld));
    if ( cellNF == 2 || cellNF == 0 )
      break;
   
    x1 = coord1(ind1, 1);
    y1 = coord1(ind1, 2);
    z1 = coord1(ind1, 3);
    if ( searchForMatchingBnd(x1, y1, z1, blk2) == true )
      cnt++;
  }
  
  // Test frontiere j = jmax
  for (E_Int i1 = 1; i1 <= im1; i1++ )
  {
    ind1 = (i1-1) + (jm1-1) * im1;
    cellNF = E_Int(cfd1( ind1, _nfld));
    if ( cellNF == 2 )
      break;
   
    x1 = coord1(ind1, 1);
    y1 = coord1(ind1, 2);
    z1 = coord1(ind1, 3);
    
    if ( searchForMatchingBnd(x1, y1, z1, blk2) == true )
      cnt++;
  }
  
  // Test frontiere i = 1
  for (E_Int j1 = 1; j1<= jm1; j1++ )
  {
    ind1 = (j1-1) * im1;
    cellNF = E_Int(cfd1( ind1, _nfld));
    if ( cellNF == 2 )
      break;
    x1 = coord1(ind1, 1);
    y1 = coord1(ind1, 2);
    z1 = coord1(ind1, 3);
    
    if ( searchForMatchingBnd(x1, y1, z1, blk2) == true )
      cnt++;
  }
  
  // Test frontiere i = imax
  for (E_Int j1 = 1; j1<= jm1; j1++ )
  {
    ind1 = (im1-1) + (j1-1) * im1;
    cellNF = E_Int(cfd1( ind1, _nfld));
    if ( cellNF == 2 )
      break;
    
    x1 = coord1(ind1, 1);
    y1 = coord1(ind1, 2);
    z1 = coord1(ind1, 3);
    if ( searchForMatchingBnd(x1, y1, z1, blk2) == true )
      cnt++;
  }
  
  if ( cnt == 0 )
    return false;
  else
    return true;
}

//=============================================================================
/* Projection du point (x,y,z) sur la frontière dir. 
   Retourne false si le projeté n'est pas situé sur cette frontiere, à
   un epsilon près.
*/
//=============================================================================
E_Boolean StructBlock::projectOrtho(E_Float x, E_Float y, E_Float z,
                                    E_Int dir)
{
  E_Float eps  = _matchTol;
 
  // changement de variables
  E_Float xx = sqrt(x*x + y*y);
  E_Float yy = z;
  E_Float xxh, yyh;
  E_Float xxm, yym, xxp, yyp;
  E_Float alphai, gammai, coefi;
  E_Float xm, ym, zm, xp, yp, zp;
  E_Float dist, distmin, dx, dy;
  E_Float xmin, xmax, ymin, ymax;
  E_Int ind, indp, imjm1;
  
  distmin = E_MAX_FLOAT;
 
  switch (dir)
  {
    case 1: // i = 1
      for (E_Int j = 0; j < _jm-1; j++)
      {
        ind = j * _im;
        indp = ind + _im;
        xm = _coord(ind,1);
        ym = _coord(ind,2);
        zm = _coord(ind,3);
        xp = _coord(indp,1);
        yp = _coord(indp,2);
        zp = _coord(indp,3);

        xxm = sqrt(xm*xm + ym*ym);
        yym = zm;
        xxp = sqrt(xp*xp + yp*yp);
        yyp = zp;   
        
        dx = E_abs(xxm-xxp);
        dy = E_abs(yym-yyp);
        xmin = E_min(xxm, xxp);
        xmax = E_max(xxm, xxp);
        ymin = E_min(yym, yyp);
        ymax = E_max(yym, yyp);
        
        if (dx > _matchTol && dy > _matchTol) //droite alphai X - Y + ci =0
        {
          alphai = (yyp-yym)/(xxp-xxm);
          gammai = yym - alphai * xxm;
          coefi = 1/(alphai*alphai + 1);
          xxh = coefi * (xx + alphai * yy - alphai * gammai);
          yyh = coefi * (alphai * xx + alphai*alphai * yy + gammai);          
        }
        else if (dx < _matchTol && dy > _matchTol) //droite X = cte
        {
          xxh = xxm;
          yyh = yy; 
        }
        else if (dx > _matchTol && dy < _matchTol) //droite Y = cte
        {
          xxh = xx;
          yyh = yym;
        }
        else
        {
          xxh = xxm;
          yyh = yym;
        }
        
        // teste si le projeté est dans le segment [xi,xi+1]
        if ( xxh > xmin - _matchTol && xxh < xmax + _matchTol &&
             yyh > ymin - _matchTol && yyh < ymax + _matchTol )
        {
          // distance entre P et son projeté
          dist = (xxh-xx)*(xxh-xx) + (yyh-yy)*(yyh-yy);
          
          if (dist < distmin) 
          {
            distmin = dist;
            if (distmin < eps) return true;
          }
        }
      }
      break;

    case 2 : // i = im
      for (E_Int j = 0; j < _jm-1; j++)
      {
        ind = _im-1 + j * _im;
        indp = ind + _im;
        xm = _coord(ind,1);
        ym = _coord(ind,2);
        zm = _coord(ind,3);
        xp = _coord(indp,1);
        yp = _coord(indp,2);
        zp = _coord(indp,3);

        xxm = sqrt(xm*xm + ym*ym);
        yym = zm;
        xxp = sqrt(xp*xp + yp*yp);
        yyp = zp;   
        
        dx = E_abs(xxm-xxp);
        dy = E_abs(yym-yyp);
        xmin = E_min(xxm,xxp);
        xmax = E_max(xxm,xxp);
        ymin = E_min(yym,yyp);
        ymax = E_max(yym,yyp);
        
        if (dx > _matchTol && dy > _matchTol)//droite alphai X - Y + ci =0
        {
          alphai = (yyp-yym)/(xxp-xxm);
          gammai = yym - alphai * xxm;
          coefi = 1/(alphai*alphai + 1);
          xxh = coefi * (xx + alphai * yy - alphai * gammai);
          yyh = coefi * (alphai * xx + alphai*alphai * yy + gammai);
        }
        else if (dx < _matchTol && dy > _matchTol)//droite X = cte
        {
          xxh = xxm;
          yyh = yy; 
        }
        else if (dx > _matchTol && dy < _matchTol)//droite Y = cte
        {
          xxh = xx;
          yyh = yym;
        }
        else
        {
          xxh = xxm;
          yyh = yym;
        }
        
        // teste si le projeté est dans le segment [xi,xi+1]
        if ( xxh > xmin - _matchTol && xxh < xmax + _matchTol &&
             yyh > ymin - _matchTol && yyh < ymax + _matchTol )
        {
          //distance entre P et son projeté
          dist = (xxh-xx)*(xxh-xx) + (yyh-yy)*(yyh-yy);
          
          if (dist < distmin) 
          {
            distmin = dist;
            if (distmin < eps) return true;
          }
        }
      }
      break;

    case 3 : // j = 1
      for (E_Int i = 0; i < _im-1; i++)
      {
        ind = i;
        indp = ind + 1;
        xm = _coord(ind,1);
        ym = _coord(ind,2);
        zm = _coord(ind,3);
        xp = _coord(indp,1);
        yp = _coord(indp,2);
        zp = _coord(indp,3);

        xxm = sqrt(xm*xm + ym*ym);
        yym = zm;
        xxp = sqrt(xp*xp + yp*yp);
        yyp = zp;   
        
        dx = E_abs(xxm-xxp);
        dy = E_abs(yym-yyp);
        xmin = E_min(xxm,xxp);
        xmax = E_max(xxm,xxp);
        ymin = E_min(yym,yyp);
        ymax = E_max(yym,yyp);

        if (dx > _matchTol && dy > _matchTol)//droite alphai X - Y + ci =0
        {
          alphai = (yyp-yym)/(xxp-xxm);
          gammai = yym - alphai * xxm;
          coefi = 1./(alphai*alphai + 1);
          xxh = coefi * (xx + alphai * yy - alphai * gammai);
          yyh = coefi * (alphai * xx + alphai*alphai * yy + gammai);
        }
        else if (dx < _matchTol && dy > _matchTol)//droite X = cte
        {
          xxh = xxm;
          yyh = yy; 
        }
        else if (dx > _matchTol && dy < _matchTol)//droite Y = cte
        {
          xxh = xx;
          yyh = yym;
        }
        else
        {
          xxh = xxm;
          yyh = yym;
        }
         
        //teste si le projeté est dans le segment [xi,xi+1]
        if ( xxh > xmin - _matchTol && xxh < xmax + _matchTol &&
             yyh > ymin - _matchTol && yyh < ymax + _matchTol )
        {
          //distance entre P et son projeté
          dist = (xxh-xx)*(xxh-xx) + (yyh-yy)*(yyh-yy);
          
          if (dist < distmin) 
          {
            distmin = dist;
            if (distmin < eps) return true;
          }
        }
      }
      break;

    case 4 : //j = jm
      imjm1 = (_jm-1)*_im;
      for (E_Int i = 0; i < _im-1; i++)
      {
        ind = i + imjm1;
        indp = ind + 1;
        xm = _coord(ind,1);
        ym = _coord(ind,2);
        zm = _coord(ind,3);
        xp = _coord(indp,1);
        yp = _coord(indp,2);
        zp = _coord(indp,3);

        xxm = sqrt(xm*xm + ym*ym);
        yym = zm;
        xxp = sqrt(xp*xp + yp*yp);
        yyp = zp;   
        
        dx = E_abs(xxm-xxp);
        dy = E_abs(yym-yyp);
        xmin = E_min(xxm,xxp);
        xmax = E_max(xxm,xxp);
        ymin = E_min(yym,yyp);
        ymax = E_max(yym,yyp);
        
        if (dx > _matchTol && dy > _matchTol)//droite alphai X - Y + ci =0
        {
          alphai = (yyp-yym)/(xxp-xxm);
          gammai = yym - alphai * xxm;
          coefi = 1/(alphai*alphai + 1);
          xxh = coefi * (xx + alphai * yy - alphai * gammai);
          yyh = coefi * (alphai * xx + alphai*alphai * yy + gammai);
        }
        else if (dx < _matchTol && dy > _matchTol)//droite X = cte
        {
          xxh = xxm;
          yyh = xx; 
        }
        else if (dx > _matchTol && dy < _matchTol)//droite Y = cte
        {
          xxh = yy;
          yyh = yym;
        }
        else
        {
          xxh = xxm;
          yyh = yym;
        }
        
        //teste si le projeté est dans le segment [xi,xi+1]
        if ( xxh > xmin - _matchTol && xxh < xmax + _matchTol &&
             yyh > ymin - _matchTol && yyh < ymax + _matchTol )
        {
          //distance entre P et son projeté
          dist = (xxh-xx)*(xxh-xx) + (yyh-yy)*(yyh-yy);
          
          if (dist < distmin) 
          {
            distmin = dist;
            if (distmin < eps) return true;
          }
        }
      }
      break;
    default:
      xxh = 0; yyh = 0;
      printf("Error: not a valid value for dir in projectOrtho.\n");
      exit(0);
  }
  return false;
}

//=============================================================================
/* Look if point x,y,z is on a match/nearmatch/nomatch bnd of blk2.
   And if the corresponding point is not blanked nor interpolated.
   Return true if it is a match/nearmatch/nomatch bnd */
//=============================================================================
E_Boolean StructBlock::searchForMatchingBnd(E_Float x, E_Float y, E_Float z,
                                            StructBlock* blk2)
{
  FldArrayF& cfd2 =  blk2->getCfdField();
  E_Int im2 = blk2->getIm();
  E_Int jm2 = blk2->getJm();
  FldArrayF& coord2 = blk2->getCoord();
  E_Float dx, dy, dz;
  E_Int ind2, cellNF;

  // TEST RACCORDS MATCH-NEARMATCH
  // Test frontiere j = 1 
  for (E_Int i2 = 1; i2 <= im2; i2++)
  { 
    ind2 = (i2-1);
    cellNF = E_Int(cfd2(ind2, _nfld));
    
    dx = coord2(ind2, 1) - x;
    dy = coord2(ind2, 2) - y;
    dz = coord2(ind2, 3) - z;
    dx = E_abs(dx);
    dy = E_abs(dy);
    dz = E_abs(dz);
    
    if ( dx < _matchTol && dy < _matchTol && dz < _matchTol && cellNF == 1)
      return true;
  }
  
  // Test frontiere j = jmax 
  for (E_Int i2 = 1; i2 <= im2; i2++)
  {
    ind2 = (i2-1) + (jm2-1) * im2;
    cellNF = E_Int(cfd2(ind2, _nfld));
    
    dx = coord2(ind2, 1) - x;
    dy = coord2(ind2, 2) - y;
    dz = coord2(ind2, 3) - z;
    dx = E_abs(dx);
    dy = E_abs(dy);
    dz = E_abs(dz);
 
    if ( dx < _matchTol && dy < _matchTol && dz < _matchTol && cellNF == 1)
      return true;
  }
    
  // Test frontiere i = 1 
  for (E_Int j2 = 1; j2 <= jm2; j2++)
  {
    ind2 = (j2-1)* im2;
    cellNF = E_Int(cfd2(ind2, _nfld));
    
    dx = coord2(ind2, 1) - x;
    dy = coord2(ind2, 2) - y;
    dz = coord2(ind2, 3) - z;
    dx = E_abs(dx);
    dy = E_abs(dy);
    dz = E_abs(dz);
    if ( dx < _matchTol && dy < _matchTol && dz < _matchTol && cellNF == 1)
      return true;
  }

  // Test frontiere i = imax 
  for (E_Int j2 = 1; j2 <= jm2; j2++)
  {
    ind2 = (im2-1) + (j2-1) * im2;
    cellNF = E_Int(cfd2(ind2, _nfld));
    
    dx = coord2(ind2, 1) - x;
    dy = coord2(ind2, 2) - y;
    dz = coord2(ind2, 3) - z;
    dx = E_abs(dx);
    dy = E_abs(dy);
    dz = E_abs(dz);
    if ( dx < _matchTol && dy < _matchTol && dz < _matchTol && cellNF == 1)
      return true;
  }

  // TEST RACCORD NOMATCH-NEARMATCH
  // Test frontiere i = 1 
  E_Int dir = 1;
  if ( blk2->projectOrtho(x, y, z, dir) == true )
    return true;

  // Test frontiere i = imax 
  dir = 2;
  if ( blk2->projectOrtho(x, y, z, dir) == true )
    return true;
 
  // Test frontiere j = 1
  dir = 3;
  if ( blk2->projectOrtho(x, y, z, dir) == true)
    return true;

  // Test frontiere j = jmax  
  dir = 4;
  if ( blk2->projectOrtho(x, y, z, dir) == true )
    return true;
  
  return false;
}
//=============================================================================
/* Look if point (x,y,z) is on a matching bnd of blk2.
   Return true if match is found */
//=============================================================================
E_Boolean StructBlock::searchForMatchingBnd2(E_Float x, E_Float y, E_Float z,
                                             StructBlock* blk2,
                                             E_Int& ind2)
{
  E_Int im2 = blk2->getIm();
  E_Int jm2 = blk2->getJm();
  FldArrayF& coord2 = blk2->getCoord();
  E_Float dx, dy, dz;

  // Test frontiere j = 1 
  for (E_Int i2 = 1; i2 <= im2; i2++)
  { 
    ind2 = (i2-1);
    
    dx = coord2(ind2, 1) - x;
    dy = coord2(ind2, 2) - y;
    dz = coord2(ind2, 3) - z;
    dx = E_abs(dx);
    dy = E_abs(dy);
    dz = E_abs(dz);
    
    if ( dx < _matchTol && dy < _matchTol && dz < _matchTol)
      return true;
  }
  
  // Test frontiere j = jmax 
  for (E_Int i2 = 1; i2 <= im2; i2++)
  {
    ind2 = (i2-1) + (jm2-1) * im2;
    
    dx = coord2(ind2, 1) - x;
    dy = coord2(ind2, 2) - y;
    dz = coord2(ind2, 3) - z;
    dx = E_abs(dx);
    dy = E_abs(dy);
    dz = E_abs(dz);
 
    if ( dx < _matchTol && dy < _matchTol && dz < _matchTol)
      return true;
  }
    
  // Test frontiere i = 1 
  for (E_Int j2 = 1; j2 <= jm2; j2++)
  {
    ind2 = (j2-1)* im2;
    
    dx = coord2(ind2, 1) - x;
    dy = coord2(ind2, 2) - y;
    dz = coord2(ind2, 3) - z;
    dx = E_abs(dx);
    dy = E_abs(dy);
    dz = E_abs(dz);
    if ( dx < _matchTol && dy < _matchTol && dz < _matchTol )
      return true;
  }

  // Test frontiere i = imax 
  for (E_Int j2 = 1; j2 <= jm2; j2++)
  {
    ind2 = (im2-1) + (j2-1) * im2;
    
    dx = coord2(ind2, 1) - x;
    dy = coord2(ind2, 2) - y;
    dz = coord2(ind2, 3) - z;
    dx = E_abs(dx);
    dy = E_abs(dy);
    dz = E_abs(dz);
    if ( dx < _matchTol && dy < _matchTol && dz < _matchTol )
      return true;
  }
  
  return false;
}
//=============================================================================
/* Check if there is a matching block between current block and interpolation
   block. In particular, pts can not be interpolable from the interpolation
   block if it's not from a much closer block 
   NOT USED */
//=============================================================================
void StructBlock::checkValidityOfInterpolationCell(
  E_Int noBlk,
  vector<StructBlock*>& vectOfBlks)
{
  E_Float eps = 1.e-10;
  E_Int nbI1 = 0;
  E_Int nbTotPrev1 = 0;
  E_Int ind1, ind2;
  E_Float dist1, dist2;
  E_Float x1, y1, z1;
  E_Boolean erased;
 
  vector<StructBlock*> testBlks = _overlappingBlks;
  vector<StructBlock*> matchingBlks = _matchingBlks;
  FldArrayI indir = _indirection;
  E_Int matchingBlksSize = matchingBlks.size();
  if ( matchingBlksSize == 0 )
    return;

  E_Int vectOfBlksSize =  vectOfBlks.size();
  E_Int testBlksSize =  testBlks.size();
  for (E_Int v = 0; v < vectOfBlksSize; v++)
  {
    erased = false;
    for (E_Int v1 = 0; v1 < testBlksSize; v1++)
    {
      if ( testBlks[v1] == vectOfBlks[v] )
      {
        nbI1 = indir[v];
        nbTotPrev1 = 0;
        
        for ( E_Int i = 0; i < v ; i++ )
          nbTotPrev1 = nbTotPrev1 + indir[i];
        
        if ( nbI1 != 0)
        {
          for (E_Int n = nbTotPrev1 ; n < nbTotPrev1 + nbI1; n++)       
          {
            ind1 = _interpolableNode[n];
            ind2 = _interpolationNode[n];
            x1 = _coord(ind1, 1);
            y1 = _coord(ind1, 2);
            z1 = _coord(ind1, 3);
            
            computeDistanceOfNodeToCell( vectOfBlks, ind1, x1, y1, z1,
                                         v, ind2, dist1);
            
            vector<StructBlock*>& matchingBlks2 =
              vectOfBlks[v]->getListOfMatchingBlks();
            E_Int matchingBlks2Size = matchingBlks2.size();
            for (E_Int vv = 0; vv < vectOfBlksSize; vv++)
            {
              for (E_Int v2 = 0; v2 < matchingBlksSize; v2++)
              {
                for (E_Int v3 = 0; v3 < matchingBlks2Size; v3++)
                {
                  if ( vectOfBlks[vv] == matchingBlks[v2] &&
                       vectOfBlks[vv] == matchingBlks2[v3] )
                  {
                    E_Int nPts3 = matchingBlks[v2]->getNbOfMeshPts();
                    for (E_Int ind3 = 0; ind3 < nPts3; ind3++)
                    {
                      computeDistanceOfNodeToCell( vectOfBlks, ind1,
                                                   x1, y1, z1,
                                                   vv, ind3, dist2);
                      if ( dist2 < dist1-eps )
                      {
                        erased = true;
                        goto nextBlk;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
    nextBlk:;
    if ( erased == true )
    { 
      for (E_Int ind = nbTotPrev1; ind < nbTotPrev1 + nbI1 ; ind++)
        _interpolableNode[ind] = -1;
      _indirection[v] = 0;
    }
  }
}

//=============================================================================
/* Given a vector of blks, eliminate double elements  */
//=============================================================================
void StructBlock::eliminateDoubleElts(vector<StructBlock*>& vectOfBlks)
{
  vector<StructBlock*> newVector;
  E_Int  vectOfBlksSize = vectOfBlks.size();
  for (E_Int i = 0; i < vectOfBlksSize; i++)
  {
    E_Boolean found = false;
    for (E_Int j = 0; j < i; j++)
    {
      if ( vectOfBlks[i] == vectOfBlks[j])
      {
        found = true;
        break;
      }
    }
    if (found == false)
    {
      newVector.push_back( vectOfBlks[i]);
    }
  }
  
  vectOfBlks.clear();
  vectOfBlks = newVector;
}

//=============================================================================
/* Search For matching bnds . Return true if match is found*/
//=============================================================================
E_Boolean StructBlock::isAMatchingBnd(E_Float x1, E_Float y1, E_Float z1,
                                      StructBlock* blk2)
{
  FldArrayF& cfd2 =  blk2->getCfdField();
  E_Int im2 = blk2->getIm();
  E_Int jm2 = blk2->getJm();
  FldArrayF& coord2 = blk2->getCoord();
  E_Float dx, dy, dz;
  E_Int ind2;
  E_Int cellNF;
  
  // test frontiere j = 1 
  for (E_Int i2 = 1; i2 <= im2; i2++)
  {
    ind2 = (i2-1);
    cellNF = E_Int(cfd2(ind2, _nfld));
    if ( cellNF == 2 && i2 != 1 && i2 != im2 )
      break;
    
    dx = coord2(ind2, 1) - x1;
    dy = coord2(ind2, 2) - y1;
    dz = coord2(ind2, 3) - z1;
    dx = E_abs(dx);
    dy = E_abs(dy);
    dz = E_abs(dz);
    if ( dx < _matchTol && dy < _matchTol && dz < _matchTol )
      return true;
  }
  
  // test frontiere j = jmax 
  for (E_Int i2 = 1; i2 <= im2; i2++)
  {
    ind2 = (i2-1) + (jm2-1) * im2;
    cellNF = E_Int(cfd2(ind2, _nfld));
    if ( cellNF == 2  && i2 != 1 && i2 != im2 )
      break;
    
    dx = coord2(ind2, 1) - x1;
    dy = coord2(ind2, 2) - y1;
    dz = coord2(ind2, 3) - z1;
    dx = E_abs(dx);
    dy = E_abs(dy);
    dz = E_abs(dz);
    if ( dx < _matchTol && dy < _matchTol && dz < _matchTol )
      return true;
  }
    
  // test frontiere i = 1 
  for (E_Int j2 = 1; j2 <= jm2; j2++)
  {
    ind2 = (j2-1)* im2;
    cellNF = E_Int(cfd2(ind2, _nfld));
    if ( cellNF == 2  && j2 != 1 && j2 != jm2)
      break;
    
    dx = coord2(ind2, 1) - x1;
    dy = coord2(ind2, 2) - y1;
    dz = coord2(ind2, 3) - z1;
    dx = E_abs(dx);
    dy = E_abs(dy);
    dz = E_abs(dz);
    if ( dx < _matchTol && dy < _matchTol && dz < _matchTol )
      return true;         
  }

  // test frontiere i = imax 
  for (E_Int j2 = 1; j2 <= jm2; j2++)
  {
    ind2 = (im2-1) + (j2-1) * im2;
    cellNF = E_Int(cfd2(ind2, _nfld));
    if ( cellNF == 2  && j2 != 1 && j2 != jm2)
      break;
    
    dx = coord2(ind2, 1) - x1;
    dy = coord2(ind2, 2) - y1;
    dz = coord2(ind2, 3) - z1;
    dx = E_abs(dx);
    dy = E_abs(dy);
    dz = E_abs(dz);
    if ( dx < _matchTol && dy < _matchTol && dz < _matchTol )
      return true;          
  }
  return false;
}
//============================================================================
/* Given a direction dir, return the array of indices of pts on the
   corresponding bnd
   dir = 1 : border i = 1  
         2 : border i = im
         3 : border j = 1
         4 : border j = jm
*/
//============================================================================
void StructBlock::createArrayOfBorders(E_Int dir, FldArrayI& indBord)
{
  E_Int im = _im;
  E_Int jm = _jm;
  indBord.malloc(im+jm);

  switch (dir)
  {
    case 1: // i = 1 
      for (E_Int j = 1; j <= jm ; j++)
        indBord[j-1] = (j-1) * im ;
      indBord.reAlloc(jm);
      break;
    case 2: // i = im 
      for (E_Int j = 1; j <= jm ; j++)
        indBord[j-1] = (im-1) + (j-1) * im ;
      indBord.reAlloc(jm);
      break;
    case 3: // j = 1
      for (E_Int i = 1; i <= im ; i++)
        indBord[i-1] = i-1;
      indBord.reAlloc(im);
      break;
    case 4:// j = jm
      for (E_Int i = 1; i <= im ; i++)
        indBord[i-1] = i-1 + (jm-1) * im;
      indBord.reAlloc(im);
      break;
    default:
      printf("Error: Not a good value for dir... Choose from 1 to 4.\n");
      exit(0);
      break;
  }
}

//=========================================================================
/* Given a segment [ind11, ind12] on a boundary of direction dir2 
   and a matching block blk2, add the segment [ind11, ind12] in the string
   list. Return false if the segment is not added in the string. 
   out : n : index of next element of the string.   
*/
//=========================================================================
E_Boolean StructBlock::addPtsInString(E_Int ind11, E_Int ind12, E_Int dir2,
                                      StructBlock* blk2, E_Int& n)
{
  E_Float x11, y11, z11;
  E_Float x12, y12, z12; 
  E_Float dx1, dy1, dz1;
  E_Float dx2, dy2, dz2;
  E_Int ind21, ind22;

  FldArrayF& coord2 = blk2->getCoord();
  FldArrayF& cfd2 = blk2->getCfdField();
  FldArrayIS& iblank2 = blk2->getIBlankArray();
  FldArrayI indBor2;

  blk2->createArrayOfBorders(dir2, indBor2);
  
  // peut etre a optimiser en le mettant en argument
  x11 = _coord(ind11, 1);
  y11 = _coord(ind11, 2);
  z11 = _coord(ind11, 3);
  x12 = _coord(ind12, 1);
  y12 = _coord(ind12, 2);
  z12 = _coord(ind12, 3);

  // WARNING test jusqu a imax-1
  for ( E_Int i2 = 0; i2 < indBor2.getSize()-1; i2++)
  {
    ind21 = indBor2[i2];
    ind22 = indBor2[i2+1];
    
    if ( cfd2(ind21, _nfld) < 2  && cfd2(ind22, _nfld) < 2 &&
         (iblank2[ind21] == 0 || iblank2[ind22] == 0 ) )
    {
      dx1 = coord2(ind21, 1) - x11;
      dy1 = coord2(ind21, 2) - y11;
      dz1 = coord2(ind21, 3) - z11;
      dx2 = coord2(ind21, 1) - x12;
      dy2 = coord2(ind21, 2) - y12;
      dz2 = coord2(ind21, 3) - z12;
      dx1 = E_abs(dx1);
      dy1 = E_abs(dy1);
      dz1 = E_abs(dz1);
      dx2 = E_abs(dx2);
      dy2 = E_abs(dy2);
      dz2 = E_abs(dz2);
      
      if ( dx1 < _matchTol && dy1 < _matchTol && dz1 < _matchTol )
      {
        dx2 = coord2(ind22, 1) - x12;
        dy2 = coord2(ind22, 2) - y12;
        dz2 = coord2(ind22, 3) - z12;
        dx2 = E_abs(dx2);
        dy2 = E_abs(dy2);
        dz2 = E_abs(dz2);
     
        if ( dx2 < _matchTol && dy2 < _matchTol && dz2 < _matchTol )
        {
          _pointsForStrings[n] = ind11;
          _pointsForStrings[n+1] = ind12;
          n = n+2;
          return true;
        }
      }
      
      else if ( dx2 < _matchTol && dy2 < _matchTol && dz2 < _matchTol )
      {
        dx1 = coord2(ind22, 1) - x11;
        dy1 = coord2(ind22, 2) - y11;
        dz1 = coord2(ind22, 3) - z11;
        dx1 = E_abs(dx1);
        dy1 = E_abs(dy1);
        dz1 = E_abs(dz1);
        if ( dx1 < _matchTol && dy1 < _matchTol && dz1 < _matchTol )
        {
          _pointsForStrings[n] = ind11;
          _pointsForStrings[n+1] = ind12;
          n = n+2;
          return true;
        }
      }
    }
  }
  return false;
}
//=============================================================================
void StructBlock::boundingBoxOfCell(
  E_Int ind,
  E_Float& xmax, E_Float& ymax, E_Float& zmax, 
  E_Float& xmin, E_Float& ymin, E_Float& zmin)
{ 
  E_Int imjm = _im*_jm;
  E_Int k = ind/imjm;
  E_Int j = ( ind - k * imjm )/_im;
  E_Int i = ind - j * _im + k * imjm;

  E_Int alpha = 1;
  E_Int beta  = 1;
  E_Int gamma = 1;
            
  if (i == _im-1) alpha = -1;
  
  if (j == _jm-1) beta = -1;
  
  if ( k == _km-1) gamma = -1;
  
  if ( _im == 1 ) alpha = 0;
  if ( _jm == 1 ) beta = 0;
  if ( _km == 1 ) gamma = 0;
  
  FldArrayI indtab(8);
  indtab[0] = ind;
  indtab[1] = (i+alpha) + j*_im + k*imjm;
  indtab[2] = (i+alpha) + (j+beta)*_im + k*imjm;
  indtab[3] = i + (j+beta)*_im + k*imjm;
  indtab[4] = i + j*_im + (k+gamma)*imjm;
  indtab[5] = (i+alpha) + j*_im + (k+gamma)*imjm;
  indtab[6] = (i+alpha) + (j+beta)*_im + (k+gamma)*imjm;  
  indtab[7] = i + (j+beta)*_im + (k+gamma)*imjm;
  
  E_Int size = _coord.getSize();
  k6boundboxofstructcell_( indtab.begin(), size,
                           _coord.begin(1), _coord.begin(2), _coord.begin(3), 
                           xmin, xmax, ymin, ymax, zmin, zmax);
}

//=============================================================================
/* Return the block Id */
//=============================================================================
E_Int StructBlock::getBlockId()
{
  return _id;
}

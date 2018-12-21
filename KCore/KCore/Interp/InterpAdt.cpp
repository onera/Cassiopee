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
# include "Interp/Interp.h"
# include <stack>
# include "CompGeom/compGeom.h"
using namespace std;
using namespace K_FLD;

struct stackData
{
    K_INTERP::IntTreeNode* current;
    E_Float xxmax, yymax, zzmax, ppmax, qqmax, rrmax;
    E_Float xxmin, yymin, zzmin, ppmin, qqmin, rrmin;
    E_Int coupure;
};

//=============================================================================
/* Destructor */
//=============================================================================
K_INTERP::InterpAdt::~InterpAdt()
{
  destroy();
}

//=============================================================================
/* Destroy data structure */
//=============================================================================
void K_INTERP::InterpAdt::destroy()
{
  // delete the ADT tree
  IntTreeNode* current = _tree;
  stack<stackData> stack;
  stackData dataForStack;

  dataForStack.current = current;
  stack.push(dataForStack);
  
  while (stack.size() != 0)
  {
    dataForStack = stack.top();
    stack.pop();
    current = dataForStack.current;
    
    if (current->_left != NULL)
    {
      dataForStack.current = current->_left;
      stack.push(dataForStack);
    }
    if (current->_right != NULL)
    {
      dataForStack.current = current->_right;
      stack.push(dataForStack);
    }
    delete current;
  }
   
  _tree = NULL;
}

//=============================================================================
/* Constructor */
//=============================================================================
K_INTERP::InterpAdt::InterpAdt(E_Int npts, 
                               E_Float* xD, E_Float* yD, E_Float* zD,
                               void* a1, void* a2, void* a3, E_Int& built) 
{
  if (a1 != NULL && a2 != NULL && a3 != NULL)// structure
  {
    E_Int ni = *(E_Int*)a1;
    E_Int nj = *(E_Int*)a2;
    E_Int nk = *(E_Int*)a3;
    built = buildStructAdt(ni, nj, nk, xD, yD, zD);
  }
  else //non structure
  {
    FldArrayI& cEV = *(FldArrayI*)a1;
    built = buildUnstrAdt(npts, cEV, xD, yD, zD);
  }
}

//=============================================================================
/* Construction de l'adt pour un maillage structure */
//=============================================================================
E_Int K_INTERP::InterpAdt::buildStructAdt(E_Int ni, E_Int nj, E_Int nk,
                                          E_Float* x, E_Float* y, E_Float* z)
{  
  _tree = NULL;
 
  K_COMPGEOM::boundingBox(ni*nj*nk, x, y, z, 
                          _xmin, _ymin, _zmin, 
                          _xmax, _ymax, _zmax);
  if (nk == 1) // traitement 2D
  {
    if (K_FUNC::E_abs(_zmax-_zmin) < K_CONST::E_GEOM_CUTOFF) //cas 1 plan en k=zmin = OK
    {
      _zmax = _zmin + 1.;
    }
    else //cas non plan !!!!!
      return 0;
  }
  E_Int nij = ni*nj;
  /* Insert all points (cells) in the tree */
  E_Float xmax,ymax,zmax,xmin,ymin,zmin;
  E_Int ind, ind2;
  
  if (nk > 1)
  { 
    for (E_Int k = 0; k < nk-1; k++)
      for (E_Int j = 0; j < nj-1; j++)
        for (E_Int i = 0; i < ni-1; i++)
        {
          ind = i+ni*j+nij*k;
          xmin = +K_CONST::E_INFINITE;
          ymin = +K_CONST::E_INFINITE;
          zmin = +K_CONST::E_INFINITE;
          xmax = -K_CONST::E_INFINITE;
          ymax = -K_CONST::E_INFINITE;
          zmax = -K_CONST::E_INFINITE;
    
          xmin = K_FUNC::E_min(xmin, x[ind]);
          ymin = K_FUNC::E_min(ymin, y[ind]);
          zmin = K_FUNC::E_min(zmin, z[ind]);
          xmax = K_FUNC::E_max(xmax, x[ind]);
          ymax = K_FUNC::E_max(ymax, y[ind]);
          zmax = K_FUNC::E_max(zmax, z[ind]);
        
          /* Compute bounding box of cell */
          ind2 = ind+1;
          xmin = K_FUNC::E_min(xmin, x[ind2]);
          ymin = K_FUNC::E_min(ymin, y[ind2]);
          zmin = K_FUNC::E_min(zmin, z[ind2]);
          xmax = K_FUNC::E_max(xmax, x[ind2]);
          ymax = K_FUNC::E_max(ymax, y[ind2]);
          zmax = K_FUNC::E_max(zmax, z[ind2]);
            
          ind2 = ind+ni;
          xmin = K_FUNC::E_min(xmin, x[ind2]);
          ymin = K_FUNC::E_min(ymin, y[ind2]);
          zmin = K_FUNC::E_min(zmin, z[ind2]);
          xmax = K_FUNC::E_max(xmax, x[ind2]);
          ymax = K_FUNC::E_max(ymax, y[ind2]);
          zmax = K_FUNC::E_max(zmax, z[ind2]);

          ind2 = ind+nij;
          xmin = K_FUNC::E_min(xmin, x[ind2]);
          ymin = K_FUNC::E_min(ymin, y[ind2]);
          zmin = K_FUNC::E_min(zmin, z[ind2]);
          xmax = K_FUNC::E_max(xmax, x[ind2]);
          ymax = K_FUNC::E_max(ymax, y[ind2]);
          zmax = K_FUNC::E_max(zmax, z[ind2]);
        
          ind2 = ind+ni+1;
          xmin = K_FUNC::E_min(xmin, x[ind2]);
          ymin = K_FUNC::E_min(ymin, y[ind2]);
          zmin = K_FUNC::E_min(zmin, z[ind2]);
          xmax = K_FUNC::E_max(xmax, x[ind2]);
          ymax = K_FUNC::E_max(ymax, y[ind2]);
          zmax = K_FUNC::E_max(zmax, z[ind2]);
        
          ind2 = ind+nij+1;
          xmin = K_FUNC::E_min(xmin, x[ind2]);
          ymin = K_FUNC::E_min(ymin, y[ind2]);
          zmin = K_FUNC::E_min(zmin, z[ind2]);
          xmax = K_FUNC::E_max(xmax, x[ind2]);
          ymax = K_FUNC::E_max(ymax, y[ind2]);
          zmax = K_FUNC::E_max(zmax, z[ind2]);
        
          ind2 = ind+ni+nij;
          xmin = K_FUNC::E_min(xmin, x[ind2]);
          ymin = K_FUNC::E_min(ymin, y[ind2]);
          zmin = K_FUNC::E_min(zmin, z[ind2]);
          xmax = K_FUNC::E_max(xmax, x[ind2]);
          ymax = K_FUNC::E_max(ymax, y[ind2]);
          zmax = K_FUNC::E_max(zmax, z[ind2]);

          ind2 = ind2+1;
          xmin = K_FUNC::E_min(xmin, x[ind2]);
          ymin = K_FUNC::E_min(ymin, y[ind2]);
          zmin = K_FUNC::E_min(zmin, z[ind2]);
          xmax = K_FUNC::E_max(xmax, x[ind2]);
          ymax = K_FUNC::E_max(ymax, y[ind2]);
          zmax = K_FUNC::E_max(zmax, z[ind2]);

          /* insert this cell in ADT tree */
          insert(ind, xmin, ymin, zmin,
                 xmax, ymax, zmax);
        }
  }
  else //2D avec un plan en z constant 
  {
    for (E_Int j = 0; j < nj-1; j++)
      for (E_Int i = 0; i < ni-1; i++)
      {
        ind = i+ni*j;
        xmin = +K_CONST::E_INFINITE;
        ymin = +K_CONST::E_INFINITE;
        zmin = +K_CONST::E_INFINITE;
        xmax = -K_CONST::E_INFINITE;
        ymax = -K_CONST::E_INFINITE;
        zmax = -K_CONST::E_INFINITE;
    
        xmin = K_FUNC::E_min(xmin, x[ind]);
        ymin = K_FUNC::E_min(ymin, y[ind]);
        xmax = K_FUNC::E_max(xmax, x[ind]);
        ymax = K_FUNC::E_max(ymax, y[ind]);
        
        /* Compute bounding box of cell */
        ind2 = ind+1;
        xmin = K_FUNC::E_min(xmin, x[ind2]);
        ymin = K_FUNC::E_min(ymin, y[ind2]);
        xmax = K_FUNC::E_max(xmax, x[ind2]);
        ymax = K_FUNC::E_max(ymax, y[ind2]);
            
        ind2 = ind+ni;
        xmin = K_FUNC::E_min(xmin, x[ind2]);
        ymin = K_FUNC::E_min(ymin, y[ind2]);
        xmax = K_FUNC::E_max(xmax, x[ind2]);
        ymax = K_FUNC::E_max(ymax, y[ind2]);
   
        ind2 = ind+ni+1;
        xmin = K_FUNC::E_min(xmin, x[ind2]);
        ymin = K_FUNC::E_min(ymin, y[ind2]);
        xmax = K_FUNC::E_max(xmax, x[ind2]);
        ymax = K_FUNC::E_max(ymax, y[ind2]);
               
        /* insert this cell in ADT tree */
        insert(ind, xmin, ymin, _zmin, xmax, ymax, _zmax);
      }
  }
  return 1;
}
//=============================================================================
/* Construction de l'adt pour un maillage non structure (TETRA) */
//=============================================================================
E_Int K_INTERP::InterpAdt::buildUnstrAdt(E_Int npts, FldArrayI& connect, 
                                         E_Float* x, E_Float* y, E_Float* z)
{  
  _tree = NULL;
  K_COMPGEOM::boundingBox(npts, x, y, z, 
                          _xmin, _ymin, _zmin, 
                          _xmax, _ymax, _zmax);
  /* Insert all points (cells) in the tree */
  E_Float xmax,ymax,zmax,xmin,ymin,zmin;
    E_Int nelts = connect.getSize();
  E_Int ind, ind2;
  E_Int* c1 = connect.begin(1);
  E_Int* c2 = connect.begin(2);
  E_Int* c3 = connect.begin(3);
  E_Int* c4 = connect.begin(4);

  for (E_Int et = 0; et < nelts; et++)
  {
    xmin = +K_CONST::E_INFINITE;
    ymin = +K_CONST::E_INFINITE;
    zmin = +K_CONST::E_INFINITE;
    xmax = -K_CONST::E_INFINITE;
    ymax = -K_CONST::E_INFINITE;
    zmax = -K_CONST::E_INFINITE;

    /* Compute bounding box of cell */
    ind = c1[et]-1;
    xmin = K_FUNC::E_min(xmin, x[ind]);
    ymin = K_FUNC::E_min(ymin, y[ind]);
    zmin = K_FUNC::E_min(zmin, z[ind]);
    xmax = K_FUNC::E_max(xmax, x[ind]);
    ymax = K_FUNC::E_max(ymax, y[ind]);
    zmax = K_FUNC::E_max(zmax, z[ind]);
    
    ind2 = c2[et]-1;
    xmin = K_FUNC::E_min(xmin, x[ind2]);
    ymin = K_FUNC::E_min(ymin, y[ind2]);
    zmin = K_FUNC::E_min(zmin, z[ind2]);
    xmax = K_FUNC::E_max(xmax, x[ind2]);
    ymax = K_FUNC::E_max(ymax, y[ind2]);
    zmax = K_FUNC::E_max(zmax, z[ind2]);
    
    ind2 = c3[et]-1;
    xmin = K_FUNC::E_min(xmin, x[ind2]);
    ymin = K_FUNC::E_min(ymin, y[ind2]);
    zmin = K_FUNC::E_min(zmin, z[ind2]);
    xmax = K_FUNC::E_max(xmax, x[ind2]);
    ymax = K_FUNC::E_max(ymax, y[ind2]);
    zmax = K_FUNC::E_max(zmax, z[ind2]);
    
    ind2 = c4[et]-1;
    xmin = K_FUNC::E_min(xmin, x[ind2]);
    ymin = K_FUNC::E_min(ymin, y[ind2]);
    zmin = K_FUNC::E_min(zmin, z[ind2]);
    xmax = K_FUNC::E_max(xmax, x[ind2]);
    ymax = K_FUNC::E_max(ymax, y[ind2]);
    zmax = K_FUNC::E_max(zmax, z[ind2]);
        
    /* insert this cell in ADT tree */
    insert(et, xmin, ymin, zmin, xmax, ymax, zmax);
  }
  return 1;
}

//=============================================================================
// Insert cell in ADT tree
//=============================================================================
void K_INTERP::InterpAdt::insert(E_Int ind,                                  
                                 E_Float xmin, E_Float ymin, E_Float zmin,
                                 E_Float xmax, E_Float ymax, E_Float zmax)
{
  E_Int coupure = 0;
  E_Float xxmax,yymax,zzmax,xxmin,yymin,zzmin;
  E_Float ppmax,qqmax,rrmax,ppmin,qqmin,rrmin;
  
  xxmin = _xmin;  // boite de recherche: initialement la bounding box du maillage
  yymin = _ymin;
  zzmin = _zmin;
  ppmin = _xmin;
  qqmin = _ymin;
  rrmin = _zmin;
  xxmax = _xmax;
  yymax = _ymax;
  zzmax = _zmax;
  ppmax = _xmax;
  qqmax = _ymax;
  rrmax = _zmax;
  
  IntTreeNode** current = &_tree;

  while (*current != NULL)
  {
    /* descent */
    switch (coupure)
    {
      // dichotomie direction 1
      case 0:

        xxmax = K_CONST::ONE_HALF*(xxmax+xxmin);
         
        if (xmin <= xxmax)
        {
          current = &(*current)->_left;
        }
        else
        {
          xxmax = 2*xxmax-xxmin;
          xxmin = K_CONST::ONE_HALF*(xxmax+xxmin);
          current = &(*current)->_right;
        }
         
        break;
        
        // dichotomie direction 2  
      case 1:
            
        yymax = K_CONST::ONE_HALF*(yymax+yymin);
            
        if (ymin <= yymax)
        {
          current = &(*current)->_left;
        }
        else
        {
          yymax = 2*yymax-yymin;
          yymin = K_CONST::ONE_HALF*(yymax+yymin);
          current = &(*current)->_right;
        }
           
        break;

        // dichotomie direction 3     
      case 2:
        zzmax = K_CONST::ONE_HALF*(zzmax+zzmin);

        if (zmin <= zzmax)
        {
          current = &(*current)->_left;
        }
        else
        {
          zzmax = 2*zzmax-zzmin;
          zzmin = K_CONST::ONE_HALF*(zzmax+zzmin);
          current = &(*current)->_right;
        }
         
        break;
        // dichotomie direction 4     
      case 3:

        ppmax = K_CONST::ONE_HALF*(ppmax+ppmin);
            
        if (xmax <= ppmax)
        {
          current = &(*current)->_left;
        }
        else
        {
          ppmax = 2*ppmax-ppmin;
          ppmin = K_CONST::ONE_HALF*(ppmax+ppmin);
          current = &(*current)->_right;
        }
           
        break;

        // dichotomie direction 5     
      case 4:
         
        qqmax = K_CONST::ONE_HALF*(qqmax+qqmin);
            
        if (ymax <= qqmax)
        {
          current = &(*current)->_left;
        }
        else
        {
          qqmax = 2*qqmax-qqmin;
          qqmin = K_CONST::ONE_HALF*(qqmax+qqmin);
          current = &(*current)->_right;
        }
          
        break;
    
        // dichotomie direction 6     
      case 5:
           
        rrmax = K_CONST::ONE_HALF*(rrmax+rrmin);
            
        if (zmax <= rrmax)
        {
          current = &(*current)->_left;
        }
        else
        {
          rrmax = 2*rrmax-rrmin;
          rrmin = K_CONST::ONE_HALF*(rrmax+rrmin);
          current = &(*current)->_right;
        }

        break;

    }
    coupure++;     if (coupure > 5) coupure = 0;  
  } /* while */
       
  *current = new IntTreeNode(ind, xmax, ymax, zmax,
                             xmin, ymin, zmin);
}
//=============================================================================
/* Recherche de la liste des cellules candidates. 
   IN: x,y,z: coord du pt a interpoler ou a extrpoler
   IN: alphaTol: tolerance pour l'extrapolation en nbre de bbox.
   Si alphaTol=0., une cellule est candidate si sa bbox contient x,y,z
   Si alphaTol=1., une cellule est candidate si sa bbox doublee a gauche et 
   a droite contient x,y,z.
   OUT: listOfCandidateCells: la liste des indices des cellules candidates.
   Retourne la taille de listOfCandidateCells.
*/
//=============================================================================
E_Int K_INTERP::InterpAdt::getListOfCandidateCells(
  E_Float x, E_Float y, E_Float z,
  list<E_Int>& listOfCandidateCells, E_Float alphaTol)
{
  E_Float dx, dy, dz;
  dx = 0.1*alphaTol*(_xmax-_xmin);
  dy = 0.1*alphaTol*(_ymax-_ymin);
  dz = 0.1*alphaTol*(_zmax-_zmin);
  if (x < _xmin - dx) return 0;
  if (y < _ymin - dy) return 0;
  if (z < _zmin - dz) return 0;
  if (x > _xmax + dx) return 0;
  if (y > _ymax + dy) return 0;
  if (z > _zmax + dz) return 0;
  if (_tree == NULL) return 0;
  E_Float a1, a2, a3, a4, a5, a6;
  E_Float b1, b2, b3, b4, b5, b6;
  E_Float xmax, ymax, zmax, xmin, ymin, zmin;
  stack<stackData> stack;
  stackData dataForStack;
  a1 = _xmin; a2 = _ymin; a3 = _zmin;
  a4 = x; a5 = y; a6 = z;
   
  b1 = x; b2 = y; b3 = z;
  b4 = _xmax; b5 = _ymax; b6 = _zmax;

  K_INTERP::IntTreeNode* current;
  current = _tree;

  /* Get down the tree */
  E_Int coupure = 0;
  E_Int coupure2;
   
  E_Float xxmax, yymax, zzmax, xxmin, yymin, zzmin;
  E_Float ppmax, qqmax, rrmax, ppmin, qqmin, rrmin;

  xxmin = _xmin;    // boite de recherche
  yymin = _ymin;
  zzmin = _zmin;
  ppmin = _xmin;
  qqmin = _ymin;
  rrmin = _zmin;
  xxmax = _xmax;
  yymax = _ymax;
  zzmax = _zmax;
  ppmax = _xmax;
  qqmax = _ymax;
  rrmax = _zmax;
  dataForStack.current = current;
  dataForStack.xxmax = xxmax;
  dataForStack.yymax = yymax;
  dataForStack.zzmax = zzmax;
  dataForStack.ppmax = ppmax;
  dataForStack.qqmax = qqmax;
  dataForStack.rrmax = rrmax;
  dataForStack.xxmin = xxmin;
  dataForStack.yymin = yymin;
  dataForStack.zzmin = zzmin;
  dataForStack.ppmin = ppmin;
  dataForStack.qqmin = qqmin;
  dataForStack.rrmin = rrmin;
  dataForStack.coupure = coupure;
  
  stack.push(dataForStack);
   
  while (stack.size() != 0)
  {
    dataForStack = stack.top();
    stack.pop();

    current = dataForStack.current;
    xxmin = dataForStack.xxmin;
    yymin = dataForStack.yymin;
    zzmin = dataForStack.zzmin;
    ppmin = dataForStack.ppmin;
    qqmin = dataForStack.qqmin;
    rrmin = dataForStack.rrmin;
    xxmax = dataForStack.xxmax;
    yymax = dataForStack.yymax;
    zzmax = dataForStack.zzmax;
    ppmax = dataForStack.ppmax;
    qqmax = dataForStack.qqmax;
    rrmax = dataForStack.rrmax;
    coupure = dataForStack.coupure;
    coupure2 = coupure+1;
    if (coupure2 > 5) coupure2 = 0;
    
    /* examine node */
    current->getCellBB(xmax, ymax, zmax, xmin, ymin, zmin);
    dx = alphaTol*(xmax-xmin)+K_CONST::E_GEOM_CUTOFF;
    dy = alphaTol*(ymax-ymin)+K_CONST::E_GEOM_CUTOFF;
    dz = alphaTol*(zmax-zmin)+K_CONST::E_GEOM_CUTOFF;
    /*
    if (xmin>a1-K_CONST::E_GEOM_CUTOFF && xmin<b1+K_CONST::E_GEOM_CUTOFF &&
        ymin>a2-K_CONST::E_GEOM_CUTOFF && ymin<b2+K_CONST::E_GEOM_CUTOFF &&
        zmin>a3-K_CONST::E_GEOM_CUTOFF && zmin<b3+K_CONST::E_GEOM_CUTOFF &&
        xmax>a4-K_CONST::E_GEOM_CUTOFF && xmax<b4+K_CONST::E_GEOM_CUTOFF &&
        ymax>a5-K_CONST::E_GEOM_CUTOFF && ymax<b5+K_CONST::E_GEOM_CUTOFF &&
        zmax>a6-K_CONST::E_GEOM_CUTOFF && zmax<b6+K_CONST::E_GEOM_CUTOFF)
    {
      listOfCandidateCells.push_back(current->_ind);
    }
    */
    if (x < xmax+dx && x > xmin-dx &&
        y < ymax+dy && y > ymin-dy &&
        z < zmax+dz && z > zmin-dz)
    {
      listOfCandidateCells.push_back(current->_ind);
    }

    switch (coupure)
    {
      // dichotomie direction 1
      case 0:
        xxmax = K_CONST::ONE_HALF*(xxmax+xxmin);
          
        if (a1 <= xxmax)
        {
          dataForStack.current = current->_left;
          dataForStack.xxmax = xxmax;
          dataForStack.coupure = coupure2;
          if (dataForStack.current != NULL) stack.push(dataForStack);
        }
            
        xxmax = 2*xxmax-xxmin;
        xxmin = K_CONST::ONE_HALF*(xxmax+xxmin);
        if (b1 >= xxmin)
        {
          dataForStack.current = current->_right;
          dataForStack.xxmax = xxmax;
          dataForStack.xxmin = xxmin;
          dataForStack.coupure = coupure2;
          if (dataForStack.current != NULL) stack.push(dataForStack);
        }
            
        break;
            
        // dichotomie direction 2  
      case 1:
          
        yymax = K_CONST::ONE_HALF*(yymax+yymin);
          
        if (a2 <= yymax)
        {
          dataForStack.current = current->_left;
          dataForStack.yymax = yymax;
           dataForStack.coupure = coupure2;
          if (dataForStack.current != NULL) stack.push(dataForStack);
        }
          
        yymax = 2*yymax-yymin;
        yymin = K_CONST::ONE_HALF*(yymax+yymin);
        if (b2 >= yymin)
        {
          dataForStack.current = current->_right;
          dataForStack.yymax = yymax;
          dataForStack.yymin = yymin;
          dataForStack.coupure = coupure2;
          if (dataForStack.current != NULL) stack.push(dataForStack);
        }
            
        break;

        // dichotomie direction 3     
      case 2:
            
        zzmax = K_CONST::ONE_HALF*(zzmax+zzmin);
          
        if (a3 <= zzmax)
        {
          dataForStack.current = current->_left;
          dataForStack.zzmax = zzmax;
           dataForStack.coupure = coupure2;
          if (dataForStack.current != NULL) stack.push(dataForStack);
        }
          
        zzmax = 2*zzmax-zzmin;
        zzmin = K_CONST::ONE_HALF*(zzmax+zzmin);
        if (b3 >= zzmin)
        {
          dataForStack.current = current->_right;
          dataForStack.zzmax = zzmax;
          dataForStack.zzmin = zzmin;
           dataForStack.coupure = coupure2;
          if (dataForStack.current != NULL) stack.push(dataForStack);
        }
          
        break;
          
        // dichotomie direction 4     
      case 3:
          
        ppmax = K_CONST::ONE_HALF*(ppmax+ppmin);
          
        if (a4 <= ppmax)
        {
          dataForStack.current = current->_left;
          dataForStack.ppmax = ppmax;
          dataForStack.coupure = coupure2;
          if (dataForStack.current != NULL) stack.push(dataForStack);
        }
          
        ppmax = 2*ppmax-ppmin;
        ppmin = K_CONST::ONE_HALF*(ppmax+ppmin);
        if (b4 >= ppmin)
        {
          dataForStack.current = current->_right;
          dataForStack.ppmax = ppmax;
          dataForStack.ppmin = ppmin;
          dataForStack.coupure = coupure2;
          if (dataForStack.current != NULL) stack.push(dataForStack);
        }
          
        break;
          
        // dichotomie direction 5     
      case 4:
          
        qqmax = K_CONST::ONE_HALF*(qqmax+qqmin);
          
        if (a5 <= qqmax)
        {
          dataForStack.current = current->_left;
          dataForStack.qqmax = qqmax;
          dataForStack.coupure = coupure2;
          if (dataForStack.current != NULL) stack.push(dataForStack);
        }
          
        qqmax = 2*qqmax-qqmin;
        qqmin = K_CONST::ONE_HALF*(qqmax+qqmin);
        if (b5 >= qqmin)
        {
          dataForStack.current = current->_right;
          dataForStack.qqmax = qqmax;
          dataForStack.qqmin = qqmin;
          dataForStack.rrmin = rrmin;
          dataForStack.coupure = coupure2;
          if (dataForStack.current != NULL) stack.push(dataForStack);
        }
          
        break;
          
        // dichotomie direction 6     
      case 5:
          
        rrmax = K_CONST::ONE_HALF*(rrmax+rrmin);
          
        if (a6 <= rrmax)
        {
          dataForStack.current = current->_left;
          dataForStack.rrmax = rrmax;
          dataForStack.coupure = coupure2;
          if (dataForStack.current != NULL) stack.push(dataForStack);
        }
          
        rrmax = 2*rrmax-rrmin;
        rrmin = K_CONST::ONE_HALF*(rrmax+rrmin);
        if (b6 >= rrmin)
        {
          dataForStack.current = current->_right;
          dataForStack.rrmax = rrmax;
          dataForStack.rrmin = rrmin;
          dataForStack.coupure = coupure2;
          if (dataForStack.current != NULL) stack.push(dataForStack);
        }
          
        break;
    } /* switch */      
  }   
  E_Int size = listOfCandidateCells.size();
  return size;
}

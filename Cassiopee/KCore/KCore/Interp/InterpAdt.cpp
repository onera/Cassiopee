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
# include <stack>
# include "CompGeom/compGeom.h"
# include "Interp/InterpAdt.h"
#include "Loc/loc.h"

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

  if (_cylCoord) { delete [] _xlc; delete [] _ylc; delete [] _zlc; }
}

//=============================================================================
/* Constructor en coord cartesienne */
//=============================================================================
K_INTERP::InterpAdt::InterpAdt(E_Int npts, 
                               E_Float* xD, E_Float* yD, E_Float* zD,
                               void* a1, void* a2, void* a3, E_Int& built):
  InterpData()
{
  _cylCoord = false; 
  _centerX = 0; _centerY = 0; _centerZ = 0;
  _axisX = -1; _axisY = -1; _axisZ = -1;
  _theta_min = K_CONST::E_MAX_FLOAT;
  _theta_max =-K_CONST::E_MAX_FLOAT;

  if (a1 != NULL && a2 != NULL && a3 != NULL) // structure
  {
    _topology = 1;
    E_Int ni = *(E_Int*)a1;
    E_Int nj = *(E_Int*)a2;
    E_Int nk = *(E_Int*)a3;
    built = buildStructAdt(ni, nj, nk, xD, yD, zD);
  }
  else //non structure
  {
    _topology = 2;
    FldArrayI& cEV = *(FldArrayI*)a1;
    built = buildUnstrAdt(npts, cEV, xD, yD, zD);
  }
}

//=============================================================================
/* Constructor en coord cylindrique */
//=============================================================================
K_INTERP::InterpAdt::InterpAdt(E_Int npts, 
                               E_Float* xD, E_Float* yD, E_Float* zD,
                               void* a1, void* a2, void* a3, 
                               E_Float centerX, E_Float centerY, E_Float centerZ,
                               E_Float axisX, E_Float axisY, E_Float axisZ,
                               E_Float thetaShift, E_Int depth,  
                               E_Int& built):
    InterpData()
{
    // keep data for cart2Cyl
    _cylCoord = true;
    _centerX = centerX; _centerY = centerY; _centerZ = centerZ;
    _axisX = axisX; _axisY = axisY; _axisZ = axisZ; _thetaShift = thetaShift;

    E_Float* coordX = new E_Float[npts];
    E_Float* coordY = new E_Float[npts];
    E_Float* coordZ = new E_Float[npts];
    
    E_Float *rt, *thetat;
    thetat = NULL; rt = NULL;
    E_Float eps = 1.e-12;

    if (axisX > eps && axisY < eps && axisZ < eps) // axe X
    {
      rt = coordY; thetat = coordZ;
      for (E_Int i = 0; i < npts; i++) coordX[i] = xD[i]; 
    }
    else if (axisY > eps && axisX < eps && axisZ < eps) // axe Y
    {
      rt = coordZ; thetat = coordX;
      for (E_Int i = 0; i < npts; i++) coordY[i] = yD[i];
    }
    else if (axisZ > eps && axisY < eps && axisX < eps) // axe Z
    {
      rt = coordX; thetat = coordY;
      for (E_Int i = 0; i < npts; i++) coordZ[i] = zD[i];
    }

    if (a1 != NULL && a2 != NULL && a3 != NULL) // structure
    {
        _topology = 1;
        E_Int ni = *(E_Int*)a1;
        E_Int nj = *(E_Int*)a2;
        E_Int nk = *(E_Int*)a3;

        // cart2Cyl coordinates
        K_LOC::cart2Cyl(npts, xD, yD, zD,
                        centerX, centerY, centerZ, 
                        axisX, axisY, axisZ, 
                        rt, thetat, ni, nj, nk, depth, thetaShift=_thetaShift);

        /*
        E_Float thetarefmin = K_CONST::E_MAX_FLOAT;
        E_Float thetarefmax =-K_CONST::E_MAX_FLOAT;
        for (E_Int i = 0; i < npts; i++)
        {
          thetarefmin = K_FUNC::E_min(thetarefmin, thetat[i]);
          thetarefmax = K_FUNC::E_max(thetarefmax, thetat[i]);
        }
        _theta_min = thetarefmin;
        _theta_max = thetarefmax;
        // printf(" thetamin = %g %g \n", _theta_min, _theta_max);

        */
        // printf("ni=%d %d %d\n",*(E_Int*)a1,*(E_Int*)a2,*(E_Int*)a3);
        // for (E_Int i = 0; i < npts; i++) printf("%g %g\n", thetat[i], coordZ[i]);
        // printf("BlkInterpAdt: axis = %g %g %g\n", axisX, axisY, axisZ); fflush(stdout);        

        built = buildStructAdt(ni, nj, nk, coordX, coordY, coordZ);
    }
    else //non structure
    {
        _topology = 2;
        FldArrayI& cEV = *(FldArrayI*)a1;

        // cart2Cyl coordinates
        K_LOC::cart2Cyl(npts, xD, yD, zD,
                        centerX, centerY, centerZ, 
                        axisX, axisY, axisZ, 
                        rt, thetat, thetaShift=_thetaShift);

        built = buildUnstrAdt(npts, cEV, coordX, coordY, coordZ);
    }
    _xlc = coordX; _ylc = coordY; _zlc = coordZ;
    //delete [] coordX; delete [] coordY; delete [] coordZ;
}

//=============================================================================
// passe le vecteur de points fourni en cylindrique
void K_INTERP::InterpAdt::cart2Cyl(E_Int npts, E_Float* x, E_Float* y, E_Float* z)
{
  E_Float* coordX = new E_Float[npts];
  E_Float* coordY = new E_Float[npts];
  E_Float* coordZ = new E_Float[npts];
    
  // tetaShift x,y,z
  E_Float* xR=NULL; E_Float* yR=NULL; E_Float* zR=NULL;
  E_Float *rt, *thetat;
  thetat = NULL; rt = NULL;
  E_Float eps = 1.e-12;
  if (_axisX > eps && _axisY < eps && _axisZ < eps) // axe X
  {
    rt = coordY; thetat = coordZ;
  }
  else if (_axisY > eps && _axisX < eps && _axisZ < eps) // axe Y
  {
    rt = coordZ; thetat = coordX;
  }
  else if (_axisZ > eps && _axisY < eps && _axisX < eps) // axe Z
  {
    rt = coordX; thetat = coordY;
  }
  // cart2Cyl coordinates
  E_Int nit = 0; E_Int njt = 0; E_Int nkt = 0;
  K_LOC::cart2Cyl(npts, x, y, z,
                  _centerX, _centerY, _centerZ, 
                  _axisX, _axisY, _axisZ, 
                  rt, thetat, nit, njt, nkt, _thetaShift);
  /*
  E_Float PI2 = 2.*K_CONST::E_PI;
#pragma omp parallel default(shared)
  {    
#pragma omp for         
    for (E_Int i = 0; i < npts; i++)
    {
      if ( thetat[i] < _theta_min) thetat[i] += PI2;
      else if (thetat[i] > _theta_max) thetat[i] -= PI2;
    }
  }
  */

  delete [] coordX; delete [] coordY; delete [] coordZ;
  if (_thetaShift != 0.) { delete [] xR; delete [] yR; delete [] zR; }
}
// Passe le pt fourni en cylindrique
void K_INTERP::InterpAdt::cart2Cyl(E_Float& x, E_Float& y, E_Float& z)
{
    E_Float eps = 1.e-12;
    E_Float* rt=NULL; E_Float* thetat=NULL;
    E_Float Xo, Yo, Zo;

    if (_axisX > eps && _axisY < eps && _axisZ < eps) // axe X
    {
      rt = &Yo; thetat = &Zo; Xo = x;
    }
    else if (_axisY > eps && _axisX < eps && _axisZ < eps) // axe Y
    {
      rt = &Zo; thetat = &Xo; Yo = y;
    }
    else if (_axisZ > eps && _axisY < eps && _axisX < eps) // axe Z
    {
      rt = &Xo; thetat = &Yo; Zo = z;
    }

    // cart2Cyl coordinates
    E_Int nit = 0; E_Int njt = 0; E_Int nkt = 0; E_Int deptht = 0;
    K_LOC::cart2Cyl(1, &x, &y, &z,
                    _centerX, _centerY, _centerZ, 
                    _axisX, _axisY, _axisZ, 
                    rt, thetat, nit, njt, nkt, deptht, _thetaShift);
    /*
    E_Float PI2 = 2.*K_CONST::E_PI;
    if ( thetat[0] < _theta_min) thetat[0] += PI2;
    else if (thetat[0] > _theta_max) thetat[0] -= PI2; */

    x = Xo; y = Yo; z = Zo;
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
   IN: x,y,z: coord du pt a interpoler ou a extrapoler
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
// ============================================================================
/* Find the interpolation cell for point (x,y,z) for a tetra kmesh
   le numero de l'elt trouve est retourne dans noelt */
// ============================================================================
short 
K_INTERP::InterpAdt::searchExtrapolationCellUnstruct(E_Float* xt, E_Float* yt, E_Float* zt,
                                                     E_Float* cellNp,
                                                     FldArrayI& connectEV,
                                                     E_Float x, E_Float y, E_Float z,
                                                     E_Int& noelt, FldArrayF& cf,
                                                     E_Int nature, E_Int extrapOrder, E_Float constraint)
{
  // cylindrical coords modification
  if (_cylCoord) cart2Cyl(x,y,z);
  
  // datas for interpolation cell (tetrahedra)
  E_Float xp, yp, zp;
  E_Float xq, yq, zq;
  E_Float xr, yr, zr;
  E_Float xs, ys, zs;
  E_Int indp, indq, indr, inds;
  E_Float xi, yi, zi;

  E_Int foundInDomain = 0;
  E_Int noeltsave = 0;
  E_Int ncf = 8;
  FldArrayF cf_sav(ncf);

  // Init cf
  cf.setAllValuesAtNull();

  // search list of candidate cells
  list<E_Int> listOfCandidateCells;
  E_Float alphaTol = K_FUNC::E_max( (2*constraint-5.)*0.1, 0.);
  E_Int found = getListOfCandidateCells(x, y, z, listOfCandidateCells, alphaTol);
  if (found == 0) return 0; // listOfCandidateCells empty

  /* Find the right cell among candidate cells */
  list<E_Int>::iterator itr;
  E_Float sum;
  E_Float sum_coef;
  E_Float saved_sum_coef = K_CONST::E_MAX_FLOAT;
  E_Float max_coef;
  E_Int* cn1 = connectEV.begin(1);
  E_Int* cn2 = connectEV.begin(2);
  E_Int* cn3 = connectEV.begin(3);
  E_Int* cn4 = connectEV.begin(4);

  // Loop on all candidate cells : computation of coefficients cf
  for (itr = listOfCandidateCells.begin();
       itr != listOfCandidateCells.end();
       itr++)
  {
    // index of candidate cell
    noelt = *itr;
    
    // compute interpolation coefficients for tetrahedra
    indp = cn1[noelt]-1;
    indq = cn2[noelt]-1;
    indr = cn3[noelt]-1;
    inds = cn4[noelt]-1;
    
    xp = xt[indp]; yp = yt[indp]; zp = zt[indp];
    xq = xt[indq]; yq = yt[indq]; zq = zt[indq];
    xr = xt[indr]; yr = yt[indr]; zr = zt[indr];
    xs = xt[inds]; ys = yt[inds]; zs = zt[inds];
    coeffInterpTetra(x, y, z, xp, yp, zp, xq, yq, zq,
                     xr, yr, zr, xs, ys, zs, xi, yi, zi);
    
    // cellN for points of interpolation cell
    E_Float cellNmp = 1.;
    E_Float cellNmq = 1.;
    E_Float cellNmr = 1.;
    E_Float cellNms = 1.;
    
    if (cellNp != NULL)
    {
      cellNmp = cellNp[indp];
      cellNmq = cellNp[indq];
      cellNmr = cellNp[indr];
      cellNms = cellNp[inds];
    }
    // extrapolation - order 0
    cf.setAllValuesAtNull();
    if (extrapOrder == 0)
    {
      for (E_Int i = 0; i < cf.getSize(); i++) cf[i]=0.;
      // test which point of the tetrahedron will be used for extrapolation
      max_coef = K_FUNC::E_max(xi,yi); max_coef = K_FUNC::E_max(max_coef,zi); 
      if (max_coef <= 0.5) cf[0] = 1.; // point P
      else if (K_FUNC::fEqualZero(max_coef,xi)) cf[1] = 1.; // point Q
      else if (K_FUNC::fEqualZero(max_coef,yi)) cf[2] = 1.; // point R
      else cf[3] = 1.; // point S         
    }
    else if (extrapOrder == 1) // extrapolation - order 1
    {
      cf[0] = 1-xi-yi-zi; 
      cf[1] = xi;
      cf[2] = yi;
      cf[3] = zi;
      sum_coef = K_FUNC::E_abs(cf[1])+K_FUNC::E_abs(cf[2])+K_FUNC::E_abs(cf[3]);
      // if sum_coef exceeds constraint, degenerate to order 0
      if (sum_coef >= constraint)
      {
        for (E_Int i = 0; i < cf.getSize(); i++) cf[i] = 0.;
        // test will point of the tetrahedra will be used for extrapolation
        max_coef = K_FUNC::E_max(xi,yi); max_coef = K_FUNC::E_max(max_coef,zi); 
        if (max_coef <= 0.5) cf[0] = 1.; // point P
        else if (K_FUNC::fEqualZero(max_coef,xi)) cf[1] = 1.; // point Q
        else if (K_FUNC::fEqualZero(max_coef,yi)) cf[2] = 1.; // point R
        else cf[3] = 1.; // point S         
      }
    }
    // Keep only valid extrapolation cell (no blanked cell in tetrahedron)
    sum = cellNmp*cellNmq*cellNmr*cellNms;
    if (sum > K_CONST::E_CUTOFF)
    {
      // Keep extrapolation cell with minimum distance (~ minimum sum of absolute coefficients) 
      // from the interpolated point
      sum_coef = K_FUNC::E_abs(cf[1])+K_FUNC::E_abs(cf[2])+K_FUNC::E_abs(cf[3]);
      if (sum_coef < saved_sum_coef)
      {
        if (sum_coef < constraint) foundInDomain = 1;
        else foundInDomain = 0;
        saved_sum_coef = sum_coef;
        noeltsave = noelt;
        cf_sav = cf;
      }
    }
  }

  noelt = noeltsave;
  cf = cf_sav;
  return foundInDomain;
}
// ============================================================================
/* Find the interpolation cell for point (x,y,z) in case of a structured donor
   zone */
// ============================================================================
short K_INTERP::InterpAdt::searchExtrapolationCellStruct(
  E_Int ni, E_Int nj, E_Int nk, 
  E_Float* xl, E_Float* yl, E_Float* zl,
  E_Float* cellNp,
  E_Float x, E_Float y, E_Float z,
  E_Int& ic, E_Int& jc, E_Int& kc,
  FldArrayF& cf,
  E_Int nature, E_Int extrapOrder, E_Float constraint)
{
  // cylindrical coords modification
  if (_cylCoord) { cart2Cyl(x,y,z); xl = _xlc; yl = _ylc; zl = _zlc; }
  
  //E_Int dim = 3;
  //if (nk == 1) dim = 2;
  E_Int i,j,k;
  E_Int icsav = 0;
  E_Int jcsav = 0;
  E_Int kcsav = 0;
  E_Int foundInDomain = 0;
  E_Int ncf = 8;
  FldArrayF cf_sav(ncf);
  E_Float* cfp = cf.begin();

  // Init cf
  cf.setAllValuesAtNull();

  // search list of candidate cells
  list<E_Int> listOfCandidateCells;
  E_Float alphaTol = K_FUNC::E_max( (2*constraint-5.)*0.1, 0.);
  E_Int found = getListOfCandidateCells(x,y,z,listOfCandidateCells, alphaTol);
  
  if (found == 0) return 0; // listOfCandidateCells empty

  /* Find the right cell among candidate cells */
  list<E_Int>::iterator itr;
  E_Float sum_coef;
  E_Float saved_sum_coef = K_CONST::E_MAX_FLOAT;
  E_Float saved_max_diff_coef = K_CONST::E_MAX_FLOAT;
  E_Float diff_coeff = K_CONST::E_MAX_FLOAT;
  E_Int ind;
  // parameters for extrapolation routine
  E_Int is, js, ks, ret;
  E_Int nij;

  // Loop on all candidate cells: computation of coefficients cf
  for (itr = listOfCandidateCells.begin();
       itr != listOfCandidateCells.end();
       itr++)
  {
    // 1D-index of interpolation cell
    ind = *itr;
 
    // (i,j,k)-indices of interpolation cell
    nij = ni*nj;
    k = ind/nij; 
    j = (ind-k*nij)/ni;
    i = ind-j*ni-k*nij;
    k++; j++; i++;

    // compute interpolation coefficients for hexahedra
    // (is, js, ks): neighbour cell (not used here)
    ret = getExtrapolationCoeffForCell(x, y, z, i, j, k, cf, ni, nj, nk, 
                                       xl, yl, zl, cellNp, is, js, ks, 
                                       nature, constraint, diff_coeff);

    // Keep extrapolation cell with minimum distance (~ minimum sum of absolute coefficients)
    // from the interpolated point
    sum_coef = K_FUNC::E_abs(cfp[0])+K_FUNC::E_abs(cfp[1])+K_FUNC::E_abs(cfp[2])+K_FUNC::E_abs(cfp[3])
    +K_FUNC::E_abs(cfp[4])+K_FUNC::E_abs(cfp[5])+K_FUNC::E_abs(cfp[6])+K_FUNC::E_abs(cfp[7]);

    if (ret == 1 && sum_coef < saved_sum_coef+1.e-6 && diff_coeff< saved_max_diff_coef)
    {
      saved_max_diff_coef = diff_coeff;
      foundInDomain = 1;
      saved_sum_coef = sum_coef;
      icsav = i; jcsav = j; kcsav = k;
      cf_sav = cf;
    }
  }
  ic = icsav; jc = jcsav; kc = kcsav;
  cf = cf_sav;

  if (extrapOrder == 0) // passage a l'ordre 0
  {
    // On reconstruit un xi,eta,zeta comme si les coeff avaient
    // ete obtenus en tri-lineaire
    //printf("xyz: %f %f %f\n", x, y, z);
    //printf("%f %f %f %f %f %f %f %f\n", cfp[0],cfp[1],cfp[2],cfp[3],cfp[4],cfp[5],cfp[6],cfp[7]);
    E_Float xi, eta, zeta;
    // eval, xi, eta, zeta
    xi = cfp[1]+cfp[3]+cfp[5]+cfp[7];
    eta = cfp[2]+cfp[3]+cfp[6]+cfp[7];
    zeta = cfp[4]+cfp[5]+cfp[6]+cfp[7];

    //printf("xi: %f %f %f\n", xi,eta,zeta);
    if (xi < 0) xi = 0.;
    if (xi > 1) xi = 1.;
    if (eta < 0) eta = 0.;
    if (eta > 1) eta = 1.;
    if (zeta < 0) zeta = 0.;
    if (zeta > 1) zeta = 1.;

    cfp[0] = (1-xi)*(1-eta)*(1-zeta);
    cfp[1] = xi*(1-eta)*(1-zeta);
    cfp[2] = (1-xi)*eta*(1-zeta);
    cfp[3] = xi*eta*(1-zeta);
    cfp[4] = (1-xi)*(1-eta)*zeta;
    cfp[5] = xi*(1-eta)*zeta;
    cfp[6] = (1-xi)*eta*zeta;
    cfp[7] = xi*eta*zeta;
  }

  return foundInDomain;
}
// ============================================================================
/* Find the interpolation cell for point (x,y,z) for a tetra kmesh
   le numero de l'elt trouve est retourne dans noelt */
// ============================================================================
short 
K_INTERP::InterpAdt::searchInterpolationCellUnstruct(E_Float* xt, E_Float* yt, E_Float* zt,
                                                     FldArrayI& connectEV,
                                                     E_Float x, E_Float y, E_Float z,
                                                     E_Int& noelt, FldArrayF& cf)
{ 
  // cylindrical coords modification
  if (_cylCoord) { cart2Cyl(x,y,z); xt = _xlc; yt = _ylc; zt = _zlc; }
  
  cf.setAllValuesAtNull();

  // recherche de la liste des cellules candidates
  list<E_Int> listOfCandidateCells; 
  E_Int found = getListOfCandidateCells(x, y, z, listOfCandidateCells);
  if (found == 0) return 0; // listOfCandidateCells vide

  const E_Float EPS = _EPS_TETRA;
  E_Float xp, yp, zp;
  E_Float xq, yq, zq;
  E_Float xr, yr, zr;
  E_Float xs, ys, zs;
  E_Int indp, indq, indr, inds;
  E_Float xi, yi, zi, sum;
  E_Int et;
  list<E_Int>::iterator itr;
  E_Int* cn1 = connectEV.begin(1);
  E_Int* cn2 = connectEV.begin(2);
  E_Int* cn3 = connectEV.begin(3);
  E_Int* cn4 = connectEV.begin(4);

  for (itr = listOfCandidateCells.begin(); 
       itr != listOfCandidateCells.end();
       itr++)
  {
    et = *itr;
    indp = cn1[et]-1;
    indq = cn2[et]-1;
    indr = cn3[et]-1;
    inds = cn4[et]-1;

    xp = xt[indp]; yp = yt[indp]; zp = zt[indp];
    xq = xt[indq]; yq = yt[indq]; zq = zt[indq];
    xr = xt[indr]; yr = yt[indr]; zr = zt[indr];
    xs = xt[inds]; ys = yt[inds]; zs = zt[inds];

    coeffInterpTetra(x, y, z, xp, yp, zp, xq, yq, zq,
                     xr, yr, zr, xs, ys, zs, xi, yi, zi);

    sum = xi+yi+zi;
    if (xi > -EPS && yi > -EPS && zi > -EPS && sum < K_CONST::ONE+3*EPS)
    {
      cf[0] = 1-sum; 
      cf[1] = xi; cf[2] = yi; cf[3] = zi;
      noelt = et;
      return 1;
    }
  }
  noelt = -1;
  return 0;
}

// ============================================================================
/* Find the interpolation cell for point (x,y,z) in case of a structured donor
   zone */
// ============================================================================
short K_INTERP::InterpAdt::searchInterpolationCellStruct(
  E_Int ni, E_Int nj, E_Int nk, 
  E_Float* xl, E_Float* yl, E_Float* zl,
  E_Float x, E_Float y, E_Float z,
  E_Int& ic, E_Int& jc, E_Int& kc,
  FldArrayF& cf)
{
  // cylindrical coords modification
  //printf("1. %g %g %g\n",x,y,z); fflush(stdout);
  if (_cylCoord) { cart2Cyl(x,y,z); xl = _xlc; yl = _ylc; zl = _zlc; }
  //printf("2. %g %g %g\n",x,y,z); fflush(stdout);

  // Sortie DBX pour orphans
  /*
  if (_cylCoord)
  {
    FILE* ptrFile = fopen("mesh.dat", "w");
    fprintf(ptrFile, "TITLE=toto\n");
    fprintf(ptrFile, "ZONE NI=%d NJ=%d NK=%d\n", ni,nj,nk);
    fprintf(ptrFile, "VARIABLES = x,y,z\n");
    
    //for (E_Int i = 0; i < ni*nj*nk; i++)
    //{
    //  fprintf(ptrFile, "%g %g %g\n", xl[i], yl[i], zl[i]);
    //}
    fclose(ptrFile);
    exit(0);
  }
  */
  // END DBX


  // recherche de la liste des cellules candidates
  list<E_Int> listOfCandidateCells; 
  E_Int found = getListOfCandidateCells(x, y, z, listOfCandidateCells);
  /* Find the right cell among candidate cells */
  
  if (found == 0) return 0; // listOfCandidateCells vide
  
  list<E_Int>::iterator itr;
  E_Bool JUMP;
  E_Int ind, i, isomm, is, js, ks;
  E_Float xi, yi, zi;
  E_Float xt[15], yt[15], zt[15];
  E_Int nij = ni*nj;
  
  /* For 8 cells, apply the technique of jump to find interpolation cell */
  JUMP = false;
  itr = listOfCandidateCells.begin();
  ind = *itr;
  for (i = 0; i < 8; i++)
  { 
    coordHexa(ind, ni, nj, nk,
              xl, yl, zl,
              ic, jc, kc,
              xt, yt, zt);
    
    if (getCellJump(x, y, z,
                    xt, yt, zt,              
                    isomm,
                    xi, yi, zi) == false)
    {
      if (i == 7)
      {
        JUMP =  false;
        break;
      }
      /* Apply a technique of jump */
      is = K_FUNC::E_min(8, E_Int((xi+K_FUNC::E_sign(xi))*K_CONST::ONE_HALF));
      is = K_FUNC::E_max(-8, is);
      js = K_FUNC::E_min(8, E_Int((yi+K_FUNC::E_sign(yi))*K_CONST::ONE_HALF));
      js = K_FUNC::E_max(-8, js);
      ks = K_FUNC::E_min(8, E_Int((zi+K_FUNC::E_sign(zi))*K_CONST::ONE_HALF));
      ks = K_FUNC::E_max(-8, ks);
      ind = ind+is+js*ni+ks*nij;

      kc = ind/nij;
      E_Int kcnij = kc*nij;
      jc = (ind-kcnij)/ni;
      E_Int jcni = jc*ni;
      ic = ind-jcni-kcnij;

      if (ic<0 || ic>ni-2 || jc<0 || jc>nj-2 || kc<0 || kc>nk-2)
      {
        JUMP =  false;
        break;
      }
    }
    else
    {        
      for (itr = listOfCandidateCells.begin();
           itr != listOfCandidateCells.end();
           itr++)
      {
        if (ind == *itr)
        {
          JUMP = true;
          goto saut;
        }
      }
      JUMP =  false;
      break;
    }
  }

  /* If the technique of jump succeeds, find the interpolation coefficients
     in the cell by cut it in 24 tetrahedras */
  saut:
  if (JUMP)  
  {
    if (getCoeffInterpHexa(x, y, z,
                          isomm,
                          xi, yi, zi, 
                          xt, yt, zt,
                          cf) == true)
    {      
      return 1;
    }
    else
    {
      JUMP =  false;
      listOfCandidateCells.erase(itr);
    }
  }
 
  /* If the technique of jump fails, we test all candidate cells */
  if (JUMP ==  false) 
  {    
    for (itr = listOfCandidateCells.begin();
         itr != listOfCandidateCells.end();
         itr++)
    {
      ind = (*itr);
          
      coordHexa(ind, ni, nj, nk,
                xl, yl, zl,
                ic, jc, kc,
                xt, yt, zt);
      if (coeffInterpHexa(x, y, z,
                          xt, yt, zt,
                          cf) == true)
      {
        return 1;
      }
    }
  }

  return 0;
}

short K_INTERP::InterpAdt::searchInterpolationCellCartO2(E_Int ni, E_Int nj, E_Int nk,
                                                         E_Float x, E_Float y, E_Float z,
                                                         E_Int& ic, E_Int& jc, E_Int& kc,
                                                         FldArrayF& cf)
{ 
  return -1;
} 

short K_INTERP::InterpAdt::searchInterpolationCellCartO3(E_Int ni, E_Int nj, E_Int nk,
                                                         E_Float x, E_Float y, E_Float z,
                                                         E_Int& icHO, E_Int& jcHO, E_Int& kcHO,
                                                         FldArrayF& cf)
{return -1;}

short K_INTERP::InterpAdt::searchInterpolationCellCartO4(E_Int ni, E_Int nj, E_Int nk,
                                                         E_Float x, E_Float y, E_Float z,
                                                         E_Int& icHO, E_Int& jcHO, E_Int& kcHO,
                                                         FldArrayF& cf)
{return -1;}


short K_INTERP::InterpAdt::searchExtrapolationCellCart(E_Int ni, E_Int nj, E_Int nk, 
                                                       E_Float* xl, E_Float* yl, E_Float* zl,
                                                       E_Float* cellNp,
                                                       E_Float x, E_Float y, E_Float z,
                                                       E_Int& ic, E_Int& jc, E_Int& kc,
                                                       FldArrayF& cf,
                                                       E_Int nature, E_Int extrapOrder, E_Float constraint)
{return -1;}  

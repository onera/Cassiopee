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

#include "BlkInterp.h"
#include <stack>
using namespace std;
using namespace K_FLD;
using namespace K_FUNC;

struct stackData
{
    K_KINTERP::BlkIntTreeNode* current;
    E_Float xxmax, yymax, zzmax, ppmax, qqmax, rrmax;
    E_Float xxmin, yymin, zzmin, ppmin, qqmin, rrmin;
    E_Int coupure;
};

//=============================================================================
/* Destructor */
//=============================================================================
K_KINTERP::BlkInterpAdt::~BlkInterpAdt()
{
  destroy();
}

//=============================================================================
/* Destroy data structure */
//=============================================================================
void K_KINTERP::BlkInterpAdt::destroy()
{
  // delete the ADT tree
  BlkIntTreeNode* current = _tree;
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
K_KINTERP::BlkInterpAdt::BlkInterpAdt(KMesh& mesh) :
  BlkInterpWithKMesh(mesh)
{
  buildAdt();
}


//=============================================================================
/* Build the ADT */
//=============================================================================
void K_KINTERP::BlkInterpAdt::buildAdt()
{
  if (_mesh.isStructured() == true) buildStructAdt();
  else buildUnstructAdt();
}
//=============================================================================
/* Construction de l adt pour un maillage structure */
//=============================================================================
void K_KINTERP::BlkInterpAdt::buildStructAdt()
{
  KMesh& mesh = _mesh;
  
  _tree = NULL;

  /* Find the bounding box of mesh */
  mesh.boundingBox(_xmax, _ymax, _zmax,
                   _xmin, _ymin, _zmin);

  // Mesh characteristics
  E_Int ni = mesh.getIm();
  E_Int nj = mesh.getJm();
  E_Int nk = mesh.getKm();
  E_Int nij = ni*nj;
  E_Float* x = mesh.getXVector();
  E_Float* y = mesh.getYVector();
  E_Float* z = mesh.getZVector();
  
  /* Insert all points (cells) in the tree */
  E_Float xmax,ymax,zmax,xmin,ymin,zmin;
  E_Int ind, ind2;
  
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
    
        xmin = E_min(xmin, x[ind]);
        ymin = E_min(ymin, y[ind]);
        zmin = E_min(zmin, z[ind]);
        xmax = E_max(xmax, x[ind]);
        ymax = E_max(ymax, y[ind]);
        zmax = E_max(zmax, z[ind]);
        
        /* Compute bounding box of cell */
        ind2 = ind+1;
        xmin = E_min(xmin, x[ind2]);
        ymin = E_min(ymin, y[ind2]);
        zmin = E_min(zmin, z[ind2]);
        xmax = E_max(xmax, x[ind2]);
        ymax = E_max(ymax, y[ind2]);
        zmax = E_max(zmax, z[ind2]);
            
        ind2 = ind+ni;
        xmin = E_min(xmin, x[ind2]);
        ymin = E_min(ymin, y[ind2]);
        zmin = E_min(zmin, z[ind2]);
        xmax = E_max(xmax, x[ind2]);
        ymax = E_max(ymax, y[ind2]);
        zmax = E_max(zmax, z[ind2]);

        ind2 = ind+nij;
        xmin = E_min(xmin, x[ind2]);
        ymin = E_min(ymin, y[ind2]);
        zmin = E_min(zmin, z[ind2]);
        xmax = E_max(xmax, x[ind2]);
        ymax = E_max(ymax, y[ind2]);
        zmax = E_max(zmax, z[ind2]);
        
        ind2 = ind+ni+1;
        xmin = E_min(xmin, x[ind2]);
        ymin = E_min(ymin, y[ind2]);
        zmin = E_min(zmin, z[ind2]);
        xmax = E_max(xmax, x[ind2]);
        ymax = E_max(ymax, y[ind2]);
        zmax = E_max(zmax, z[ind2]);
        
        ind2 = ind+nij+1;
        xmin = E_min(xmin, x[ind2]);
        ymin = E_min(ymin, y[ind2]);
        zmin = E_min(zmin, z[ind2]);
        xmax = E_max(xmax, x[ind2]);
        ymax = E_max(ymax, y[ind2]);
        zmax = E_max(zmax, z[ind2]);
        
        ind2 = ind+ni+nij;
        xmin = E_min(xmin, x[ind2]);
        ymin = E_min(ymin, y[ind2]);
        zmin = E_min(zmin, z[ind2]);
        xmax = E_max(xmax, x[ind2]);
        ymax = E_max(ymax, y[ind2]);
        zmax = E_max(zmax, z[ind2]);

        ind2 = ind2+1;
        xmin = E_min(xmin, x[ind2]);
        ymin = E_min(ymin, y[ind2]);
        zmin = E_min(zmin, z[ind2]);
        xmax = E_max(xmax, x[ind2]);
        ymax = E_max(ymax, y[ind2]);
        zmax = E_max(zmax, z[ind2]);

        /* insert this cell in ADT tree */
        insert(ind, xmin, ymin, zmin,
               xmax, ymax, zmax);
      }
}
//=============================================================================
/* Construction de l'adt pour un maillage non structure (TETRA) */
//=============================================================================
void K_KINTERP::BlkInterpAdt::buildUnstructAdt()
{
  KMesh& mesh = const_cast<KMesh&>(_mesh);
  
  _tree = NULL;

  /* Find the bounding box of mesh */
  mesh.boundingBox(_xmax, _ymax, _zmax,
                   _xmin, _ymin, _zmin);

  // Mesh characteristics
  E_Float* x = mesh.getXVector();
  E_Float* y = mesh.getYVector();
  E_Float* z = mesh.getZVector();
  
  /* Insert all points (cells) in the tree */
  E_Float xmax,ymax,zmax,xmin,ymin,zmin;
  
  FldArrayI& connect = mesh.getConnectivity();
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
    xmin = E_min(xmin, x[ind]);
    ymin = E_min(ymin, y[ind]);
    zmin = E_min(zmin, z[ind]);
    xmax = E_max(xmax, x[ind]);
    ymax = E_max(ymax, y[ind]);
    zmax = E_max(zmax, z[ind]);
    
    ind2 = c2[et]-1;
    xmin = E_min(xmin, x[ind2]);
    ymin = E_min(ymin, y[ind2]);
    zmin = E_min(zmin, z[ind2]);
    xmax = E_max(xmax, x[ind2]);
    ymax = E_max(ymax, y[ind2]);
    zmax = E_max(zmax, z[ind2]);
    
    ind2 = c3[et]-1;
    xmin = E_min(xmin, x[ind2]);
    ymin = E_min(ymin, y[ind2]);
    zmin = E_min(zmin, z[ind2]);
    xmax = E_max(xmax, x[ind2]);
    ymax = E_max(ymax, y[ind2]);
    zmax = E_max(zmax, z[ind2]);
    
    ind2 = c4[et]-1;
    xmin = E_min(xmin, x[ind2]);
    ymin = E_min(ymin, y[ind2]);
    zmin = E_min(zmin, z[ind2]);
    xmax = E_max(xmax, x[ind2]);
    ymax = E_max(ymax, y[ind2]);
    zmax = E_max(zmax, z[ind2]);
        
    /* insert this cell in ADT tree */
    insert(et, xmin, ymin, zmin, xmax, ymax, zmax);
  }
}

//=============================================================================
// Insert cell in ADT tree
//=============================================================================
void K_KINTERP::BlkInterpAdt::insert(E_Int ind, 
                                     E_Float xmin, E_Float ymin, E_Float zmin,
                                     E_Float xmax, E_Float ymax, E_Float zmax)
{
  E_Int coupure = 0;

  E_Float xxmax,yymax,zzmax,xxmin,yymin,zzmin;
  E_Float ppmax,qqmax,rrmax,ppmin,qqmin,rrmin;
  
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
  
  BlkIntTreeNode** current = &_tree;

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
  
  *current = new BlkIntTreeNode(ind, xmax, ymax, zmax,
                                xmin, ymin, zmin);
}

// ==================== Interp/BlkInterpAdt.cpp ===========================
